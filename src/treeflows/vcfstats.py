#!/bin/python3

#import vcf2est as est
import gzip
import cyvcf2
from itertools import compress
import numpy as np
from statistics import mean, median
from treeflows import vcf_core as _vcf
from pathlib import Path
from scipy.stats import norm, poisson, nbinom
from concurrent.futures import ProcessPoolExecutor, as_completed
from treeflows.gff import GenomeMask


def get_depth_global(infix):
    """Print per-site depth (DP) for the first few variants in a VCF prefix."""
    vcf = cyvcf2.VCF(infix + '.vcf.gz')

    for n, variant in enumerate(vcf):
        print(f'{variant.CHROM}\t{variant.start}\t{variant.INFO["DP"]}')
        if n > 10:
            break

def get_depth_global_window_fixed(infix, outfix, window=100, missthresh = 3): # process is curently extremely naive!!!!! 
    """Compute sliding-window mean depth/missingness and write to `<outfix>.depths`."""
    # can possibly make a bp slider... that needs to be in the middle of the window? 
    vcf = cyvcf2.VCF(infix + '.vcf.gz')
    nsamps = len(vcf.samples)
    depths = np.zeros((window, nsamps))
    missness = np.zeros(window)
    slider_bp = np.zeros(window)
    with open(outfix + '.depths', "w")as ofi:
        ofi.write("position\tmean_depth_individual\tnonmiss_prop\tsnp_density_100\n")
    #print(slider)

    #minpos & maxpos. 

    variter = iter(vcf)
    variant = next(variter)

    dpt = variant.gt_depths
    depths = np.reshape(dpt, newshape=(1, -1))
    slider_bp = np.full(1, variant.start)

    missness = np.full(1, np.mean(np.greater_equal(dpt, missthresh)))


    csite = variant.start
    cmin = csite - window
    cmid = csite - (window / 2)
    reset = False

    while True: # i.e. the entire loop... (switch back and forth)... 

        # load up a new site
        try: 
            variant = next(variter)
            vsite = variant.start
            #slider_bp = np.concatenate((slider_bp, np.full(1, vsite))) # extend the snp window to encompass the new value... 
            # if variant.start > minpos and variant.start < maxpos:
            #     print(variant.start)
        except StopIteration:
            variant = None
            break

        # now check in site numbers
        while csite < vsite:
            #print(len(slider_bp))

            csite += 1
            cmin += 1
            cmid += 1
            if slider_bp[0] < cmin and len(slider_bp) > 1:
                slider_bp = slider_bp[1:] # drop values at the far end... 
                depths = depths[1:]
                missness = missness[1:]
                with open(outfix + '.depths', 'a') as ofi:
                    ofi.write(f'{median(slider_bp)}\t{np.mean(depths)}\t{np.mean(missness)}\t{len(slider_bp) / window * 100}\n')
            elif slider_bp[0] < cmin and len(slider_bp) <= 1:
                reset = True
                csite = vsite
                cmin = csite - window
                cmid = csite - (window / 2)
                break

        if reset:
            slider_bp = np.full(1, variant.start)
            dpt = variant.gt_depths
            depths = np.reshape(dpt, newshape=(1, -1))
            missness = np.full(1, np.mean(np.greater_equal(dpt, missthresh)))
        if not reset:
            slider_bp = np.concatenate((slider_bp, np.full(1, variant.start)))
            dpt = variant.gt_depths
            depths = np.concatenate((depths, np.reshape(dpt, newshape=(1, -1))))
            missness = np.concatenate((missness, np.full(1, np.mean(np.greater_equal(dpt, missthresh)))))
        reset = False
        with open(outfix + '.depths', 'a') as ofi:
            ofi.write(f'{median(slider_bp)}\t{np.mean(depths)}\t{np.mean(missness)}\t{len(slider_bp) / window * 100}\n')


def init_depth_global_window_fixed(infix, window, missthresh): # process is curently extremely naive!!!!! 
    """Initialize arrays/state for windowed depth/missingness tracking."""
    # can possibly make a bp slider... that needs to be in the middle of the window? 
    vcf = cyvcf2.VCF(infix + '.vcf.gz')
    nsamps = len(vcf.samples)
    depths = np.zeros((window, nsamps))
    missness = np.zeros(window)
    slider_bp = np.zeros(window)

    variter = iter(vcf)
    variant = next(variter)

    dpt = variant.gt_depths
    depths = np.reshape(dpt, newshape=(1, -1))
    slider_bp = np.full(1, variant.start)

    missness = np.full(1, np.mean(np.greater_equal(dpt, missthresh)))
    snps = np.full(1, variant.is_snp)

    csite = 0 # set to 1... 
    cmin = csite - window
    cmid = csite - (window / 2)
    counter = 0

    return slider_bp, depths, snps, missness, csite, cmin, cmid, counter

def build_depth_global_window_fixed(infix, outfix, window, missthresh, slider_bp, depths, missness, snps, csite, cmin, cmid, thin, counter):
    """Update windowed depth/missingness stats over a VCF and append results."""
    vcf = cyvcf2.VCF(infix + '.vcf.gz')
    variter = iter(vcf)
    reset = False

    while True: # i.e. the entire loop... (switch back and forth)... 

        # load up a new site
        try: 
            variant = next(variter)
            vsite = variant.start
            #slider_bp = np.concatenate((slider_bp, np.full(1, vsite))) # extend the snp window to encompass the new value... 
            # if variant.start > minpos and variant.start < maxpos:
            #     print(variant.start)
        except StopIteration:
            variant = None
            break

        # now check in site numbers
        while csite < vsite:
            #print(len(slider_bp))

            csite += 1
            cmin += 1
            cmid += 1
            counter += 1
            if slider_bp[0] < cmin and len(slider_bp) > 1:
                slider_bp = slider_bp[1:] # drop values at the far end... 
                depths = depths[1:]
                missness = missness[1:]
                snps = snps[1:]
                # if counter % thin == 0:
                #     with open(outfix + '.depths', 'a') as ofi:
                #         ofi.write(f'{cmid}\t{np.mean(depths)}\t{np.mean(missness)}\t{np.mean(snps)}\n')
            elif slider_bp[0] < cmin and len(slider_bp) <= 1:
                reset = True
                csite = vsite
                cmin = csite - window
                cmid = csite - (window / 2)
                break

        if reset:
            slider_bp = np.full(1, variant.start)
            dpt = variant.gt_depths
            depths = np.reshape(dpt, newshape=(1, -1))
            missness = np.full(1, np.mean(np.greater_equal(dpt, missthresh)))
            snps = np.full(1, variant.is_snp)
        if not reset:
            slider_bp = np.concatenate((slider_bp, np.full(1, variant.start)))
            dpt = variant.gt_depths
            dpt[dpt<0] = 0
            depths = np.concatenate((depths, np.reshape(dpt, newshape=(1, -1))))
            missness = np.concatenate((missness, np.full(1, np.mean(np.greater_equal(dpt, missthresh)))))
            snps = np.concatenate((snps, np.full(1, variant.is_snp)))
        reset = False
        if counter % thin == 0:
            with open(outfix + '.depths', 'a') as ofi:
                ofi.write(f'{cmid}\t{np.mean(depths)}\t{np.mean(missness)}\t{np.mean(snps)}\n')

    return slider_bp, depths, missness, snps, csite, cmin, cmid, counter


def build_depth_global_all(infix, outfix, window, missthresh, slider_bp, depths, missness, csite, cmin, cmid, thin, counter):
    """Alternative traversal updating windowed stats across all sites."""
    vcf = cyvcf2.VCF(infix + '.vcf.gz')
    variter = iter(vcf)
    reset = False

    for variant in vcf:
        vsite = variant.start
        csite += 1
        cmin += 1



        # now check in site numbers
        while csite < vsite:
            #print(len(slider_bp))

            csite += 1
            cmin += 1
            cmid += 1
            counter += 1
            if slider_bp[0] < cmin and len(slider_bp) > 1:
                slider_bp = slider_bp[1:] # drop values at the far end... 
                depths = depths[1:]
                missness = missness[1:]
                if counter % thin == 0:
                    with open(outfix + '.depths', 'a') as ofi:
                        ofi.write(f'{cmid}\t{np.mean(depths)}\t{np.mean(missness)}\t{len(slider_bp) / window * 100}\n')
            elif slider_bp[0] < cmin and len(slider_bp) <= 1:
                reset = True
                csite = vsite
                cmin = csite - window
                cmid = csite - (window / 2)
                break

        if reset:
            slider_bp = np.full(1, variant.start)
            dpt = variant.gt_depths
            depths = np.reshape(dpt, newshape=(1, -1))
            missness = np.full(1, np.mean(np.greater_equal(dpt, missthresh)))
        if not reset:
            slider_bp = np.concatenate((slider_bp, np.full(1, variant.start)))
            dpt = variant.gt_depths
            depths = np.concatenate((depths, np.reshape(dpt, newshape=(1, -1))))
            missness = np.concatenate((missness, np.full(1, np.mean(np.greater_equal(dpt, missthresh)))))
        reset = False
        if counter % thin == 0:
            with open(outfix + '.depths', 'a') as ofi:
                ofi.write(f'{cmid}\t{np.mean(depths)}\t{np.mean(missness)}\t{len(slider_bp) / window * 100}\n')

    return slider_bp, depths, missness, csite, cmin, cmid, counter

def get_chrom_stats_list(fixlist, outfix, window=100, missthresh=3, thin=100):
    """Compute combined depth/missingness/SNP density stats across VCF prefixes."""
    with open(outfix + '.depths', "w")as ofi:
        ofi.write("position\tmean_depth_individual\tnonmiss_prop\tsnp_density\n")
    slider_bp, depths, missness, snps, csite, cmin, cmid, counter = init_depth_global_window_fixed(infix = fixlist[0], window=window, missthresh=missthresh)
    for infix in fixlist:
        print(f'Processing {infix}.vcf.gz')
        slider_bp, depths, missness, csite, cmin, cmid, counter = build_depth_global_window_fixed(
            infix=infix, outfix=outfix, window=window, missthresh=missthresh, slider_bp=slider_bp, 
            depths=depths, missness=missness, snps=snps, csite=csite, cmin=cmin, cmid=cmid, thin=thin, counter=counter
        )
    print("operation finished")


#get_depth_global_window_fixed(fname, fname, window=1000000, missthresh = 10)

# make_freqfile(fname, 'tester', popfile)


def _win_idx(pos: int, window: int) -> int:
    return (pos - 1) // window

def _safe_depths(var):
    dpt = var.depths
    dpt[dpt < 0] = 0
    return dpt

def _process_one_vcf_file(
    vcf_path: str,
    window: int,
    min_miss_frac: float,
    interpolate: bool,
    update_interval: int = 0,
    maskfile=None,
    genint = None 
):
    """
    Worker: process a single .vcf.gz and return (samples, chrom_dict, repcount, skipcount).
    chrom_dict[chrom] = {"depths": (W,N), "sitecounts": (W,), "window_indices": (W,)}
    """
    import numpy as np
    from treeflows import vcf_core  # ensure import inside worker if needed by multiprocessing

    vcf = vcf_core.VCFReader(vcf_path)
    samples = vcf.samples
    ll = len(samples)

    if genint is not None:
        _chrom, _start, _end = parse_genint(genint)
        #vcf.set_region(chrom, start, end)
    else:
        _chrom, _start, _end = None, None, None

    if maskfile is not None:
        mask = GenomeMask.from_gff3(maskfile, chromosomes=_chrom, pos_min=_start, pos_max=_end)

    chrom_dict = {}
    repcount = 0
    skipcount = 0

    vv = iter(vcf)

    try:
        var = next(vv)
        while var.missing_frac > min_miss_frac or (maskfile is not None and mask.pos_in_mask(var.CHROM, var.POS)):
            var = next(vv)
            skipcount += 1
    except StopIteration:
        # empty file
        return samples, chrom_dict, repcount, skipcount

    chrom = var.CHROM
    idx = _win_idx(var.POS, window)

    windowlist = []
    if interpolate and idx > 0:
        # fill initial empty windows (optional)
        for ii in range(idx):
            windowlist.append((np.zeros(ll), 0, ii))

    depths_window = np.zeros(ll, dtype=float)
    depths_window += _safe_depths(var)
    sites_window = 1
    repcount += 1

    def _finalize_chrom(_chrom, _windowlist, _depths_window, _sites_window, _idx):
        _windowlist.append((_depths_window, _sites_window, _idx))
        depths = np.vstack([wl[0] for wl in _windowlist]) if _windowlist else np.zeros((0, ll))
        sitecounts = np.array([wl[1] for wl in _windowlist], dtype=int) if _windowlist else np.array([], dtype=int)
        win_indices = np.array([wl[2] for wl in _windowlist], dtype=int) if _windowlist else np.array([], dtype=int)
        return {"depths": depths, "sitecounts": sitecounts, "window_indices": win_indices}

    while True:
        try:
            var = next(vv)
            while var.missing_frac > min_miss_frac or (maskfile is not None and mask.pos_in_mask(var.CHROM, var.POS)):
                var = next(vv)
                skipcount += 1
            repcount += 1
        except StopIteration:
            # finalize last chromosome
            chrom_dict[chrom] = _finalize_chrom(chrom, windowlist, depths_window, sites_window, idx)
            break

        # chromosome switch
        if var.CHROM != chrom:
            chrom_dict[chrom] = _finalize_chrom(chrom, windowlist, depths_window, sites_window, idx)

            # reset for new chrom
            chrom = var.CHROM
            idx = _win_idx(var.POS, window)
            windowlist = []
            if interpolate and idx > 0:
                for ii in range(idx):
                    windowlist.append((np.zeros(ll), 0, ii))
            depths_window = np.zeros(ll, dtype=float)
            sites_window = 0

        # window switch / catch up
        new_idx = _win_idx(var.POS, window)
        while new_idx != idx:
            if new_idx < idx:
                raise RuntimeError("new index less than current index (unsorted VCF?)")

            if interpolate:
                idx += 1
            else:
                idx = new_idx

            windowlist.append((depths_window, sites_window, idx))
            depths_window = np.zeros(ll, dtype=float)
            sites_window = 0

        # accumulate
        depths_window += _safe_depths(var)
        sites_window += 1

        if update_interval and repcount % update_interval == 0:
            # keep worker prints minimal or you’ll get messy output
            print(f"Processed {repcount:,} variants in {vcf_path}", flush=True)

    return samples, chrom_dict, repcount, skipcount


def _merge_chrom_dicts(chrom_dicts_in_order, ll):
    """
    Merge list of chrom_dict (one per file) preserving file order.
    """
    merged = {}
    for cd in chrom_dicts_in_order:
        for chrom, wd in cd.items():
            if chrom not in merged:
                merged[chrom] = {
                    "depths": wd["depths"],
                    "sitecounts": wd["sitecounts"],
                    "window_indices": wd["window_indices"],
                }
            else:
                merged[chrom]["depths"] = np.vstack([merged[chrom]["depths"], wd["depths"]])
                merged[chrom]["sitecounts"] = np.concatenate([merged[chrom]["sitecounts"], wd["sitecounts"]])
                merged[chrom]["window_indices"] = np.concatenate([merged[chrom]["window_indices"], wd["window_indices"]])
    return merged

def get_genint_str(filefix, pref, suff):
    reg = filefix.remove_prefix(pref).remove_suffix(suff)
    return reg

def parse_genint(genint):
    chrom, coords = genint.split(':')
    start, end = map(int, coords.split('-'))
    return chrom, start, end

def get_depths_dict_parallel(
    fixlist: list,
    window=10_000,
    n_jobs=8,
    min_miss_frac=0.5,
    interpolate=False,
    update_interval: int = 0, 
    maskfile=None,
    gen_intvs=None,
    use_ints=False): # silence worker progress

    fixfiles = [f + ".vcf.gz" for f in fixlist]
    genints = gen_intvs if use_ints else [None] * len(fixfiles) # if not using genints, just pass None to workers

    # run workers
    results = [None] * len(fixfiles)
    sample_ref = None
    ll = None
    total_rep = 0
    total_skip = 0

    with ProcessPoolExecutor(max_workers=n_jobs) as ex:
        futs = {
            ex.submit(
                _process_one_vcf_file,
                vcf_path,
                window,
                min_miss_frac,
                interpolate,
                update_interval, 
                maskfile,
                genints[i]
            ): i
            for i, vcf_path in enumerate(fixfiles)
        }

        for fut in as_completed(futs):
            i = futs[fut]
            samples, chrom_dict, repcount, skipcount = fut.result()
            results[i] = chrom_dict
            total_rep += repcount
            total_skip += skipcount

            if sample_ref is None:
                sample_ref = list(samples)
                ll = len(sample_ref)
            else:
                # hard requirement: sample order identical across files
                if list(samples) != sample_ref:
                    raise ValueError(
                        f"Sample list/order mismatch in {fixfiles[i]}.\n"
                        f"First few: {samples[:5]} vs {sample_ref[:5]}"
                    )

    print(f"Processed {len(fixfiles)} files. Variants kept: {total_rep:,}. Variants skipped: {total_skip:,}.")

    # merge in file order
    chrom_dict = _merge_chrom_dicts(results, ll)

    # from here, continue exactly as you already do:
    # - build glob_depths / glob_sitecounts from chrom_dict
    # - apply window_depth_thresh and your NB/Poisson/OptionB logic
    # - map back to chrom/positions and write outputs

    return chrom_dict, sample_ref

def nb_tail_pvalues_with_overdisp(k, lam, phi, overdisp_factor=1.0):
    """
    k: (W,N) observed counts (int or float, will be rounded)
    lam: (W,N) expected mean counts
    phi: scalar (global NB2 dispersion)
    overdisp_factor: scalar multiplier on phi

    returns: p_low, p_high (W,N)
    """
    eps = 1e-12

    # inflate dispersion
    phi_eff = float(phi) * float(overdisp_factor)
    if not np.isfinite(phi_eff) or phi_eff <= 0:
        phi_eff = 1e-6

    r = 1.0 / (phi_eff + eps)              # "size" parameter
    p = r / (r + lam + eps)                # success prob; mean=lam, var=lam+lam^2/r

    k_int = np.rint(k).astype(np.int64)

    p_low = nbinom.cdf(k_int, r, p)
    p_high = nbinom.sf(k_int - 1, r, p)
    return p_low, p_high, phi_eff


def get_depths_outliers_parallel(fixlist: list, outfix, window=10_000, write=True, threshold_frac=0.05, alpha=1e-6, overdisp_factor=1.0,
                        min_sites=20, min_miss_frac = 0.5, window_depth_thresh=5, n_jobs=8,
                        distrib = "nbinom", interpolate=False, q=False, update_interval = 100_000, use_median=True, 
                        maskfile=None, gen_intvs=None, use_ints=False):
    """Compute depth outlier windows across VCF prefixes in parallel."""

    if distrib not in ("poisson", "nbinom", "normal", "ragweed"):
        raise ValueError("distrib must be one of 'poisson', 'nbinom', 'normal', or 'ragweed'") # in future,  update ragweed to broader... 
    
    chrom_dict, samples = get_depths_dict_parallel(
        fixlist,
        window,
        n_jobs,
        min_miss_frac,
        interpolate,
        update_interval=update_interval, 
        maskfile=maskfile,
        gen_intvs=gen_intvs,
        use_ints=use_ints
    )
    fail_window_positions = get_outliers_from_chrom_dict(
        chrom_dict=chrom_dict,
        samples=samples,
        window=window,
        threshold_frac=threshold_frac,
        alpha=alpha,
        min_sites=min_sites,
        window_depth_thresh=window_depth_thresh,
        distrib=distrib,
        overdisp_factor=overdisp_factor, 
        use_median=use_median
    )

    if write:
        write_outlier_windows(fail_window_positions, outfix, samples, window, use_median=use_median)
    return fail_window_positions, chrom_dict, samples


def get_outliers_from_chrom_dict(
    chrom_dict,
    samples,
    window=10_000,
    threshold_frac=0.05,
    alpha=1e-6,
    min_sites=20,
    window_depth_thresh=5,
    distrib="nbinom", 
    overdisp_factor=1.0, 
    use_median=True):
    """Given chrom_dict and sample list, compute outlier windows."""

    if distrib not in ("poisson", "nbinom", "normal", "ragweed"):
        raise ValueError("distrib must be one of 'poisson', 'nbinom', 'normal', or 'ragweed'")
    ll = len(samples)
    idvs = samples
    # global summary: add rows from all chromosomes into one array
    glob_depths = np.vstack([chrom_dict[chrom]["depths"] for chrom in chrom_dict])
    glob_sitecounts = np.concatenate([chrom_dict[chrom]["sitecounts"] for chrom in chrom_dict])
    # get idv mean depth in array shaped such that we can normalize by idv average depth
    #idv_mean_depths = np.sum(glob_depths, axis=0) / np.sum(glob_sitecounts)

    ### Poisson outlier model ----
    k = glob_depths.astype(float)       # (n_windows, n_idv)
    E = glob_sitecounts.astype(float)       # (n_windows,)

    # prevent divide-by-zero in arithmetic
    E_safe = E.copy()
    E_safe[E_safe == 0] = np.nan

    # 1) per-sample global mean depth per site (size factor basis), normalized to mean 1
    ### NORMALIZED MEAN INDIVIDUAL DEPTH ###
    #idv_mean_depths = np.nansum(k, axis=0) / np.nansum(E_safe)   # mean DP per site, per sample (n_idv,)
    #s = idv_mean_depths / np.nanmean(idv_mean_depths)  # (n_idv,) - mean 1 across idvs (normalized mean depth, individuals)
    depth_per_site = k / (E_safe[:, None]) # per site depth
    usable_w = (E >= min_sites) & np.isfinite(E_safe)
    idv_scale = np.nanmedian(depth_per_site[usable_w, :], axis=0) if use_median else np.nanmean(depth_per_site[usable_w, :], axis=0) # median depth per site, per sample (n_idv,) across usable windows
    s = idv_scale / np.nanmedian(idv_scale) if use_median else idv_scale / np.nanmean(idv_scale) # (n_idv, ) median 1 across idvs (normalized median depths, indivdiuals)

    # 2) window baseline mu_w: median normalized depth per site across samples
    ### NORMED PER SITE & INDIVIDUAL (mean-normed) # TODO: figure out if this should be median??? could be affected by outliers??? 
    per_site_norm = k / (E_safe[:, None] * s[None, :])  # (n_windows, n_idv) # Y_{w, i} = k/(E*s)
    Y = per_site_norm # NORMED PER SITE/INDIVIDUAL

    # MEDIAN VALUE PER WINDOW (normed by mean previously)
    mu = np.nanmedian(per_site_norm, axis=1) if use_median else np.nanmean(per_site_norm, axis=1)  # (n_windows,)

    # expected lambda (mean counts) for each (window, individual)
    # EXPECTED COUNTS FOR WINDOW/INDIVIDUAL
    lam = (mu[:, None] * E_safe[:, None]) * s[None, :]  # (n_windows, n_idv)
    # SELECTION MASK FOR DATA CALCULATION
    usable = (E >= min_sites) & np.isfinite(mu) & (mu > 0) & (mu >= window_depth_thresh) # windows that aren't deep enough skipped

    # 3) estimate global NB dispersion phi (NB2)
    #     Var = mean + phi * mean^2
    # -------------------
    eps = 1e-12
    
    # across individuals, on the normalized per-site scale Y = k/(E*s)
    m_w = np.nanmean(Y, axis=1)  # (n_windows,)             # WINDOW MEAN # median??? 
    v_w = np.nanvar(Y, axis=1, ddof=1)   # (n_windows,)     # WINDOW VARIANCE

    phi_w = (v_w - m_w) / (m_w**2 + eps)  # (n_windows,)    # PHI (?) FOR NBINOM? 
    phi_w = np.where(phi_w > 0, phi_w, 0.0)  # set negative dispersions to 0

    # robust global phi (ignore non-usable windows)
    phi = np.nanmedian(phi_w[usable]) if use_median else np.nanmean(phi_w[usable])                       # NOW SET TO MEDIAN (global phi median?)
    # if phi collapses to 0 or NaN, fall back to small dispersion
    if not np.isfinite(phi) or phi <= 0:
        phi = 1e-6

    def alpha_from_sd(val):
        """Convert SD threshold to two-sided alpha."""
        return 2 * norm.sf(val)

    # -------------------
    # NB model
    # 4) NB tail p-values (one-sided)
    # SciPy nbinom parametrization
    # n = r (number of successes), p = success prob
    # mean = r*(1-p)/p = lam
    # var = lam + lam^2 / r
    # so r = 1/phi, p = r / (r + lam)
    if distrib == "nbinom":
        #r = 1.0 / phi
        #p = r / (r + lam)  # (n_windows, n_idv)

        # nbinom needs integer k
        #k_int = np.rint(k).astype(np.int64)

        # low tail: P(X <= k )
        #p_low = nbinom.cdf(k_int, r, p)

        # high tail: P(X >= k ) = sf(k-1)
        #p_high = nbinom.sf(k_int-1, r, p)
        p_low, p_high, phi_eff = nb_tail_pvalues_with_overdisp(k, lam, phi, overdisp_factor=overdisp_factor) # CALCULATED FROM DEPTH ARRAY, EXPECTED, PHI & OVERDISPERSE
        r = 1.0 / phi_eff
        p_nb = r / (r + lam)

        k_lower_cut = nbinom.ppf(alpha, r, p_nb) # NBINOM QUANTILE FOR ALPHA
        k_upper_cut = nbinom.isf(alpha, r, p_nb) # NB

    # -----------------------
    # POISSON model
    # 4) Poisson tail p-values (one-sided)
    # low tail: P(X <= k )
    elif distrib == "poisson":
        p_low = poisson.cdf(k, lam) # DEPTH ARRAY VS EXPECTED
        
        # largest k such that P(X <= k) <= alpha
        k_lower_cut = poisson.ppf(alpha, lam) # POISSON QUANTILE FOR ALPHA

        # high tail: P(X >= k ) = sf(k-1)
        p_high = poisson.sf(k-1, lam) 
        # smallest k such that P(X >= k) <= alpha
        k_upper_cut = poisson.isf(alpha, lam) # POISSON INVERSE SF FOR ALPHA

    elif distrib == "normal":
        # normal approximation
        sd = np.sqrt(lam + phi * lam**2) # CONVERT TO NORMAL APPROX
        z_scores = (k - lam) / sd

        p_low = norm.cdf(z_scores)
        p_high = norm.sf(z_scores)

        k_lower_cut = norm.ppf(alpha, loc=lam, scale=sd)
        k_upper_cut = norm.isf(alpha, loc=lam, scale=sd)

    elif distrib == "ragweed": # Try to implement something approximating the ragweed paper (2sd for each window)
        # mean & sd from each individual window... for best results, use mean mode... 
        # m_w & v_w here... now for sd_w
        # make dispersion factor here 2x SD
        sd_w = np.sqrt(v_w) # SD PER WINDOW
        z_scores = (k - lam) / sd_w[:, None]

        p_low = norm.cdf(z_scores)
        p_high = norm.sf(z_scores)

        # ragweed update: TODO evaluate this later
        alpha = alpha_from_sd(2) # 2 SD threshold
        k_lower_cut = norm.ppf(alpha, loc=lam, scale=sd_w[:, None]) # should be windowed... in future can generalize this... 
        k_upper_cut = norm.isf(alpha, loc=lam, scale=sd_w[:, None])

    # 5) IDENTIFY OUTLIER INDIDIVUALS IN WINDOWS based on alpha threshold and usable windows
    fail_lower = (p_low <= alpha) & usable[:, None]
    fail_upper = (p_high <= alpha) & usable[:, None]
    fail_array = fail_lower | fail_upper

    # 6) Cohort-level filter: window must have >= threhold_frac outliers
    fail_frac_upper = fail_upper.mean(axis=1)
    fail_frac_lower = fail_lower.mean(axis=1)
    fail_frac = fail_array.mean(axis=1) # FAIL ARRAY (MEAN PER WINDOW) # TODO: update this so that only 1 direction counts... 
    #fail_windows = np.where(fail_frac >= threshold_frac)[0] # SELECT BASED ON THRESHOLD
    upper_fail = fail_frac_upper >= threshold_frac
    lower_fail = fail_frac_lower >= threshold_frac
    fail_windows = np.where(upper_fail | lower_fail)[0] # SELECT ON THRESHOLD - MUST BE IN AT LEAST ONE

    depth_lower_cut = k_lower_cut / E_safe[:, None] # CONVERT TO DEPTH CUTS
    depth_upper_cut = k_upper_cut / E_safe[:, None]

    if distrib == "nbinom":
        print(f"NB dispersion phi_eff={phi_eff:.3g} (x {overdisp_factor}); "
              f"flagged {len(fail_windows)} windows at alpha={alpha}, frac>={threshold_frac}, min_sites={min_sites}", flush=True)
    ### PREPARE FINAL OUTPUT ###

    fail_window_positions = []
    cum_window_counts = []
    cum_count = 0
    for chrom in chrom_dict:
        win_dict = chrom_dict[chrom]
        nwin = len(win_dict["sitecounts"])
        start_idx = cum_count
        cum_count += nwin
        cum_window_counts.append((chrom, cum_count, start_idx))
    for fw in fail_windows:
        for cc in cum_window_counts:
            chrom, cum_count, start_idx = cc
            if fw < cum_count:
                # this is the chromosome
                win_dict = chrom_dict[chrom]
                wstart = win_dict["window_indices"][fw - start_idx] * window + 1
                depth_meds_window = mu[fw]
                #depth_sd_window = std_depths[fw]
                fail_upper_list = np.where(fail_upper[fw])[0]
                fail_lower_list = np.where(fail_lower[fw])[0]
                depth_lower_cut_window = np.nanmedian(depth_lower_cut[fw])
                depth_upper_cut_window = np.nanmedian(depth_upper_cut[fw])

                #print(f"{chrom}:{wstart}-{wstart+window-1}, frac: {fail_frac[fw]:.3f}, high depth: {', '.join([idvs[i] for i in fail_upper_list])}, low depth: {', '.join([idvs[i] for i in fail_lower_list])}")
                fail_window_positions.append((chrom, wstart, fail_lower_list, fail_upper_list, fail_frac[fw], 
                                              depth_meds_window, depth_lower_cut_window, depth_upper_cut_window))
                break
                # else continue
    print(f"{'Median' if use_median else 'Mean'} outlier fraction per window:", np.median(fail_frac) if use_median else np.mean(fail_frac))
    print("95th percentile outlier fraction per window:", np.percentile(fail_frac, 95))
    return fail_window_positions


def write_outlier_windows(fail_window_positions, outfix, samples, window=10_000, use_median=True):
    """Write outlier windows to file."""
    idvs = samples
    with open(outfix + '.windows_final', 'w') as ofi, open(outfix + '.windows_low', 'w') as ofil, open(outfix + '.windows_high', 'w') as ofih:
        ofi.write(f'chrom_position\tlow_depths\thigh_depths\tfrac\tdepth_{ "median" if use_median else "mean" }\tdepth_lower_cut\tdepth_upper_cut\n')
        ofil.write(f'chrom_position\twindow_start\tlow_depth_ids\tfrac\tdepth_{ "median" if use_median else "mean" }\tdepth_lower_cut\tdepth_upper_cut\n')
        ofih.write(f'chrom_position\twindow_start\thigh_depth_ids\tfrac\tdepth_{ "median" if use_median else "mean" }\tdepth_lower_cut\tdepth_upper_cut\n')
        for fwp in fail_window_positions:
            chrom, wstart, fail_lower_list, fail_upper_list, fail_frac, depth_meds_window, depth_lower_cut_window, depth_upper_cut_window = fwp
            lower_ids = ','.join([idvs[i] for i in fail_lower_list])
            upper_ids = ','.join([idvs[i] for i in fail_upper_list])
            depth_lower_cut_window = float(depth_lower_cut_window)
            depth_upper_cut_window = float(depth_upper_cut_window)
            ofi.write(f'{chrom}:{wstart}-{wstart+window-1}\t{wstart}\t{lower_ids}\t{upper_ids}\t{fail_frac:.3f}\t{depth_meds_window}\t{depth_lower_cut_window:.2f}\t{depth_upper_cut_window:.2f}\n')
            if len(lower_ids) > 0:
                ofil.write(f'{chrom}:{wstart}-{wstart+window-1}\t{wstart}\t{lower_ids}\t{fail_frac:.3f}\t{depth_meds_window}\t{depth_lower_cut_window:.2f}\t{depth_upper_cut_window:.2f}\n')
            if len(upper_ids) > 0:
                    ofih.write(f'{chrom}:{wstart}-{wstart+window-1}\t{wstart}\t{upper_ids}\t{fail_frac:.3f}\t{depth_meds_window}\t{depth_lower_cut_window:.2f}\t{depth_upper_cut_window:.2f}\n')

# for tract in tractlist:
#     print(batchid + '_' + tract + '_' + postfix + '.vcf.gz')

def win_idx(p, window):
        # function that works out which window index a snp should be in (chromosome agnostic)
        return (p - 1) // window

def safe_depths(v):
    # function that ensures no strange depths read from the vcf file (due to e.g. encoding issues)
    dpt = v.depths
    dpt[dpt < 0] = 0
    return dpt

def _get_depths_file(filename, window, min_miss_frac, update_interval):
    """Function will build core elements for multi-thread. returns dict"""

    ## INIT
    skipcount = 0
    repcount = 0

    vcf = _vcf.VCFReader(filename)
    ll = len(vcf.samples)
    vv = iter(vcf)
    var = next(vv)
    pos = var.POS
    chrom = var.CHROM
    idx = win_idx(pos, window)

    windowlist = [] # try here. no interpolation as other blocks exist

    depths = safe_depths(var)                   # np array of depths for site
    depths_window = np.zeros(ll)                # instantiate window depth
    depths_window += depths                     # update window depths with site depth
    sites_window = 1                            # set initial site window (???)



    while True:

        ### VARIANT SWITCHING BLOCK
        try:
            var = next(vv)
            while var.missing_frac > min_miss_frac:
                var=next(vv)
                skipcount += 1
            repcount += 1
        except StopIteration:
            windowlist.append((depths_window, sites_window, idx, chrom))
            break

        if var.CHROM != chrom:
            windowlist.append((depths_window, sites_window, idx, chrom))
            chrom = var.CHROM
            idx = win_idx(var.POS)
            depths_window = np.zeros(ll)
            sites_window = 0

        new_idx = win_idx(var.POS)
        while new_idx != idx:
            # this while loop enables us to catch up if there are empty windows
            if new_idx < idx:
                print("Error: new index less than current index!")
                exit(1)
            idx = new_idx
            windowlist.append((depths_window, sites_window, idx, chrom))
            depths_window = np.zeros(ll)
            sites_window = 0

        # DEFAULT BLOCK - PROCESS SITE
        dpt = safe_depths(var)
        depths_window += dpt
        sites_window += 1
        if repcount % update_interval == 0:
            print(f"Processed {repcount/1_000_000}M variants... current pos {var.POS:,} on {var.CHROM}", flush=True)

    return windowlist, skipcount, repcount

def _get_depths_all(fixlist: list, window, min_miss_frac, update_interval):
    """Aggregate data across all files/chromosomes"""

    # going to be multivar, but will start simple... 
    winlist = []
    skipcount = 0
    repcount = 0

    fixfiles = [f + '.vcf.gz' for f in fixlist] # (order) file list for processing
    for filename in fixfiles:
        w, s, r = _get_depths_file(filename, window, min_miss_frac, update_interval)
        winlist.extend(w)
        skipcount += s
        repcount += r

    return winlist, skipcount, repcount



def get_depth_global_list(fixlist, outfix, window = 100, missthresh = 3, thin = 10): # assumes 
    """Compute global (INFO/DP-based) depth statistics across VCF prefixes."""
    with open(outfix + '.gdepth', 'w') as ofi:
        ofi.write('position\tgdepth\n')

    vcf = cyvcf2.VCF(fixlist[0] + '.vcf.gz')
    nidv = len(vcf.samples)
    vv = iter(vcf)
    var = next(vv)
    pos = var.start
    pos_lower = pos - window
    pmid = int(pos - (window / 2))
    gdepth = np.full(1, var.INFO["DP"])
    init = True

    for infix in fixlist:
        vcf = cyvcf2.VCF(infix + '.vcf.gz')
        for variant in vcf:
            if init:
                #print('init')
                gdepth = np.concatenate((gdepth, np.full(1, variant.INFO["DP"])))

            else:
                try:
                    gdepth = np.concatenate((gdepth[1:], np.full(1, variant.INFO["DP"])))
                except KeyError:
                    #print('key!')
                    continue
            
            pos = variant.start
            pos_lower = pos - window
            pmid = int(pos - (window / 2))
            if pos_lower > 0 and init:
                init = False
            if not init and pmid % thin == 0:
                with open(outfix + '.gdepth', 'a') as ofi:
                    ofi.write(f'{pmid}\t{round(np.mean(gdepth) / nidv, 3)}\n')


def get_miss_global_list(fixlist, outfix, window = 100, missthresh = 3, thin = 10): # assumes 
    """Compute global callable fraction statistics across VCF prefixes."""
    with open(outfix + '.gmiss', 'w') as ofi:
        ofi.write('position\tnonmissing_frac\n')

    vcf = cyvcf2.VCF(fixlist[0] + '.vcf.gz')
    nidv = len(vcf.samples)
    vv = iter(vcf)
    var = next(vv)
    pos = var.start
    pos_lower = pos - window
    pmid = int(pos - (window / 2))
    dpt = var.gt_depths
    missness = np.full(1, np.mean(np.greater_equal(dpt, missthresh)))
    init = True

    for infix in fixlist:
        vcf = cyvcf2.VCF(infix + '.vcf.gz')
        for variant in vcf:
            dpt = variant.gt_depths
            if init:
                missness = np.concatenate((missness, np.full(1, np.mean(np.greater_equal(dpt, missthresh)))))
                #gdepth = np.concatenate((gdepth, np.full(1, variant.INFO["DP"])))

            else:
                missness = np.concatenate((missness[1:], np.full(1, np.mean(np.greater_equal(dpt, missthresh)))))
                # try:
                #     gdepth = np.concatenate((gdepth[1:], np.full(1, variant.INFO["DP"])))
                # except KeyError:
                #     #print('key!')
                #     continue
            
            pos = variant.start
            pos_lower = pos - window
            pmid = int(pos - (window / 2))
            if pos_lower > 0 and init:
                init = False
            if not init and pmid % thin == 0:
                with open(outfix + '.gmiss', 'a') as ofi:
                    ofi.write(f'{pmid}\t{round(np.mean(missness), 3)}\n')


def get_all_global_list(fixlist, outfix, window = 100, missthresh = 3, thin = 10): # assumes 
    """Compute combined global depth/missingness/SNP stats across VCF prefixes."""
    with open(outfix + '.gstats', 'w') as ofi:
        ofi.write('position\tgdepth\tnonmissing_frac\tsnp_density\n')

    vcf = cyvcf2.VCF(fixlist[0] + '.vcf.gz')
    nidv = len(vcf.samples)
    vv = iter(vcf)
    var = next(vv)
    pos = var.start
    pos_lower = pos - window
    pmid = int(pos - (window / 2))
    dpt = var.gt_depths
    snps = np.full(1, var.is_snp)
    missness = np.full(1, np.mean(np.greater_equal(dpt, missthresh)))
    gdepth = np.full(1, var.INFO["DP"])
    init = True

    for infix in fixlist:
        vcf = cyvcf2.VCF(infix + '.vcf.gz')
        for variant in vcf:
            dpt = variant.gt_depths
            if init:
                snps = np.concatenate((snps, np.full(1, variant.is_snp)))
                missness = np.concatenate((missness, np.full(1, np.mean(np.greater_equal(dpt, missthresh)))))
                gdepth = np.concatenate((gdepth, np.full(1, variant.INFO["DP"])))

            else:
                try:
                    gdepth = np.concatenate((gdepth[1:], np.full(1, variant.INFO["DP"])))
                except KeyError:
                    #print('key!')
                    continue
                snps = np.concatenate((snps[1:], np.full(1, variant.is_snp)))
                missness = np.concatenate((missness[1:], np.full(1, np.mean(np.greater_equal(dpt, missthresh)))))

            
            pos = variant.start
            pos_lower = pos - window
            pmid = int(pos - (window / 2))
            if pos_lower > 0 and init:
                init = False
            if not init and pmid % thin == 0:
                with open(outfix + '.gstats', 'a') as ofi:
                    ofi.write(f'{pmid}\t{round(np.mean(gdepth) / nidv, 3)}\t{round(np.mean(missness), 3)}\t{round(np.mean(snps), 3)}\n')

def get_snps_global_list(fixlist, outfix, window = 100, missthresh = 3, thin = 10): # assumes 
    """Compute global SNP density statistics across VCF prefixes."""
    with open(outfix + '.gsnps', 'w') as ofi:
        ofi.write('position\tsnp_density\n')

    vcf = cyvcf2.VCF(fixlist[0] + '.vcf.gz')
    nidv = len(vcf.samples)
    vv = iter(vcf)
    var = next(vv)
    pos = var.start
    pos_lower = pos - window
    pmid = int(pos - (window / 2))
    dpt = var.gt_depths
    snps = np.full(1, var.is_snp)
    init = True

    for infix in fixlist:
        vcf = cyvcf2.VCF(infix + '.vcf.gz')
        for variant in vcf:
            #dpt = variant.gt_depths
            if init:
                snps = np.concatenate((snps, np.full(1, variant.is_snp)))
                #missness = np.concatenate((missness, np.full(1, np.mean(np.greater_equal(dpt, missthresh)))))
                #gdepth = np.concatenate((gdepth, np.full(1, variant.INFO["DP"])))

            else:
                snps = np.concatenate((snps[1:], np.full(1, variant.is_snp)))
                #missness = np.concatenate((missness[1:], np.full(1, np.mean(np.greater_equal(dpt, missthresh)))))
                # try:
                #     gdepth = np.concatenate((gdepth[1:], np.full(1, variant.INFO["DP"])))
                # except KeyError:
                #     #print('key!')
                #     continue
            
            pos = variant.start
            pos_lower = pos - window
            pmid = int(pos - (window / 2))
            if pos_lower > 0 and init:
                init = False
            if not init and pmid % thin == 0:
                with open(outfix + '.gsnps', 'a') as ofi:
                    ofi.write(f'{pmid}\t{round(np.mean(snps), 3)}\n')


# tracts = ['../../rawvcfs/' + batchid + '_' + t + '_' + postfix for t in -]
# tracta = tracts[0]
# get_all_global_list(tracts, 'tnew', window=10000, missthresh=3, thin = 1000)


def get_nonrepeats(infix, outfix, window=1000, interval = 1, chrom=1): # this is the new function
    """Derive non-repeat intervals by scanning SNP density (workflow-specific)."""
    chromdict = {1: 'NC_035107.1', 2: 'NC_035108.1', 3: 'NC_035109.1'}
    clendict = {1: 310827022, 2: 474425716, 3: 409777670}
    chromname = chromdict[chrom]
    chromlength = clendict[chrom]
    vstarts = []
    vends = []
    vlens = []
    count = 0

    # set up the window... these three will slide together... but not write until winstart is in the frame... 
    focus = int(0 - window / 2)
    winstart = 0 - window
    winend = 0
    with open(outfix, 'w') as ofi:
        ofi.write("chrom\tpos\tnonrepeat_frac\n")

    with gzip.open(infix) as ifi:
        iiter = iter(ifi)
        while (True):
            try:
                i = next(iiter)
            except StopIteration:
                i = None
                break
            ii = i.decode()
            if not ii.startswith(chromname):
                continue

            # now we have isolated the core files... 
            # First, we add the new chromosomal feature to the file
            parsed = ii.strip().split(sep='\t')
            vstart_new = int(parsed[3])
            vend_new = int(parsed[4])
            vlen_new = vend_new - vstart_new

            # now we attempt to move our window... we will iterate to the end of the next feature (using winend)
            while winend < vend_new: # while the window hasn't caught up to the new feature end... the == case is handled in this loop... 
                winstart += 1
                focus += 1
                winend += 1
                # iterator time
                count += 1

                # get rid of features no longer in the window
                if vends and winstart > vends[0]: # if the end of sequence no longer contained
                    vstarts = vstarts[1:]
                    vends = vends[1:]
                    vlens = vlens[1:]

                # now compute features... 
                win_replen = sum(vlens) # raw, and now need to deal with the edge cases
                # first: overlap at start
                if vstarts and winstart > vstarts[0]: # if there is an overlap... other condition should be dealt with... 
                    win_replen -= (winstart - vstarts[0])

                # next: overlap at end
                if winend >= vstart_new: # we have the new overlap
                    overlap = (winend - vstart_new + 1)
                    if overlap > window:
                        overlap = window
                    win_replen += overlap

                # now: add new end to sequence if necessary... 
                if winend == vend_new: # we have achieved complete overlap!
                    vstarts.append(vstart_new)
                    vends.append(vend_new)
                    vlens.append(vlen_new)

                # loop is updated. now for the writing step... 
                # first we check the interval variable. 
                if winstart > 0 and count % interval == 0: # if the windows are efficient 
                    # perform calculation... 
                    nonrepeats_raw = window - win_replen # the bp that are not repetitive sequences
                    nonrepeats_fraction = nonrepeats_raw / window # the fraction
                    with open(outfix, 'a') as ofi:
                        ofi.write(f"{chromname}\t{focus}\t{round(nonrepeats_fraction, 3)}\n")
        
        ## Clean up at the end of the file. 
        while winend < chromlength:
            winstart += 1
            focus += 1
            winend += 1
            count += 1

             # get rid of features no longer in the window
            if vends and winstart > vends[0]: # if the end of sequence no longer contained
                vstarts = vstarts[1:]
                vends = vends[1:]
                vlens = vlens[1:]

            # now compute features... 
            win_replen = sum(vlens) # raw

            # first: overlap at start (no overlap at end during cleanup)
            if vstarts and winstart > vstarts[0]: # if there is an overlap... other condition should be dealt with... 
                win_replen -= (winstart - vstarts[0])

            # WRITING STEP 
            # first we check the interval variable. 
            if n % interval == 0: # if the windows are efficient . 
                # perform calculation... 
                nonrepeats_raw = window - win_replen # the bp that are not repetitive sequences
                nonrepeats_fraction = nonrepeats_raw / window # the fraction
                with open(outfix, 'a') as ofi:
                    ofi.write(f"{chromname}\t{focus}\t{round(nonrepeats_fraction, 3)}\n")


def get_depth_global_window_idv(infix, outfix, window=100, missthresh = 3, interval = 1): # process is curently extremely naive!!!!! 
    """Compute windowed depth/missingness per individual and write to file."""
    # can possibly make a bp slider... that needs to be in the middle of the window? 
    vcf = cyvcf2.VCF(infix + '.vcf.gz')
    nsamps = len(vcf.samples)
    depths = np.zeros((window, nsamps))
    missness = np.zeros(window)
    slider_bp = np.zeros(window)
    with open(outfix + '.depths', "w")as ofi:
        ofi.write("position\tmean_depth_individual\tnonmiss_prop\n")
    #print(slider)

    for n, variant in enumerate(vcf):
        if not variant.start == variant.end - 1:
            print("multisite warning")
        dpt = variant.gt_depths
        if n < window:
            depths[n] = dpt
            slider_bp[n] = variant.start # note that this assumes no non-SNP data... might need to sanity check this... 
            missness[n] = np.mean(np.greater_equal(dpt, missthresh))
        #print(f'{variant.CHROM}\t{variant.start}\t{variant.INFO["DP"]}')
        if n >= window:
            dpt = variant.gt_depths
            depths = np.concatenate((depths[1:], np.reshape(dpt, newshape=(1, -1))))
            slider_bp = np.concatenate((slider_bp[1:], np.full(1, variant.start)))
            missness = np.concatenate((missness[1:], np.full(1, np.mean(np.greater_equal(dpt, missthresh)))))
            #print(slider)
            #slider_bp = np.concatenate(slider_bp[1:], variant.start)
            #print(f'{median(slider_bp)}\t{mean(slider)}')
            if n % interval == 0:
                with open(outfix + '.depths', 'a') as ofi:
                    ofi.write(f'{median(slider_bp)}\t{np.mean(depths)}\t{np.mean(missness)}\n')
    #print(slider)

#get_depth_global_window_idv(tracta, 'old_test', window=100000, missthresh=3, interval=100)

def vcf_trio_counts(vcffile, focus, out1, out2):
    ''''''

    vcf = _vcf.VCFReader(vcffile)
    vcf.set_samples_ordered([focus, out1, out2]) # need to be in order... (order is not going to be possible)

    discrim = 0
    m1 = 0
    m2 = 0
    mid = 0

    for site in vcf:
        if not site.is_biallelic:
            continue

        # now for biallelic sites... 
        geno = np.array(site.genotypes)[site._sortkey]
        geno = np.array([sum(g[:2]) for g in geno])
        if any(geno < 0):
            continue
        if any(geno[1:] == 1):
            continue
        if geno[1] == geno[2]:
            continue
        discrim += 1
        if geno[0] == geno[1]:
            m1 += 1
        elif geno[0] == geno[2]:
            m2 += 1
        elif geno[0] == 1:
            mid += 1
        else:
            print('bad code')

    print(f"Discriminatory sites: {discrim}. O1 match: {m1}. O2 match: {m2}. Intermediate: {mid}") 
    return [discrim, m1, m2, mid]

def vcf_trio_counts_in(vcffile, focus, out1, out2, in1):
    ''''''

    vcf = _vcf.VCFReader(vcffile)
    vcf.set_samples_ordered([focus, out1, out2, in1]) # need to be in order... (order is not going to be possible)

    discrim = 0
    m1 = 0
    m2 = 0
    mid = 0

    for site in vcf:
        if not site.is_biallelic:
            continue

        # now for biallelic sites... 
        geno = np.array(site.genotypes)[site._sortkey]
        geno = np.array([sum(g[:2]) for g in geno])
        if any(geno < 0):
            continue
        if any(geno[1:] == 1):
            continue
        if geno[1] != geno[3]: #ingroup constraint - going to use Saudi... 
            continue
        if geno[1] == geno[2]:
            continue
        discrim += 1
        if geno[0] == geno[1]:
            m1 += 1
        elif geno[0] == geno[2]:
            m2 += 1
        elif geno[0] == 1:
            mid += 1
        else:
            print('bad code')

    print(f"Ingroup discriminatory sites: {discrim}. O1 match: {m1}. O2 match: {m2}. Intermediate: {mid}") 
    return [discrim, m1, m2, mid]

def vcf_trio_2_counts_in(vcffile, focus, out1, out2, in1):
    ''''''

    vcf = _vcf.VCFReader(vcffile)
    vcf.set_samples_ordered([focus, out1, out2, in1]) # need to be in order... (order is not going to be possible)

    discrim = 0
    m1 = 0
    m2 = 0
    mid = 0

    for site in vcf:
        if not site.is_biallelic:
            continue

        # now for biallelic sites... 
        geno = np.array(site.genotypes)[site._sortkey]
        geno = np.array([sum(g[:2]) for g in geno])
        if any(geno < 0):
            continue
        if any(geno[1:] == 1):
            continue
        if geno[1] != geno[3]: #ingroup constraint - going to use Saudi... 
            continue
        if geno[1] == geno[2]:
            continue
        discrim += 1
        if geno[0] == geno[1]:
            m1 += 1
        elif geno[0] == geno[2]:
            m2 += 1
        elif geno[0] == 1:
            mid += 1
        else:
            print('bad code')

    print(f"Ingroup discriminatory sites: {discrim}. O1 match: {m1}. O2 match: {m2}. Intermediate: {mid}") 
    return [discrim, m1, m2, mid]



def vcf_quicksum(vcffile, sitestart=None, siteend=None):
    """Calculates simple site-based stats for a vcf file. 
    Specifically: Sites (over a region), % coverage/missing, 
    Counts of invariant vs SNP, biallelic vs multi, singleton vs other"""
    vcf = _vcf.VCFReader(vcffile)
    ## summary vars
    n_sites = 0
    n_monos = 0
    n_snps = 0
    n_bi = 0
    n_single = 0
    n_double = 0
    n_triple = 0
    n_multi = 0
    n_other = 0
    sitepos = 0

    for site in vcf:
        if not n_sites % 500_000:
            print(site.POS)
        # window logic
        sitepos = site.POS
        if (siteend is not None) and sitepos > siteend:
            break
        if (sitestart is not None) and sitepos < sitestart:
            continue
        if sitestart is None:
            sitestart = site.POS

        if not site.passes_filter(monoallelic=True, multiallelic=True): 
            n_other += 1
            continue
        n_sites += 1
        if site.is_monoallelic:
            n_monos += 1
            continue
        if site.is_nonsnp:
            n_other == 1
            continue
        n_snps += 1
        if site.is_multiallelic:
            n_multi += 1
            continue
        n_bi += 1
        if site.has_singleton:
            n_single += 1
        elif site.has_doubleton:
            n_double += 1
        elif site.has_tripleton:
            n_triple += 1
    
    # close the window
    if siteend is None:
        siteend = sitepos

    # calculations
    segment_length = siteend - sitestart + 1 # taking 1-based vcf into account
    frac_coverage = n_sites / segment_length

    ret = f"Segment length: {segment_length} bp ({segment_length/1_000_000:.2n})\n"
    ret += f"Site coverage: {frac_coverage:.2%} ({1-frac_coverage:.2%} missing)\n"
    ret += f"Total sites: {n_sites:n}\n"
    ret += f"\tInvariant: {n_monos}.\tSNPs: {n_snps}\n"
    ret += f"SNPS: Multiallelic: {n_multi:n}.\tBiallelic: {n_bi:n}\n"
    ret += f" (biallelic) Singletons: {n_single:n}\tDoubles: {n_double:n}\tTriples: {n_triple:n}\tRest: {n_snps - n_single - n_double - n_triple:n}\n"
    ret += f" (other): {n_other:n}"
    print(ret)
    return ret
        

def vcf_idv_misstats(vcf_prefix: str|Path, outprefix: str|Path = None, filterfile: str|Path = None, finalize=False, finalprefix=None):
    """Compute per-individual depth and callable-site counts for a VCF prefix.

    Args:
        vcf_prefix: Prefix/path without `.vcf.gz` suffix.
        outprefix: Output prefix for intermediate stats.
        filterfile: Optional VCF to intersect against (uses `VCFCompare`).
        finalize: If True, write a final output at `finalprefix`.
        finalprefix: Output prefix for the final aggregated file.

    Returns:
        None.
    """
    vcf_file = str(vcf_prefix) + ".vcf.gz"
    outprefix = str(vcf_file) if outprefix is None else str(outprefix)
    usefilter = False if filterfile is None else True
    if finalize and finalprefix is None:
        finalprefix = outprefix
    sitecount = 0
    if usefilter:
        vcf = _vcf.VCFCompare(vcf_file, filterfile)
    else:
        vcf = _vcf.VCFReader(vcf_file)

    samples = vcf.samples

    depthcounts = np.zeros(len(samples), dtype=int)
    nonmisscounts = np.zeros(len(samples), dtype=int)

    if usefilter:
        for site1, _ in vcf:
            sitecount += 1
            depths = np.array(site1.depths, dtype=int)
            depths[depths < 0] = 0
            nonmiss = (depths >=3).astype(int)
            depthcounts += depths
            nonmisscounts += nonmiss

    else:
        for site in vcf:
            sitecount += 1
            depths = np.array(site.depths, dtype=int)
            depths[depths < 0] = 0
            nonmiss = (depths >= 3).astype(int)
            depthcounts += depths
            nonmisscounts += nonmiss

    print(f"Processed {sitecount} sites...")
    #print(depthcounts)
    #print(nonmisscounts)

    with open(outprefix + ".idvstat", "w") as ofi: 

        ofi.write(f"{sitecount}\n")
        np.savetxt(ofi, samples[np.newaxis, :], fmt='%s', delimiter=" ")
        np.savetxt(ofi, depthcounts[np.newaxis, :], fmt='%d', delimiter=" ")
        np.savetxt(ofi, nonmisscounts[np.newaxis, :], fmt='%i', delimiter=" ")

    if finalize:

        with open(finalprefix + ".miss", "w") as ofi:
            nonmiss_frac = nonmisscounts / sitecount
            depth_avg = depthcounts / sitecount

            for idv, nonmiss, nonmiss_frac, depth in zip(samples, nonmisscounts, nonmiss_frac, depth_avg):
                ofi.write(f"{idv}\t{nonmiss}\t{nonmiss_frac}\t{depth}\n")

            


def vcf_idv_misstats_gather(stat_list: list, outfile, template_vcf = None):
    """Gather multiple per-individual stats files into a single table."""
    # vcf_list here coudl be any files... but not really
    if template_vcf is None:
        template_vcf = stat_list[0]

    vcf = _vcf.VCFReader(template_vcf) # a terrible dependency
    samples_start = [str(samp) for samp in vcf.samples]

    sitecount = 0
    depthcounts = np.zeros(len(samples_start), dtype=int)
    nonmisscounts = np.zeros(len(samples_start), dtype=int)

    for stat_file in stat_list:

        with open(str(stat_file) + ".idvstat") as ifi:
            sitecount += int(next(ifi).strip())
            next(ifi)
            dc = np.fromstring(next(ifi), dtype=int, sep=" ")
            nmc  = np.fromstring(next(ifi), dtype=int, sep=" ")
            depthcounts += dc
            nonmisscounts += nmc

    nonmiss_frac = nonmisscounts / sitecount
    depth_avg = depthcounts / sitecount

    with open(outfile, "w") as ofi: 
        ofi.write("idv\tnonmiss\tnonmiss_frac\tdepth_avg\n")
        
        for idv, nonmiss, nonmiss_frac, depth in zip(samples_start, nonmisscounts, nonmiss_frac, depth_avg):

            ofi.write(f"{idv}\t{nonmiss}\t{nonmiss_frac}\t{depth}\n")

