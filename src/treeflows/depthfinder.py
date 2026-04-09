#!/bin/python3

import numpy as np
from treeflows import vcf_core as _vcf
from scipy.stats import norm
from concurrent.futures import ProcessPoolExecutor, as_completed
from treeflows.gff import GenomeMask
from treeflows.Spath import Mpath
from pathlib import Path


def get_genint_str(filefix, pref, suff):
    reg = filefix.remove_prefix(pref).remove_suffix(suff)
    return reg

def parse_genint(genint):
    chrom, coords = genint.split(':')
    start, end = map(int, coords.split('-'))
    return chrom, start, end

def _win_idx(pos: int, window: int) -> int:
    return (pos - 1) // window

def _safe_depths(var):
    dpt = var.depths
    dpt[dpt < 0] = 0
    return dpt

def alpha_from_sd(val):
        """Convert SD threshold to two-sided alpha."""
        return float(2 * norm.sf(val))

def _transform_depth(Y, transform="none", log_pseudocount=0.5):
    """Transform normalized depth matrix."""
    if transform == "none":
        return Y
    elif transform == "sqrt":
        return np.sqrt(Y)
    elif transform == "log":
        return np.log(Y + log_pseudocount)
    else:
        raise ValueError("transform must be 'none', 'log', or 'sqrt'")

def _inverse_transform_depth(Yt, transform="none", log_pseudocount=0.5):
    if transform == "none":
        return Yt
    elif transform == "sqrt":
        return Yt ** 2
    elif transform == "log":
        return np.exp(Yt) - log_pseudocount
    else:
        raise ValueError("transform must be 'none', 'log', or 'sqrt'")

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

    vcf = _vcf.VCFReader(vcf_path)
    samples = vcf.samples
    ll = len(samples)

    if genint is not None:
        _chrom, _start, _end = parse_genint(genint)
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


def _merge_chrom_dicts(chrom_dicts_in_order):
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


def get_depths_dict_parallel(
    fixlist: list|str,
    window=10_000,
    n_jobs=8,
    min_miss_frac=0.5,
    interpolate=False,
    update_interval: int = 0, 
    maskfile=None,
    gen_intvs=None,
    use_ints=False): # silence worker progress

    if isinstance(fixlist, str):
        fixlist = [fixlist,]

    # WORK OUT THE SUFFIXES.NO .BAMS YET
    f1 = fixlist[0]
    if f1.endswith('.vcf.gz') or f1.endswith('.vcf'):
        fixfiles = fixlist
    elif Path(f1 + '.vcf.gz').exists():
        fixfiles = [f + '.vcf.gz' for f in fixlist]
    elif Path(f1 + '.vcf').exists():
        fixfiles = [f + '.vcf' for f in fixlist]
    else:
        raise ValueError(f"Unknown file format for {f1}. Needs to be '.vcf.gz' or '.vcf'") # not yet implementing .bcf as issues with vcf_core

    genints = gen_intvs if use_ints else [None] * len(fixfiles) # if not using genints, just pass None to workers

    # run workers
    results = [None] * len(fixfiles)
    sample_ref = None
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
            else:
                # hard requirement: sample order identical across files
                if list(samples) != sample_ref:
                    raise ValueError(
                        f"Sample list/order mismatch in {fixfiles[i]}.\n"
                        f"First few: {samples[:5]} vs {sample_ref[:5]}"
                    )

    print(f"Processed {len(fixfiles)} files. Variants kept: {total_rep:,}. Variants skipped: {total_skip:,}.")

    # merge in file order
    chrom_dict = _merge_chrom_dicts(results)

    return chrom_dict, sample_ref

def get_outliers_from_chrom_dict_normalized(
    chrom_dict,
    window=10_000,
    threshold_frac=0.05,
    alpha=1e-6,
    min_sites=20,
    window_depth_thresh=5.0,
    use_median=True,
    min_spread=1e-6,
    transform="none",
    log_pseudocount=0.5
):
    """
    Detect depth outlier windows using a fully normalized depth scale.

    Model scale:
        Y[w, i] = k[w, i] / (E[w] * s[i])

    where:
        k = raw summed depth in a window
        E = usable site count in a window
        s = per-individual depth scale factor

    All comparisons, z-scores, and cutoffs are made on Y.
    """

    # -----------------------------
    # 1) stack arrays across chromosomes
    # -----------------------------
    glob_depths = np.vstack([chrom_dict[chrom]["depths"] for chrom in chrom_dict]).astype(float)
    glob_sitecounts = np.concatenate(
        [chrom_dict[chrom]["sitecounts"] for chrom in chrom_dict]
    ).astype(float)

    k = glob_depths                      # (n_windows, n_idv)
    E = glob_sitecounts                  # (n_windows,)

    E_safe = E.copy()
    E_safe[E_safe <= 0] = np.nan

    # -----------------------------
    # 2) estimate per-individual scale factors
    # -----------------------------
    depth_per_site = k / E_safe[:, None]   # raw per-site depth
    usable_for_scale = (E >= min_sites) & np.isfinite(E_safe)

    if use_median:
        idv_depth = np.nanmedian(depth_per_site[usable_for_scale, :], axis=0)
        scale_ref = np.nanmedian(idv_depth)
    else:
        idv_depth = np.nanmean(depth_per_site[usable_for_scale, :], axis=0)
        scale_ref = np.nanmean(idv_depth)

    idv_scale = idv_depth / scale_ref

    # protect against bad individuals
    bad_scale = ~np.isfinite(idv_scale) | (idv_scale <= 0)
    idv_scale[bad_scale] = np.nan

    # -----------------------------
    # 3) normalized depth matrix
    # -----------------------------
    Y = k / (E_safe[:, None] * idv_scale[None, :])


    window_center_raw = np.nanmedian(Y, axis=1) if use_median else np.nanmean(Y, axis=1)
    usable = (
        (E >= min_sites)
        & np.isfinite(window_center_raw)
        & (window_center_raw >= window_depth_thresh)
    )
    # TRANSFORM
    Yt = _transform_depth(Y, transform=transform, log_pseudocount=log_pseudocount)

    # -----------------------------
    # 4) per-window center
    # -----------------------------
    if use_median:
        mu = np.nanmedian(Yt, axis=1)
    else:
        mu = np.nanmean(Yt, axis=1)

    # -----------------------------
    # 5) per-window spread
    # -----------------------------
    if use_median: # we will use Median Absolute Deviation in this case
        # robust sigma estimate: 1.4826 * MAD
        abs_dev = np.abs(Yt - mu[:, None])
        sigma = 1.4826 * np.nanmedian(abs_dev, axis=1)
    else: # SD used if mean is used
        sigma = np.nanstd(Yt, axis=1, ddof=1) # ddof 1 for sample stddev

    sigma = np.where(np.isfinite(sigma) & (sigma > min_spread), sigma, min_spread)

    # -----------------------------
    # 6) z-scores and p-values, all on normalized scale
    # -----------------------------
    z = (Yt - mu[:, None]) / sigma[:, None]

    p_low = norm.cdf(z)
    p_high = norm.sf(z)

    # normalized cutoffs
    y_lower_cut = mu + norm.ppf(alpha) * sigma
    y_upper_cut = mu + norm.isf(alpha) * sigma

    # -----------------------------
    # 7) individual-level failures
    # -----------------------------
    fail_lower = (p_low <= alpha) & usable[:, None]
    fail_upper = (p_high <= alpha) & usable[:, None]
    fail_array = fail_lower | fail_upper

    # cohort-level fractions
    fail_frac_lower = np.nanmean(fail_lower, axis=1)
    fail_frac_upper = np.nanmean(fail_upper, axis=1)
    fail_frac = np.nanmean(fail_array, axis=1)

    upper_fail = fail_frac_upper >= threshold_frac
    lower_fail = fail_frac_lower >= threshold_frac
    fail_windows = np.where(upper_fail | lower_fail)[0]

    # -----------------------------
    # 8) convert failing windows back to genomic positions
    # -----------------------------
    fail_window_positions = []

    cum_window_counts = []
    cum_count = 0
    for chrom in chrom_dict:
        nwin = len(chrom_dict[chrom]["sitecounts"])
        start_idx = cum_count
        cum_count += nwin
        cum_window_counts.append((chrom, cum_count, start_idx))

    for fw in fail_windows:
        for chrom, cum_end, start_idx in cum_window_counts:
            if fw < cum_end:
                win_dict = chrom_dict[chrom]
                local_idx = fw - start_idx
                wstart = win_dict["window_indices"][local_idx] * window + 1

                fail_upper_list = np.where(fail_upper[fw])[0]
                fail_lower_list = np.where(fail_lower[fw])[0]

                fail_window_positions.append(
                    {
                        "chrom": chrom,
                        "start": wstart,
                        "end": wstart + window - 1,
                        "window_index_global": int(fw),
                        "window_index_chrom": int(local_idx),
                        "fail_lower_list": fail_lower_list,
                        "fail_upper_list": fail_upper_list,
                        "fail_frac": float(fail_frac[fw]),
                        "fail_frac_lower": float(fail_frac_lower[fw]),
                        "fail_frac_upper": float(fail_frac_upper[fw]),
                        "mu_trans": float(mu[fw]),
                        "mu_norm": float(_inverse_transform_depth(mu[fw], transform=transform, log_pseudocount=log_pseudocount)),
                        "sigma_trans": float(sigma[fw]),
                        "lower_cut_trans": float(y_lower_cut[fw]),
                        "upper_cut_trans": float(y_upper_cut[fw]),
                        "lower_cut_norm": float(_inverse_transform_depth(y_lower_cut[fw], transform=transform, log_pseudocount=log_pseudocount)),
                        "upper_cut_norm": float(_inverse_transform_depth(y_upper_cut[fw], transform=transform, log_pseudocount=log_pseudocount)),
                    }
                )
                break

    print(
        f"{'Median' if use_median else 'Mean'} outlier fraction per window:",
        np.nanmedian(fail_frac) if use_median else np.nanmean(fail_frac)
    )
    print("95th percentile outlier fraction per window:", np.nanpercentile(fail_frac, 95))

    return {
        "fail_window_positions": fail_window_positions,
        "Y": Y,
        "mu": mu,
        "sigma": sigma,
        "z": z,
        "p_low": p_low,
        "p_high": p_high,
        "y_lower_cut": y_lower_cut,
        "y_upper_cut": y_upper_cut,
        "idv_scale": idv_scale,
        "usable": usable,
        "fail_frac": fail_frac,
        "fail_frac_lower": fail_frac_lower,
        "fail_frac_upper": fail_frac_upper,
    }


def write_outlier_windows(fail_window_positions, outfix, samples, window=10_000, use_median=True):
    """Write outlier windows to file.

    Compatible with the normalized-scale outlier function, where each entry in
    fail_window_positions is a dict with keys such as:
        chrom, start, end,
        fail_lower_list, fail_upper_list,
        fail_frac, fail_frac_lower, fail_frac_upper,
        mu_norm, sigma_norm,
        lower_cut_norm, upper_cut_norm
    """
    idvs = list(samples)

    with open(outfix + '.windows_final', 'w') as ofi, \
         open(outfix + '.windows_low', 'w') as ofil, \
         open(outfix + '.windows_high', 'w') as ofih:

        center_label = "median" if use_median else "mean"

        ofi.write(
            f'chrom_position\twindow_start\tfrac\tlow_frac\thigh_frac\t'
            f'norm_depth_{center_label}\tnorm_depth_lower_cut\tnorm_depth_upper_cut\t'
            f'low_depth_ids\thigh_depth_ids\n'
        )
        ofil.write(
            f'chrom_position\twindow_start\tfrac\tlow_frac\t'
            f'norm_depth_{center_label}\tnorm_depth_lower_cut\tnorm_depth_upper_cut\tlow_depth_ids\n'
        )
        ofih.write(
            f'chrom_position\twindow_start\tfrac\thigh_frac\t'
            f'norm_depth_{center_label}\tnorm_depth_lower_cut\tnorm_depth_upper_cut\thigh_depth_ids\n'
        )

        for fwp in fail_window_positions:
            chrom = fwp["chrom"]
            wstart = int(fwp["start"])
            wend = int(fwp.get("end", wstart + window - 1))

            fail_lower_list = fwp["fail_lower_list"]
            fail_upper_list = fwp["fail_upper_list"]

            fail_frac = float(fwp["fail_frac"])
            fail_frac_lower = float(fwp["fail_frac_lower"])
            fail_frac_upper = float(fwp["fail_frac_upper"])
            depth_center_window = float(fwp["mu_norm"])
            depth_lower_cut_window = float(fwp["lower_cut_norm"])
            depth_upper_cut_window = float(fwp["upper_cut_norm"])

            lower_ids = ','.join(idvs[i] for i in fail_lower_list)
            upper_ids = ','.join(idvs[i] for i in fail_upper_list)

            chrom_pos = f'{chrom}:{wstart}-{wend}'

            ofi.write(
                f'{chrom_pos}\t{wstart}\t{fail_frac:.2f}\t'
                f'{fail_frac_lower:.2f}\t{fail_frac_upper:.2f}\t{depth_center_window:.4f}\t'
                f'{depth_lower_cut_window:.4f}\t{depth_upper_cut_window:.4f}\t{lower_ids}\t{upper_ids}\n'
            )

            if len(fail_lower_list) > 0:
                ofil.write(
                    f'{chrom_pos}\t{wstart}\t{fail_frac:.2f}\t'
                    f'{fail_frac_lower:.2f}\t{depth_center_window:.4f}\t'
                    f'{depth_lower_cut_window:.4f}\t{depth_upper_cut_window:.4f}\t{lower_ids}\n'
                )

            if len(fail_upper_list) > 0:
                ofih.write(
                    f'{chrom_pos}\t{wstart}\t{fail_frac:.2f}\t'
                    f'{fail_frac_upper:.2f}\t{depth_center_window:.4f}\t'
                    f'{depth_lower_cut_window:.4f}\t{depth_upper_cut_window:.4f}\t{upper_ids}\n'
                )

def get_depths_outliers_parallel(
    fixlist: list,
    outfix,
    window=10_000,
    n_jobs=8,
    min_miss_frac=0.5,
    interpolate=False,
    update_interval=100_000,
    maskfile=None,
    gen_intvs=None,
    use_ints=False,
    threshold_frac=0.05,
    alpha=1e-6,
    min_sites=20,
    window_depth_thresh=5,
    use_median=True,
    min_spread=0.01,
    transform="none",
    log_pseudocount=0.5,
    write=True,
):
    """
    Compute normalized depth outlier windows across VCF prefixes in parallel.

    Workflow order
    --------------
    1. Aggregate per-window depth information from VCF files.
    2. Convert depths to a normalized per-site, per-individual scale.
    3. Detect outlier windows based on normalized depth deviations.
    4. Optionally write outlier windows to output files.

    Parameters
    ----------
    fixlist : list
        List of VCF prefixes (without `.vcf.gz`).

    outfix : str
        Output prefix. If `write=True`, output files will be written using this prefix.

    window : int, default 10_000
        Size (bp) of non-overlapping windows used for depth aggregation.

    n_jobs : int, default 8
        Number of worker processes to use for parallel VCF processing.

    min_miss_frac : float, default 0.5
        Sites with missingness above this threshold are skipped during depth aggregation.

    interpolate : bool, default False
        Whether to interpolate missing windows during the initial sweep.

    update_interval : int, default 100_000
        Interval for progress updates during depth extraction.

    maskfile : str or None, default None
        GFF3 file containing genomic regions to mask.

    gen_intvs : list or None, default None
        List of genomic intervals of the form `CHROM:START-END` to restrict analysis to.
        Only relevant if `maskfile` is provided and `use_ints=True`.

    use_ints : bool, default False
        If True, restrict analysis to `gen_intvs`. If False, process the full genome
        apart from masked regions.

    threshold_frac : float, default 0.05
        Minimum fraction of outlier individuals required to call a window an outlier.

    alpha : float, default 1e-6
        Per-individual tail probability cutoff used for outlier calling.

    min_sites : int, default 20
        Minimum number of usable sites required for a window to be tested.

    window_depth_thresh : float, default 5
        Minimum normalized window depth required for a window to be tested.

    use_median : bool, default True
        If True, use median-based summaries; otherwise use mean-based summaries.

    min_spread : float, default 0.01
        Minimum allowed spread on the normalized depth scale. Prevents unrealistically
        tight cutoffs when depth is extremely consistent across individuals.

    transform : str, default "none"
        Optional transformation to apply to normalized depths before outlier detection.
        One of "none", "log", or "sqrt".

    log_pseudocount : float, default 0.5
        Pseudocount added to normalized depths before log transformation to prevent log(0).

    write : bool, default True
        Whether to write output files.

    Returns
    -------
    fail_outliers : dict
        Output from `get_outliers_from_chrom_dict_normalized()`.

    chrom_dict : dict
        Per-chromosome aggregated depth information.

    samples : list
        Sample IDs in the same order used throughout the analysis.
    """
    chrom_dict, samples = get_depths_dict_parallel(
        fixlist=fixlist,
        window=window,
        n_jobs=n_jobs,
        min_miss_frac=min_miss_frac,
        interpolate=interpolate,
        update_interval=update_interval,
        maskfile=maskfile,
        gen_intvs=gen_intvs,
        use_ints=use_ints,
    )

    fail_outliers = get_outliers_from_chrom_dict_normalized(
        chrom_dict=chrom_dict,
        window=window,
        threshold_frac=threshold_frac,
        alpha=alpha,
        min_sites=min_sites,
        window_depth_thresh=window_depth_thresh,
        use_median=use_median,
        min_spread=min_spread,
        transform=transform,
        log_pseudocount=log_pseudocount,
    )

    if write:
        write_outlier_windows(
            fail_window_positions=fail_outliers["fail_window_positions"],
            outfix=outfix,
            samples=samples,
            window=window,
            use_median=use_median,
        )

    return fail_outliers, chrom_dict, samples

    
