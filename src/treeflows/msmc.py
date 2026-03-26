from treeflows.vcf_core import VCFReader
from treeflows import config, refdata, _msmc2tools_functions as _mtools
import random
import subprocess
from pathlib import Path

_default_segmentation = "1*2+25*1+1*2+1*3"
_near_forced_segmentation = "1*2+25*1+1*5+1*5+1*5"
_short_segmentation = "1*2+15*1+1*2+1*3"
_ultranear_segmentation = "1*2+30*1+1*5+1*5+1*5+1*10"

_segdict = {
    "default": _default_segmentation,
    "near_forced": _near_forced_segmentation,
    "short": _short_segmentation,
    "ultranear": _ultranear_segmentation
}

def vcf_to_msmc(vcffile, out_prefix):
    """Convert a VCF into MSMC2 input format.

    Writes tab-delimited lines of the form `chrom  pos  sites  hap_string` where
    `hap_string` contains allele characters for each haplotype (or `?` for
    missing).

    Args:
        vcffile: Input VCF path.
        out_prefix: Output filename to write MSMC2 input to.

    Returns:
        None.
    """
    vcf = VCFReader(vcffile)

    sites = 0

    with open(out_prefix, "w") as ofi:

        for site in vcf:
            sites += 1
            if not site.is_biallelic or site.is_nonallelic:
                continue

            genotypes = [h for g in site.genotypes for h in g[:2]]
            genos = "".join([str(site.vcf_alleles[g]) if g >= 0 else '?' for g in genotypes])
            outline = f"{site.CHROM}\t{site.POS}\t{sites}\t{genos}\n"
            ofi.write(outline)
            sites = 0


def msmc2_internal(pop, infile=None, out_prefix=None, n_idv=4, threads=1, short=True, test=False, ref=None, vcf_ref=None, conf=None, segmentation="default"):
    """Run MSMC2 within a population by subsampling individuals.

    Args:
        pop: Population label (as used in the reference metadata).
        infile: MSMC2 input file path.
        out_prefix: Output prefix for MSMC2 results.
        n_idv: Number of individuals to sample.
        threads: Number of threads passed to MSMC2.
        short: If True, use the "short" population label accessor.
        test: If True, return the `-I` argument string without running MSMC2.
        ref: Reference metadata CSV filename.
        vcf_ref: VCF reference file.
        conf: FileConfig configuration object.
        segmentation: Segmentation string to use for MSMC2.

    Returns:
        If `test=True`, returns the MSMC2 `-I` argument string; otherwise None.
    """
    vcf_ref = "oct_aeg_apacout.vcf" if vcf_ref is None else vcf_ref
    vcf_idvs = VCFReader(vcf_ref).samples # will correspond to pop_short
    ref = "oct_aeg_apacout.csv" if ref is None else ref
    conf = config.FileConfig.from_yaml() if conf is None else conf
    aeg = refdata.load_aegdata(conf.refdir / ref)
    # reorder to only include individuals in the VCF and ensure the same order(via pop_short column)
    aeg = aeg[aeg["id"].isin(vcf_idvs)].reset_index(drop=True)
    aeg = aeg.set_index("id").loc[vcf_idvs].reset_index()
    print(aeg)
    popids = aeg.get_idx_pop_short(pop) if short else aeg.get_idx_pop(pop) 
    popids = sorted(random.sample(popids, min(len(popids), n_idv)))
    pophaps = aeg.get_nodes_idx(popids)
    popstring = ",".join(str(i) for i in pophaps)
    if test:
        return popstring
    subprocess.run(f"msmc2 -r 1 -t {threads} -I {popstring} -p {_segdict[segmentation]} -o {out_prefix} {infile}", shell=True)
    #return popstring

def msmc2_internal_unphased(pop, infile=None, out_prefix=None, n_idv=4, threads=1, short=True, test=False, ref=None, conf=None, segmentation="default"):
    """Run MSMC2 within a population using unphased individual pair strings.

    Args:
        pop: Population label.
        infile: MSMC2 input file path.
        out_prefix: Output prefix for MSMC2 results.
        n_idv: Number of individuals to sample.
        threads: Number of threads passed to MSMC2.
        short: If True, use the "short" population label accessor.
        test: If True, return the `-I` argument string without running MSMC2.
        ref: Reference metadata CSV filename.
        conf: FileConfig configuration object.
        segmentation: Segmentation string to use for MSMC2.

    Returns:
        If `test=True`, returns the MSMC2 `-I` argument string; otherwise None.
    """
    ref = "oct_aeg_apacout.csv" if ref is None else ref
    conf = config.FileConfig.from_yaml() if conf is None else conf
    aeg = refdata.load_aegdata(conf.refdir / ref)
    popids = aeg.get_idx_pop_short(pop) if short else aeg.get_idx_pop(pop) 
    popids = sorted(random.sample(popids, min(len(popids), n_idv)))
    pophaps = ["-".join(str(x) for x in aeg.get_nodes_idx([popid])) for popid in popids]
    
    popstring = ",".join(str(i) for i in pophaps) ### CURRENTLY EMPTY... 
    popstring = ",".join(str(i) for i in popids) ### CURRENTLY EMPTY... 
    if test:
        return popstring
    subprocess.run(f"msmc2 -r 1 -t {threads} -I {popstring} -p {_segdict[segmentation]} -o {out_prefix} {infile}", shell=True)
    #return popstring

def msmc2_between(pop1, pop2, infile=None, out_prefix=None, n_idv=4, threads=1, short=True, test=False, ref=None, vcf_ref=None, conf=None, segmentation="default"):
    """Run MSMC2 between two populations by pairing sampled haplotypes.

    Args:
        pop1: First population label.
        pop2: Second population label.
        infile: MSMC2 input file path.
        out_prefix: Output prefix for MSMC2 results.
        n_idv: Max individuals to sample from each population.
        threads: Number of threads passed to MSMC2.
        short: If True, use the "short" population label accessor.
        test: If True, return the `-I` argument string without running MSMC2.
        ref: Reference metadata CSV filename.
        vcf_ref: VCF reference file.
        conf: FileConfig configuration object.
        segmentation: Segmentation string to use for MSMC2.


    Returns:
        If `test=True`, returns the MSMC2 `-I` argument string; otherwise None.
    """
    vcf_ref = "oct_aeg_apacout.vcf" if vcf_ref is None else vcf_ref
    vcf_idvs = VCFReader(vcf_ref).samples # will correspond to pop
    ref = "oct_aeg_apacout.csv" if ref is None else ref
    conf = config.FileConfig.from_yaml() if conf is None else conf
    aeg = refdata.load_aegdata(conf.refdir / ref)
    # reorder to only include individuals in the VCF and ensure the same order(via pop column)
    aeg = aeg[aeg["id"].isin(vcf_idvs)].reset_index(drop=True)
    aeg = aeg.set_index("id").loc[vcf_idvs].reset_index()
    ids1 = aeg.get_idx_pop_short(pop1) if short else aeg.get_idx_pop(pop1)
    ids2 = aeg.get_idx_pop_short(pop2) if short else aeg.get_idx_pop(pop2)
    ids1 = sorted(random.sample(ids1, min(len(ids1), n_idv)))
    ids2 = sorted(random.sample(ids2, min(len(ids2), n_idv)))
    haps1 = aeg.get_nodes_idx(ids1)
    haps2 = aeg.get_nodes_idx(ids2)
    pairings = [f"{h1}-{h2}" for h1 in haps1 for h2 in haps2]
    #pairings = [f"{i1}-{i2}" for i1 in ids1 for i2 in ids2]  # paired by individual
    popstring = ",".join(pair for pair in pairings) # CURRENTLY EMPTY
    print(aeg)
    print("ids1: ", ids1)
    print("ids2: ", ids2)
    #print("haps 1: ", haps1)
    #print("Haps 2: ", haps2)
    print("Pairings: ", pairings)
    print("Threads: ", threads)
    print("idvs: ", popstring)
    print("out: ", out_prefix)
    print("in: ", infile)
    print(f"msmc2 -r 1 -t {threads} -I {popstring} -p {_segdict[segmentation]} -o {out_prefix} {infile}")
    if test:
        return popstring
    process = subprocess.run(f"msmc2 -r 1 -t {threads} -I {popstring} -p {_segdict[segmentation]} -o {out_prefix} {infile}", shell=True, 
                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    print("OUT: ", process.stdout)
    print("ERR: ", process.stderr)
    #return popstring

def msmc2_combinepops(poplist, in_dir=None, in_prefix="ms_", sep="-", outfile=None, append=True):
    """Combine MSMC2 output files across multiple populations.

    Args:
        poplist: List of population labels.
        in_dir: Input directory for MSMC2 result files.
        in_prefix: Input prefix for MSMC2 result files.
        sep: Separator between population labels in MSMC2 filenames.
        outfile: Output filename for combined results.
        append: If True, append to existing output file.

    Returns:
        Combined MSMC2 cross-coalescence string.
    """ 
    combolist = []
    if outfile is not None:
        mode = "a" if append else "w"
        ofi = open(outfile, mode)
        if not append:
            ofi.write("time_index\tleft_time_boundary\tright_time_boundary\tlambda_00\tlambda_01\tlambda_11\tpop0\tpop1\n")

    in_dir = Path(in_dir) if in_dir is not None else Path.cwd()
    for pop1 in poplist:
        in_file1 = Path(f"{in_dir}/{in_prefix}{pop1}.final.txt")
        if not in_file1.exists():
            print(f"Missing MSMC2 output file for population {pop1}: {in_file1}")
            continue
        for pop2 in poplist:
            if pop1 == pop2:
                continue
            in_file2 = Path(f"{in_dir}/{in_prefix}{pop2}.final.txt")
            joint_file = Path(f"{in_dir}/{in_prefix}{pop1}{sep}{pop2}.final.txt")
            if not joint_file.exists():
                joint_file = Path(f"{in_dir}/{in_prefix}{pop2}{sep}{pop1}.final.txt")
            if not in_file2.exists() or not joint_file.exists():
                print(f"Missing an MSMC2 output file combination for population {pop1} ({pop2}).")
                continue
            ccc = _mtools.get_ccc(joint_file, in_file1, in_file2, pop1, pop2)
            if outfile is not None:
                ofi.write(ccc + "\n")
                combolist.append(f"Written combination {pop1}-{pop2}")
            else:
                combolist.append(ccc)
    
    if outfile is not None:
        ofi.close()
    return "\n".join(combolist)

def msmc2_combinepops_all(in_dir=None, in_prefix="ms_", sep="-", outfile=None, short=True, ref=None, conf=None):
    """Combine MSMC2 output files across all populations found in input prefix.

    Args:
        in_dir: Input directory for MSMC2 result files (default: current working directory).
        in_prefix: Input prefix for MSMC2 result files (default: "ms_").
        sep: Separator between population labels in MSMC2 filenames (default: "-").
        outfile: Output filename for combined results (default: "aggregated.crosscoal.txt").
        short: If True, use the "short" population label accessor (default: True).
        ref: Reference metadata CSV filename (default: "oct_aeg_apacout.csv").
        conf: FileConfig configuration object (default: from local yaml config).
    """ 
    in_dir = Path.cwd() if in_dir is None else Path(in_dir)
    ref = "oct_aeg_apacout.csv" if ref is None else ref
    conf = config.FileConfig.from_yaml() if conf is None else conf
    aeg = refdata.load_aegdata(conf.refdir / ref)
    combolist = []
    outfile = in_dir / "aggregated.crosscoal.txt" if outfile is None else Path(outfile)
    popnames = aeg["pop_short"] if short else aeg["pop"]
    popnames = sorted(list(set([p for p in popnames])))
    print("Initial pops: ", popnames)
    msmc2_combinepops(popnames, in_dir, in_prefix, sep, outfile, append=False)
    
    return "\n".join(combolist)