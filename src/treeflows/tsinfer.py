##### Package Imports #####

# general
from pathlib import Path
from treeflows.Spath import Mpath
import subprocess
import shutil
from cyvcf2 import VCF

from bio2zarr.vcf import explode as v_explode, encode as v_encode # updated version
from bio2zarr.cli import check_overwrite_dir as v_check_overwrite_dir

#NOTE icf module being renamed to vcf module

import tsinfer, tskit, tsdate

"""
Implementation notes: 

FORMAT: 
Need to get into the .vcz zarr format. THis has turned out to be tricky because I need to balance file number on the 
one hand with not crashing zarr via numcodecs chunk takes too much memory error. 
At the moment we are being prudent (ish) with 10k idv chunk size, but it is still gobbling up 120k files, meaning we 
need to be careful before we do anything drastic. I can deal with this via a very very careful zip after conversion, 
but can't have too many of these running at once or we'll crash everything. It would perhaps be good to pre-calculate this? 

A future possibiity might be to add to zar in stages, but I suspect this won't work so well... zip archives can technically be 
edited in stages, but it is unclear what would happen after that... would need a way of writing zar directly to zip

# note that my typical problem is that files are zipped and have indexes added at the end. These are the things that we need to deal with in the long run... there is normally only
one 'core' suffix. the others are ancilliary. 

NOTE: if a site mask is est, the ancestral states array should only specify alleles for the unmasked sites... 

"""

def strip_biosuffixes(p, getlist=False) -> Path:
    """Strip common bioinformatics suffixes from a path.

    Args:
        p: Input path-like.
        getlist: If True, return `(stripped_path, removed_suffixes)`.

    Returns:
        If `getlist=True`, `(Path, list[str])`; otherwise the stripped `Path`.
    """

    p = Path(p)
    suffixlist = ['.gz', '.vcf', '.bcf', '.icf', '.vcz', '.bam', '.sam', '.fa', '.fasta', '.bai', '.fai', '.tbi', '.csi', '.zip', '.trees'] # add more later if necessary
    striplist = []
    while True:
        if p.suffix in suffixlist:
            striplist.append(p.suffix)
            p = p.with_suffix('')
        else:
            break
    if getlist:
        return p, striplist
    return p



#### VCF2ZARR - a program, so we will try and run it from the commanline, meaning subprocess is the way to go... 

    # we will work with n(individuals) - which we have to set manually - and N(sites) - which we can make 100_000

def vcf_to_bcf_idx(vcf, rem=False):
    '''Converts .vcf(.gz) file to .bcf file and indexes it (with bcftools) '''
    vcf_r = strip_biosuffixes(vcf) # a path
    vcf_o = vcf_r.with_suffix('.bcf')
    
    subprocess.run(f'bcftools view -Ob -o {vcf_o} {vcf}', shell=True)
    subprocess.run(f'bcftools index {vcf_o}', shell=True)
    if rem:
        Path(vcf).unlink()
    
 
def run_vcf2zarr(invcf, outicf = None, outvcz=None, varchunks=100_000, zip=False, rem_icf = True, temp = '/tmp/', threads=1, show_progress=False, force=True): # use temp directory? 
    '''
    Runs `vcf2zarr explode & encode` processes on .bcf (or .vcf) file'''

    invcf = Path(invcf) # maybe go and use Mpath? (though I think there were vulnerabilities)... Spath??? 
    vcfs = [invcf,]
    vcf = VCF(invcf)
    n_idv = len(vcf.samples)
    
    core = strip_biosuffixes(invcf)
    icf = core.with_suffix('icf') if outicf is None else Path(outicf)

    
    # explode
    v_check_overwrite_dir(icf, force)
    v_explode(icf, vcfs, worker_processes=threads, show_progress=True) # may work... 
    if outvcz is None:
        outvcz = Path(temp) / core.with_suffix('.vcz')

    if zip:
        outvcz = strip_biosuffixes(outvcz, False).with_suffix('.vcz.zip')
        tvcz = Path(temp) / ( core.with_suffix('.vcz').name) 
        print(outvcz)
        v_check_overwrite_dir(tvcz, force)
        v_encode(icf, tvcz, variants_chunk_size=varchunks, samples_chunk_size=n_idv, worker_processes = threads, show_progress=show_progress)
        print('zipping to final loc')
        subprocess.run(f'(cd {tvcz} && pwd && zip -rq {outvcz.resolve()} .)', shell=True)
    else:
        v_check_overwrite_dir(outvcz, force)
        v_encode(icf, outvcz, variants_chunk_size=varchunks, samples_chunk_size=n_idv, worker_processes = threads, show_progress=show_progress)

    if rem_icf:
        shutil.rmtree(icf)
    print('done')

def infer_tree_basic(vcz, treename=None, threads=0): # try and run in one wild dash!!!
    """Infer a tree sequence from a bio2zarr `.vcz` archive and write `.trees`."""
    vcz=Path(vcz)
    root = strip_biosuffixes(vcz, False)
    if treename is None: treename = root.with_suffix('.trees')

    # ts infer load... 
    vdata = tsinfer.VariantData(vcz, ancestral_state = 'variant_AA') # eventually load dates here... 
    inferred_ts = tsinfer.infer(vdata, num_threads=threads) # fancy settings later... 
    simp_ts = inferred_ts.simplify()
    simp_ts.dump(treename)
    return simp_ts

def infer_tree_params(vcz, treename=None, recombination_rate=4e-7, mismatch=0.01, threads=0):
    """Infer a tree sequence from `.vcz` with explicit recombination/mismatch params."""
    vcz = Path(vcz)
    root = strip_biosuffixes(vcz, False)
    if treename is None: treename = root.with_suffix('.trees')
    vdata = tsinfer.VariantData(vcz, ancestral_state = 'variant_AA')
    inferred_ts = tsinfer.infer(vdata, recombination_rate=recombination_rate, mismatch_ratio=mismatch, num_threads=threads)# gonna just simplify
    simp_ts = inferred_ts.simplify()
    simp_ts.dump(treename)
    return simp_ts


def infer_tree_multibatch():
    """Placeholder for a batched tsinfer pipeline (not yet implemented)."""
    # Ancestor Generation. Lots of cores. Limiting factor is gentype array in RAM... high water mark for memory
    # Ancestor Matching... 
        # parallelism limited to ancestors whose inheritors are already matched... (all inheritors must be matched in earlier group)... 
    #tsinfer.match_ancestors_batch_init() # -> metadata.json to work_dir... (hmm)... + ancestor groupinng... 
    pass

def date_tree(ts: Path|str|tskit.TreeSequence, mutation_rate, out: Path|str = None, preprocess = False, match_segs=True, eps=1e-10) -> tskit.TreeSequence:
    """Dates tree sequence. If `match segs`, runs the method appropriate for estimating coalescence rates
        Default eps is 1e-10. """
    if isinstance(ts, Path) or isinstance(ts, str):
        ts = tskit.load(ts)
    if preprocess:
        ts = tsdate.preprocess_ts(ts)
    ts_dated = tsdate.variational_gamma(ts, mutation_rate = mutation_rate, match_segregating_sites=match_segs, eps=eps)
    if out is not None:
        ts_dated.dump(out)
    return ts_dated
    




if __name__ == '__main__':

    cwd = Path('__file__').parent.resolve() # not an absolute path!!!
    invcf = cwd / 'mj_init_may24_apacest_i10_chr2a_merged.vcf.gz'
    inbcf = cwd / 'mj_init_may24_apacest_i10_chr2a_merged.bcf'
    dir1 = cwd / 'nodate' / 'tsinfer'
    outicf = dir1 / 'chrom2a_apacest_i10.icf'
    outvcz = dir1 / 'chrom2a_apacest_i10.vcz.zip'
    outtree = dir1 / 'inf' / 'chrom2a_apacest_i10_inf.trees'
    print(outvcz)
    vcf_to_bcf_idx(invcf)
    run_vcf2zarr(inbcf, outicf, outvcz, threads=40, zip=True, show_progress=True, varchunks=10_000)
    infer_tree_basic(outvcz, treename=outtree, threads=40)
