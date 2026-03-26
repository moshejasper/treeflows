##### Package Imports #####

from pathlib import Path
from cyvcf2 import VCF, Writer
import subprocess
import numpy as np

path2singer = 'x'

'''
IMPLEMENTATION NOTES.... 

We basically cant do anything with missing data here, so we need to impute the whole arm. This is a step in the prior process... 

'1. impute the missing data following current protocols (import from prior, & make sure it all works)

Have been no updates for 8 months, so can comfortably use the current versions. 

NOTES FROM PAGE... 
parallel breaks into 1MB chunks, runs on them all, then reconstructs (different from others)... 
IT doesn't enable appropriate breaks at ends, but otherwise looks good.... 

Possible to scale the mutation rate to reflect underlying data presence? (a good idea?).... (and recombination map also...)... 
need to set -polar to 0.99 (or thereabouts)... possibly need to remove non-ancestral sites... 

parallele singer is an interesting way to go... 

Use Ne, m, ratio... pi = 4 * Ne * m. (pi here is nucleotide diversity)... --> need to calculate nuc divers... (how)... (possibly windows?)... 
(get rid of regions with no sequenced sites) -> characterise missingness... (pairwise diversity)... 

Probalby use PIXY to calcualte pi... --> PRAGMATIC: get the pipeline running, then worry about this!... (a general issue)... 

"only support phased, high quality genomes. polymorphic sites with missingness will be excluded... "
probably needs chromosome with only one name in header... but will see. 

# has a default 'compute_Ne' function (which runs if you don't specify it) which basically looks at diversity over window (assuming all sites present) and adds up
'diversity' in non-mono-allelic sites. It works by assuming all sites present in data (possibly this is far more tractible than alternatives?). There is a suggestion that I need to filter my 
data to just variant stuff to help these programs run easier... it is very unclear to me if this is true or not. But at least with Singer, it seems pointless to include sites that aren't variant
in the end... (makes the files much more tractible)... probably pipeline looks like: (1) strip invariant/multivar/ancestry-indeterminate; (2) run in chunks; (3) connect at the end... (some kind of boots approach)
The big risk is that the 'chunkiness' will lead to poor results. (probably a risk in all analyses...)... It can ultimately be dealt with via a recombination rate adjuster... (across windows)... 

Idea would be you have 'default rate'... then you adjust it according to a 'missing data mask'... or something. (we could look into this with the existing recombination map...)... (probably a combined approach)
So then you basically rescale the recombination map in blocks of e.g. 1000... (with or without 'base map') to guide things... this would then sit in the background with the singer map... a more 'fine-tuned' 
approach to the missing data question. Still effectively worried about elimination of local signals of selection... but this is probably more of an issue with mutation... 
 -> so need to keep the mutation parameter scaled 
 If mut map provided, basically Ne will be computed from average rate across that section of map, then used to calculate Ne across the same section... (unclear how this works with pop structure)... 

 FOR NOW: will use default rate (somewhat rescaled) -> parallele SINGER -> recombine... 

 PARALLEL SINGER (index)... note that this doesn't currently run wiht maps... (mut/rec ratio also exists)... hmm... NOTE: abandon non-parallel singer... 

 QUOTES FROM PAPER. 
 As for MCMC sampling, we uniformized the number of iterations and thinning scheme across all methods. We drew 100 samples with the thinning interval set to 20 for ARGweaver, Relate, 
  and SINGER, and used 1,000 iterations for burn-in.

For all simulation benchmarks involving ARGweaver and SINGER, posterior averages were taken for the statistics of interest, such as pairwise coalescence time, allele age, etc. 
 Since Relate outputs averages of MCMC iterations while tsdate results are averages from a probability table, we simply use their results from a single output as they are 
  effectively posterior averages.

  This means that we need to let thin*sample = 4000 (the SINGER burning), then take 100 samples, thinning every 60 its (i.e. 6000 its + 4000 its)... Definitely need to break the chomosome up... 

At this point, running is better than not. So we are going to try and just implement parallel singer on imputed data... VCF format... (as we probably can't get the rest to run without failing... )
Try and replicate vcf index... 

NOTE: THE FILE SEEM TO REQUIRE OPEN, UNZIPPED VCF FILES!!!... 

NOTE2... as it doesn't talk aobut this anywhere, it is likely that SINGER assumes ref (0) is ancestral allele. Will need to recode in some setting for ths to be true. 
Can run the 'random polarity' one first... 

'''

def flip_biallelic_genotype(gt):
    phased = gt[-1]
    rest = gt[:2]
    flipped = [1 if a==0 else 0 if a==1 else -1 for a in rest] # 
    return flipped + [phased]

def vcf_reorder_ancestral_biallelic(invcf, outvcf):
    vcf = VCF(invcf)
    writer = Writer(outvcf, vcf)

    total = 0
    flipped = 0
    unflipped = 0
    skipped = 0
    invar = 0

    for variant in vcf:
        total += 1
        aa = variant.INFO.get('AA')
        if len(variant.ALT) == 0:
            writer.write_record(variant)
            invar += 1
            continue
        ref, alt = variant.REF, variant.ALT[0] # assumes biallelic
        if len(alt) == 0:
             writer.write_record(variant)
             invar += 1
             continue
        if aa is None or aa not in (ref, alt):
            skipped += 1
            continue
        
        if aa == alt:
            flipped += 1
            variant.REF, variant.ALT = alt, [ref]
            new_gts = [flip_biallelic_genotype(gt) for gt in variant.genotypes]
            variant.set_format('GT', np.array(new_gts, dtype=int)) # an excellent line
        else:
            unflipped += 1

        writer.write_record(variant)
    writer.close()
    print(f"Total sites: {total}\tSkipped: {skipped}\tInvariant: {invar}\tFlipped: {flipped}\tUnflipped: {unflipped}")





# this is a nice time...NOTE: 
def strip_biosuffixes(p, getlist=True):

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



def vcf_to_unzip(vcf, out=None, rem=False):
    vcf, suffs = strip_biosuffixes(vcf)
    if out is None and vcf.with_suffix('.vcf').exists():
        return
    elif out is None:
        subprocess.run(f"bcftools view -Ov -o {vcf}.vcf {vcf}.vcf.gz", shell=True)
    else:
        out = strip_biosuffixes(out, False)
        subprocess.run(f"bcftools view -Ov -o {out}.vcf {vcf}.vcf.gz", shell=True)
    if rem:
        vcf.with_suffix('.vcf.gz').unlink()

# 1mb chunk parallel.... 

def run_parallel_singer(vcf, Ne, m, L=1e7, output=None, n=100, thin=20, polar=0.99, cores=1, execmode=""): # Ne is an issue... 
    if len(execmode) > 0: # if it has a mode... 
        execmode += "-"
    vcf = strip_biosuffixes(vcf, False)
    if not vcf.with_suffix('.vcf').exists():
        if vcf.with_suffix('.vcf.gz').exists():
            vcf_to_unzip(vcf, rem=True)

    if output is None:
        output = vcf
    else:
        output = strip_biosuffixes(output, False)

    subprocess.run(f"{execmode}parallel_singer -Ne {Ne} -m {m} -L {int(L)} -vcf {vcf} -output {output} -n {n} -thin {thin} -polar {polar} -num_cores {cores}", shell=True)  
    # all in one job... 40 cores??? 


def run_basic_singer(vcf, outfix, start, end, Ne=50_000, m=5e-7, n=100, thin=20, polar=0.99, execmode=""):
    if len(execmode) > 0: # if it has a mode... 
        execmode += "-"
    subprocess.run(f"{execmode}singer_master -Ne {Ne} -m {m} -vcf {vcf} -output {outfix} -start {start} -end {end} -n {n} -thin {thin} -polar {polar} ", shell=True)


## output

if __name__ == "__main__":
    wdir = Path(__file__).parent # local directory... 
    in_vcf = wdir / "mj_init_may24_apacest_i10_chr2a_merged.vcf.gz" # 7GB file
    dir1 = wdir / 'nodate' / 'singer'
    reordered_vcf = dir1 / "chr2a_apacest_i10_reorder.vcf.gz"
    midvcf = dir1 / "chr2a_apacest_i10.vcf"
    out = dir1 / "chr2a_apacest_i10_inf.trees" # should work???
    #vcf_reorder_ancestral_biallelic(in_vcf, reordered_vcf)
    #vcf_to_unzip(reordered_vcf, midvcf, rem=True)
    run_parallel_singer(midvcf, 50000, 4.89e-9, output=out, polar=0.99, cores=40) # note
    # 4.89*10-9... (but: all sites includes... so mutation rate needs to be adjusted to reflect missing data... 1/10 sites? we will scalte. 5*10-8)
