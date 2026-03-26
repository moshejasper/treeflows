import numpy as np
# try and use proper vcf files... 
import sys
from cyvcf2 import VCF
import gzip
from pathlib import Path
from treeflows.vcf_core import VCFReader
from treeflows.slurm import SlurmJob
import subprocess
import argparse
import pandas as pd
import random

### attempts to build new program... 

# it sounds like I can only use one position in the file... (?)... perhaps we stick to AA? ... 
# I am going to try and write a temp fasta file that can be used to reconstruct a major fasta file at a later point... default 'N'... 


def vcf_to_like_bgl(infile, outfile): 
    """Converts  (or extracts) genotype likelihoods from gatk-derived vcf and creates beagle 3 (gzipped) likelihood file."""
    vcf = VCFReader(infile)
    with gzip.open(outfile, 'wb') as ofi:
        hbegin = 'marker\talleleA\talleleB\t'
        sampnames = '\t'.join([x for y in zip(vcf.samples, vcf.samples, vcf.samples) for x in y])
        header = hbegin + sampnames + '\n'
        #print(header)
        ofi.write(header.encode())

        counter = 0

        for site in vcf:
            # if counter >= 10_000:
            #     break
            # for now we are going to trust the simple info-based information, though we probably need to update this at some point... 
            if not site.passes_filter(): # should make it biallelic
                continue
            # now we are going to be more nuanced. 
            if len(site.ALT) > 1: # just getting rid of bad alleles... this shoudl be run after recalibration... 
                continue
            # at this point we will just trust things... 
            counter += 1

            # main code
            pref = site.var.gt_phred_ll_homref
            phet = site.var.gt_phred_ll_het
            palt = site.var.gt_phred_ll_homalt                  # (n_samples, 1)

            locid = f'{site.CHROM}_{site.POS}'
            vinf = '\t'.join([locid, site.REF, site.ALT[0]])

            phreds = np.dstack((pref, phet, palt))              # (n_samples, 1, 3)
            likes = phred_to_likelihood(phreds)                 # (n_samples, 1, 3)
            likes = likes.reshape(len(vcf.samples), 3)          # (n_samples, 3)
            likes = likes.round(decimals=4) # limiting step
            joined = '\t'.join('\t'.join(map(str, row)) for row in likes)
            joined = vinf + '\t' + joined + '\n'

            #likes = [phred_to_likelihood(phred) for phred in phreds]
            #print(likes)
            ofi.write(joined.encode())
        

def phred_to_likelihood(val):
    '''
    val: array of phred-scaled likelihoods; last axis = 3 (AA, AB, BB)
        e.g. shape (n_samples, 1, 3) or (n_samples, 3)
    returns: same shape, but converted to probabilities summing to 1 
        along the last axis'''
    val = np.asarray(val, dtype=float)

    # bad if ANY of 3 genotype phreds is nan
    bad = np.isnan(val).any(axis=-1)

    # convert Phred to probs: p = 10^(-Q/10)
    arr1 = np.power(10.0, -val / 10.0)
    # sum across genotype dimensions, keep that dimension for broadcasting
    arr2=arr1.sum(axis=-1, keepdims=True)
    out =  arr1 / arr2

    bad2 = ~np.isfinite(out).all(axis=-1)

    # merge conditions
    bad = bad | bad2
    # uniform replacement vector
    uniform = np.full(val.shape[-1], 1.0/val.shape[-1])
    out[bad] = uniform
    return out

def run_NGSadmix(infile, outfile, K, maf=0.025, miss=0.5, threads=1):
    """Runs NGSadmix. 
    Parameters:
        infile - full path of gzipped likelihoods file
        outfile - path prefix of output files
        K - number of clusters to infer 
        maf - minor allele frequency cutoff
        miss - maximum level of missingness before site is exlcuded
        threads - number of threads to proceed on"""
    result=subprocess.run(f"NGSadmix -likes {infile} -K {K} -o {outfile} -minMaf {maf} -misTol {miss} -P {threads}", shell=True, text=True, stderr=subprocess.PIPE)
    print(result.stderr)
    if result.returncode != 0:
        exit(f"ADMIX failed")
    # allows ancestral population frequencies... admixture proporitions, outfile-prefix... + snp stats
    # seed (reproduce, test convergence)... 
    # missing tolerance (0.05)... default 0.05 (5% presence = pass)
    # has a minMaf filter (5%). Need to work out what this is doing? - min informative idvs. (essentially a MAC reading...)
    # [probably run this at two missingnesses and at mac 3 and maf 0.05]

def run_NGSadmix_multiN(infile, outfile, kmin=1, kmax=10, maf=0.025, miss=0.5, threads=1):
    """ Runs NGSadmix across multiple K"""
    for k in range(kmin, kmax+1):
        run_NGSadmix(infile, str(outfile) + f'_{k}', k, maf, miss, threads)


# future script will scrape logfiles to extract likelihoods, as well as scrape other files (perhaps with pop markers?) and extract things for graphing

# def conversion_2(infile, outfile):
#     """Experimental. Super quick, but currently writes Phred"""
#     command = f"bcftools view -i 'GT=\"alt\"' {infile} "
#     command += "| bcftools query -f '%CHROM:%POS\t%REF\t%ALT{0}[\t%PL{0}\t%PL{1}\t%PL{2}]\n' - "
#     subprocess.run(f"{command} > {outfile}", shell=True)

#### SECTION FOR POSTPROCESSING #####

def adx_get_idxs(adx_file, ancestry, thresh=0.75):
    adx = pd.read_csv(adx_file, sep=" ", header=None).dropna(axis=1) # assume in order
    return [i for i in adx.iloc[adx[ancestry] >=thresh, ancestry]] # gives me indexes... 

def adx_get_idxs_all(adx_file, thresh=0.75):
    df = pd.read_csv(adx_file, sep=" ", header=None).dropna(axis=1)
    return [[x for x in df[df[anc] >= thresh].index] for anc in range(df.shape[1])]

def adx_idx_permute(adxlist):
    # we want to preserve list shapes... 
    grouplens = [len(a) for a in adxlist]
    adx_flat = [a for group in adxlist for a in group]
    random.shuffle(adx_flat) # modifies in place... 
    newlist = [[adx_flat.pop() for _ in range(g)] for g in grouplens]
    return newlist

def adx_idx_downsample(adxlist, max=5):
    return [sorted(random.sample(x, k=max if len(x) > max else len(x))) for x in adxlist]

if __name__ == '__main__':
    myjob = SlurmJob(threads = 4, threadmem=20, jobtime="10:00:00")
    myjob.submit_slurm(jobname="vcf_to_like")
    # infile = Path('mj_init_may24_apacest_chr2a_merged.vcf.gz')
    # likefile = Path('chr2a_apac_m10_likes.gz')
    # outprefix = Path('ngs_test')
    # #conversion_2("ptest_estfilter.vcf.gz", "likes.txt")
    # vcf_to_like_bgl('ptest_estfilter.vcf.gz', 'chr2a_apac_m10_likes.gz')
    # run_NGSadmix_multiN(likefile, outprefix)
    # #test against vcftools query


if __name__ == "__slurm__":
    infile = Path('mj_init_may24_apacmaj_m10_chr2a_merged.vcf.gz').resolve()
    likefile = Path('chr2a_apac_m10_likes.gz').resolve()
    outprefix = Path('chr2a_apac_m10_').resolve()
    #vcf_to_like_bgl(infile, likefile)

    # set up job... 
    my_array = SlurmJob(threads=10, threadmem=15, jobtime="2-00:00:00", runblock="__array__") # in future make it its own class... 
    for k in range(1, 10+1):
        my_array.submit_slurm(jobname = f"NGSadmix-{k}", argstring = f"--infile {likefile} --outprefix {outprefix} --K {k} --maf {0.025} --miss {0.5} --threads {10}")

if __name__ == "__array__":
    # we are going to be an argument-based server
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', type=str)
    parser.add_argument('--outprefix', type=str)
    parser.add_argument('--K', type=int)
    parser.add_argument('--maf', type=float)
    parser.add_argument('--miss', type=float)
    parser.add_argument('--threads', type=int)
    args = parser.parse_args()

    run_NGSadmix(Path(args.infile).resolve(), Path(args.outprefix + str(args.K)).resolve(), args.K, args.maf, args.miss, args.threads)

