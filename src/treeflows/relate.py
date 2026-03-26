
##### Package Imports #####

# general
from pathlib import Path
import shutil
#from cyvcf2 import VCF
#from pyfaidx import Fasta
from treeflows.vcf_core import VCFReader
from treeflows import refdata
import gzip
import subprocess
import os

path2relate = "/data/gpfs/projects/punim1778/programs/relate/relate/"
path2relate_parallel = path2relate + "scripts/RelateParallel/RelateParallel.sh"
path2relate_popsize = path2relate + "scripts/EstimatePopulationSize/EstimatePopulationSize.sh"
path2relate_tskit = path2relate + "bin/RelateFileFormats"
path2relate_mutrate = path2relate + "bin/RelateMutationRate"
path2relate_extract = path2relate + "bin/RelateExtract"


# custom

##### Implementation Notes #####

"""
RELATE uses the haps/sample file format (which is the output of SHAPEIT2) as input
Data must be phased (& separated by chromosome)
Relate provides function to convert to haps/sample:
    :   "RelateFileFormats --mode ConvertFromVCF --haps example.haps --sample example.sample -i example [vcf without extension
    : there is also a more extensive shell script, 'PrepareInputFIles.sh' which handles possibly important info... 

The order of events is: 
*    1. Convert to haps/sample 
*    2. Remove non-biallelic SNPs
        a. custom script exists, but I could always do this prior
*    3. Determine ancestral allele: flip if necessary
*        a. needs ancestral genome as FASTA file [JUST DID IT MANUALLY]
            - generate this using custom pyfaidx-based function (one for each)
        b. all of these are to some extent replicable... but I shouldn't re-invent the wheel
*    4. Remove samples from dataset
        a. I will do in advance, but could possibly be useful here, depending (might save space)
*    5. Filter SNPs using genomic mask [via .dist file - generated in custom manner]
        a. FASTA file - N = get rid, P = pass. 
            - deletes SNPs that are non-passing
*            - outputs .dist file -> [BUILT MANUALLY]
                > contains distances between SNPs
                > adjusted for regions that are not passing in the mask
                * i.e. distances are shortened where samples have been deleted
        b. A single mask is used
            - needs to be adjusted for a variety of use cases... all combined into one master fasta file
            - main benefit is the .dist file (but this can potentially be produced via another method? )
    6. Generate SNP annotations (mostly pop-related)

Relate has scripts for parallelizing at various scales (as does TINFER) - they seem quite complex. 
PARRALELIZATION RROCESS

1. Divide Chromosomes into Chunks 
[
    2. Paint -> distance matrices as temporary files
    3. Build Topology -> Estimate Tree topologies (tskit equiv)
    4. Find Equivalent branches (in parallel trees) (tskit equiv)
    5. Infer Branch Lenghts (tsdate equiv)
    6. Combine Sections (concatanate files within chunk)
]
7. Finalize (combine chunks, generate output)

Note that in this process, relate uses 20,000 bp overlaps along the chunks, then (owing to its internal consistency) erases the end
 section in one with the start section in the next. This cannot be directly replicated with 'parallel chunks' without forfeiting what 
 Relate is doing, so I might need to consider this. Basically, my chunks can't be too small... if I am running a chunk-based analysis, then I need to 
 only be interested in aggregates across the genome, and I need to interpret them correctly. Jury is still out. Also, how to combine tskit chunks more 
 effectively

]


"""


##### Helper Functions #####

# Generate ancestral genome FASTA... (?)... 
    # Need: 
    # (1) template genome
    # (2) vcf with AA tag... 

def make_ancestral_fasta(invcf, infasta):
    """Currently not in use"""
    #infa = Fasta(infasta)
    vcf = VCFReader(invcf) # get the struture right


def make_hap_and_samps_file(vcf, out): 
    """Single-threaded workhorse script
    Takes input vcf (with AA ancestral allele data & all relevant (nonvariant) sites present). Produces
     .haps, .sample, & .dist files as output, ready for Relate inference. 
    Filtering of individuals is currently assumed to have happened elsewhere in .vcf land. 

    Parameters:
        vcf     -   (Path/string)   -   full path to .VCF input file
        out     -   (Path/string)   -   path to relate file prefix. (endings added later)
    
    Requirements:
        fasta2sfs (+ pyfaidx, cyvcf2)
        numpy > 2 (args conda env has a numpy 1.x version dependency) -> needed for fasta2sfs

    Details:
        1. Writes all individuals in VCF file to .samples file format (with missing set to '0')
        2. Iterate through vcf outputting biallelic sites with known ancestry to .haps, and
            counting known (invariant) allele states between them. [probably adjust formula for unknown biallelic]. 
            Write these latter distances to .dist file
    """
    vcf = VCFReader(vcf)
    out = Path(out)

    with open(out.with_suffix('.sample'), 'w') as samp:
        samp.write("ID_1 ID_2 missing\n0 0 0\n")
        for sample in vcf.samples:
            samp.write(f"{sample} {sample} 0\n")
    
    with open(out.with_suffix('.haps'), 'w') as hapfile, open(out.with_suffix('.dist'), 'w') as distfile:
        ambigs = 0
        prevline = "#pos dist\n"
        distfile.write(prevline)
        prev_pos = -1 # default... 
        poscount = 0
        # CHROM NUM :: SNP ID :: POS :: ANC ALLELE :: ALT ALLELE :: (space delineated)
        for site in vcf:
            poscount += 1
            if not site.is_biallelic:
                continue

            ancestral = site.INFO['AA']
            if ancestral == 'N':
                ambigs += 1
                continue
            derived = site.geno_alleles[1 - site.geno_alleles.index(ancestral)] # super reliant on biallelic!!!

            # we need to flip alleles... 
            genos = [h for g in site.genotypes for h in g[:2]] # should do... 
            if ancestral != site.REF: #i.e. a flip event... 
                genos = [1 if g == 0 else 0 if g == 1 else -1 for g in genos] # flip... full form

            hapfile.write(f"0 . {site.POS} {ancestral} {derived} {" ".join(str(g) for g in genos)}\n")
            # distfile... 
            if prev_pos == -1:
                prev_pos = site.POS
            else:
                distfile.write(f"{prev_pos} {poscount}\n") # only using sites I have access to... 
                prev_pos = site.POS
            poscount = 0 # will only be activated if biallelic site present, meaning it should always be reset after writing. 
        distfile.write(f"{prev_pos} 1\n") # final line of file... 

        
def make_map_relate(inmap, mapfile, chrval='X2'):
    """Takes *very specific* file format (emailed to me from study) strip out chromosome 
     corresponding to chrval, and return Relate file format with recombination rate calculated. 
    
     Parameters: 
        inmap      -   (Path/string)   -   Path to original map file. which has space-separated format: 
                    - set   (not important)
                    - map   (chromosome column)
                    - mkr   (marker - not used)
                    - phys  (physical genome position)
                    - gen   (genetic map position) [format???]
        mapfile     -   (Path/string)   -   Path to output map file, with space-separated format:
                    - pos   (physical genome position)
                    - COMBINED_rate (recombination rate between POS & POS_next)
                    - Genetic_Map (genetic map position)
                    , where rate is per Mb (?). 
        chrval      -   (string)        -   selection value to match 'map' chromosomal name to be extracted from inmap
    
    Details:
        Transforms genetic map (type1) into Relate-compatible, including calculating recombination rates. 
         For recombination rates, calculated between current position [i] and next position [i+1] using the formula:
            RATE = [(PHYS_POS[i+1] - PHYS_POS[i]) / (GEN_POS[i+1] - GEN_POS[i])] * 1,000,000
        The beginning and end of the chromosome are set to have the same rates as the first and last calculated 
         regions respectively. 
     """
    genprev = 0
    physprev = 0
    rate=0
    with open(inmap) as imap, open(mapfile, 'w') as omap: 
        omap.write(f"pos COMBINED_rate Genetic_Map\n")
        for x in imap:
            try:
                set, map, mkr, phys, gen = x.strip().split()
            except ValueError:
                break
            if map != chrval:
                continue
            rate = (float(gen) - genprev) / (int(phys) - physprev) * 1e6
            omap.write(f"{physprev} {rate} {genprev}\n")
            physprev, genprev = int(phys), float(gen)

    pass

def make_map_relate_simple(mapfile, rate):
    """creates simple constant-rate map, where input rate is multiplied by 1,000,000 to give Relate-compatibel rate"""
    with open(mapfile, 'w') as omap:
        omap.write(f"pos COMBINED_rate Genetic_Map\n")
        omap.write(f"0 {rate * 1e6} 0\n")
        

##### Main Functions #####

def run_relate_parallel(infix, outfix, m, N, map=None, mem=10, threads=1):
    """Runs parallel relate with specific input parameters (should be run on multithread job)"""

    subprocess.run(f"bash {path2relate_parallel} -m {m} -N {N} --haps {infix}.haps --sample {infix}.sample --dist {infix}.dist -o {outfix} --map {map} --memory {mem} --threads {threads}", shell=True)


def relate_est_pop_size(in_prefix, out_prefix, poplabels, mut_rate, gen_years, epochs="1.5,5,0.1", thresh=0.5, threads=1):
    """Note: epochs here are converted to powers of 10, and divided by years_per_gen. a 0 is appended at start. gen_years is currently 1/15, or 0.067. 
    Gen_Years will need to be used a second time in the graphing. There, popsize must be scaled (/2), and time must be scaled (/gen_years) to make year.
    Thresh is set to mark default. Can potentially increase it to e.g 0.9"""
    # try running without 'pop of interest'... if this fails, need more complicated script (like msmc - adapt from that)
    subprocess.run(f"bash {path2relate_popsize} -i {in_prefix} -o {out_prefix} --poplabels {poplabels} -m {mut_rate} --years_per_gen {gen_years} --bins {epochs} --threshold {thresh} --threads {threads}", shell=True)

def make_relate_popfile(fileprefix, ref_file):
    """ Generate correctly formated Relate poplabels file from reference metadata file and Relate .sample file."""
    ref = refdata.load_aegdata(ref_file)
    fileprefix = Path(fileprefix)
    with open(fileprefix.with_suffix(".sample"), 'r') as samp, open(fileprefix.with_suffix(".poplabels"), 'w') as poplab:
        next(samp) # skip header
        next(samp) # skip header2
        poplab.write("ID POP GROUP SEX\n")
        for line in samp:
            idv = line.strip().split()[0]
            pop = ref.get_pop_short_idv(idv)
            poplab.write(f"{idv} {pop} {pop} NA\n")

def rel_annotate(fileprefix, poplabels):
    """Generate Relate SNP annotations for a run prefix.

    Args:
        fileprefix: Prefix for `.haps/.sample/.mut` files.
        poplabels: Poplabels file path passed to Relate.

    Returns:
        None. Writes Relate annotation outputs to `fileprefix.*`.
    """
    subprocess.run(f"{path2relate_tskit} --mode GenerateSNPAnnotations --haps {fileprefix}.haps --sample {fileprefix}.sample --mut {fileprefix}.mut --poplabels {poplabels} -o {fileprefix}", shell=True)
    

def relate_to_treeseq(relfix, treefix=None):
    """Converts to tree sequence using internal relate script"""
    treefix = relfix if treefix is None else treefix

    subprocess.run(f"{path2relate_tskit} --mode ConvertToTreeSequence -i {relfix} -o {treefix}", shell=True)

def relate_est_mut_rate(in_prefix, out_prefix=None, gentime = 1/15, dist=None):
    """Estimate average mutation rate using RelateMutationRate.

    Args:
        in_prefix: Input Relate prefix (expects `.anc/.mut` etc.).
        out_prefix: Output prefix (defaults to `in_prefix`).
        gentime: Years per generation.
        dist: Distance file prefix/path (defaults to `in_prefix`).

    Returns:
        None.
    """

    if out_prefix is None: out_prefix = in_prefix
    if dist is None:
        dist = in_prefix # + .dist? 
    subprocess.run(f"{path2relate_mutrate} --mode Avg -i {in_prefix} -o {out_prefix} --years_per_gen {gentime} --dist {dist}", shell=True)

def relate_extract_subpops(in_prefix, poplabels, pops: list, out_prefix): # handles things on disk rather than tree sequences... (?)... -> unclear results... 
    """Extract subpopulation subtrees from a Relate run.

    Args:
        in_prefix: Input Relate prefix.
        poplabels: Poplabels file path.
        pops: List of population labels to include.
        out_prefix: Output prefix.

    Returns:
        None.
    """
    anc = str(in_prefix) + ".anc"
    mut = str(in_prefix) + ".mut"
    subpops = ",".join(pops)
    subprocess.run(f"{path2relate_extract} --mode SubTreesForSubpopulation --anc {anc} --mut {mut} --poplabels {poplabels} --pop_of_interest {subpops} -o {out_prefix}", shell=True)

def relate_lowmut_filter(in_prefix, out_prefix, thresh=0.5):
    """Removes trees & associated SNPs from .anc/.mut files onto which less than a given number of SNPs are mapped
        Note that this will (possibly?) tend to increase remaining tree lengths (based on edges)"""
    anc = str(in_prefix) + ".anc"
    mut = str(in_prefix) + ".mut"
    subprocess.run(f"{path2relate_extract} --mode RemoveTreesWithFewMutations --anc {anc} --mut {mut} --threshold {thresh} -o {out_prefix}", shell=True)


##### Operating system #####

if __name__ == '__main__':
    cwd = Path(__file__).parent.resolve()
    invcf = cwd / "mj_init_may24_apacest_i10_chr2a_merged.vcf.gz"
    outdir = cwd / "nodate" / "relate"
    relfile = outdir / "chr2a_apacest_i10" # file endings left off... 
    inmap = cwd / "Linkage_Phys_map.txt"
    outmap = cwd / "relate_map_chr2.txt"
    outmap_simp = outdir / "relate_map_simple.txt"
    make_hap_and_samps_file(invcf, relfile)
    make_map_relate(inmap, outmap)
    os.chdir(outdir)
    outfile = Path("chr2a_apacest_i10")
    print(outfile)
    run_relate_parallel(relfile, outfile, 4.89e-9, 50000, map=outmap_simp, mem=20, threads=40)
