#!/bin/python3

import subprocess
from pathlib import Path
from cyvcf2 import VCF, Writer
from pyfaidx import Fasta
from .vcf_core import VCFReader, VCFWriter, VCFCompare, VarParsed
from .vcf2est import vcf_polish_all
from .Spath import Mpath
# import vcf2est as est
# import vcf2tsinfer as v2t
import numpy as np
from math import log
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

logger = logging.getLogger(__name__)
logger.setLevel('INFO')


# strategy will be to match chromosome, etc. (in fasta) via having a '1 fasta to many vcfs' ratio and using the vcf to work out what the fasta is. 
# we will generate sfs input files from individual vcfs, run the inference, then write new vcfs and fastas downstream (gather many fastas to make 'ancestral' file)
# in future this will need some kind of pileup from the formosus genome (and its ancestral instincts?) also. still to work this out... 


gen_intvs = ("NC_035108.1",)# "NC_035109.1")#"NC_035107.1", "NC_035108.1",)# "NC_035109.1") #, "NC_035159.1") #NC_035107.1
gen_lengths = (474425716,)# 409777670)#310827022,  474425716,)# 409777670) #, 16790) #310827022
int_sublength = 10000000 # was historically using 10,000,000 (this will match the various files)


def vcf_num_sites(invcf):
    vcf = VCF(invcf + '.vcf.gz')
    vcount = 0
    for var in vcf:
        vcount += 1
    print(f"Total sites: {vcount}")

# going to write a new file... (a vcf aligner and header restorer)... (hmm)... will be tricky... 

def vcf_format_test(invcf, outvcf):
    vcf = VCF(invcf)
    print(vcf.raw_header)

def vcf_gather_preinfer(outvcf, chromname="NC_035108.1", chromlength=474425716, sublength=10000000, workdir = '/data/gpfs/projects/punim1778/Projects/aegypti/2023/data/previous_studies/core_files/',
                        fprefix='mj_init_may24_', fpostfix='_raw_aegpac_bgla'):
    
    tractlist = []
    snpcount = 0
    nsplits = chromlength //sublength + 1
    for n in range(nsplits):
        start = n*sublength + 1
        end = (n + 1) * int_sublength
        if end > chromlength:
            end = chromlength
        gentract = f"{chromname}:{start}-{end}"
        tractlist.append(gentract)

    # now to run the main file... 
    #1. generate template vcf... # going to make a hacky thing... 
    # tvcf = VCF(tempvcf)
    ftemplate0 = f"{workdir}{fprefix}{tractlist[0]}{fpostfix}.vcf.gz"
    vcf0 = VCF(ftemplate0)
    vwriter = Writer(outvcf, vcf0, 'wz')
    vwriter.write_header()
    # with gzip.open(outvcf, 'w') as ov:
    #     ov.write(tvcf.raw_header.encode())
    # hwriter = Writer(outvcf, tvcf, 'wz')
    # hwriter.write_header()
    # hwriter.close
    # hopefully a clean file with the vwriter still going strong? (unless we hit a race condition)

    for tract in tractlist:
        print(tract)
        infile = f"{workdir}{fprefix}{tract}{fpostfix}.vcf.gz"
        vcf = VCF(infile) # load in new vcf in sequence... this is primarily a gather operation... we are just summing up vcf locs under
        for var in vcf:
            vwriter.write_record(var)
            snpcount += 1
    vwriter.close()
    print(f"Total SNPS in final sequence: {snpcount}")


def vcf_post_beagle_rehead(invcf, outvcf, tempvcf): # write header from tempvcf
    # we will probably try and rework them... 
    vcf = VCF(invcf)
    vtemp = VCF(tempvcf)
    vwrite = Writer(outvcf, vtemp, 'wz')
    vwrite.write_header()
    # try and fix the headers... 
    for hline in vtemp.raw_header.split(sep='\n'):
        if hline.startswith('##contig=<ID=NC'):
            vcf.add_to_header(str(hline).strip())
    vcf.add_info_to_header({'ID': 'AA', 'Description': 'inferred ancestral allele (est-sfs)', 'Type': 'Character', 'Number': '1'})
    for n, (var, ancvar) in enumerate(zip(vcf, vtemp)): 

        if var.POS != ancvar.POS:
            exit(f'Positions not compatible: :{var.POS} & {ancvar.POS}')

        ancval = ancvar.INFO['AA']
        var.INFO['AA'] = str(ancval)
        
        vwrite.write_record(vwrite.variant_from_string(str(var).strip())) # a terrible workaround to stop some kind of terrible issue!!!


def est_postprocess_group(invcf, outvcf, tempfolder, filter=True, pthresh=0.05): # I will need to write some sort of config file... 
    '''Postprocesses a jointly prepaired .vcf .est and .anc files after running EST-SFS. 
This function assumes that the .vcf and .est files were processed with the preprocess function 
 and exactly match. It takes the .anc output files (with ancestral probabilities) and recombines them
 according to the template order of the .vcf file. Currently there are no inbuilt checks to ensure these 
 exactly match, so only files produced in this specific way should be used.
    Arguments:
        invcf       -   VCF file (produced via est_preprocess) to use as template for ancestry assignment
        tempfolder  -   ancestry folder (produced from est_preprocess) where config and EST-SFS in/output files are located
        outvcf      -   VCF file name to write ancestry-adjusted VCF to
        filter      -   currently not implemented
        pthresh     -   Probability threshold for ancestry filtering. Sets upper (& lower) bounds for ancestral allele assignment
        '''
    

    ## load config files
    tempfolder = Path(tempfolder)
    fbins = np.load(tempfolder / '.fbins.npy')
    fnames = np.load(tempfolder / '.fnames.npy', allow_pickle=True)
    fidx = np.load(tempfolder / '.fidx.npy')
    idx = np.load(tempfolder / '.idx.npy')
    bins = np.load(tempfolder / '.bins.npy')
    bincount = np.load(tempfolder / '.bincount.npy')
    binwidth, cmin, bincountmin = np.load(tempfolder / '.pri.npy')
    monoallelic, multiallelic, keep_nonsnp, keep_ambigs, keep_singles = np.load(tempfolder / '.prb.npy')
    # print(fnames)
    # print(bins)
    # print(binwidth)

    ## set up files
    vcf = VCF(Path(invcf))
    vcf.add_info_to_header({'ID': 'AA', 'Description': 'inferred ancestral allele (est-sfs)', 'Type': 'Character', 'Number': '1'})
    vcf.add_info_to_header({'ID': 'AAP', 'Description': 'Pr(ancestral allele known)', 'Type': 'Float', 'Number': '1'})
    vcf_out = Writer(Path(outvcf), vcf, 'wz')
    vcf_out.write_header()

    ancreaders = [open(fname.with_suffix('.anc')) for fname in fnames]
    estreaders = [open(fname) for fname in fnames]



    sitecount = 0
    flipcount = 0 

    for site in vcf:
        sitecount += 1
        # if sitecount > thresh:
        #     break
        # assume filtering exists... 
        v = VarParsed(site)

        alleles = v.geno_alleles # the proper kind

        an = v.an
        if an < cmin: exit('Data is invalid!')
        aix = an // binwidth # index value... 
        bincount[idx[aix]] += 1 # here is the key filter step. we use the index to sort everything out. 
        fid = fidx[aix] # file index

        ancline = ancreaders[fid].readline()
        while ancline.startswith('0'):
            ancline = ancreaders[fid].readline()
        estline = estreaders[fid].readline()
        # now can skip things, etc. 
        # site.INFO['AN'] = v.an
        # site.INFO['DP'] = str(v.dp).encode()

        ancline_parsed = ancline.strip().split()
        pval = float(ancline_parsed[2]) # pr(reference allele is ancestral (assume two alleles))

        if v.is_monoallelic:
            site.INFO['AAP'] = pval
            site.INFO['AA'] = v.geno_maj
            vcf_out.write_record(site)
            continue

        #est_parsed = estline.strip().split(sep='\t')
        if pval > 1 - pthresh:
            # majority allele is go. assumes biallelic
            site.INFO['AAP'] = pval
            site.INFO['AA'] = v.geno_maj
            

        elif pval < pthresh: # flipped - we won't worry about minor alleles just yet
            if v.is_biallelic:
                site.INFO['AAP'] = 1 - pval
                site.INFO['AA'] = v.geno_2nd
            else:
                #continue # will fix in future...
                site.INFO['AAP'] = 1 - pval
                site.INFO['AA'] = v.geno_2nd # assuming things get filtered out later... (e.g. by to_biallelic singleton removers... )

        else: # these are the hard cases... we'll rely on previous files 
            if not keep_ambigs:
                continue
            site.INFO['AAP'] = pval
            site.INFO['AA'] = 'N' # parental mismatches still to be dealt with... 

        # print(site.POS)
        # print(est_parsed)
        # print(alleles)
        # print(pval)
        # print('\n')

        vcf_out.write_record(site)
    
    for reader in ancreaders:
        reader.close()
    for reader in estreaders:
        reader.close()







def est_vcf_postprocess(invcf, outvcf, inanc, filter=True, pthresh=0.05):

    monoalleles = 0
    multialleles = 0
    minority = 0
    flipped = 0
    unflipped = 0
    ambiguous_sites = 0
    allsites = 0

    
    vcf = VCF(Path(invcf))
    alleledict = {"A": 0, "C": 1, "G": 2, "T": 3}
    alleledict_rev = {0: "A", 1: "C", 2: "G", 3: "T"}

    vcf.add_info_to_header({'ID': 'AA', 'Description': 'inferred ancestral allele (est-sfs)', 'Type': 'Character', 'Number': '1'})

    vcf_out = Writer(Path(outvcf), vcf, "wz") # compressed vcf
    vcf_out.write_header()
    with open(inanc, 'r') as ianc:
        ancline = next(ianc)
        while ancline.startswith('0'):
            ancline = next(ianc) # this should walk to the start of the sequence... 

        for n, var in enumerate(vcf): # iterate ancline at ends as it is already active
            # currently dealing with snps and monoallelic sites, where some snps are not in the sequence... 
            allsites += 1

            if len(var.ALT) == 0: # (get rid of monoallelic sites)
                ancline = next(ianc)
                monoalleles += 1
                continue

            if 2 * var.num_hom_alt + var.num_het < 3: # mac 3
                ancline = next(ianc)
                monoalleles += 1 # a little dubious a definition here... this affects our ability to run e.g. heterozygosity analyses. 
                continue

            if len(var.ALT) > 1: # (get rid of hypothetical multiallelic sites)
                ancline = next(ianc)
                multialleles += 1
                continue

            # need to work out whether referece is major allele!!!
            refcount = 2 * var.num_hom_ref + var.num_het
            altcount = 2 * var.num_hom_alt + var.num_het
            ref_is_maj = refcount > altcount # NECESSARY AS EST-SFS IS COMPLETELY UNAWARE AS TO WHAT THE MAJOR ALLELE IS, BUT IT IS DERIVED FROM THIS VCF ORIGINALLY!!!

            ancline_parsed = ancline.strip().split() # hopefully this will be enough... # need to makes this about the ... 
            pval = float(ancline_parsed[2]) # pr(reference allele is ancestral (assume two alleles))
            
            # upper filter
            if pval > 1 - pthresh: #unflipped allele
                if ref_is_maj:
                    ancbp = var.REF
                    unflipped += 1
                else:
                    ancbp = var.ALT[0]
                    flipped += 1
            elif pval < pthresh: # flipped allele
                if ref_is_maj:
                    ancbp = var.ALT[0] # probably ok? need to test... 
                    flipped += 1
                    minority += 1
                else:
                    ancbp = var.REF
                    unflipped += 1
                    minority += 1
            else: # ambiguous
                ambiguous_sites += 1
                ancline = next(ianc)
                continue
            
            ### CODE USING ANCESTRY OF OUTROUPS (IGNORE FOR NOW)
            # # 20/80 principle... #probably encode this
            # if float(ancline_parsed[2]) < 0.5: # i.e. if the main site needs to be flipped. # we are going to take a window here... 
            #     # new window approach::: confidence of 0.9 to 0.1... 
            #     ancvals = [float(t) for t in ancline_parsed[3:]] # these should typically be float values... 
            #     ancidx = int(max(range(len(ancvals)), key=ancvals.__getitem__)) # gives the index (0-3)... 
            #     ancbp = alleledict_rev[ancidx] # should give the value... 
            #     if not ancbp in alleles: # dismiss... 
            #         ancline = next(ianc)
            #         orthogonals += 1
            #         continue
            #     flipped += 1
            # else:
            #     ancbp = var.REF
            #     unflipped += 1
            
            # have now identified ancestral site... ready to code... 
            var.INFO['AA'] = str(ancbp)
            vcf_out.write_record(var) # does this all work? 
        
        # should probably check file integrity here
    vcf_out.close()
    retained = flipped + unflipped + ambiguous_sites
    print(f"Total sites: {allsites}. Monoallelic: {monoalleles}. Multiallellic: {multialleles}. Retained: {retained}")
    print(f"Biallelic: {retained} ({retained/allsites * 100 :.2f}). -  Flipped: {flipped} ({flipped / retained * 100 :.2f}). Unflipped: {unflipped} ({unflipped / retained * 100 :.2f}). Ambiguous: {ambiguous_sites} ({ambiguous_sites / retained * 100 :.2f}). MAJ_overrides: {minority}")

# TO DO: WRITE A RULE FLIP WHEN HOME... (general flip)... (etc.)... (ok)... 


def assign_ancestrals_majority(invcf, outvcf, multiallelic=False, monoallelic=False, keep_nonsnp=False, keep_ambigs=False, keep_singles=False, repair=False):
    '''Assign ancestral alleles (INFO:AA) to be the most common (present) allele in the current data'''

    if monoallelic: 
        keep_singles = True

    sitecount = 0
    misscount = 0
    multisitecount = 0
    monositecount = 0
    orthogcount = 0
    ambiguous_count = 0
    other_count = 0
    nonsnp_variants = 0
    singletons = 0
    repaircount = 0
    good=0

    vcf = VCF(Path(invcf))
    vcf.add_info_to_header({'ID': 'AA', 'Description': 'inferred ancestral allele (est-sfs)', 'Type': 'Character', 'Number': '1'})
    vcf.add_info_to_header({'ID': 'RRP', 'Description': 'flags when ref needs repair', 'Type': 'Character', 'Number': '1'})
    vcf_out = Writer(Path(outvcf), vcf, "wz") # compressed vcf
    vcf_out.write_header()

    flipcount = 0

    for var in vcf:
        # if sitecount > thresh:
        #     break
        refrepair = False
        #altid = 0
        sitecount += 1
        v = VarParsed(var)

        ### test bit
        if v.is_nonsnp:
            nonsnp_variants += 1
            if not keep_nonsnp:
                continue
        #alleles = v.vcf_alleles
        if not v.has_data:
            misscount += 1
            continue

        #v_a, v_c = v.get_genocounts() # better structure, but probably not good enough. Perhaps I need a shadow variant? 

        # allele count work... 
        if v.is_nonallelic:
            misscount += 1
            continue
        
        if v.is_monoallelic:
            monositecount += 1
            if monoallelic:
                if not v.has_ref:
                    if repair:
                        var.INFO['RRP'] = v.get_geno_ref()
                var.INFO['AA'] = v.get_geno_ref()
                vcf_out.write_record(var)
                good += 1
            continue

        if v.is_multiallelic:
            multisitecount += 1
            if not multiallelic:
                continue

        if not v.has_ref:
            other_count += 1
            if not repair:
                continue
            refrepair = True

        #altid = v.get_geno_alt_idx() # changed definition so this slices from 1: but is indexed from 0:

        #refcount = v.v_c[0]
        #altcount = v.v_c[altid]
        #alleles = v.geno_alleles

        #majorfreq = max(refcount, altcount) / (refcount + altcount) # should be equal and opposite... 

        # possibly future ambiguous site threshold... 

        blanks = ['*', 'N']

        #######################

        if not keep_singles and sum(v.v_c) - max(v.v_c) < 3: # mac 3...
            monositecount += 1
            continue 

        majallele = v.geno_maj
        if majallele in blanks:
            if not keep_ambigs:
                continue
            var.INFO['AA'] = 'N'
        else: var.INFO['AA'] = majallele

        if not majallele == var.REF:
            flipcount += 1
        
        if refrepair and var.INFO['AA'] != 'N': # set to ancestral if ref is ambiguous
            repaircount += 1
            var.INFO['RRP'] = var.INFO['AA']
        vcf_out.write_record(var)
        good += 1
    vcf_out.close()

    with open(Path(outvcf).with_suffix('.stats'), 'w') as ostat:
        rem1 = sitecount - misscount - monositecount - multisitecount - orthogcount - singletons - ambiguous_count
        if keep_singles: rem1 += singletons
        if keep_ambigs: rem1 += ambiguous_count
        ostat.write(f"Infile: {invcf}. Outfile: {outvcf}\n")
        ostat.write(f"Total sites: {sitecount}. Missing ancestral information: {misscount}. Monosites: {monositecount}. Singletons (& dbs): {singletons} Multisites: {multisitecount}. Orthogs: {orthogcount}. Candidates: {rem1}\n")
        maj1 = rem1 - other_count - flipcount
        if keep_ambigs: maj1 -= ambiguous_count
        ostat.write(f"Candidates: {rem1}. Flipped: {flipcount}. Ambiguous: {ambiguous_count}. Other: {other_count}. Majority: {maj1}. RefRepair: {repaircount}\n")
        kept_alleles = flipcount + maj1
        if kept_alleles == 0: kept_alleles = 1
        ostat.write(f"Final alleles: {kept_alleles}. Reference: {maj1} ({maj1 / kept_alleles * 100 :.2f}%). Flipped: {flipcount} ({flipcount / kept_alleles * 100 :.2f}%)\n")


# to do: write tool to impute this into a .fasta rather than alternative means... 


def estimate_ancestrals_basic(invcf, outvcf, fasta1, fasta2, multiallelic=False, monoallelic=False, keep_nonsnp=False, keep_ambigs=False, keep_singles=False, repair=False):
    '''Assign ancestral alles (INFO:AA) to a VCF in  a rule-based manner using two fasta outgroup sequences. 
    The first outgroup (fasta1) is assumed to be more closely related than the second (fasta2). The process is in order and as follows:
        1. If both outgroups match an allele, assign it to be ancestral
        2. If one allele exceeds 95% threshold within spp., assign it to be ancestral
        3. If only one outgroup contains a relevant allele, assign it to be ancestral
        4. If both outgroups contradict, assign ancestral state to be unkown
    
    This was initially developed for biallelic SNPs. Extension to monoallelic sites is trivial; mileage may vary with other options.
    It is currently being modified to enable writing of ancestral fasta file along with (or instead of ) vcf modification'''

    # key change: assign alleles to unknown rather than skipping them (not for other site kinds)

    if monoallelic:
        keep_singles = True

    sitecount = 0
    misscount = 0
    multisitecount = 0
    monositecount = 0
    orthogcount = 0
    ambiguous_count = 0
    other_count = 0
    nonsnp_variants = 0
    singletons = 0
    repaircount = 0

    vcf = VCF(Path(invcf))
    fasta1 = Fasta(Path(fasta1), sequence_always_upper=True) # can't remember what this last bit means
    fasta2 = Fasta(Path(fasta2), sequence_always_upper=True)
    vcf.add_info_to_header({'ID': 'AA', 'Description': 'inferred ancestral allele (est-sfs)', 'Type': 'Character', 'Number': '1'})
    vcf.add_info_to_header({'ID': 'RRP', 'Description': 'flags when ref needs repair', 'Type': 'Character', 'Number': '1'})

    vcf_out = Writer(Path(outvcf), vcf, "wz") # compressed vcf
    vcf_out.write_header()
    flipcount = 0

    # now going to open basic files... (1)

    for var in vcf:
        refrepair = False # we will always repair to the majority allele
        altid = 0
        sitecount += 1
        fastavar1 = fasta1[var.CHROM][var.POS-1]
        fastavar2 = fasta2[var.CHROM][var.POS-1]

        v = VarParsed(var)

        if v.is_nonsnp:
            nonsnp_variants += 1
            if not keep_nonsnp:
                continue

        alleles = v.vcf_alleles
        if not v.has_data:
            misscount += 1
            continue

        # allele count work... 
        if v.is_nonallelic:
            misscount += 1
            continue
        
        if v.is_monoallelic:
            monositecount += 1
            if monoallelic:
                if not v.has_ref:
                    if repair:
                        var.INFO['RRP'] = v.get_geno_ref()
                var.INFO['AA'] = v.get_geno_ref()
                vcf_out.write_record(var)
            continue

        if v.is_multiallelic:
            multisitecount += 1
            if not multiallelic:
                continue

        if not v.has_ref:
            other_count += 1
            if not repair:
                continue
            refrepair = True

        altid = v.get_geno_alt_idx() # changed definition so this slices from 1: but is indexed from 0:

        refcount = v.v_c[0]
        altcount = v.v_c[altid]
        alleles = v.geno_alleles

        majorfreq = max(refcount, altcount) / (refcount + altcount) # should be equal and opposite... 

        blanks = ['*', 'N']

        if fastavar1 == 'N' and fastavar2 == 'N' and majorfreq <= 0.95:
            if keep_ambigs:
                ambiguous_count += 1
                var.INFO['AA'] = 'N'
                vcf_out.write_record(var)
            else: misscount += 1
            continue

        if majorfreq > 0.98: # 
            singletons += 1
            if keep_singles:
                majallele = v.geno_maj
                if majallele in blanks:
                    if not keep_ambigs:
                        continue
                    var.INFO['AA'] = 'N'
                else: var.INFO['AA'] = majallele
        
                if not majallele == var.REF:
                    flipcount += 1
                if refrepair and not var.INFO['AA'] == 'N':
                    var.INFO['RRP'] = var.INFO['AA']
                vcf_out.write_record(var)
            continue
                
        # possibly future ambiguous site threshold... 

        if not fastavar1 in alleles and not fastavar2 in alleles and majorfreq <= 0.95: 
            orthogcount += 1
            if keep_ambigs:
                var.INFO['AA'] = 'N'
                vcf_out.write_record(var)
            continue


        # now apply simple rules... 
        #1. If unanimity of other sequences, force them to be ancestral. 
        if fastavar1 in alleles and fastavar2 in alleles and fastavar1 == fastavar2: # make sure these are appropriate alelles
            var.INFO['AA'] = str(fastavar1)
            if str(fastavar1) != var.REF:
                flipcount += 1
        elif majorfreq > 0.95: # we are going to rectify based entirely on frequency at this threshold... 
            majallele = v.geno_maj
            if majallele in blanks:
                if not keep_ambigs:
                    continue
                var.INFO['AA'] = 'N'
            else: var.INFO['AA'] = majallele
    
            if not majallele == var.REF:
                flipcount += 1
            if refrepair and not var.INFO['AA'] == 'N':
                var.INFO['RRP'] = var.INFO['AA']
            vcf_out.write_record(var)
            continue

        elif fastavar1 in alleles and fastavar2 not in alleles: # these are getting a little fine-grained
            var.INFO['AA'] = str(fastavar1)
            if str(fastavar1) != var.REF:
                flipcount += 1
        elif fastavar2 in alleles and fastavar1 not in alleles: # these are getting a little fine-grained
            var.INFO['AA'] = str(fastavar2)
            if str(fastavar2) != var.REF:
                flipcount += 1
        elif fastavar1 in alleles and fastavar2 in alleles and fastavar1 != fastavar2:
            # these are going to be hard to decide... 
            ambiguous_count += 1
            if keep_ambigs:
                var.INFO['AA'] = 'N'
                vcf_out.write_record(var)
            continue
        else:
            other_count += 1
            print(f"Ref: {var.REF}. Alt: {var.ALT}. Fasta1: {fastavar1}. Fasta2: {fastavar2}. refcount: {refcount}. altcount: {altcount}. MajAF: {majorfreq}")
            continue

        # now we can write the variable... 
        if refrepair and var.INFO['AA'] != 'N': # set to ancestral if ref is ambiguous
            repaircount += 1
            var.INFO['RRP'] = var.INFO['AA']
        vcf_out.write_record(var)
    vcf_out.close()
    with open(Path(outvcf).with_suffix('.stats'), 'w') as ostat:
        rem1 = sitecount - misscount - monositecount - multisitecount - orthogcount - singletons - ambiguous_count
        if keep_singles: rem1 += singletons
        if keep_ambigs: rem1 += ambiguous_count
        ostat.write(f"Infile: {invcf}. Outfile: {outvcf}\n")
        ostat.write(f"Total sites: {sitecount}. Missing ancestral information: {misscount}. Monosites: {monositecount}. Singletons (& dbs): {singletons} Multisites: {multisitecount}. Orthogs: {orthogcount}. Candidates: {rem1}\n")
        maj1 = rem1 - other_count - flipcount
        if keep_ambigs: maj1 -= ambiguous_count
        ostat.write(f"Candidates: {rem1}. Flipped: {flipcount}. Ambiguous: {ambiguous_count}. Other: {other_count}. Majority: {maj1}. RefRepair: {repaircount}\n")
        kept_alleles = flipcount + maj1
        if kept_alleles == 0: kept_alleles = 1
        ostat.write(f"Final alleles: {kept_alleles}. Reference: {maj1} ({maj1 / kept_alleles * 100 :.2f}%). Flipped: {flipcount} ({flipcount / kept_alleles * 100 :.2f}%)\n")





def fasta2est(invcf, fasta1, fasta2, outfile='test.est', outvcf = 'test.vcf.gz', writevcf=True, filtermiss=True): # these last two are to enable some control... I want to create a vcf alongside and also 


    #thresh = 1000000
    misscount = 0
    multisitecount = 0
    monositecount = 0
    orthogcount = 0
    minorcount1 = 0
    majorcount1 = 0
    minorcount2 = 0
    majorcount2 = 0
    alleledict = {"A": 0, "C": 1, "G": 2, "T": 3}
    # if fasta2:
    #     exit('No usage yet for second fasta file')

    vcf = VCF(Path(invcf))
    fasta1 = Fasta(Path(fasta1), sequence_always_upper=True) # can't remember what this last bit means
    fasta2 = Fasta(Path(fasta2), sequence_always_upper=True)
    vcf_out = Writer(Path(outvcf), vcf, "wz") # compressed vcf
    vcf_out.write_header()

    with open(outfile, 'w') as ofi: # this is our data-file... 

        #
        for n, var in enumerate(vcf): # walk through vcf as main strategy /// note: -- needs to be included above... 
            # if n > thresh:
            #     break # catchall

            fastavar1 = fasta1[var.CHROM][var.POS-1] 
            fastavar2 = fasta2[var.CHROM][var.POS-1]
            # at this point we will assume that vcf filtering has occurred previously... (perhaps we want to write a vcf here?) # include vcf writing in the future? 

            if fastavar1 == 'N' and fastavar2 == 'N': # extremely limited... # note... these sites would not be passed on to downstream fasta... which is a problem... 
                misscount += 1
                continue

            # at this point, should only be working with values present in appropriate file. Now want to get the vcf variants (from the file)... (ref & alt)... 
            # get rid of irrelevant sites... 
            if var.is_mnp or var.is_indel or var.is_sv: # flying blind... 
                continue # and not counting these... hopefully this leaves us in the ballpark of appropriate variation
            alleles = [var.REF] + var.ALT
            if len(alleles) > 2:
                multisitecount += 1
                continue
            # have to use allele counts... 
            if var.aaf == 0: # probably this means monosite
                monositecount += 1 # these are now allowed... 
                #continue
            if not fastavar1 in alleles and not fastavar2 in alleles: # i.e. if outgroup doesn't match the sample data... [in future we will include checks for formosus, but make them more gentle]
                orthogcount += 1
                continue

            # now we can assume that there is a match between ancestrals and alleles, so can get up to more comprehensive processing
            anc_template1 = ["0", "0", "0", "0"]
            if not fastavar1 in ['N', '*']:
                anc_template1[alleledict[str(fastavar1)]] = "1" # assumes fasta and only one allele
            anc1 = ",".join(anc_template1)
            anc_template2 = ["0", "0", "0", "0"]
            if not fastavar2 in ['N', '*']: # other solution would be to increase the lookup, but would have some issues... 
                anc_template2[alleledict[str(fastavar2)]] = "1"
            anc2 = ",".join(anc_template2)
            # now to try and add up the raw files... 
            ## we will assume biallelic snps!... 
            pop_template = ["0", "0", "0", "0"]
            pop_template[alleledict[var.REF]] = str(int(2 * var.num_hom_ref + var.num_het)) # 

            if not len(var.ALT) == 0 and not var.ALT[0] == '*': # get rid of singletons...
                pop_template[alleledict[var.ALT[0]]] = str(int(2 * var.num_hom_alt + var.num_het)) # this REQUIRES biallelic sites!!! WHICH WE WILL NOT HAVE... MIGHT BE FINE THOUGH
            pop1 = ",".join(pop_template)
            if not fastavar1 == var.REF:
                minorcount1 += 1
            else:
                majorcount1 += 1
            if not fastavar2 == var.REF:
                minorcount2 += 1
            else:
                majorcount2 += 1

            # now write the line... 
            ofi.write("\t".join([pop1, anc1, anc2]) + '\n')
            vcf_out.write_record(var)
    print('done')
    vcf_out.close()
    print(f"Missing: {misscount}\tMonosite: {monositecount}\tMultisite: {multisitecount}\tOrthogonal: {orthogcount}\t Minority1: {minorcount1}\tMajority1: {majorcount1}\tMinority2: {minorcount2}\tMajority2: {majorcount2}")
    total_sites = misscount + monositecount + multisitecount + orthogcount + minorcount1 + majorcount1
    print(f"Total sites: {total_sites}\tCoded percentage: {(majorcount1 + minorcount1) / total_sites * 100 :.2f}%\tMinority1 (as % of coded): {minorcount1 / (majorcount1 + minorcount1) * 100:.2f}%")


# writing meta function for est-sfs inference. we need to essentially pre-process a file according to a few statistics... let's try a single file... 

def est_preprocess(in_vcf, out_vcf, fasta1, fasta2, tempfolder, binwidth=5, cmin=20, cmax=False, bincountmin = 10000, 
                   monoallelic=True, multiallelic=False, keep_nonsnp=False, keep_ambigs=False, keep_singles=False):
    """This utility function gets key statistics with which to pre-process a vcf prior to passing it to broader processing files"""

    vcf = VCF(Path(in_vcf))
    fasta1 = Fasta(Path(fasta1), sequence_always_upper=True) # can't remember what this last bit means
    fasta2 = Fasta(Path(fasta2), sequence_always_upper=True)

    if not cmax:
        cmax = len(vcf.samples) * 2

    sitecount = 0

    # set up binning... 

    bins = np.array([n for n in range(0, cmax+1, binwidth)])
    idx = [n for n in range(len(bins))]
    bincount = np.zeros(len(bins), dtype=np.uint32)

    # step 1. get the arrays. we can probably make this faster by just scrapping lowqual arrays, but it doesn't really matter. 

    for site in vcf:
        sitecount += 1
        # if sitecount > thresh:
        #     break

        v = VarParsed(site)
        # now we want to try and get a site count... 
        if not v.passes_filter(monoallelic, multiallelic, keep_nonsnp): # missing biallelic, ambigs, & singletons at the moment... 
            continue
        if v.has_del: continue
        ac = v.an
        if ac < cmin: continue
        if ac > cmax: ac = cmax
        aix = ac // binwidth

        fastavar1 = fasta1[site.CHROM][site.POS-1] 
        fastavar2 = fasta2[site.CHROM][site.POS-1]

        if fastavar1 == 'N' and fastavar2 == 'N': # extremely limited... # note... these sites would not be passed on to downstream fasta... which is a problem... 
            continue

        alleles = v.geno_alleles # the proper kind

        if not fastavar1 in alleles and not fastavar2 in alleles: # i.e. 
            continue

        bincount[aix] += 1

    # Step 2. Rework the data to produce an index object that enables us to map alleles in the second pass to the appropriate place. 
    bincount2 = np.zeros(len(bincount), dtype=np.uint32)
    rmax = max(idx) # last good r value... 
    rprev = 0 # last bad r value... 
    rprevprev = max(idx)
    for r in reversed(idx):
        if bins[r] < cmin:
            idx[:r+1] = (0 for n in range(r+1))
            break
        rcount = bincount2[r] + bincount[r]
        if rcount < bincountmin:
            if bins[r-1] >= cmin:
                bincount2[r-1] += rcount
                bincount2[r] = 0
                if not rprev:
                    rprev = r
            else: # we are at the lower limit. basically we either discard or round down... 
                rcount += bincount2[rmax]
                bincount2[r] = rcount
                bincount2[rmax] = 0
                for n in range(rprev, r, -1):
                    idx[n] = r
                idx[rmax] = r
                if idx[rprevprev] == rmax:
                    for n in range(rprevprev, rmax, -1):
                        idx[n] = r
        else:
            bincount2[r] = rcount
            rmax = r # now this is the upper limit... 
            if rprev:
                for n in range(rprev, r, -1):
                    idx[n] = r
                rprevprev = rprev
            rprev = 0 # reset this variable
    # print(idx)
    # print(bincount)
    # print(sum(bincount))
    # print(bincount2)
    # print('\n\n')

    
    # step 3. go over the vcf again and generate appropriate mappings (virtually)

    # work on mappings... 
    fmapping = np.unique(idx, return_inverse=True)
    fbins = bins[fmapping[0]]
    tempfolder = Path(tempfolder)
    fnames = [tempfolder / f"est_{fbin}" for fbin in fbins]
    writers = [open(fname, 'w') for fname in fnames]
    # print(fnames)
    fidx = fmapping[1]
    # print(fbins)
    # print(fidx)
    # the scheme is that the bin maps onto file index (fidx) which maps on to the writer (writers) whcih writes the entry

    rng = np.random.default_rng() # fix later? doing things once, so probably not that bad an option... in long run want a seed... 
    bincount = np.zeros(len(bins), dtype=np.uint32)
    sitecount = 0
    misscount = 0
    multisitecount = 0
    monositecount = 0
    orthogcount = 0
    minorcount1 = 0
    majorcount1 = 0
    minorcount2 = 0
    majorcount2 = 0

    vcf.close()
    vcf = VCF(Path(in_vcf)) # have to load it twice as the file fails otherwise
    vcf_out = Writer(Path(out_vcf), vcf, "wz") # compressed vcf
    vcf_out.write_header()
    alleledict = {"A": 0, "C": 1, "G": 2, "T": 3}
    for site in vcf:
        sitecount += 1
        # if sitecount < 50000:
        #     continue
        # if sitecount > thresh:
        #     break

        fastavar1 = fasta1[site.CHROM][site.POS-1] 
        fastavar2 = fasta2[site.CHROM][site.POS-1]
        # at this point we will assume that vcf filtering has occurred previously... (perhaps we want to write a vcf here?) # include vcf writing in the future? 

        if fastavar1 == 'N' and fastavar2 == 'N': # extremely limited... # note... these sites would not be passed on to downstream fasta... which is a problem... 
            misscount += 1
            continue

        v = VarParsed(site)
        if not v.passes_filter(monoallelic, multiallelic, keep_nonsnp):
            continue

        if v.has_del:
            continue

        alleles = v.geno_alleles # the proper kind

        if not fastavar1 in alleles and not fastavar2 in alleles: # i.e. 
            orthogcount += 1
            continue

        # now we can assume that there is a match between ancestrals and alleles, so can get up to more comprehensive processing
        anc_template1 = ["0", "0", "0", "0"]
        if not fastavar1 in ['N', '*']:
            anc_template1[alleledict[str(fastavar1)]] = "1" # assumes fasta and only one allele
        anc1 = ",".join(anc_template1)
        anc_template2 = ["0", "0", "0", "0"]
        if not fastavar2 in ['N', '*']: # other solution would be to increase the lookup, but would have some issues... 
            anc_template2[alleledict[str(fastavar2)]] = "1"
        anc2 = ",".join(anc_template2)
        # now to try and add up the raw files... 
        ## we will assume biallelic snps!... 
        pop_template = ["0", "0", "0", "0"]

        # now we want to try and get a site count... 
        gts = v.gts
        ac = v.an
        if ac < cmin: continue
        if ac > cmax: ac = cmax
        aix = ac // binwidth # index value... 
        bincount[idx[aix]] += 1 # here is the key filter step. we use the index to sort everything out. 
        gts = rng.choice(gts, size=bins[idx[aix]], replace=False) # set shuffle to True? 

        gc = np.unique_counts(gts)
        v_a = gc[0]
        v_c = gc[1]

        alleles = [v.vcf_alleles[u] for u in v_a] # makes equivalent to geno_alleles... 

        for i, g in enumerate(alleles):
            pop_template[alleledict[g]] = str(v_c[i]) # should be allele count at at right spot... and multiallelic... [we have gotten rid of random alleles... ]

        pop1 = ",".join(pop_template)
        if not fastavar1 == site.REF:
            minorcount1 += 1
        else:
            majorcount1 += 1
        if not fastavar2 == site.REF:
            minorcount2 += 1
        else:
            majorcount2 += 1

        # write the line
        writers[fidx[aix]].write('\t'.join([pop1, anc1, anc2]) + '\n')
        vcf_out.write_record(site)

    for writer in writers:
        writer.close()

    np.save(tempfolder / '.fbins.npy', fbins)
    np.save(tempfolder / '.fnames.npy', fnames)
    np.save(tempfolder / '.fidx.npy', fidx)
    np.save(tempfolder / '.idx.npy', idx)
    np.save(tempfolder / '.bins.npy', bins)
    np.save(tempfolder / '.bincount.npy', bincount)
    paramvals_int = np.array([binwidth, cmin, bincountmin], dtype=np.int32)
    paramvals_bool = np.array([monoallelic, multiallelic, keep_nonsnp, keep_ambigs, keep_singles], dtype=np.bool)
    np.save(tempfolder / '.pri.npy', paramvals_int)
    np.save(tempfolder / '.prb.npy', paramvals_bool)
    (tempfolder / '._preprocessed').touch()


def run_est(tempfolder, model=1, nrandom=10, estprogram=None, threads=1, repair=False):
    '''Runs EST-SFS given a tempfolder prepared using prior steps in the pipeline. 
    Parameters:
        tempfolder      -   Path object showing temporary folder where est-sfs data is stored for the segment to be processed
        model           -   [Integer] Number of model to be applied (default 1): as follows: 
                                Jukes-Cantor:           0
                                Kimura 2-parameter:     1   *default
                                Rate-6 model:           2
        nrandom         -   [Integer] Numer of random starting runs to use (default 10) [rec at least 10 for rate-6]
        estprogram      -   [String] either None or path object showing location of EST-SFS program to run
        threads         -   [Integer] Number of threads to use (default 1)
        repair          =   [Bool] whether to attempt to complete an incomplete run [TODO: fix this with a new results object... ]
'''
    # load files & set up values
    if estprogram is None: 
        estprogram = 'est-sfs'
    else:
        estprogram = Path(estprogram)
    tempfolder = Path(tempfolder)
    fnames = np.load(tempfolder / '.fnames.npy', allow_pickle=True)
    rng = np.random.default_rng()
    config = tempfolder / '.estconfig'

    # write shared config file once
    with open(config, 'w') as ecf:
        ecf.write(f'n_outgroup 2\n')
        ecf.write(f'model {model}\n')
        ecf.write(f'nrandom {nrandom}')

    # Pre-generate unique seeds for each file
    seeds = rng.integers(1, 999_999_999, size=len(fnames), endpoint=True)

    def run_one(fname, seed):
        fname=str(fname) # in case path or scalar
        print(fname)

        seed_file = f"{fname}.seed"
        sfs_file = f"{fname}.sfs"
        anc_file = f"{fname}.anc"
        check_file = Path(f"{fname}.done")
        if repair and check_file.exists(): 
            return False

        # Write seed file
        with open(seed_file, "w") as sfi:
            sfi.write(str(int(seed)))

        # Run est-sfs
        cmd = [
            estprogram, 
            str(config), 
            fname, 
            seed_file, 
            sfs_file, 
            anc_file,
        ]
        # Raises CalledProcessError if non-zero exit
        subprocess.run(cmd, check=True)
        check_file.touch()

    # Run in parallel, up to `threads` subprocesses at once
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(run_one, fname, seed): (fname, seed)
            for fname, seed in zip(fnames, seeds)
        }

    # Optional: iterate to catch/print errors early
    for fut in as_completed(futures):
        fname, seed = futures[fut]
        try:
            fut.result()
        except Exception as e:
            # error handling customize
            raise RuntimeError(f"run_est failed for {fname} (seed {seed})") from e
        
    (tempfolder / "._ran").touch()
        
    # for fname in fnames:
    #     print(fname)
    #     with open(f'{fname}.seed', 'w') as sfi:
    #         sfi.write(str(rng.integers(1, 999_999_999, endpoint=True)))
    #     subprocess.run(f'{estprogram} {config} {fname} {fname}.seed {fname}.sfs {fname}.anc', shell=True)
        

def run_est_inference_full(in_vcf, out_vcf, outfasta1, outfasta2, tempfolder=None, binwidth=5, cmin=20, bincountmin=10000,
                           monoallelic=True, multiallelic=True, keep_ambigs=False, keep_nonsnp=False, keep_singles=True, 
                           model=1, nrandom=10, estprogram=None, pthresh=0.05, threads=1, force=True, cleanup=True, repair=False, override=False):
    '''Runs the entirety of the EST-SFS pipeline. Note that on spartan a module has to be loaded for this (GSL)
    PARAMETERS:
                        FILES AND FOLDERS
        in_vcf      -   VCF file (from upstream pipeline). Expect .vcf.gz
        out_vcf     -   VCF file to write final results to. Will include AA info column
        outfasta1   -   FASTA file containing first outgroup consensus sequence (mapped to main genome) - should be indexed
        outfasta2   -   FASTA file containing second outgroup consensus sequence (mapped to main genome) - should be indexed
        tempfolder  -   Folder to contain processing files. SHould be unique. Will be prefixed with `.tmp_`. If None` will be constructed from output VCF stem.

                        PREPROCESSING
        binwidth    -   [Integer] Allele count binwidth for splitting up EST-SFS run. 
        cmin        -   [Integer] Minimimum allele count for passing through EST-SFS. Counts below this will be dropped from downstream VCF files currently. 
        bincountmin -   [Integer] Minimum number of sites that a site bin must have to be valid for EST-SFS processing. Effects how subset chunking occurs. 
        monoallelic -   [Bool] whether to keep (real) monoallelic sites for the analysis. Defaults to True as EST-SFS requires invariant sites to work properly
        multiallelic-   [Bool] whether to keep (real) multiallelic sites. Defaults to False as EST-SFS only supplies info about majority allele, so minor alles are ambig.
        keep_ambigs -   [Bool] whether to preserve sites of ambiguous ancestry in downstream analysis. Defaults to False. Currently not implemented. 
        keep_nonsnp -   [Bool] whether INDELs, etc. should be allowed. Defaults to False as EST-SFS can't use them (I think). 
        keep_singles-   [Bool] whether to keep (true) singletons. Defaults to true as EST-SFS relies on this for site inference. Currently not implemented

                        EST-SFS PROGRAM
        model       -   [Integer] Which mutation model to apply [0: JC, 1: Kimura 2-P, 2: Rate-6] (default 1)
        nrandom     -   [Integer] Number of random starting points to use for likelihood search (default 10) rec at least 10
        estprogram  -   [None or filepath] - location of est-sfs program to run analysis with. Defaults if None. 

                        POSTPROCESSING
        pthresh     -   [Float] Probability threshold for retaining ancestral alleles (if keep_ambigs is False). Alleles kept if p < pthresh or p > (1 - pthresh)

                        MULTITHREAD
        threads     -   [Integer] How many threads to use for the main step (default 1)

                        TEMPDATA MANAGEMENT
        override       -   [Bool] whether to override (i.e. delete) files in `tempfolder` before running the analysis. If false, will throw error if prior files present
        cleanup     -   [Bool] whether to cleanup (delete) the tempdirectory and files after analysis. Default True. 
        repair      -   [Bool] whether to attempt to repair an incomplete run. Default False
        force       -   [Bool] whether to complete analysis at all costs. Will switch `repair` to `override` if inappropriate intermediate files are present. Default True. 
        '''

    # set up folders and files. 
    active_repair = False
    ran = False
    in_vcf, out_vcf, outfasta1, outfasta2 = map(Path, (in_vcf, out_vcf, outfasta1, outfasta2))
    if tempfolder is None:
        tvcf = out_vcf.stem
        if tvcf.endswith('.vcf'):
            tvcf = tvcf[:-4]
        tempfolder = out_vcf.parent / ('.tmp_' + tvcf)
    else: 
        tempfolder = Path(tempfolder)
        tempfolder = tempfolder.parent / ('.tmp_' + tempfolder.name)
    if tempfolder.is_dir():
        if override:
            for f in tempfolder.iterdir():
                f.unlink()
        else:
            if len([f for f in tempfolder.iterdir()]) > 0:
                if repair:
                    if (tempfolder / '._preprocessed').exists():
                        active_repair = True
                        if (tempfolder / '._ran').exists():
                            ran = True
                    else:
                        if force:
                            print(f"Preprocess step not completed. As `force` is True, overwriting...")
                            for f in tempfolder.iterdir():
                                f.unlink()
                        else:
                            exit(f"Junk files in folder and preprocessing step is not completed. Consider setting `force` to True to override")
                else:
                    if force:
                        for f in tempfolder.iterdir():
                            f.unlink()
                    else:
                        if (tempfolder / '._preprocessed').exists():
                            exit(f'Directory {tempfolder} already exists, and preprocessing is complete. Consider setting `force` to True to override or `repair` to true to continue')
                        else:
                            exit(f'Directory {tempfolder} already exists! Consider setting `force` to True to override')
    else:
        tempfolder.mkdir()

    if estprogram is None:
        estprogram='est-sfs'

    # intermediate files
    anc_vcf = tempfolder / 'est_anc.vcf.gz' # not that specific at this point - others handled by intermediate... 

        # now run the script. 
    if not active_repair:
        est_preprocess(in_vcf, anc_vcf, outfasta1, outfasta2, tempfolder, binwidth=binwidth, cmin=cmin, bincountmin=bincountmin, 
                    monoallelic=monoallelic, multiallelic=multiallelic, keep_nonsnp=keep_nonsnp, keep_ambigs=keep_ambigs, 
                    keep_singles=keep_singles)
    if not ran:
        run_est(tempfolder, model=model, nrandom=nrandom, estprogram=estprogram, threads=threads, repair=repair)
    est_postprocess_group(anc_vcf, out_vcf, tempfolder, filter=False, pthresh=pthresh)

    if cleanup:
        for f in tempfolder.iterdir():
            f.unlink()
        tempfolder.rmdir()
    
    

def load_samplefile(fname):
    '''Returns list (str) of sample names to include in filtering event'''
    with open(fname) as fi:
        return [s.strip() for s in fi if len(s.strip()) > 0]

def quick_filter(infile, outfile, samplefile=None, max_miss=0.1, remove_multi_singletons = True, depth=False, depth_thresh=10):

    invcf = VCFReader(infile)
    if samplefile is not None:
        invcf.set_samples(load_samplefile(samplefile))
    lowcount = int(remove_multi_singletons)

    if remove_multi_singletons:
        infodict = {'ID': 'MTS','Number': 0, 'Type': 'Flag', 'Description': f'{lowcount}-tons converted to missing'}
        invcf.add_info_to_header(infodict)
    logger.info(invcf.samples)
    outvcf = VCFWriter(outfile, invcf)
    outvcf.write_header()

    #print(invcf.samples[invcf.idx])
    #print(outvcf.samples)

    count = 0
    nonallelic =0
    rem_single = 0
    missing = 0
    nonancestral=0
    noref = 0
    dels = 0
    good = 0
    ins = 0

    sitestart = 0
    site_curr = 0

    for site in invcf:
        count += 1
        site_curr = int(site.var.POS)
        if not sitestart:
            sitestart = int(site.var.POS)
        if not depth and site.missing_frac > max_miss:
            missing += 1
            continue
        elif depth and site.low_depth_frac(depth_thresh) > max_miss:
            missing += 1
            continue
        if not site.is_biallelic and not site.is_monoallelic: # we're going to try 
            if remove_multi_singletons and site.has_lowcount(lowcount):
                v1, v2 = site.geno_alleles, site.v_c
                site = site.lowcount_to_miss(lowcount)
                if site:
                    site.INFO['MTS'] = True # we'll see if this works... 
                    if not site.is_biallelic and not site.is_monoallelic:
                        nonallelic += 1
                        continue
                    else:
                        if not depth and site.missing_frac > max_miss: # reckeck after running operation
                            missing += 1
                            continue
                        elif depth and site.low_depth_frac_alleles(depth_thresh) > max_miss:
                            missing += 1
                            continue
                        rem_single += 1
                        
                else:
                    nonallelic += 1
                    continue
            else:
                nonallelic += 1
                continue
        if 'AA' in dict(site.var.INFO).keys() and not site.var.INFO['AA'] in site.geno_alleles:
            nonancestral += 1
            continue
        if not site.has_ref:
            noref += 1
            continue
        if site.has_del:
            dels += 1
            continue
        if site.is_nonsnp:
            ins += 1
            continue
        # if site.recode: # we are going to ride close to the wind here... 
        #     if not site.is_ordered:
        #         continue
        good += 1
        outvcf.write_record(site)

    outvcf.close()
    span_total = site_curr - sitestart
    logger.info(f"Total sites: {count}\tmultiallelic: {nonallelic}\tmulti_restored: {rem_single}\tmissing: {missing}\tno_ref: {noref}\tdel: {dels}\tins: {ins}\tno_ancestral: {nonancestral}")
    logger.info(f"Good: {good}\t{good/count * 100:.2f}%\tActual: {good/span_total * 100:.2f}%")



def vcf_make_age_basic(invcf, outvcf): # next project... 
    vcf = VCFReader(invcf)
    # if not vcf.vcf.contains('AA'):
    #     exit('alleles must be polarised, but no ancestral allele tag present!')
    infodict = {'ID': 'AGE','Number': 1, 'Type': 'Float', 'Description': 'allele age'}
    vcf.add_info_to_header(infodict)
    outvcf = VCFWriter(outvcf, vcf)
    for site in vcf:
        if site.allelicity > 1:
            altfreq = 1 - max(site.v_c) / sum(site.v_c) # freq if not ref... 
            aage = ((-2 * altfreq)/(1 - altfreq)* log(altfreq)) * 10_000 // 1 / 10000 # maybe make these more specific? 
            site.INFO['AGE'] = aage
            site.INFO['AF'] = altfreq * 1000 // 1 / 1000
        outvcf.write_record(site)

def vcf_make_age_rule(invcf, outvcf, mask01, mask02, mask03):
    pass



def vcf_substitute(vcf, template_vcf, out_vcf):
    '''Subsitute missing alleles from template site... '''

    vcfs = VCFCompare(vcf, template_vcf)

    # vcfs.subsample(5)
    # print(vcfs.samples)
    mismatch = 0
    snps = 0
    mts = 0
    count = 0
    good = 0

    outvcf = VCFWriter(out_vcf, vcfs, vcfs.vcf1.header) # opens the file
    outvcf.write_header()

    for n, (source, template) in enumerate(vcfs):

        if not source.genotypes_replace(template): # need sanity check here... 
            mismatch += 1
            continue

        if source.is_biallelic:
            snps += 1

        good += 1
        outvcf.write_record(source)
        count = n
    
    outvcf.close()
    logger.info(f"Alleles total: {count}\tAlleles written: {good}\tIncompatible alleles: {mismatch}\tRemoved MTS alleles: {mts}\tSNPs: {snps}")


def vcf_diagnostic(invcf, fa1, fa2, thresh=1000):

    invcf = VCFReader(invcf)
    fasta1 = Fasta(Path(fa1), sequence_always_upper=True)
    fasta2 = Fasta(Path(fa2), sequence_always_upper=True)
    count = 0
    refcount=0
    altcount = 0
    majcount = 0
    mincount = 0
    outmatch = 0
    outmismatch = 0
    monomismatch = 0
    missing = 0

    for site in invcf:
        # if site.INFO['AAP'] < 0.99:
        #     print(site.INFO['AAP'], '\t', site.INFO['AA'], '\t', site.var.REF, '\t', site.geno_maj, '\t', site.var.ALT)

        if not site.is_multiallelic:
            if site.is_monoallelic and site.INFO['AA'] != site.var.REF:
                monomismatch += 1
                print(site.has_ref)
            continue


        if site.INFO['AA'] == 'N':
            missing += 1
            continue
        print(site.POS)

        fv1 = fasta1[site.CHROM][site.POS-1]
        fv2 = fasta2[site.CHROM][site.POS-1]

        aa = site.var.INFO['AA']
        if aa == site.var.REF: refcount += 1
        else: altcount += 1
        if aa == site.geno_maj: 
            majcount += 1
        else: 
            mincount += 1
            print(site.geno_maj, aa, site.geno_2nd, '|', fv1, fv2, [v for v in site.geno_alleles], [int(c) for c in site.v_c], aa == fv1 or fv1=='N' and aa == fv2)
        if aa == fv1 or aa == fv2: 
            outmatch += 1
        else: 
            outmismatch += 1


    print(f'REF ancs: {refcount}\tALT ancs: {altcount}\tMAJ ancs: {majcount}\tMIN ancs: {mincount}\tOUT: {outmatch}\tNOT: {outmismatch}\tMONOMIS: {monomismatch}\tAMBIG: {missing}')

def vcf_diagnostic2(in1, in2, fasta1, fasta2, thresh=1000):

    vcfs = VCFCompare(in1, in2)
    fasta1 = Fasta(Path(fasta1), sequence_always_upper=True)
    fasta2 = Fasta(Path(fasta2), sequence_always_upper=True)
    count = 0

    for site1, site2 in vcfs:
        if site1.allelicity > 2:# and (site1.INFO['AA'] != site1.var.REF or site2.INFO['AA'] != site2.geno_maj):

            count += 1
            fastavar1 = fasta1[site1.CHROM][site1.POS-1] # masc
            fastavar2 = fasta2[site1.CHROM][site1.POS-1] # brom
            print(f'{site1.POS} DEFAULT\t{site1.var.REF}\t{site1.var.ALT}\tMAJ:{site1.INFO['AA']} {site1.INFO['AA'] == site1.var.REF}\tEST:{site2.INFO['AA']} {site2.INFO['AA']==site2.geno_maj}\t{site2.INFO['AAP']}\tMASC: {fastavar1}\tBROM: {fastavar2}\t{site1.geno_alleles} {[int(c) for c in site1.v_c]}')

            if count >= thresh:
                break


if __name__ == "__main__":
    base = Path('/data/gpfs/projects/punim1778/')
    #fasta_in = base / "genomes/GCF_002204515.2_AaegL5.0_genomic.fa"
    workdir = base / "Projects/aegypti/2023/data/previous_studies/core_files/"
    outdir = base / "Projects/aegypti/2023/scripts/"
    fasta_in1 = workdir / "mascloc_consensus.fa"
    fasta_in2 = workdir / "bromloc_consensus.fa"
    #vcf_in = workdir / "chr2loc2_NC_035108.1:9451401-10451400_d3m50mac2_bgl_rehead.vcf.gz"
    vcf_in = workdir / "mj_init_may24_NC_035108.1:340000001-350000000_aegallraw_bglm.vcf.gz"
    vcf_in = workdir / "mj_init_may24_NC_035108.1:470000001-474425716_aegallraw_bglm.vcf.gz"
    anc_in = workdir / "mj_init_may24_NC_035108.1:470000001-474425716_raw_aegpac_anc_pvalues.txt"

    vcf_bgl = workdir / "mj_init_may24_NC_035108.1:470000001-474425716_aegallraw_bgl.vcf.gz"
    #temp_in = workdir / "mj_init_may24_NC_035108.1:470000001-474425716_raw_aegpac_ancfilt.vcf.gz"
    vcf_in = 'ptest_bglm.vcf.gz'
    vcf_out = outdir / "postprocess_test_ancfilter.vcf.gz"
    vcf_out = outdir / "ptest_est_majfilter.vcf.gz"
    vcf_out1 = outdir / "ptest_majfilter.vcf.gz"
    vcf_out2 = outdir / "ptest_estfilter.vcf.gz"
    vcf_out3 = outdir / "ptest_rulefilter.vcf.gz"
    #vcf_post_beagle_rehead(invcf=vcf_in, tempvcf=temp_in, outvcf=vcf_out)

    invcf = 'basetemp.vcf'
    outfix = 'chr2_maj_gather'
    stripvcf = 'chr2_maj_stripped'
    #vcf_gather_preinfer(outvcf=outfix + '.vcf.gz', fpostfix = '_aegpac_maj_bgla')
    #est.strip_header_chrom_info(infile=outfix, ofile = stripvcf)
    #vcf_num_sites(stripvcf)
    # #v2t.vcf2ts_sample(stripvcf)
    assign_ancestrals_majority(vcf_in, vcf_out1, keep_ambigs=True, monoallelic=True, multiallelic=True, repair=True)
    # estimate_ancestrals_basic(vcf_in, vcf_out3, fasta_in1, fasta_in2, keep_ambigs=True, monoallelic=True, multiallelic=True, repair=True)
    # #est_vcf_postprocess(vcf_in, anc_in, vcf_out)
    # #fasta2est(vcf_in, fasta_in1, fasta_in2)
    # #est_preprocess(vcf_in, vcf_out, fasta_in1, fasta_in2, out='.tester', monoallelic=True, multiallelic=False, bincountmin=5000)
    # # have to import i.e. module load GSL/2.7
    # #est_postprocess_group(vcf_out, 'estpost.vcf.gz', '.tester')
    run_est_inference_full(vcf_in, vcf_out2, fasta_in1, fasta_in2, force=True, bincountmin=10000, keep_ambigs=True, monoallelic=True, multiallelic=True)
    # #tester = VCFReader(vcf_in)
    # #writ = VCFWriter('tester.vcf.gz', tester.header)
    #writ.rectest()
    #tester.write_header('tester.vcf.gz')
    #quick_filter(vcf_out, 'tester.vcf.gz', 'aeg_asiapac_new.txt', max_miss=0.05)
    if not Path(str(vcf_bgl) + '.tbi').exists():
        subprocess.run(f'tabix {vcf_bgl}', shell=True)
    # quick_filter('dummy.vcf.gz', 't4.vcf.gz', 'aeg_asiapac_new.txt', max_miss=0.9)
    #subprocess.run(f'bcftools annotate -x ^FORMAT/GT,FORMAT/PL,FORMAT/AD,FORMAT/RGQ,FORMAT/GQ,^INFO/AC,INFO/AF,INFO/MQ,INFO/AA,INFO/RRP,INFO/MTS {'t4.vcf.gz'} | bcftools view -a -Oz -o t5.vcf.gz', shell=True)
    # vcf_substitute('t5.vcf.gz', vcf_bgl, 't6.vcf.gz')
    #vcf_make_age_basic('t6.vcf.gz', 'd2.vcf.gz')
    dcf_maj = workdir / "mj_init_may24_NC_035108.1:220000001-230000000_aegallraw_ancrule.vcf.gz"
    dcf_est = workdir / "mj_init_may24_NC_035108.1:220000001-230000000_aegallraw_bglm.vcf.gz"
    #vcf_diagnostic(vcf_out3, fasta_in1, fasta_in2, thresh=1000)
    vcf_diagnostic2(vcf_out1, vcf_out2, fasta_in1, fasta_in2, thresh=100)



# note: cyvcf2 adds an index via using the VCF.set_index(idxfile) command. This is very helpful!!! (I probably need to adopt this in future)... 