#!/bin/python3

import gzip
#from Bio import bgzf
from itertools import compress
import subprocess
from pathlib import Path
import cyvcf2
from treeflows import vcf_core
import numpy as np
#import pysam
#from pysam import VariantFile

spartan = True

if spartan:
    bin = '/data/gpfs/projects/punim1778/bin'
else:
    bin = '~/bin'

alleledict = {"A": 0, "C": 1, "G": 2, "T": 3}

mascids = ["SRR11006707", "SRR11006708", "SRR11006709", "SRR11006710"]
est_chrom_length_dict = {1: 310827022, 2: 474425716, 3: 409777670}
est_chrom_name_dict = {1: 'NC_035107.1', 2: 'NC_035108.1', 3: 'NC_035109.1'}
est_chrom_alias_dict = {1: 'X1', 2: 'X2', 3: 'X3'}

path2singer = '/data/gpfs/projects/punim1778/programs/singer/SINGER/releases'
# load vcf file...

def vcf_get_fraction_missing(infile: str) -> float:
    """Return site density (sites per bp) over the observed VCF interval."""
    vcf = cyvcf2.VCF(infile)
    sites  = 0
    pos_start = 0
    pos_last = 0
    for site in vcf:
        if pos_start == 0:
            pos_start = site.POS
        pos_last = site.POS
        sites += 1
    return sites / (pos_last - pos_start) # authentic interval... 

def vcf_get_fraction_missing_range(infile: str, start: int, end: int) -> float:
    """Return site density (sites per bp) for a specified coordinate range."""
    vcf = cyvcf2.VCF(infile)
    sites  = 0
    for site in vcf:
        sites += 1
    return sites / (end - start) # authentic interval... 

def vcf_get_site_count(infile: str) -> int:
    """Count and print the number of non-header records in `infile + '.vcf.gz'`."""
    with gzip.open(infile + '.vcf.gz', "r") as ifi:
        sitecount = 0

        for ll in ifi:
            jj = ll.decode()
            if jj.startswith("#"):
                continue
            sitecount += 1
        print(f"VCF Sitecount: {sitecount}")
        return(sitecount)

vcf_get_site_count.__doc__ = "Returns how many sites (lines) a .vcf.gz file contains. Assumes each line corresponds to a unique site"

def state_get_site_count(infile: str) -> int:
    """Count and print the number of sites in a `.min4.phy.state` file."""
    with open(infile + '.min4.phy.state') as ifi:
        sitecount = 0

        for ll in ifi:
            if not ll.startswith('Node1\t'):
                continue
            sitecount += 1
        print(f"State Sitecount: {sitecount}")
        return(sitecount)


(310827022, 474425716, 409777670)

def poslimit(pos, clen):
    """Clamp a position to at most `clen`."""
    if pos > clen:
        return(clen)
    return(pos)

def range_finder_multi(startpos, endpos, chrom, prefix, flen, postfix):
    """Return a list of interval-VCF filenames covering a genomic range."""
    chromdict = {1: 'NC_035107.1', 2: 'NC_035108.1', 3: 'NC_035109.1'}
    clendict = {1: 310827022, 2: 474425716, 3: 409777670}
    clen = clendict[chrom]
    start_base = startpos // flen 
    end_base = (endpos // flen)
    if start_base == end_base: # if everything is contained in the same sequence
        return([f'{prefix}{chromdict[chrom]}:{startpos // flen * flen + 1}-{poslimit((startpos // flen +1)*flen, clen)}{postfix}.vcf.gz', ])
    else:
        flist = []
        for base in range(start_base, end_base + 1):
            bfile = f'{prefix}{chromdict[chrom]}:{base * flen + 1}-{poslimit((base +1)*flen, clen)}{postfix}.vcf.gz'
            flist.append(bfile)
    return(flist)

def vcf_filter_range_multifile(startpos: int, endpos: int, chrom: int, prefix: str = 'mj_init_may24_', flen = 10_000_000, postfix = '_subsnps', outpost='_subsnps', prefix_new=None) -> int:
    """Write a ranged subset VCF by gathering records across interval-partitioned VCFs."""
    if prefix_new is None:
        prefix_new = prefix
    chromdict = {1: 'NC_035107.1', 2: 'NC_035108.1', 3: 'NC_035109.1'}
    ### This file only works with a max of two files... 
    clendict = {1: 310827022, 2: 474425716, 3: 409777670}
    clen = clendict[chrom]
    file1 = f'{prefix}{chromdict[chrom]}:{startpos // flen * flen + 1}-{poslimit((startpos // flen +1)*flen, clen)}{postfix}.vcf.gz' # not safe for final chromosome length... (?)... 
    file2 = f'{prefix}{chromdict[chrom]}:{endpos // flen * flen + 1}-{poslimit((endpos // flen + 1)*flen, clen)}{postfix}.vcf.gz'
    print(file1)
    print(file2)
    if file1 == file2:
        print('files are the same!')
        same = True
    else: 
        print('files are different!')
        same = False
    flist = range_finder_multi(startpos=startpos, endpos=endpos, chrom=chrom, prefix=prefix, flen=flen, postfix=postfix)

    start_vcf = cyvcf2.VCF(file1)
    # with pysam.VariantFile(file1, 'rb') as sv:
    #     for n, var in enumerate(sv.fetch()):
    #         if n < 20:
    #             print(var)
    #         else:
    #             break
    owrite = cyvcf2.Writer(f'{prefix_new}{chromdict[chrom]}:{startpos}-{endpos}{outpost}.vcf.gz', start_vcf, mode='wz')
    owrite.write_header()
    for ff in flist:
        vcf = cyvcf2.VCF(ff)
        for variant in vcf:
            if variant.POS >= startpos and variant.POS <= endpos:
                owrite.write_record(variant)
            elif variant.POS > endpos:
                break
    owrite.close()
    print('done!')



def get_out_miss(infile: str, truncate: int =0, outids: list = mascids) -> tuple:
    """Return (missing_count, missing_fraction) for outgroup calls in a VCF."""

    with gzip.open(infile, "r") as ifi:

        linecount = 0
        misscount = 0
        nlines=truncate
        for ll in ifi:
            jj = ll.decode()
            if jj.startswith("#CHROM"):
                outfilter = getoutfilter(jj, outids)
            if jj.startswith("#"):
                continue
            if truncate:
                if linecount > nlines:
                    missfrac = misscount / linecount
                    print(f"Missing ancestral sites: {misscount}\nFraction: {missfrac}")
                    return(misscount, missfrac)
            linecount += 1
            genotypes = parse_genos(jj)
            outgroup = getgroup(genotypes, outfilter)
            refcount, altcount = getcounts(outgroup)
            if refcount + altcount <= 0:
                misscount += 1
        missfrac = misscount / linecount
        print(f"Missing ancestral sites: {misscount}\nFraction: {missfrac}")
        return(misscount, missfrac)

get_out_miss.__doc__ = """Returns a tuple (misscount, missfrac) of the number and fraction of sites 
within a .vcf.gz file (infile) that have missing data for outgroups (outids)"""

def get_allele_mean(infile: str, truncate: int=0, outids: list = mascids):
    """Compute mean number of called alleles per site for the ingroup."""

    with gzip.open(infile, "r") as ifi:

        sitecount = 0
        allelecount = 0
        nsites =truncate
        for ll in ifi:
            jj = ll.decode()
            if jj.startswith('#CHROM'):
                infilter = getinfilter(jj, outids)
            if jj.startswith('#'):
                continue
            if truncate:
                if sitecount > nsites:
                    print(f"Mean alleles per site: {allelecount / sitecount}")
                    return(allelecount / sitecount)
            genotypes = parse_genos(jj)
            ingroup = getgroup(genotypes, infilter)
            refcount, altcount = getcounts(ingroup)
            if refcount + altcount <= 0:
                continue
            allelecount += refcount + altcount
            sitecount += 1
        mean_alleles = allelecount / sitecount
        print(f"Mean alleles per site: {mean_alleles}")
        return(mean_alleles)
    
def fasta2est(infile: str, outfile: str, truncate: int=0, fixed: int=0, outids: list=mascids):
    """Placeholder for creating EST-SFS inputs from FASTA (not implemented)."""
    # need to compare the files... 
    pass

def vcf2est(infile: str, outfile: str, truncate: int=0, fixed: int=0, outids: list = mascids):
    """Convert a VCF (with outgroup samples) to an EST-SFS input file."""

    with open(outfile, "w") as ofi:
        ofi.write("")

    with gzip.open(infile, "r") as ifi:

        linecount = 0
        nlines=truncate
        for ll in ifi:
            # decoding as we assume it is gzipped
            jj = ll.decode()
            # skip over comments

            if jj.startswith("#CHROM"):
                infilter, outfilter = getgroupfilters(jj, outids)


            if jj.startswith("#"):
                continue         
            ### temporary linecount code
            if truncate:
                linecount += 1
                if linecount > nlines:
                    return(0)

            # parse lines presume snp, but not biallelic... 
            genotypes, reff, altt = parse(jj)
            ingroup, outgroup = splitgroups(genotypes, infilter, outfilter)

            # pass to next function (but could easily do more here...)
            outline = make_estline(ingroup, outgroup, reff, altt, fixed)

            with open(outfile, "a") as ofi:
                ofi.write(outline)
        return(0)

def vcf_make_dummy_chrom(infix, outfix, startpos, endpos, name, chrom_true, rem=False):
    """Rewrite a VCF interval onto a dummy contig with remapped coordinates."""
    chrlength = endpos - startpos + 1
    with gzip.open(infix + '.vcf.gz', 'r') as ifi, gzip.open(outfix + '.vcf.gz', 'w') as ofi:
        writeswitch = False
        for ii in ifi: 
            jj = ii.decode()
            if jj.startswith('##contig'):
                writeswitch = True
                continue
            if writeswitch:
                ofi.write(f"##contig=<ID={name},length={chrlength},assembly=AaegL5>\n".encode())
                writeswitch = False
            if jj.startswith('#'):
                ofi.write(ii)
                continue
            
            # now to handle sequence recoding... 
            parsedlist = jj.strip().split(sep='\t')
            old_pos = int(parsedlist[1])
            if parsedlist[0] != chrom_true: continue # filter out irrelevant chromosomes
            if old_pos < startpos or old_pos > endpos: continue # filter out snps outside of range... 
            new_pos = old_pos - startpos + 1
            newlinelist = [name, str(new_pos)]
            newlinelist.extend(parsedlist[2:])
            outline = '\t'.join(newlinelist) + '\n'
            ofi.write(outline.encode())
    if rem:
        Path(infix + '.vcf.gz').unlink()
    return(0)







def zip_ests(target: str, ancestral: str, outfile: str):
    """Combine (zip) target and ancestral EST files line-by-line into one file."""
    with open(target) as targ, open(ancestral) as anc, open(outfile, "w") as out:
        for t, a in zip(targ, anc):
            tt = t.strip().split()[0]
            aa = a.strip().split()[1]
            out.write(f'{tt} {aa}\n')

def vcf_filter_ancestrals(infile: str, outfile: str, outids: list=mascids, rem=False):
    """Filter sites where the outgroup is entirely missing and write a new VCF."""
    with gzip.open(infile + '.vcf.gz', "r") as ifi, gzip.open(outfile + '.vcf.gz', "wb") as ofi:
        misscount = 0
        sitecount = 0

        for ll in ifi:
            jj = ll.decode()
            if jj.startswith('#CHROM'):
                outfilter = getoutfilter(jj, outids)
            if jj.startswith("#"):
                ofi.write(ll)
                continue
            sitecount += 1
            genotypes = parse_genos(jj)
            outgroup = getgroup(genotypes, outfilter)
            refcount, altcount = getcounts(outgroup)
            if refcount + altcount <= 0:
                misscount += 1
                continue
            ofi.write(ll)
        print(f"Filtered {misscount} sites with missing ancestral data. {sitecount-misscount} remaining.")
    if rem:
        Path(infile + '.vcf.gz').unlink()
    return(0)
    
def vcf_filter_ancestrals_template(focusfile: str, ancfile: str, outfile: str, outids: list=mascids, rem=False):
    """Filter `focusfile` sites using outgroup-missingness evaluated in `ancfile`."""
    # filters vcf based on ancestry presence in another vcf, but otherwise leaves the vcf the same. (needs to be combined later with a 
    # zip operation on the two .est files produced from both vcfs... )
    with gzip.open(focusfile + '.vcf.gz', "r") as focfi, gzip.open(ancfile + '.vcf.gz', "r") as ancfi, gzip.open(outfile + '.vcf.gz', "wb") as ofi:
        misscount = 0
        sitecount = 0

        for fl, al in zip(focfi, ancfi):
            aj = al.decode()
            if aj.startswith('#CHROM'):
                fj = fl.decode()
                if not fj.startswith('#CHROM'):
                    exit("Files do not have even headers!!!")
                outfilter = getoutfilter(aj, outids)
                if not outfilter == getoutfilter(fj, outids):
                    exit("Files do not have same arrangement of individuals!")
            if aj.startswith("#"):
                ofi.write(fl)
                continue
            sitecount += 1
            if not checkpos(aj, fl.decode()):
                exit("Files do not have same arrangement of sites!")
            genotypes = parse_genos(aj)
            outgroup = getgroup(genotypes, outfilter)
            refcount, altcount = getcounts(outgroup)
            if refcount + altcount <= 0:
                misscount += 1
                continue
            ofi.write(fl)
        print(f"Filtered {misscount} sites with missing ancestral data in template. {sitecount-misscount} remaining.")
    if rem:
        Path(focusfile + '.vcf.gz').unlink()
    return(0)

def strip_header_chrom(infile, ofile, chrom='NC_035109.1'):
    """Remove contig header lines except for the specified contig."""
    with gzip.open(infile + '.vcf.gz', "r") as ifi, gzip.open(ofile + '.vcf.gz', "w") as ofi:
        for iline in ifi:
            if not iline.decode().startswith('##contig'):
                ofi.write(iline)
                continue
            elif not iline.decode().startswith(f'##contig=<ID={chrom}'):
                continue
            ofi.write(iline)
    return(0)

def strip_header_info(infile, ofile, info = 'AA'):
    """Remove INFO header lines except for the specified INFO ID."""
    with gzip.open(infile + '.vcf.gz', "r") as ifi, gzip.open(ofile + '.vcf.gz', "w") as ofi:
        for iline in ifi:
            if not iline.decode().startswith('##INFO'):
                ofi.write(iline)
                continue
            elif not iline.decode().startswith(f'##INFO=<ID={info}'):
                continue
            ofi.write(iline)
    return(0)

def strip_header_chrom_info(infile, ofile, chrom='NC_035108.1', info = 'AA'):
    """Filter header to keep only a specific contig and INFO field definitions."""
    with gzip.open(infile + '.vcf.gz', "r") as ifi, gzip.open(ofile + '.vcf.gz', "w") as ofi:
        header = True
        for iline in ifi:
            if header and iline.decode().startswith('#'): # we are in header territory
                hline = iline.decode()
                if hline.startswith('#CHROM'):
                    header = False
                    ofi.write(iline)
                    continue
                elif hline.startswith('##INFO') and not hline.startswith(f'##INFO=<ID={info}'):
                    continue
                elif hline.startswith('##contig') and not hline.startswith(f'##contig=<ID={chrom}'):
                    continue
                ofi.write(iline)
            else:
                ofi.write(iline)
    return(0)

def vcf_reheader(infile, ofile, tempfile, rem=False): # going to modify so that we can get rid of the extra chromosomes... 
    """Reheader a VCF by copying the template header (with NC contigs) then data."""
    with gzip.open(infile + '.vcf.gz', "r") as ifi, gzip.open(tempfile + '.vcf.gz', "r") as tfi, gzip.open(ofile + '.vcf.gz', "w") as ofi:

        tempiter = iter(tfi)
        initer = iter(ifi)
        while True:
            
            # go through template file and copy lines... 
            while True:
                try:
                    templine = next(tempiter)
                except StopIteration:
                    templine = None
                    break
                tt = templine.decode()
                if not tt.startswith('##'):
                    break
                if tt.startswith('##contig') and not tt.startswith('##contig=<ID=NC'):
                    continue
                ofi.write(templine)

            # go through other file
            while True:
                try:
                    inline = next(initer)
                except StopIteration:
                    inline = None
                    break
                if not inline.decode().startswith('##'): # we lose some info here... 
                    ofi.write(inline)
                    break
            while True:
                try:
                    inline = next(initer)
                except StopIteration:
                    inline = None
                    break
                ofi.write(inline)
            if templine is None and inline is None:
                break
        if rem:
            Path(infile + '.vcf.gz').unlink()
        return(0)

def vcf_flip_ancestrals_majority(invcf, outvcf, phased=True, rem=False):
    """Add an `AA` INFO field based on majority allele, without changing genotypes."""

    vcf = cyvcf2.VCF(invcf)
    adict = {'ID': 'AA', 'Number': '1', 'Type': 'String', 'Description': "Ancestral allele (majority-derived)"}
    vcf.add_info_to_header(adict)
    vwrite = cyvcf2.Writer(outvcf, vcf, 'wz')

    vwrite.write_header()
    for variant in vcf:
        # want to work out whether to flip it... 
        if variant.num_hom_alt > variant.num_hom_ref: # if we need to flip
            variant.INFO['AA'] = variant.ALT[0]
        else: 
            pass
            #variant.INFO['AA'] = variant.REF
        
        vwrite.write_record(variant)





def vcf_flip_ancestrals_iqtree(invcf, outvcf, statefile, node='Node1', phased=True, complex=False, postfix='.min4.phy', rem=False): # assuming Node1 corresponds to route... somewhat probably? NEEDS REFINEMENT
    """Flip VCF alleles based on ancestral states inferred by IQ-TREE `.state` output."""

    with gzip.open(invcf + '.vcf.gz', 'r') as ifi, gzip.open(outvcf + '.vcf.gz', 'w') as out, open(statefile + f'{postfix}.state') as sfi:
        flipcount = 0
        linecount = 0
        ancount = 0
        if complex:
            coalcount = 0
        initer = iter(ifi)
        anciter = iter(sfi)

        while True:
            # read the two files in parallel
            # 
            # move thorugh the infile. write header and keep going if line is commented
            while True:
                try:
                    inline = next(initer)
                except StopIteration:
                    inline = None
                    break
                if not inline.decode().startswith('#'):
                    break
                out.write(inline)
            
            # move through the reffile. keep going if line commented or opener... 
            while True:
                try:
                    refline = next(anciter)
                except StopIteration:
                    refline = None
                    break
                if (refline.startswith(f'{node}\t')): # '#' is the comment line... then nodes are numbered according ot a tricky scheme. I assume this the route
                        break
            
            # check that files are same length
            if (inline is None or refline is None) and not (inline is None and refline is None):
                exit("Non-commented sections of files are of different lengths!")
            
            # terminate loop if file ends have been reached
            if inline is None and refline is None:
                break

            # main recoding task
            # check flipped status of refline
            linecount += 1
            ii = inline.decode()
            iref, ialt = parse_ref(ii)
            if  not complex and checkmatch_ancestral(refline, iref, ialt): # test if reference matches ancestral
                out.write(inline)
                continue

            if complex:
                cmatch = checkmatch_ancestral_complex(refline, iref, ialt) # assign value: 0 = ref/anc don't match, 1 = they do, 2 = unresolved (-> remove)... 
                if cmatch == 2:
                    coalcount += 1
                    continue
                elif cmatch == 1:
                    out.write(inline)
                    continue

            # now flip the line
            flipcount += 1
            outline = make_flippedline(ii, phased) + '\n'
            out.write(outline.encode())
    print(f'flipcount: {flipcount}. linecount: {linecount}. ratio: {flipcount/linecount}')
    if complex:
        print(f'no_coalescent_count: {coalcount}. ratio: {flipcount / coalcount}')
    if rem:
        Path(invcf + '.vcf.gz').unlink()
    return(0)

            
def vcf_flip_ancestrals_est(invcf, outvcf, pfile, phased=True, rem=True): # only works  if on same chromosomal interval
    """Flip VCF alleles based on EST-SFS probability file output."""
    # presumptions: pfile relates to invcf; 1/1 correspondence between alleles... 
    with gzip.open(invcf + '.vcf.gz', "r") as ifi, gzip.open(outvcf + '.vcf.gz', "w") as out, open(pfile) as pfi:
        flipcount = 0
        linecount = 0
        initer = iter(ifi)
        anciter = iter(pfi)

        while True:
            # we are in the main loop... now we want to read our two files in parallel

            # move through the infile. write header and keep going if line is commented. 
            while True:
                try:
                    inline = next(initer)
                except StopIteration:
                    inline = None
                    break
                if not inline.decode().startswith('#'):
                    break
                out.write(inline)

            # move through the reffile. keep going if line starts with 0
            while True:
                try:
                    refline = next(anciter)
                except StopIteration:
                    refline = None
                    break
                if not refline.startswith('0'):
                    break
            
            # check that files are same length
            if (inline is None or refline is None) and not (inline is None and refline is None):
                exit("Non-commented sections of files are of different lengths!")
            
            # terminate loop if file ends have been reached
            if inline is None and refline is None:
                break

            # main recoding task... 

            # check flip status of refline. 
            linecount += 1
            if not checkflipped(refline):
                out.write(inline)
                continue

            # now flip the line... 
            flipcount += 1
            outline = make_flippedline(inline.decode(), phased) + '\n'
            out.write(outline.encode())
    if rem:
        Path(invcf + '.vcf.gz').unlink()
    print(f'flipcount: {flipcount}. linecount: {linecount}. ratio: {flipcount/linecount}')
    return(0)


def make_flippedline(rawline, phased=True): # assumes line has been decoded
    """Return a VCF record line with REF/ALT swapped and genotypes flipped."""
    metadata, genotypes = parse_full(rawline)
    reff = metadata[3]
    altt = metadata[4]

    if phased:
        flipped_genotypes = [flipidv_phased(i) for i in genotypes]
    else: 
        flipped_genotypes = [flipidv_unphased(i) for i in genotypes]
    final = list(metadata)
    final[3] = altt
    final[4] = reff
    final.extend(flipped_genotypes)
    return('\t'.join(final)) 

def flipidv_phased(idv: str):
    """Flip a phased genotype string by swapping alleles (0<->1)."""
    return(flip_gt(idv))

def flipidv_unphased(idv: str):
    """Flip an unphased sample field, including allele depths (AD) order."""
    idv_parsed = idv.strip().split(sep=':')
    idv_latter = ':'.join(idv_parsed[2:])
    gt = idv_parsed[0]
    ad = idv_parsed[1]    
    gt_flipped = flip_gt(gt)
    ad_parsed = ad.split(sep=',')
    ad_flipped = ','.join(reversed(ad_parsed))
    return(':'.join([gt_flipped, ad_flipped, idv_latter]))

def idv_checkmiss(idv: str)->bool:
    """Return True if the sample genotype is missing."""
    idv_parsed = idv.strip().split(sep=':')
    gt = idv_parsed[0]
    return(gt_checkmiss(gt))

def gt_checkmiss(gt: str, phased=False):
    """Return True if a GT string is missing (`./.` or similar)."""
    if phased:
        exit('should not be phased')
    g = gt.split(sep='/')[0]
    return(g == '.')

def idvs_checkmiss_all(idvs: list)->bool: # for working with selection masks -> return TRUE if all idvs missing
    """Return True if all genotypes in a list of sample fields are missing."""
    genos = [idv.strip().split(sep=':')[0].split(sep='/')[0] for idv in idvs] # possibly add | separator at some point... though may be imputed? 
    return(all(geno == '.' for geno in genos))

def idv_makemiss(idv: str, phased: bool = True)->str:
    """Return a missing GT string in phased form (`.|.`)."""
    if not phased:
        exit('should be phased')
    return('.|.')



def flip_gt(gt: str):
    """Flip a diploid GT string by reversing alleles and applying `flip`."""
    dv = gt[1]
    gt_parsed = gt.split(sep=dv)
    gt_reversed = list(reversed(gt_parsed))
    gt_flipped = [flip(c) for c in gt_reversed]
    return(dv.join(gt_flipped))

def flip(c: str):
    """Flip an allele code ('0' <-> '1'); leave '.' unchanged."""
    if c == '.': return(c)
    elif c == '1': return('0')
    elif c == '0': return('1')
    exit('Invalid genotype!')

def getgroupfilters(titleline, outidvs):
    """Return boolean masks selecting ingroup/outgroup samples from a VCF header line."""
    ids = parse_genos(titleline)

    infilter = [x not in outidvs for x in ids]
    outfilter = [x in outidvs for x in ids]
    return(infilter, outfilter) 

def getoutfilter(titleline, outidvs):
    """Return a boolean mask selecting outgroup samples from a VCF header line."""
    ids = parse_genos(titleline)
    return([x in outidvs for x in ids])

def getinfilter(titleline, outidvs):
    """Return a boolean mask selecting ingroup samples from a VCF header line."""
    ids = parse_genos(titleline)
    return([x not in outidvs for x in ids])



def format_alleles_raw(refcount, altcount, reff, altt):
    """Format allele counts as a 4-vector string (A,C,G,T) for EST-SFS."""
    outvals = ["0", "0", "0", "0"]
    outvals[alleledict[reff]] = str(int(refcount))
    outvals[alleledict[altt]] = str(int(altcount))
    return(",".join(outvals))

def format_alleles_ancestral(refcount, altcount, reff, altt):
    """Format a one-hot ancestral allele vector (majority allele)."""
    outvals = ["0", "0", "0", "0"] # workaround
    if refcount + altcount > 0:
        if refcount >= altcount:
            outvals[alleledict[reff]] = "1"
        else:
            outvals[alleledict[altt]] = "1"
    return(",".join(outvals))

def format_alleles_fixed(refcount, altcount, reff, altt, total):
    """Rescale counts to a fixed total while preserving proportions (approx)."""
    if refcount == 0 and altcount == 0:
        return(format_alleles_raw(refcount, altcount, reff, altt))
    if refcount == 0:
        return(format_alleles_raw(0, total, reff, altt))
    if altcount == 0:
        return(format_alleles_raw(total, 0, reff, altt))
    reffrac = (refcount + altcount) / refcount
    altfrac = (refcount + altcount) / altcount
    rcnew = total // reffrac # floor division
    acnew = total // altfrac
    if not rcnew + acnew == total:
        if (total % reffrac) >= (total % altfrac):
            rcnew += 1
        else:
            acnew += 1 # possibly not fool-proof
    return(format_alleles_raw(rcnew, acnew, reff, altt))

def make_estline(ingroup, outgroup, reff, altt, fixed=0):
    """Build a single EST-SFS input line from ingroup/outgroup genotypes."""
    refcount, altcount = getcounts(ingroup)
    refcounto, altcounto = getcounts(outgroup)

    if fixed:
        instr = format_alleles_fixed(refcount, altcount, reff, altt, fixed)
    else:
        instr = format_alleles_raw(refcount, altcount, reff, altt)
    outstr = format_alleles_ancestral(refcounto, altcounto, reff, altt)

    return(" ".join([instr, outstr]) + '\n')



    
def getcounts(genotypelist): # this is a list of raw genotype/other format markers for individuals... 
    """Count ref/alt alleles (0/1) from a list of VCF sample fields."""
    refcount = 0
    altcount = 0

    for i in range(len(genotypelist)):   # iterate through individuals [currently all]
        gtraw = genotypelist[i].split(sep=':')[0] # extract the genotype value
        if gtraw.startswith('.'):
            continue
        gtsplit = gtraw.split(sep='/')
        for j in gtsplit: # increase counters
            if j == '0':
                refcount += 1
            if j == '1':
                altcount += 1 # implicit error here - won't work with multi-allelic sites
    return(refcount, altcount)

def parse(rawline):
    """Parse a VCF record line into (sample_fields, REF, ALT)."""
    parsedlist = rawline.strip().split(sep='\t')
    reff = parsedlist[3]
    altt = parsedlist[4]
    genotypes = parsedlist[9:]
    if len(altt) > 1:
                if altt != "ALT":
                    print("Alt is > 1!!!")
    return(genotypes, reff, altt)

def parse_ref(rawline):
    """Parse a VCF record line into (REF, ALT)."""
    parsedlist = rawline.strip().split(sep='\t')
    reff = parsedlist[3]
    altt = parsedlist[4]
    if len(altt) > 1:
                if altt != "ALT":
                    print("Alt is > 1!!!")
    return(reff, altt)

def parse_genos(rawline):
    """Extract sample fields from a VCF line (assumes 9 fixed columns)."""
    return(rawline.strip().split(sep='\t')[9:])

def parse_full(rawline, biallelic=True):
    """Parse a VCF record line into (metadata_fields, sample_fields)."""
    parsedlist = rawline.strip().split(sep='\t')
    metadata = parsedlist[:9]
    genotypes = parsedlist[9:]
    altt = parsedlist[4]
    if len(altt) > 1 and biallelic:
                if altt != "ALT":
                    print("Alt is > 1!!!")
    return(metadata, genotypes)

def checkpos(rawline1, rawline2):
    """Return True if two VCF record lines share the same POS value."""
    pos1 = rawline1.strip().split(sep='\t')[1]
    pos2 = rawline2.strip().split(sep='\t')[1]
    return(pos1 == pos2)

def splitgroups(genlist, infilter, outfilter):
    """Split a sample field list into ingroup/outgroup using boolean masks."""
    ingroup = list(compress(genlist, infilter))
    outgroup = list(compress(genlist, outfilter))
    return(ingroup, outgroup)

def getgroup(genlist, filter):
    """Select a subset of sample fields using a boolean mask."""
    return(list(compress(genlist, filter)))


def get_flipcount(pfile: str)-> int:
    """Count how many sites in an EST-SFS probability file imply a flip (<0.5)."""
    with open(pfile) as pf:
        flipcount=0
        sitecount=0
        for l in pf:
            if l.startswith("0"):
                continue
            sitecount += 1
            if checkflipped(l):
                flipcount += 1
        print(f"{flipcount} sites flipped out of {sitecount}total; fraction flipped: {flipcount/sitecount}")
        return(flipcount)
    
def checkmatch_ancestral(rawline, reff, altt):
    """Return True if inferred ancestral allele equals REF; validate membership."""
    anc = rawline.strip().split()[2]
    if not anc.upper() in {'T', 'A', 'G', 'C'}:
        exit(f"Incorrect genotype in inferred line: {anc}")
    if not anc.upper() in {reff.upper(), altt.upper()}:
        exit(f"Ancestral allele does not match vcf options: {anc} vs ({reff}, {altt}")
    if anc.upper() == reff.upper():
        return(True)
    return(False) # should be fine as we have previously checked set membership

def checkmatch_ancestral_complex(rawline, reff, altt): # this is an out of date script, misunderstanding iqtree's output... 
    """Complex ancestral match checker supporting ambiguous codes (legacy)."""
    anc = rawline.strip().split()[2]
    if not anc.upper() in {'A', 'G', 'R', 'N', 'D', 'C', 'Q', 'E', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'}:
        exit(f"Incorrect genotype in inferred line: {anc}")
    if not anc.upper() in {'A', 'C', 'G', 'T'}:
        return(2) # 2 = skip case (presumption is sites haven't coalesced)
    if not anc.upper() in {reff.upper(), altt.upper()}:
        exit(f"Ancestral allele does not match vcf options: {anc} vs ({reff}, {altt})")
    if anc.upper() == reff.upper():
        return(1) # the 'True case'
    return(0) # the False case

def checkflipped(rawline):
    """Return True if the reference-ancestral probability indicates flipping (<0.5)."""
    refprob = rawline.strip().split()[2]
    if float(refprob) < 0.5:
        return(True)
    else:
        return(False)

def make_map_beagle(inmap, outmap, chrom_num=1): # can't expose the deeper structs anymore... 
    """Convert a linkage map into Beagle map format for a given chromosome."""
    chromcode = est_chrom_name_dict[chrom_num]
    alias = est_chrom_alias_dict[chrom_num]
    with open(inmap, "r") as imap, open(outmap, "w") as omap:

        for iline in imap:
            iparse = iline.strip().split(sep='\t')
            if not iparse[1] == alias:
                continue
            varid = iparse[2]
            genpos = iparse[4]
            physpos = iparse[3]
            omap.write(f'{chromcode}\t{varid}\t{genpos:.5f}\t{physpos}\n')
        return(0)
    
def make_map_shapeit(inmap, outmap, chrom_num=1, ratecalc=False, minrate=0.001, maxrate=25, defaultrate=0.1): # now essentially alias with dif default
    """Alias for `make_map_relate` with SHAPEIT-oriented defaults."""
    make_map_relate(inmap, outmap, chrom_num, ratecalc, minrate, maxrate, defaultrate)
    return 0

def make_map_relate(inmap, outmap, chrom_num=1, ratecalc=True, minrate=0.001, maxrate=25, defaultrate=0.1): # going to remove alias now (and sub with dict function) # implemented at top... 
    """Convert a linkage map into Relate/Shapeit-style recombination map format."""
    # define a simple ratelimit function, rr. 

    def rr(rate_val, rmin = minrate, rmax=maxrate):
        if rmin is not None and rate_val < rmin:
            return rmin
        if rmax is not None and rate_val > rmax:
            return rmax
        return rate_val
    alias = est_chrom_alias_dict[chrom_num] # constants (???)
    with open(inmap, 'r') as imap, open(outmap, 'w') as omap: # let's try and get the main field... 

        physnext = 0 # i.e. the start of the genome (?) (I think this is right)
        gennext = 0 # we'll see

        for iline in imap:
            iparse = iline.strip().split(sep='\t')
            if not iparse[1] == alias:
                continue
            physpos = physnext # load prev
            genpos = gennext # load the previous values
            physnext = int(iparse[3])
            gennext = float(iparse[4])
            rate = (gennext - genpos)/(physnext - physpos) * 1e6 if ratecalc else defaultrate
            omap.write(f'{physpos}\t{rr(rate):.3f}\t{genpos:.5f}\n') # current frequencies... 
    # now round out the map... 
    # key assumption: the rate is equal to the last rate measure... (rmin)
    with open(outmap, 'a') as omap:
        omap.write(f'{physnext}\t{rr(rate):.3f}\t{gennext:.5f}\n') # should last to the end of the file... 
    return 0

def run_beagle(infile, outfile, chrom, map=0, nthreads=8, burnin=15, iters=30): # inputs complete filename, outputs extension on prefix
    """Run Beagle 5.4 genotype imputation/phasing on a VCF prefix."""

    beaglestring = f'java -jar {bin}/beagle_5.4.jar '
    beaglestring += f'gt={infile}.vcf.gz '
    beaglestring += f'out={outfile} '
    beaglestring += f'burnin={burnin} '
    beaglestring += f'iterations={iters} '
    if map:
        beaglestring += f'map={map} '
    beaglestring += f'chrom={chrom} '
    beaglestring += f'nthreads={nthreads}'
    subprocess.run(beaglestring, shell=True, check=True, stderr=subprocess.PIPE)

def run_shapeit(infile, outfile, region, map=0, nthreads=8, log='shapeit.log'): # files are prefixes only
    """Run SHAPEIT on a region and write a bgzipped/indexed VCF."""
    shapestring = f'phase_common_static --thread {nthreads} '
    shapestring += f'--input {infile} --output {outfile}.bcf '
    if map:
        shapestring += f'--map {map} '
    shapestring += f'--log {log} --region {region}'
    print('running shapeit')
    subprocess.run(shapestring, shell=True, check=True, stderr=subprocess.PIPE)
    print('converting to .vcf.gz')
    subprocess.run(f'bcftools view {outfile}.bcf -Oz -o {outfile}.vcf.gz', shell=True, check=True, stderr=subprocess.PIPE)
    print('indexing .vcf.gz')
    subprocess.run(f'bcftools index {outfile}.vcf.gz', shell=True, check=True, stderr=subprocess.PIPE)

def strip_ancestrals(invcf, outvcf, outids=mascids):
    """Drop outgroup samples from a VCF using bcftools sample exclusion."""
    stripstring = f'bcftools view {invcf} -Oz '
    stripstring += f'-o {outvcf} '
    stripstring += f'--samples ^{",".join(outids)}'
    subprocess.run(stripstring, shell=True, check=True, stderr=subprocess.PIPE)
# program update... names will be defined in a more rigourous scheme... 


def run_parallel_singer(infile, outfile, n=20, Ne='1e4', m='3e-9', start=90_000_001, end=91_000_000,
                        chrom="NC_035109.1", L=100_000, ratio=0.1, thin=20, ncores=10): 
    """Run SINGER in parallel mode over a genomic interval."""
    singstring = f'{path2singer}/parallel_singer_new -Ne {Ne} -m {m} -L {L} -n {n} -thin {thin} '
    singstring += f'-chrom {chrom} -start {start} -end {end} -ratio {ratio} '
    singstring += f'-vcf {infile} -output {outfile} -num_cores {ncores}'
    subprocess.run(singstring, shell=True, check=True, stderr=subprocess.PIPE)

def run_singer(infile, outfile, start, end, n=20, Ne='1e4', m='3e-9', ratio=0.1, thin=20): # in and out here are prefixes
    """Run SINGER over an interval, then convert output to `.trees`."""
    print(f'unzipping file {infile}.vcf.gz')
    subprocess.run(f'gunzip {infile}.vcf.gz', shell=True, check=True, stderr=subprocess.PIPE) # fails if input is .vcf, which it is out of beagle
    print(f'running singer master over sequence {start}-{end}')
    singstring = f'{path2singer}/singer_master -Ne {Ne} -m {m} -ratio {ratio} -n {n} -thin {thin} '
    singstring += f'-vcf {infile} -output {outfile} -start {start} -end {end}'
    subprocess.run(singstring, shell=True, check=True, stderr=subprocess.PIPE)
    print('converting to .trees format')
    get_trees_singer(outfile, n)
    print(f'rezipping vcf {infile}.vcf') # possibly. TODO check if these need removing... 
    subprocess.run(f'gzip {infile}.vcf', shell=True, check=True, stderr=subprocess.PIPE)
    return(0)

def get_trees_singer(infile, n=20): # simple wrapper. 
    """Convert SINGER output to tskit `.trees` using `convert_to_tskit`."""
    treestring = f'{path2singer}/convert_to_tskit -input {infile} -output {infile} -start 0 -end {n} -step 1'
    subprocess.run(treestring, shell=True, check=True, stderr=subprocess.PIPE)

def vcf_depth_to_miss(infix, outfix=None, depth = 10, rem=False):
    """Converts all sites with depth < 'depth' to missing value (.)"""
    if outfix is None: 
        outfix = infix
        rem = False
    filterstring = f'bcftools +setGT {infix}.vcf.gz -Oz -o {outfix}.vcf.gz - -- -t q -n . '
    filterstring += f'-i "FORMAT/DP<{depth} | SMPL_MAX(FORMAT/PL)=0"'
    subprocess.run(filterstring, shell=True, check=True, stderr=subprocess.PIPE)
    if rem:
        Path(infix + '.vcf.gz').unlink()

def vcf_thin(infix, outfix, thin=1000, rem=False, idx=False, seed=None, thinmax=np.inf, thinmin=0, method="window"):
    """ Method can be either `exp` or `window` IF `thin`, window size is same as thin... """

    rng = np.random.default_rng(seed=seed) # will leave it at this... unreplicable  

    invcf = Path(str(infix) + ".vcf.gz")
    outvcf = Path(str(outfix) + ".vcf.gz")
    vcf = vcf_core.VCFReader(invcf)
    vout = vcf_core.VCFWriter(outvcf, vcf)
    allcount = 0
    trimcount = 0
    # now the filter step... need a function to trigger residual tracking of variant in thousands... perhaps exponential decay?
    if method == "exp":
        newpos = None
        for site in vcf:
            allcount += 1
            if newpos is None:
                newpos = site.POS-1 + max(min(rng.exponential(thin), thinmax), thinmin)
            if site.POS > newpos:
                trimcount += 1
                vout.write_record(site)
                newpos = None

    elif method == "window":
        window = None
        snp = None
        pcount = 0
        for site in vcf:
            allcount += 1
            if window is None:
                window = site.POS + thin
            if site.POS >= window:
                vout.write_record(snp)
                trimcount += 1
                pcount = 0
                while site.POS >= window:
                    window += thin # jump over long gaps... (assumes single chromosome)
            pcount += 1 # probability
            if rng.uniform() < 1/pcount:
                snp=site
        # clean up at end... (last window will be aborted)... 
        vout.write_record(snp)
        trimcount += 1
    vout.close()

    vcf_polish(outfix)
    if idx:
        vcf_index_tbi(outfix)
    if rem:
        vcf_fullrem(str(infix))
    print(f"Finished thinning: original snps: {allcount}. Thinned snps: {trimcount} ({trimcount/allcount*100:.2f}%)")

def vcf_filter_standard(infix, outfix=None, rem=False):
    """Placeholder for a standard bcftools filtering recipe (not implemented)."""
    if outfix is None:
        outfix = infix
        rem = False
    print("Not yet functional!")
    subprocess.run(f'bcftools filter --help', shell=True, check=True, stderr=subprocess.PIPE)
    if rem:
        Path(infix + '.vcf.gz').unlink()

def vcf_filter_nonsnp(infix, outfix=None, rem=False, all=False, idx=False):
    """Filter a VCF to SNP sites using `treeflows.vcf_core` parsing."""
    if outfix is None:
        outfix = infix
        rem = False

    vcf = vcf_core.VCFReader(f"{infix}.vcf.gz")
    vout = vcf_core.VCFWriter(f"{outfix}.vcf.gz", vcf)
    vout.write_header()
    for site in vcf:
        if site.passes_filter(monoallelic=all, multiallelic=True):
            vout.write_record(site)

    vout.close()

    vcf_polish(f"{outfix}")
 
    if idx: 
        vcf_index_tbi(outfix)
    if rem:
        vcf_fullrem(infix)

def vcf_filter_biallelic(infix, outfix=None, rem=False, all=False, idx=False):
    """Filter a VCF to (bi)allelic sites using `treeflows.vcf_core` parsing."""
    if outfix is None:
        outfix = infix
        rem = False

    vcf = vcf_core.VCFReader(f"{infix}.vcf.gz")
    vout = vcf_core.VCFWriter(f"{outfix}.vcf.gz", vcf)
    vout.write_header()
    for site in vcf:
        if site.passes_filter(monoallelic=all, biallelic=True):
            vout.write_record(site)

    vout.close()

    vcf_polish(f"{outfix}")
 
    if idx: 
        vcf_index_tbi(outfix)
    if rem:
        vcf_fullrem(infix)


def vcf_filter_mac(infix, outfix, mac, rem=False, idx=False):
    """Filter a VCF by minor allele count (MAC) using bcftools."""
    subprocess.run(f'bcftools view -Oz -o {outfix}.vcf.gz --min-ac {mac}:minor {infix}.vcf.gz', shell=True, check=True, stderr=subprocess.PIPE)
    if idx: 
        vcf_index_tbi(outfix)
    if rem:
        vcf_fullrem(infix)

def vcf_filter_missing(infix, outfix, miss, rem=False, idx=False): # note, missingness here works oppositely to vcftools (so low values = few missing sites left, not few filtered)
    """Filter a VCF by site missingness using bcftools `F_MISSING`."""
    res = subprocess.run(f'bcftools view -Oz -o {outfix}.vcf.gz -i "F_MISSING < {miss}" {infix}.vcf.gz', shell=True, check=True, stderr=subprocess.PIPE)
    if res.returncode != 0:
        print(res.stderr)
        exit(f"Error in bcftools missingness filtering: {res.returncode}")
    if idx: 
        vcf_index_tbi(outfix)
    if rem:
        vcf_fullrem(infix)

def vcf_filter_depth(infix, outfix, depth_frac, depth_threshold, rem=False, idx=False):
    """Filter by frequency of samples below depth (replaces broader allele filtering)"""
    subprocess.run(f"""bcftools view -Oz -o {outfix}.vcf.gz -i 'F_PASS(FORMAT/DP >= {depth_threshold}) > {depth_frac}' {infix}.vcf.gz""", shell=True, check=True, stderr=subprocess.PIPE)
    if idx: 
        vcf_index_tbi(outfix)
    if rem:
        vcf_fullrem(infix)

def vcf_fullrem(filefix):
    """Removes standard _.vcf.gz file and its two indexes if they exist"""
    Path(filefix + ".vcf.gz").unlink()
    if Path(filefix + ".vcf.gz.tbi").exists():
        Path(filefix + ".vcf.gz.tbi").unlink()
    if Path(filefix + ".vcf.gz.csi").exists():
        Path(filefix + ".vcf.gz.csi").unlink()

def get_popdict(popfile):
    """Parse a two-column `id<TAB>pop` file into a `pop -> [ids]` mapping."""
    # assume format is id\tloc with no header
    with open(popfile, "r") as popf:
        locdict = dict()
        for ll in popf.read().strip().split(sep='\n'):
            parsed = ll.strip().split(sep='\t')
            id, pop = parsed[0],  parsed[1]
            if not pop in locdict:
                locdict[pop] = [id,]
            else:
                locdict[pop].append(id)
    return(locdict)

def get_popmasks(vcffix, popdict): # uses filename. produces dict of locmappings (for good order)... 
    """Build per-population sample masks for a VCF prefix."""
    vcf = cyvcf2.VCF(vcffix + '.vcf.gz')
    return(get_popmasks_internal(vcf, popdict))

def get_popmasks_internal(vcf, popdict): # uses vyvcf2.VCF structure... 
    """Build per-population sample masks for a cyvcf2 VCF object."""
    idvs = vcf.samples
    popmasks = dict()
    for pop in popdict.keys():
        popidvs = popdict[pop]
        popmask = [x in popidvs for x in idvs]
        popmasks[pop] = popmask
    return(popmasks) 

def vcf_filter_missing_pops(infix, outfix, popfile, popthresh=3, rem=False):
    """Remove sites where any sufficiently-sized population has all samples missing."""
    ### This one will extract a list of location masks keyed to the vcf involved. needs a tsv/csv and a vcf? (or just the former and we take it on faith)
    # currently only eliminating sites for which ALL individuals are missing data. 

    ## phase 1: create the locdict to use to query vcf
    locdict = get_popdict(popfile)

    with gzip.open(infix + '.vcf.gz', "r") as vcff, gzip.open(outfix + '.vcf.gz', "w") as outf:
        locmasks = list()
        cutcount = 0
        linecount = 0
        for vline in vcff:
            if vline.decode().startswith("##"):
                outf.write(vline)
                continue
            if vline.decode().startswith("#CHROM"):
                vv = vline.decode()
                for loc in locdict.keys():
                    if len(locdict[loc]) >= popthresh:
                        locmask = getoutfilter(vv, locdict[loc]) # this function returns a mask of all genos keyed to loc... 
                        locmasks.append(locmask)
                outf.write(vline)
                continue
            linecount += 1
            vv = vline.decode()
            genos = parse_genos(vv)
            to_write = True
            for locmask in locmasks:
                idvs = list(compress(genos, locmask))
                if idvs_checkmiss_all(idvs):
                    to_write = False
                    cutcount += 1
                    break
            if to_write:
                outf.write(vline)
    print(f"Total SNPS: {linecount}; SNPs removed: {cutcount}; SNPs remaining: {linecount - cutcount}")
    if rem:
        Path(infix + '.vcf.gz').unlink()
    return(0)


def vcf_polish(filefix):
    """Takes a gzipped vcf (probably not in proper format) and uses bcftools to package it into proper vcf. workaround as .bgz format is very specific"""
    proc = subprocess.run(f'bcftools view -Oz -o {filefix}_temp.vcf.gz {filefix}.vcf.gz', shell=True, stderr=subprocess.PIPE, text=True, check=True) # properly gzipped... (?)
    print(proc.stderr)
    Path(f'{filefix}_temp.vcf.gz').rename(f'{filefix}.vcf.gz') # gets rid of the worst effects? 

def vcf_index_all(filefix, tbi=True, csi=True):
    """Indexes formates in both major modalities"""
    if tbi: vcf_index_tbi(filefix)
    if csi: vcf_index_csi(filefix)

def vcf_polish_all(filefix, tbi=True, csi=False):
    """Polishes and indexes vcf to ensure full compatibility with other programs"""
    vcf_polish(filefix)
    vcf_index_all(filefix, tbi=tbi, csi=csi)

def vcf_unzip(filefix): # basic utility
    """Uncompress `<filefix>.vcf.gz` to `<filefix>.vcf` using `gunzip`."""
    subprocess.run(f'gunzip {filefix}.vcf.gz', shell=True, check=True, stderr=subprocess.PIPE)

def vcf_zip(filefix, rem=False):
    """Compress `<filefix>.vcf` to `<filefix>.vcf.gz` using bgzip."""
    subprocess.run(f'bgzip -f -c {filefix}.vcf > {filefix}.vcf.gz', shell=True, check=True, stderr=subprocess.PIPE)
    if rem:
        Path(filefix + '.vcf').unlink()

def vcf_index_csi(filefix):
    """Index `<filefix>.vcf.gz` using bcftools (CSI/TBI depending on file)."""
    subprocess.run(f'bcftools index -f {filefix}.vcf.gz', shell=True, check=True, stderr=subprocess.PIPE)

def vcf_index_tbi(filefix):
    """Index `<filefix>.vcf.gz` using tabix (TBI)."""
    subprocess.run(f'tabix -f -p vcf {filefix}.vcf.gz', shell=True, check=True, stderr=subprocess.PIPE)

def vcf_filter_idvs(filefix, outfix, idvfile, rem=False): # list of individuals to include in the file... 
    """Filter a VCF to a provided sample list using bcftools `-S`."""
    subprocess.run(f'bcftools view -Oz -o {outfix}.vcf.gz -S {idvfile}.txt {filefix}.vcf.gz', shell=True, check=True, stderr=subprocess.PIPE)
    if rem:
        qq = Path(filefix + '.vcf.gz')
        qq.unlink()

def vcf_to_phylip(filefix):
    """Convert a VCF to PHYLIP format using `vcf2phylip`."""
    subprocess.run(f'vcf2phylip -i {filefix}.vcf.gz -o {filefix}.phy', shell=True, check=True, stderr=subprocess.PIPE)

def phylip_infer_ancestrals(filefix, outfix=None, nthreads=1, rootidv = 'SRR11006662', spartan=False): # outputs .state file. The outgroup is a highish coverage Kenya idv.
    """Infer ancestral states with IQ-TREE from a `.phy` alignment (writes `.state`)."""
    if outfix is None:
        outfix == filefix
    if spartan:
        execute = 'iqtree'
    else:
        execute = 'iqtree2'

    modelset = 'PMB+F+R7'
    subprocess.run(f'{execute} -s {filefix}.min4.phy -nt {nthreads} -m MFP -o {rootidv} --prefix {outfix} -asr --seqtype DNA', shell=True, check=True, stderr=subprocess.PIPE)


if __name__=="__main__":
    vcf2est("loc1m_biallelic-miss01.vcf.gz", "loc1m_biallelic-miss01_raw.est", fixed = False, truncate = 0, outids = mascids)
    #getoutmiss("loc1m_biallelic-miss01.vcf.gz", truncate = 0, outids=mascids)


# possiblility... define a python object class that takes us through the rest of the process??? 
# the idea would be that we create it with a constructor and it keeps track of the various intermediate files for us... 
# but it can also 'populate' its file reserves with pre-existing versions of the various sub-filess... 
# it could manage e.g. a branching hierarchy of file-names (& even folders!) and keep track of the variety of processing 
# steps connecting them all... (it could also make branching choices thorugh the differing domains and execute the intentions... 

# today's primary goals

#   2. Implement initial inferences with tsinfer... (sans missingness) -> probably find out how to get rid of 2 lowquals...
#       ii. figure out how to undo the interpolation process by beagle/shapeit (i.e. recode vcfs) and run this also (missingness)
#   3. Implement an object-based logic for managing the large variety of files that we are dealing with at this stage (database?)
#   4. Treeseq/TreeInfer tutorials
#   5. Sketch out how to implement global missingness (possibly by infering sequences in small clusters??? - or by a post-process filter?'
