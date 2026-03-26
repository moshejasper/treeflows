import subprocess
from cyvcf2 import VCF, Writer
from pathlib import Path
import gzip
from . import vcf2est as est # fix soon



def vcf_post_beagle_missrestore(invcf, outvcf, tempvcf): #inbuilt contig selector... 
    """postprocess beagle output to (a) restore missingness, (b) clean up headers, and (c) select NC_ chromosomes (not good longterm strategy)"""
    #
    if not Path(invcf + '.tbi').exists():
        subprocess.run(f'tabix {invcf}', shell=True)
    if not Path(tempvcf + '.tbi').exists():
        subprocess.run(f'tabix {tempvcf}', shell=True)
    vcf = VCF(invcf)
    vtemp = VCF(tempvcf)
    if not vcf.samples == vtemp.samples:
        exit("Samples are not equivalent!")
    vwrite = Writer(outvcf, vtemp, 'wz')
    #vwrite.write_header() # in long run, need to combine things a little better... 

    for n, (var, pvar) in enumerate(zip(vcf, vtemp)):

        if n > 100:
            break

        if var.POS != pvar.POS:
            exit(f'Positions not compatible: {var.POS} & {pvar.POS}')

        vwrite.write_record(var)
        vwrite.write_record(pvar)

def vcf_beagle_postprocess(beagvcf, refvcf, outvcf):
    
    with gzip.open(beagvcf) as ifi, gzip.open(outvcf, 'w') as out, gzip.open(refvcf) as tfi:

        linecount = 0
        hcount = 0

        hdict = make_headerdict()

        initer = iter(ifi)
        tempiter = iter(tfi)

            # move through the infile. write header (eventually); keep going if line is commented
        while True:
            try:
                inline = next(initer)
                # print(inline.decode())
            except StopIteration:
                inline = None
                break
            ll = inline.decode()
            if ll.startswith('##'): # standard format
                hdict = build_headerdict(hdict, ll, hcount)
                hcount += 1
                
            if not ll.startswith('#'):
                break
            pass # out.write(inline)

        # move through the template file. keep going if line commented or opener
        while True:
            try:
                refline = next(tempiter)
                #print(refline.decode())
                # print(refline.decode())
            except StopIteration:
                refline = None
                break

            rr = refline.decode()
            if rr.startswith('##'):
                hdict = build_headerdict(hdict, rr, hcount)
                hcount += 1
            if rr.startswith('#CHROM'):
                templine=rr
            if not refline.decode().startswith('#'):
                break

        # write headers... 
        #print(hdict)
        write_header(out, hdict)
        out.write(templine.encode()) # kind of poor code, but gets the line in... 

        # check files are same length
        while True:
            try:
                inline = next(initer)
            except StopIteration:
                inline = None
            try:
                refline = next(tempiter)
            except StopIteration:
                refline = None

            if (inline is None or refline is None) and not (inline is None and refline is None):
                exit("Non-commented sections of files are of different lengths!") # assumes same underlying file... 

            # terminate loop if file ends have been reached

            if inline is None and refline is None:
                break
            # if linemax and linecount > linemax:
            #     break

            # main task... first a simple test. 
            linecount += 1
            ii = inline.decode()
            rr = refline.decode()

            bgl_meta, bgl_genotypes = est.parse_full(ii, False) # not biallelic
            ref_meta, ref_genotypes = est.parse_full(rr, False)

            if not bgl_meta[:2] == ref_meta[:2]: exit("Site mismatch between BGL & REF")

            # going to work on alleles... 
            genos = [make_phased_null(geno) if nullgeno(geno) else merge_geno(geno, beag) for (geno, beag) in zip(ref_genotypes, bgl_genotypes)]
            newline = '\t'.join(ref_meta + genos) + '\n'
            out.write(newline.encode())

# need helper functions

def nullgeno(geno):
    return geno.startswith('.')

def get_gt(geno):
    return geno.split(sep=':')[0] # in future, should make more use of the format field!... 

def isflipped(bgl, ref): # will operate on the level of the geno
    bgl = get_gt(bgl)
    ref = get_gt(ref) # we are going to assume diploid and cheat a little... 
    return bgl[0] != ref[0]

def make_phased_null(geno):
    geno = geno.split(sep=":")
    gtnew = ['|'.join(geno[0].split(sep='/')), ] # should make list
    gtnew.extend(geno[1:])
    return ":".join(gtnew)

def merge_geno(ref, bgl):
    geno = ref.split(sep=":")[1:]
    beg = [get_gt(bgl), ]
    beg.extend(geno)
    return ":".join(beg)

def make_headerdict():
    return {'fileformat': dict(), 'filedate': dict(), 'source': dict(), 'reference': dict(),
        'INFO': dict(), 'FILTER': dict(), 'FORMAT':dict(), 'ALT': dict(), 'assembly': dict(), 'contig': dict(), 
        'SAMPLE': dict(), 'PEDIGREE': dict(), 'misc': dict()}

def build_headerdict(hdict, ll, n):
    lll = ll[2:-1] # trimming - this whole code needs a much better long-term solution. 
    lparse = lll.split('=<')
    if len(lparse) == 1: # i.e. the split failed & the format is irregular
        lparse = lll.split('=')
        cat = lparse[0]
        id = lparse[1]
        if not cat in hdict.keys():
            cat = 'misc'
        if not id in hdict[cat]:
            valdict = {'id': id, 'num': n, 'line': ll}
            hdict[cat][id] = valdict
            return hdict
        else:
            return hdict

    cat = lparse[0]
    if not cat in hdict.keys():
        cat = 'misc'
    val = lparse[1]
    id = val.split(',')[0].split('=')[1]
    if not id in hdict[cat]: # 
        valdict = {'id': id, 'num': n, 'line': ll}
        hdict[cat][id] = valdict
    return hdict

def lineprint(f, hdict, cat):
    for val in hdict[cat]:
        f.write(hdict[cat][val]['line'].encode())

def contigprint(f, hdict, starts='NC_'):
    for val in hdict['contig']:
        if hdict['contig'][val]['id'].startswith(starts):
            f.write(hdict['contig'][val]['line'].encode())

def write_header(f, hdict): # we will do this in stages... 
    lineprint(f, hdict, 'fileformat')
    lineprint(f, hdict, 'filedate')
    lineprint(f, hdict, 'source')
    lineprint(f, hdict, 'reference')
    lineprint(f, hdict, 'INFO')
    lineprint(f, hdict, 'FILTER')
    lineprint(f, hdict, 'FORMAT')
    lineprint(f, hdict, 'ALT')
    lineprint(f, hdict, 'assembly')
    contigprint(f, hdict)
    lineprint(f, hdict, 'SAMPLE')
    lineprint(f, hdict, 'PEDIGREE')
    lineprint(f, hdict, 'misc')



# to do later: 

# fileformat, filedate, source, INFO, FORMAT, ALT[NONREF], FILTER, FORMAT, contig (program lines)


