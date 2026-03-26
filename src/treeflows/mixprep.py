#!/bin/python3

from . import vcf2est as est
import gzip
import cyvcf2
from itertools import compress
#import os
#import multiprocessing
#import numpy as np

# we are going to use cyvcf2 here... 

# really need to fix up this whole ecosystem in the long run... 

def test_vcf(filename):
    vcf = cyvcf2.VCF(filename)
    for n, gens in enumerate(vcf):
        if n > 10:
            break
        print(gens.genotypes)


def simplify_genos(genos):
    genos_new = [x[0] for x in genos] # only good for crude sums... 
    genos_new.extend([x[1] for x in genos])
    return(genos_new)

# now need to produce the full file

def make_freqfile(vcffix, outfix, popfile, thresh=3):
    """ New version includes threshold!!! """
    vcf = cyvcf2.VCF(vcffix + '.vcf.gz')
    popdict = est.get_popdict(popfile)
    pops1 = popdict.keys()
    pops = [p for p in pops1 if len(popdict[p]) >= thresh] # should quietly drop smaller populations... (I think)... 
    popmasks = est.get_popmasks_internal(vcf, popdict)
    with gzip.open(outfix + '.frq.gz', 'wb') as ofi:
        headline = ' '.join(pops) + '\n'
        ofi.write(headline.encode())
        for n, site in enumerate(vcf):
            # if n >= 10:  
            #     break
            gens = site.genotypes
            freqlist = list()
            for pop in pops:
                popidvs = list(compress(gens, popmasks[pop])) # select the subset
                popgenos = simplify_genos(popidvs) # flatten for rest of stats... 
                maj = len([x for x in popgenos if x == 0])
                min = len([x for x in popgenos if x == 1])
                freq = ','.join([str(maj), str(min)])
                freqlist.append(freq)
            freqline = ' '.join(freqlist) + '\n'
            ofi.write(freqline.encode())
    print(f'{n} sites frequencies added for {len(pops)} pops')
            

# fname = 'mj_init_may24_NC_035109.1:100000001-110000000_missfilter'
# popfile = 'aeg_aeg.pops'

# make_freqfile(fname, 'tester', popfile)