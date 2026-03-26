import numpy as np
import subprocess 
from . import vcf_core as _vcf


# utils

def load_snplist(snpfile): 
    """Loads .snplist format (untitled, one integer line per cell, commented lines with #)
        Returns sorted unique int list of snps"""
    snps = set()
    with open(snpfile) as sfile:
        for ll in (l.strip() for l in sfile):
            if not ll or ll.startswith("#"):
                continue
            snps.add(int(ll))
    
    return(sorted(snps))

def load_samplelist(samplefile):
    """ Loads .samples file. Preserve order"""
    samplelist = []
    with open(samplefile) as sfile:
        for ll in (l.strip() for l in sfile):
            if not ll or ll.startswith("#"):
                continue

            sample = "_".join(ll.split('-')) # deals with idiosyncratic sample issue in my dataset
            samplelist += [sample]
    return samplelist


def vcf_filter_range(vcf_file, outvcf, min, max, maf=0.0, bi=False, low2miss=0):
    """Extract vcf sites from min to max (inclusive) with real maf of >= `maf` and optionally retaining biallelic sites only. 
    Currently only works for one chromosome in VCF"""

    vcf = _vcf.VCFReader(vcf_file)
    with _vcf.VCFWriter(outvcf, vcf) as vo:
        vo.write_header()

        counter = 0

        for site in vcf:
            if site.POS < min:
                continue
            if site.POS > max:
                break
            if low2miss and site.has_lowcount(low2miss):
                site = site.lowcount_to_miss(low2miss) # possible issues with underlying frame??? 

            # this effected as above is properly dealt with... 
            if site.allelicity > 1 and site.geno_2nd_frac > maf:
                if bi and site.allelicity > 2:
                    continue
                counter += 1
                vo.write_record(site)
    print(counter)

def vcf_filter_snp(vcf_file, outvcf, snplist: list):
    """Take vcf file (or list) & list of snp locations (int) and use to filter a vcf for matching sites"""
    if not type(vcf_file) == list:
        vcf_file = [vcf_file,]
    vcfi = _vcf.VCFReader(vcf_file[0])
    snps = set(snplist)
    with _vcf.VCFWriter(outvcf, vcfi) as vo:
        vo.write_header()
        counter = 0
        snps_total = len(snplist)
        for vv in vcf_file:
            snpmin = min(snps)
            snpmax = max(snps)
            scan = True
            print(f"CHecking {vv}")
            vcf = _vcf.VCFReader(vv)

            for site in vcf:
                if scan and site.POS < snpmin:
                    continue
                elif scan:
                    scan = False
                if site.POS > snpmax:
                    break
                if not site.POS in snps:
                    continue

                vo.write_record(site)
                print(site.vcf_alleles)
                conter += 1
                snps.remove(site.POS)

    print(f"missing snps: {snps_total - counter}")
    print(f"found snps:  {counter}")

def vcf_extract_genotype(vcf_file, loc, outfix, depth=3, chrom=None):
    vcf = _vcf.VCFReader(vcf_file)
    samples = vcf.samples
    total = len(samples)
    misscount, wildcount, hetcount, homcount = 0, 0, 0, 0
    with open(outfix + "_wt.txt", "w") as wt, open(outfix + "_het.txt", "w") as het, open(outfix + "_hom.txt", "w") as hom, open(outfix + "_fact.txt", "w") as fact:
        fact.write("sample\tlocweight\n")
        for site in (vcf if chrom is None else vcf(f"{chrom}:{loc}")):
            if site.POS < loc:
                continue
            if site.POS > loc:
                break
            genos = site.genotypes
            depths = site.depths
            for geno, sitedepth, sample in zip(genos, depths, samples):
                gs = sum(geno[:2])
                fact.write(f"{sample}\t{gs}\n") if sitedepth >= depth else fact.write(f"{sample}\t-2\n")
                if sitedepth < depth:
                    print(f"low depth geno: {sample}: {sitedepth}")
                    misscount += 1
                    continue
                if gs == 0: 
                    wt.write(sample + "\n")
                    wildcount += 1
                elif gs == 1:
                    het.write(sample + "\n")
                    hetcount += 1
                elif gs == 2:
                    hom.write(sample + "\n")
                    homcount += 1
                else:
                    print(f"bad geno: {sample}: {geno}")

    print(f"   WT: {wildcount}\n  HET: {hetcount}\n  HOM: {homcount}\t  MISS: {misscount}\nTOTAL: {total}")