from treeflows import vcf_core as vcf
from pathlib import Path
import numpy as np


def vcf_to_geno(in_vcf: str|Path|vcf.VCFReader, out_geno: str|Path, snps = True):
    """Replicates the geno conversion method required for the LEA package"""

    if not type(in_vcf) == vcf.VCFReader:
        in_vcf = vcf.VCFReader(in_vcf)

    # just going to go to town... assume that they are snps... 

    with open(Path(out_geno).with_suffix(".geno"), "w") as ofi:

        for site in in_vcf:

            if snps and site.is_nonsnp:
                continue
            
            # try to write the format... 

            genos = np.array(site.genotypes)[:,:2].sum(axis=1)
            genos[genos < 0] = 9
            np.savetxt(ofi, genos.reshape(1, -1), fmt='%d', delimiter="")