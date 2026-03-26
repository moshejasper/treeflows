import pyfaidx

def extract_region(infasta, outfasta, seqname, start, end, outwidth=40):
    '''Extracts region to new fasta file (inclusive). Note bp will be adjusted. need correction later. 
    
    VCF REFS HERE ARE 1-BASED AND INCLUSIVE'''

    # note. selection coords are 0-based, but then sequence names are 1-based (but not slicing)... 
    # start and stop will be 1-based here (per fasta type). And inclusive!!! Means we need conversion. 

    genome = pyfaidx.Fasta(infasta)

    with open(outfasta, 'w') as ofi:

        rawregion = genome[seqname][start-1:end]
        ofi.write(f'>{seqname}:{rawregion.start}-{rawregion.end}\n')
        chunks = len(rawregion) // outwidth # NOTE: CODE WAS BADLY WRITTEN AND IS NOW BENG REWRITTEN...
        for i in range(chunks): # this slicing is based on a string, not original sequence. 
            ofi.write(f'{rawregion[i*outwidth:(i+1)*outwidth]}\n')
        ofi.write(f'{rawregion[chunks*outwidth:]}')


def exclude_region(infasta, outfasta, blacklist, outwidth=40): # must have set this to 80 somewhere
    '''Extracts all regions but one to a new fasta file'''

    genome = pyfaidx.Fasta(infasta)

    with open(outfasta, 'w') as ofi: 

        for scaffold in genome.keys():
            if not scaffold in blacklist:
                scaf = genome[scaffold]
                ofi.write(f'>{scaf.name}\n')
                seglength = len(scaf)
                chunks = seglength // outwidth # NOTE: CODE HERE WAS BADLY WRITTEN AND HAS JUST BEEN FIXED... 
                for i in range(chunks):
                    ofi.write(f'{scaf[i*outwidth:(i+1)*outwidth]}\n')
                ofi.write(f'{scaf[chunks*outwidth:]}\n')



