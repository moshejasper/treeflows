from .vcf_core import VCFReader, VCFCompare
from pathlib import Path
import subprocess

""" Specification section of file

    Allele frequency file: tab-delineated & contains a header.
    Row: allele frequency for position in genome (rows ordered by increasging position along chromosome)
    Should be separate allele freq file for each chromosome
    ~~~
    position    x   n   folded [POS on chromosome (1-based?), allele count(0-n), folded? i.e. polarized i.e. known whether derived or ancestral... ]

    # Polarized? = 'known that the allele is derived or ancestral'... if polarized (opp of folded) 'folded' should be '0' (i.e. not folded) and 'x' should be 'derived count'... 
    # complex and will need to derive from rule-based dataset??? 
    *** if polarize and x=0, site is monomorphic and allele is identical to outgroup used to polarize site NOT RECOMMENDED TO USE x=0 ON POLARISED (Huber 2015). 
    *** if polarized and x=n, site is a substitution (monomorphic, different from outgroup used to polarize)... 0therwise, a polymorphis.m... 
    'only use polymorphisms & substitutions to scan for sweeps'... (all integers)... (can switch back and forth between folded and unfolded depending on certanity)... (e.g. fusion)... 
    
    Going to need to sort this out (comparatively?) for subgroups... (hmm)... 


    RECOMBINATION FILE - tab-delimited.. contains header... 
    Row: recombination rate (cM) between position in genome and prev pos in file... (orderd by increased position) - sep for each chrom... 
    Rate is rate up until this position... 

    position    rate
    2000    0.0
    2020    1.2

    BVALUE FILE (won't have)... (???) unsure how to generate... 

    USER-DEFINED GRID (positions at which test is calculated) - should be spanned by positions in freq input file... (just column of numbers, with no header)

"""

def vcf_to_freq_sweep_folded(vcf_file: Path|str, outname):
    """
    Function to convert vcf file into multi-file format. Will eventually adapt to somethign more extensible
    """

    vcf = VCFReader(vcf_file)

    with open(outname, "w") as ofi:
        ofi.write("position\tx\tn\tfolded\n")
        for site in vcf:

            # need site to be biallelic or alternatively monoallelic for derived allele (does this falsify my prior scripts???)... hmm... 
            # at present, no such thing as monallelic derived. probably this doesn't matter as we are comparative within species
            #   but could possibly have 'monoallelic derived' within sub-group. Need to consider this VERY CAREFULLY... (i.e. not prematurely eradicating derived)

            if site.is_biallelic: # big issues with calculation of dac... 
                dac = site.v_c[1] # assumes that ref is main, this is alt... possibly poor assumption

            elif site.is_monoallelic:
                dac = 0 if site.has_ref else site.an

            ofi.write(f"{site.POS}\t{dac}\t{site.an}\t1\n")

def vcf_to_freq_sweep_unfolded(vcf_file: Path|str, outname, use_folded=False):
    pass

def make_map_and_freq_interp(vcf_file: Path|str, inmap: Path|str, out_freq: Path|str, out_map: Path|str, chrval_map='X2', chrval_vcf='NC_035108.1'):
    """ Going to try and interpolate values for other file as well as run main... """
    vcf = VCFReader(vcf_file)
    vcfi = iter(vcf)

    with open(inmap) as imap, open(out_map, 'w') as omap, open(out_freq, 'w') as ofreq:
        ofreq.write("position\tx\tn\tfolded\n")
        omap.write("position\trate\n")
        genprev=0
        physprev=0
        phys=0
        posprev=0
        distrate_prev = 0
        site = vcfi.__next__() # setting up so can be referenced... means we will have an initial position... 
        rec_remainder = 0
        map_complete = False
        vcf_complete = False
        discard = imap.readline()
        while True:
            # main loop... we are going to get primary value of first chromosome... 

            while True:
                physprev = float(phys) if not phys=="phys" else 0
                try:
                    set, map, mkr, phys, gen = imap.readline().strip().split() # initial read... abort loop if end of file... BIG ISSUE PROBABLY... 
                except ValueError or EOFError:
                    map_complete = True
                    break

                if map != chrval_map: # rotates through irrelevant chromosomes
                    continue

                rate = (float(gen) - genprev) # probably zero for start... 
                dist = (float(phys) - physprev)
                distrate = rate / dist # should be multipliable by distance... 

                genprev = float(gen)

                if float(phys) > site.POS: # at this point, we are ready to pass to the other loop
                    break

                else: 
                    rec_remainder += rate
                    distrate_prev = distrate
                    continue
            
            if vcf_complete:
                break

            while True:
                print(site.POS)
                partrate = 0
                if posprev < physprev: # if prior vcf position is earlier than prior physmap position
                    # we need to capture the overlap
                    partrate = (physprev - posprev) * distrate_prev # the recombination amount needed
                    posprev = physprev
                
                if site.POS > float(phys): # genetic map has been overtaken.. in this case all bets are off... 
                    if not map_complete:
                        break
                    # once map complete, we pretty much just use its final value over and over... 

                recval = rec_remainder + partrate + ((site.POS - posprev) * distrate) if not map_complete else 0 # outro... 
                if site.is_biallelic: 
                    dac = site.v_c[1]
                elif site.is_monoallelic: 
                    dac = 0 if site.has_ref else site.an
                ofreq.write(f"{site.POS}\t{dac}\t{site.an}\t1\n")
                omap.write(f"{site.POS}\t{recval}\n")
                # reset values
                rec_remainder = 0
                posprev = site.POS
                try: 
                    site = vcfi.__next__()
                except StopIteration:
                    vcf_complete = True
                    break
            
            if map_complete and vcf_complete: 
                break
        
        print("finished program")




def make_recmap_sweep(inmap, mapfile, chrval='X2'):
    """Takes *very specific* file format (emailed to me from study) strip out chromosome 
     corresponding to chrval, and return sweep recomb. file format with recombination rate calculated. 
    
     Parameters: 
        inmap      -   (Path/string)   -   Path to original map file. which has space-separated format: 
                    - set   (not important)
                    - map   (chromosome column)
                    - mkr   (marker - not used)
                    - phys  (physical genome position)
                    - gen   (genetic map position) [format???]
        mapfile     -   (Path/string)   -   Path to output map file, with space-separated format:
                    - position   (physical genome position)
                    - rate (rate UP UNTIL this position (from previous) - in cM (of actual recombination, I think))
        chrval      -   (string)        -   selection value to match 'map' chromosomal name to be extracted from inmap
    
    Details:
        Transforms genetic map (type1) into Relate-compatible, including calculating recombination rates. 
         For recombination rates, calculated between current position [i] and next position [i+1] using the formula:
            RATE = [((GEN_POS[i+1] - GEN_POS[i])] 
        The beginning and end of the chromosome are set to have the same rates as the first and last calculated 
         regions respectively. 
     """
    genprev = 0
    physprev = 0
    rate=0
    with open(inmap) as imap, open(mapfile, 'w') as omap: 
        omap.write(f"position\trate")
        for x in imap:
            try:
                set, map, mkr, phys, gen = x.strip().split()
            except ValueError:
                break
            if map != chrval:
                continue
            rate = (float(gen) - genprev)
            omap.write(f"\n{phys}\t{rate}") # very important ordering to prevent fileread error (I think)
            physprev, genprev = int(phys), float(gen)

def sweep_get_spectrum(infile, outfile):
    subprocess.run(f"SweepFinder2 -f {infile} {outfile}", shell=True)

def run_sweepfinder2(freqfile, spectfile, recfile, outfile, windowsize=10_000):
    subprocess.run(f"SweepFinder2 -lrg {windowsize} {freqfile} {spectfile} {recfile} {outfile}", shell=True)

def run_sweepfinder2_nomap(freqfile, spectfile, outfile, windowsize=10_000):
    subprocess.run(f"SweepFinder2 -lg {windowsize} {freqfile} {spectfile} {outfile}", shell=True)


# to get a map running properly, need to basically interpolate the map values then lay them all down in a file. it will not be ideal at all. 