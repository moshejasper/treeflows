import subprocess
from pathlib import Path

path_to_gaia_r = Path("/data/gpfs/projects/punim1778/Projects/aegypti/2023/results/")

def run_gaia_nospace(treefile, popfile, filterthresh, outname):
    """Basic implementation of non-spatial location-based gaia (with 1:all distance matrix, etc.)"""

    subprocess.run(f"Rscript {path_to_gaia_r / "run_gaia_nospace.R"} {outname} {treefile} {popfile} {filterthresh}", shell=True)


def run_gaia_spatial(treefile, popfile, res, outname):
    """Basic implementation of spatial location-based gaia (stepping stone model).)"""
    
    subprocess.run(f"Rscript {path_to_gaia_r / "run_gaia_spatial.R"} {outname} {treefile} {popfile} {res}", shell=True)