from pathlib import Path
import yaml
import inspect
class FileConfig:

    def __init__(self, base=None, genomes=None, project_base=None, vcf_dir=None, tree_dir=None, adx_dir=None, binaries=None, refdir=None, gendbpath=None, 
                 genome_fa=None, map_raw=None, rec_rate=None, rec_rate_high=None, mut_rate=None, mismatch_high=None, mismatch_low=None, gen_time_years=None):
        """Centralize filesystem paths and default rates for the project."""
        self.base = Path("/data/gpfs/projects/punim1778") if base is None else Path(base)
        self.genomes = self.base / "genomes" if genomes is None else self.base / genomes
        self.project_base = self.base / "Projects/aegypti/2023" if project_base is None else self.base / project_base
        self.vcf_dir = self.project_base / "data/previous_studies/core_files" if vcf_dir is None else self.project_base / vcf_dir
        self.tree_dir = self.project_base / "results/tree_inference" if tree_dir is None else self.project_base / tree_dir
        self.adx_dir = self.project_base / "results/ngsadmix" if adx_dir is None else self.project_base / adx_dir
        self.binaries = self.base / "bin" if binaries is None else self.base / binaries
        self.refdir = self.base / "refdata" if refdir is None else self.base / refdir
        self.gendbpath = self.vcf_dir / "frag_gendb/aeg_gendb1_cn" if gendbpath is None else self.vcf_dir / gendbpath

        #resources
        self.genome_fa = self.genomes / "GCF_002204515.2_AaegL5.0_genomic.fa" if genome_fa is None else self.genomes / genome_fa

        #maps
        self.map_raw = self.genomes / "aegypti_linkage_physical_map.txt" if map_raw is None else self.genomes / map_raw

        # rates
        self.rec_rate = 4e-7 if rec_rate is None else rec_rate # based on weighted averages... 
        self.rec_rate_high = 4e-5 if rec_rate_high is None else rec_rate_high # based on non-weighted average
        self.mut_rate = 4.85e-09 if mut_rate is None else mut_rate
        self.mismatch_high = 1e-2 if mismatch_high is None else mismatch_high
        self.mismatch_low = 1e-4 if mismatch_low is None else mismatch_low
        self.gen_time_years = 0.067 if gen_time_years is None else gen_time_years # based on Rose 2023 Dating origin & spread # basically 15 gens per year

    @staticmethod
    def _normalize_yaml_path(path):
        """Ensure yaml/yml suffix; add .yml if missing; error on other suffixes."""
        cfg_path = Path(path)
        if cfg_path.suffix in {".yaml", ".yml"}:
            return cfg_path
        if cfg_path.suffix == "":
            return cfg_path.with_suffix(".yml")
        raise ValueError("Config file must have a .yaml or .yml extension")

    @classmethod
    def from_yaml(cls, path=None):
        """Create a FileConfig from a YAML file whose keys match __init__ args."""
        path = Path.cwd() / ".fileconfig.yml" if path is None else Path(path)
        if path.is_dir():
            path = path / ".fileconfig.yml"
        cfg_path = cls._normalize_yaml_path(path)
        if not cfg_path.exists():
            raise FileNotFoundError(cfg_path)
        with cfg_path.open() as fh:
            data = yaml.safe_load(fh) or {}
        if not isinstance(data, dict):
            raise ValueError("Config file must contain a mapping at the top level")

        # Only pass recognised __init__ parameters through to avoid surprises.
        params = inspect.signature(cls.__init__).parameters
        valid_keys = {name for name in params if name != "self"}
        kwargs = {k: v for k, v in data.items() if k in valid_keys}
        return cls(**kwargs)

    def to_yaml(self, path=None):
        """Write current __init__ parameters to a YAML file for later reloading."""
        path = Path.cwd() / ".fileconfig.yml" if path is None else path
        cfg_path = self._normalize_yaml_path(path)
        params = inspect.signature(self.__init__).parameters
        keys = [name for name in params if name != "self"]

        def _serialize(val):
            return str(val) if isinstance(val, Path) else val

        data = {k: _serialize(getattr(self, k)) for k in keys if hasattr(self, k)}
        with cfg_path.open("w") as fh:
            yaml.safe_dump(data, fh)

def get_recombination_rate_average(inmap, chr_name="X2"): # assume input format I expect. # very poor... 
    """Takes specific (non-relate) input format & a specifi chromosome (pre-determined), and 
     outputs average (non-unit-adjusted) recombination rate (input to make map relate simple)"""
    with open(inmap) as imap: 
        genprev =0
        physprev = 0
        rate=0
        count=0
        for line in imap: 
            try:
                _, map, _, phys, gen = line.strip().split()
            except ValueError:
                break
            if not map == chr_name:
                continue
            rate += (float(gen) - genprev) / (int(phys) - physprev)
            count += 1
            physprev, genprev = int(phys), float(gen)
    print(rate/count)
    return rate / count

def get_recombination_rate_weighted(inmap, chr_name="X2"):
    """ Takes specific (non-relate) input format and speficic chromosome and outputs average 
     recombination rate, but weighted by sequence length"""
    with open(inmap) as imap: 
        genprev =0
        physprev = 0
        rate=0
        count=0
        for line in imap: 
            try:
                _, map, _, phys, gen = line.strip().split()
            except ValueError:
                break
            if not map == chr_name:
                continue
            physdiff = int(phys) - physprev
            gendiff = float(gen) - genprev
            rate += gendiff / physdiff
            count += 1
            physprev, genprev = int(phys), float(gen)
    print(f"Non-weighted rate: {rate/count}")
    print(f"Weighted rate: {genprev / physprev}")
    return genprev / physprev

def get_recombination_rate_literature():
    """Currently just convenience function to hold value from other paper"""
    pass

def get_mutation_rate_literature():
    """Current mutation rate sourced from Rose et al. 2023, to be 4.89 *  10(-9)"""
    return 4.85e-09

def adjust_mutation_rate_missing(mut_rate, missing_frac):
    """Takes mutation rate and fraction missing data over interval, and returns adjusted mutation rate. Designed so that 
     diversity present between sites is appropriately accounted for. Works globally, but will not work locally where there 
     are clusters of missing data combined with clusters of present. For this, a future function to create a missing data-derived 
     mutation map will be applied
     Params: 
        mut_rate        (float) - non-adjusted mutation rate
        missing_frac    (float) - fraction data missing (0 = all present, 1 = all missing)
        
    Returns (float) adjusted mutation rate. 
    Note that this is designed to give appropriate performance for a program which assumes that all sites are present in the data. 
     if there have been other adjustments (e.g. of distance, as in Relate) then this measure will be less necessary.
     Forumula is current_rate * fraction_data_present. So 50% data present leads to 50% rate. 
     """
    new_rate = mut_rate * ( 1 - missing_frac )
    return new_rate

def get_chr_arm(chr_num, arm="left", rounding=1_000_000): # in future, adjust to deal with telomeres also
    """Return genomic coordinates for a chromosome arm.

    This is currently hard-coded for a small set of chromosomes and uses
    approximate centromere/arm boundaries, with optional rounding to convenient
    bin sizes.

    Args:
        chr_num: Integer index into the hard-coded chromosome list.
        arm: One of `"left"`, `"right"`, `"centromere"`.
        rounding: Round returned coordinates to the nearest multiple of this
            value (bp).

    Returns:
        Tuple `(start, end)` in bp coordinates.
    """
    arm_left = [(0, 150_000_000), (0, 227_000_000), (0, 196_000_000)]
    arm_right = [(154_000_000, 310_827_022), (232_000_000, 474_425_716), (201_000_000, 409_777_670)]
    if not arm in ["left", "right", "centromere"]:
        ValueError(f"Arm not one of 'left, right, centromere'!")
    if arm == "left": 
        rval = tuple((x // rounding * rounding for x in arm_left[chr_num]))
        if rval[0] == rval[1]:
            ValueError(f"Rounding value is too high!")
        return rval
    elif arm == "left":
        r_init = arm_right[chr_num]
        armstart = r_init[0] // rounding * rounding
        if armstart < r_init[0]: 
            armstart = ( r_init[0] // rounding + 1 ) * rounding
        armend = r_init[1] // rounding * rounding
        if armstart == armend:
            ValueError(f"Rounding value is too high!")
        return (armstart, armend)
    elif arm == "centromere":
        armstart_i = arm_left[chr_num][1]
        armend = arm_right[chr_num][0]
        armstart = armstart_i // rounding * rounding
        if armstart < armstart_i:
            armstart = ( armstart_i // rounding + 1 ) * rounding
        armend = armend // rounding * rounding
        if armstart == armend:
            ValueError(f"Rounding value is too high!")
        return (armstart, armend)

# it is a possibility that I could use a relate-style distance compression for tsinfer? Lots to look into... 

class SpartanConfig:

    def __init__(self, REFFILE, POPFILE, SUBSEQ, TRACKFILE, FILETAG):
        """Legacy config container for Spartan jobs and Aegypti pipeline filenames."""
        self.base = "/data/gpfs/projects/punim1778/"
        self.genome = self.base + "genomes/AaegL5mt" # need to move this across!!!
        self.genome_fa = self.base + "genomes/GCF_002204515.2_AaegL5.0_genomic.fa"
        self.binaries = self.base + "bin/"
        self.workdir = self.base + "Projects/aegypti/2023/data/previous_studies/core_files/" # for gendb access... (all currently on open filesystem) # previously core files... 
        self.datadir = self.workdir
        self.batchid =  "mj_init_may24" #"mj_bwa_feb25" #"mj_init_may24" #  "mj_aeg_masc_may30" # old mj_init_may24
        self.outdir = self.workdir # not currently used.
        self.touchdir = self.workdir + "touchdir/"
        self.reffile = self.datadir + REFFILE #  "aegmasc.txt" #formerly aeg_cn.txt
        self.popfile = self.datadir + POPFILE #"aeg_asiapac.pops"
        self.subreffile = self.datadir + SUBSEQ
        self.trackfile = self.datadir + TRACKFILE #  "aegmasc.track" #formerly aeg_cn3.track
        self.trackfile_alt = self.datadir + "bootreps_new.track"
        self.logname = self.datadir +  "aeg_wgs_bwa" #  'aegmasc_wgs_cn' # formerly aeg_wgs_cn
        self.failfile = self.datadir + f"{FILETAG}.fail" #  "aegmasc.fail"
        self.remotedir = "moshe@115.146.86.169:~/Projects/aegypti/2023/data/new_sequences/sparfiles/"
        self.gendbpath = self.datadir + "frag_gendb/aeg_gendb1_cn"
        self.gen_intvs = ("NC_035108.1", "NC_035109.1", "NC_035107.1")#, "NC_035108.1",)# "NC_035109.1") #, "NC_035159.1") #NC_035107.1
        self.gen_lengths = (474425716, 409777670, 310827022)#,  474425716,)# 409777670) #, 16790) #310827022
        self.gen_subintvs = ("NC_035108.1",) # CHR2
        self.gen_sublengths = (220_000_000,) # LEFT HALF (centromere 227-232 Mbp)
        self.int_sublength = 10_000_000 # was historically using 10,000,000 (this will match the various files)
        self.tempdir = self.base + "temp/"
        self.global_temp = "/tmp/"
        self.batchsize = 20
        self.multinode=False
        self.self_ref = self.base + "Projects/aegypti/2023/scripts/" + "aaa_spartan.py" #BIG TWEAK TO KEEP THINGS CLEAR CHANGE BACK LATER
        self.selff = "aaa_spartan"
        self.errorfile = self.datadir + "slurm_mj_0018.txt"
        self.exfolder = self.datadir + "temp_execs/"
        self.logfolder = self.datadir + "logs/"
        self.sl_temp = self.datadir + "temp_slurm.txt"
        self.py_temp = self.exfolder + "temp_py_"
        self.randpart = True # THIS OVERRIDES BELOW COMMANDS AND OPENS UP THREE SUBMISSION SOURCES!!!
        self.jobpartition = 'fos' # alt is cascade - and sapphire ... will fix this to enable more randomization in future... [based on profiling]
        # qosswitch = True
        # qos = 'fos'
        self.partitiondict = {'fos': 'fos', 'cascade': 'normal', 'sapphire': 'normal'}
        self.masc_fa = self.base + "genomes/mascloc_consensus.fa"
