from pathlib import Path
import json
from treeflows.treestats import get_gnn_alignment
#from collections import OrderedDict

# concept is of a hierarchical grouping of factor-based file permutations, where things apply to them at levels... 
def load_factor_hierarchy(jsonpath=None, rootdir=None):
    """Load a `FactorHierarchy` from a `.factorhierarchy` JSON file.

    Args:
        jsonpath: Optional explicit path to the JSON hierarchy file. If None,
            defaults to `<rootdir>/../.factorhierarchy` when `rootdir` is given,
            otherwise `treeflows/.factorhierarchy` adjacent to this module.
        rootdir: Optional root directory of the hierarchy tree on disk.

    Returns:
        A `FactorHierarchy` instance.
    """
    if jsonpath is None:
        if rootdir is not None:
            jsonpath = Path(rootdir).parent / ".factorhierarchy"
        else:
            jsonpath = Path.cwd() / ".factorhierarchy"
    with open(jsonpath) as jfile:
        hierarchy = dict(json.load(jfile))
    return FactorHierarchy(rootdir=rootdir, hierarchy=hierarchy)

class FactorHierarchy:

    def __init__(self, rootdir=None, levels: list=None, hierarchy: dict=None):
        """Create a factor-based directory hierarchy manager.

        Args:
            rootdir: Root directory for the hierarchy tree on disk.
            levels: Ordered list of hierarchy level names.
            hierarchy: Mapping level name -> list of factor names.
        """
        self.rootdir = Path.cwd() if rootdir is None else Path(rootdir)
        self.rootdir = self.rootdir.resolve()
        self.hierarchy = {} if hierarchy is None else dict(hierarchy)
        self.levels = list(self.hierarchy.keys()) if levels is None else list(levels)
        self.depth = len(self.levels)

        # include functions to store a metadata object here, perhaps in JSON... (?)... 

    @classmethod
    def from_jsonfile(cls, jsonpath="./.factorhierarchy", rootdir: Path=None):
        """Load a `FactorHierarchy` from a JSON file."""
        jsonpath = Path(jsonpath)
        # one filestructure in place at a time. 
        with open(jsonpath) as jfile:
            hierarchy = dict(json.load(jfile))
        return FactorHierarchy(rootdir=rootdir, hierarchy = hierarchy)

    def initialize(self):
        """Write the current hierarchy to `<rootdir>/.factorhierarchy`."""
        with open(self.rootdir / ".factorhierarchy", "w") as jsonfile:
            json.dump(self.hierarchy, jsonfile)

    def save(self):
        """Alias for `initialize()` (persist hierarchy JSON)."""
        with open(self.rootdir / ".factorhierarchy", "w") as jsonfile:
            json.dump(self.hierarchy, jsonfile)

    def generate_hierarchy_full(self)->None:
        """Create the full directory tree for all factor combinations."""
        self.initialize()
        endpaths = self.get_paths_depth(self.depth)
        for path in endpaths:
            path.mkdir(parents=True) # throw error on fail... 

    def update_hierarchy_full(self)-> None: 
        """Ensure the full directory tree exists (mkdir with exist_ok)."""
        self.initialize()
        endpaths = self.get_paths_depth(self.depth)
        for path in endpaths:
            path.mkdir(parents=True, exist_ok = True)

    def get_paths_depth(self, depth: int)->list: # numeric
        """Return all directory paths at a given depth of the hierarchy."""
        paths = [self.rootdir]
        for n in range(depth):
            lfactors = self.hierarchy[self.levels[n]]
            newpaths = [path / lfactor for path in paths for lfactor in lfactors]
            paths = newpaths
        return paths
    
    def get_pathdirs(self, factordict: dict, truncate=False)->list:
        """Return directory paths matching a partial factor selection.

        Args:
            factordict: Mapping level -> selected factor(s) (str or list).
            truncate: If True, only build paths down to the deepest specified level.
        """
        levels, factors = list(factordict.keys()), list(factordict.values())
        if not any(isinstance(el, list) for el in factors):
            factors = [[el] for el in factors]
        depths = [self.levels.index(level) for level in levels]
        levels = [x for _, x in sorted(zip(depths, levels))]
        factors = [x for _, x in sorted(zip(depths, factors))]
        depths = sorted(depths)
        if len(depths) == 0: depths = [0]
        maxdepth = depths[-1]+1 if truncate else self.depth

        # big timewaste done!!! 
        paths = [self.rootdir]
        for n in range(maxdepth):
            lfactors = self.hierarchy[self.levels[n]]
            if n in depths: # this needs filtering
                filtfactors = factors[depths.index(n)]
                if not type(filtfactors) is list:
                    filtfactors = list(filtfactors)
                if not all([factor in lfactors for factor in filtfactors]):
                    print(filtfactors)
                    print(lfactors)
                    raise Exception("Improper Selection Factor!")
                lfactors = [factor for factor in lfactors if factor in filtfactors] # get rid of missing on either side. probably too tolerant... 
            newpaths = [path / lfactor for path in paths for lfactor in lfactors]
            paths = newpaths
        return paths
    
    def find_pattern(self, *args, **kwargs):
        """Find files under the hierarchy root matching patterns and factor constraints."""
        # args are patterns, kwargs are semantics. 
        if len(kwargs) == 0:
            return [f for arg in args for f in self.rootdir.rglob(f"*{arg}")]
        pathdirs = self.get_pathdirs(factordict = kwargs, truncate=True) # meaning we have higher level dicts. 
        return [f for path in pathdirs for arg in args for f in path.rglob(f"*{arg}")]
    
    def write_stats(self, statfile, StatClass: str = "FileStat", *args, **kwargs):
        """Finds pattern and extracts Stat class for all valid entries"""
        pathlist = self.find_pattern(*args, **kwargs)
        statfile = Path(statfile).resolve()
        with open(statfile, "w") as stat:
            globals()[StatClass](pathlist.pop(0), self).write(stat)
            for path in pathlist:
                globals()[StatClass](path, self).write(stat, header=False)
    
class FileStat:
    """Base class for extracting info from Hierarchy system"""
    def __init__(self, filename: str | Path, hierarchy: FactorHierarchy):
        """Parse a filename within a hierarchy and compute its factor dictionary."""
        self.treepath = hierarchy.rootdir
        self.hierarchy = hierarchy
        self.filepath = Path(filename).resolve().relative_to(self.treepath)
        self.filepath_full = Path(filename).resolve()
        self.factors = self.filepath.parent.parts
        self.factordict = {}
        for i, val in enumerate(self.filepath.parent.parts):
            key = self.hierarchy.levels[i]
            if val in self.hierarchy.hierarchy[key]:
                self.factordict[key] = val
            else:
                self.factordict[key] = "NA"
        for key in self.hierarchy.hierarchy.keys():
            if not key in self.factordict.keys():
                self.factordict[key] = "NA"

    def __str__(self):
        """Tab-delimited representation of factor keys and values."""
        keys = [key for key in self.factordict.keys()]
        vals = [val for val in self.factordict.values()]
        out = "\t".join(keys) + '\n' + "\t".join(vals)
        return out
    
    def write(self, filehandle, header=True):
        """Write factor values (and optionally header keys) to an open filehandle."""
        if header:
            keys = "\t".join([str(key) for key in self.factordict.keys()])
            filehandle.write(keys + "\n")
        vals = "\t".join([str(val) for val in self.factordict.values()])
        filehandle.write(vals + "\n")

    
class GnnStat(FileStat):

    def __init__(self, *args, **kwargs):
        """Extend `FileStat` by parsing GNN filename-encoded parameters."""
        super().__init__(*args, **kwargs)
        stem = self.filepath.stem
        stemparts = stem.split("_") # now we need a meaningfull stemparts dict... (non-exclusive)... i.e. semantic rules... 
        self.factordict["gnn_type"] = stemparts[0]
        valid = True if self.factordict["gnn_type"] == "adx" else False 
        self.factordict["gnn_type"] = stemparts[0] if valid else "NA"
        self.factordict["gnn_k"] = int(stemparts[1][1:]) if valid else "NA"
        self.factordict["gnn_thresh"] = float("." + stemparts[2][1:]) if valid else "NA"
        self.factordict["gnn_nodes"] = stemparts[3][1:] if valid else "NA"
        self.factordict["gnn_has_ind_pairs"] = stemparts[4] if valid else "NA"
        self.factordict["gnn_has_pop_pairs"] = stemparts[5] if valid else "NA"
        self.factordict["gnn_is_randomized"] = stemparts[6] if valid else "NA"
        self.factordict["gnn_discrim"] = stemparts[7] == "dd" if valid else "NA"
        self.factordict["gnn_num_neighbours"] = int(stemparts[8][1:]) if valid else "NA"
        self.factordict["gnn_adx_alignment"] = float(get_gnn_alignment(self.filepath_full, self.factordict["gnn_k"])) if valid else "NA"

        """
        Rules for the string classes. 
        Things we will vary:
            0.  adx Whether this is admixture-based GNN
            1.  t9  Admixture threshold (0.9? 0.75?)
            2.  a5  Num individuals chosen per ancestral group (max) (5? all? 1?)
            3.  xi  Whether to exclude individaul 
            4.  xp  Whether to exclude idvs from same subpopulation
            5.  xr  Whether we are randomizing individuals between pop groupings (that were preselected)
            6.  rr  Whether sequence is relate or tsinfer
            7.  k5  How many K (determined by admixture)
            8.  dd  Whether discriminant filter used
            9.  n1  How many neighbouring nodes to match
        
        E.g. adx_k5_t9_a5_xi_xp_xr_dd_n1.gnn
        """
    # def _calculate_alignment(self, )

class TreeStat(FileStat):
    pass



if __name__ == "__main__":
    # test time... 

    mylevs = ["misslevel", "imputation", "program"]
    hierarchy = {
        "ancestry": ["est", "maj", "rule"], 
        "subpop": ["apac", "Aaa", "all"], 
        "misslevel": ["m10", "m25"], # could all be JSON
        "imputation": ["miss", "impute"],
        "dating": ["nodate", "date_rule", "rep_date"], 
        "program": ["tsinfer", "relate", "singer"], 
        "param": ["inf", "mme2", "mme4"],
    }
    my_hr = FactorHierarchy(hierarchy = hierarchy)
    #my_hr = load_factor_hierarchy()
    my_hr.initialize()
    #my_hr.generate_hierarchy_full()
    paths = my_hr.find_pattern(".trees", program="tsinfer")
    print("This is a test")
    print(paths)
    for path in paths:
        print(type(path))
        print(path)
        print(path.exists())
