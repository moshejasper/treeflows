"""Treestats will assume tree sequences as inputs and work to implement the python half of general-
 purpose tree programs. It will need to take in some basic things like population clustering, etc.
 
 Will often interface with a set of R packages that will run the other half of the functionality. This 
  is a split-language situation... but oh well... Next goal is to get more python-integrated"""


import tskit
import pandas as pd
import numpy as np
from random import Random
from pathlib import Path
from treeflows.refdata import AegData, load_aegdata
from treeflows.config import FileConfig
import treeflows.admix as adx
import scipy
import matplotlib.pyplot as plt
import seaborn as sns
#from textwrap import dedent

def run_gnn(ts: tskit.TreeSequence, threads) -> np.ndarray:
    """Will need feeder functions to identify my various lists. Probably re-calculate these things 
    based on other programs (e.g. meaning that I have same focal, different sample sets, or whatever)"""
    focal = [int(s) for s in ts.samples()]
    sample_sets = [[s] for s in focal] # need good functionality for getting these and mapping to tree forms... 
    gnn = ts.genealogical_nearest_neighbours(focal, sample_sets, num_threads=threads)
    return gnn


def run_gnn_idv_idv(ts: tskit.TreeSequence, threads=1) -> np.ndarray:
    """every individual vs every individual... most basic"""
    focal = [int(s) for s in ts.samples()]
    sample_sets = [[int(s)] for s in ts.samples()]
    gnn = ts.genealogical_nearest_neighbours(focal, sample_sets, num_threads=threads)
    return gnn

def run_gnn_idv_site(ts: tskit.TreeSequence, data: AegData, focal: str|list=None, sitelist: list|str=None, threads: int=1, discrim=False, num_neighbours=1) -> np.ndarray:
    """every individual vs every site (home inclusive)... 
        Note: can make trios via picking only two contrast sites... & setting discrim & group settings
    """
    focal = data.get_nodes_idv(focal)
    sample_sets = data.get_group_nodes_pop(sitelist)
    if discrim:
        gnn = ts.genealogical_nearest_neighbours_advanced(focal, sample_sets, num_threads=threads, discrim=discrim, num_neighbours=num_neighbours)
    else:
        gnn = ts.genealogical_nearest_neighbours(focal, sample_sets, num_threads=threads)
    return gnn

def run_gnn_idv_country(ts: tskit.TreeSequence, data: AegData, focal: str|list=None, countrylist: list|str=None, threads: int=1) -> np.ndarray:
    """every individual vs every site (home inclusive)... """
    focal = data.get_nodes_idv(focal)
    sample_sets = data.get_group_nodes_country(countrylist)
    gnn = ts.genealogical_nearest_neighbours(focal, sample_sets, num_threads=threads)
    return gnn

def run_gnn_idv_continent(ts: tskit.TreeSequence, data: AegData, focal: str|list=None, continentlist: list|str=None, threads: int=1) -> np.ndarray:
    """every individual vs every site (home inclusive)... """
    focal = data.get_nodes_idv(focal)
    sample_sets = data.get_group_nodes_continent(continentlist)
    gnn = ts.genealogical_nearest_neighbours(focal, sample_sets, num_threads=threads)
    return gnn

def run_gnn_site_site(ts: tskit.TreeSequence) -> np.ndarray:
    """one idv from site vs every other site (one idv)... """
    pass

def run_gnn_idv_ancestry(ts: tskit.TreeSequence, data: AegData, focal: str|list=None, ancestrylist: list=None, threads: int=1, discrim=False) -> np.ndarray:
    """every individual vs panama , america, Aaa... """
    focal = data.get_nodes_idv(focal)
    sample_sets = data.get_group_nodes_idx(ancestrylist)
    if discrim:
        gnn = ts.genealogical_nearest_neighbours_discrim(focal, sample_sets, num_threads=threads)
    else:
        gnn = ts.genealogical_nearest_neighbours(focal, sample_sets, num_threads=threads)
    return gnn


FACTORDICT = {

}

def get_tree_bits(tree: tskit.Tree):
    nsamp=tree.num_samples()

    n_edges = tree.num_edges
    #n_edges_max = (tree.num_samples - 1) * 2 # max for bifurcating... I believe... 
    info_bits = n_edges - nsamp
    maxbits = nsamp - 2
    bitratio = info_bits / maxbits # 
    return info_bits, bitratio

def get_ts_bits_quantile(ts: tskit.TreeSequence):
    # probably have to traverse this in a specific order? ... 
    l_infobits = []
    l_bitratio = []
    for tree in ts.trees():
        infobits,  bitratio = get_tree_bits(tree)
        l_infobits.append(infobits)
        l_bitratio.append(bitratio)
    return np.quantile(l_infobits, [0.025, 0.1, 0.5, 0.9, 0.975])

def get_ts_bits_all(ts: tskit.TreeSequence):
    # probably have to traverse this in a specific order? ... 
    l_infobits = []
    #l_bitratio = []
    for tree in ts.trees():
        infobits,  bitratio = get_tree_bits(tree)
        l_infobits.append(infobits)
        #l_bitratio.append(bitratio)
    return np.array(l_infobits)

def get_ts_mask_quant(ts: tskit.TreeSequence, quant, rev=False):
    '''The quantile (0<q<1) below which results are discarded'''
    bits = get_ts_bits_all(ts)
    thresh = np.quantile(bits, quant)
    if rev:
        return np.array(bits <= thresh, dtype=int)
    return np.array(bits >= thresh, dtype=int)

def get_ts_mask_thresh(ts: tskit.TreeSequence, thresh, rev=False):
    '''The quantile (0<q<1) below which results are discarded'''
    bits = get_ts_bits_all(ts)
    if rev:
        return np.array(bits <= thresh, dtype=int)
    return np.array(bits >= thresh, dtype=int)

def get_ts_nodes_treemask(ts: tskit.TreeSequence, treemask, rev=False):
    '''Index of all nodes for tree'''
    # from hyan wong June 24, 2025
    used_nodes = np.zeros(ts.num_nodes)
    for i, tree in ts.trees():
        if treemask[i] ^ rev:
            used_nodes[tree.preorder()] = True
    return used_nodes

def ts_simplify_pops(ts: tskit.TreeSequence|Path|str, sitelist, dataname: AegData|Path|str, tsout = None):
    config = FileConfig()
    if not type(dataname) == AegData:
        data = load_aegdata(config.refdir / dataname)
    filternodes = data.get_nodes_pop(sitelist)
    if not type(ts) == tskit.TreeSequence:
        ts=tskit.load(str(ts))
    ts_simp = ts.simplify(filternodes)
    if tsout is not None:
        ts_simp.dump(str(tsout))
    return ts_simp

def ts_simplify_levnames(ts: tskit.TreeSequence|Path|str, namelist, dataname: AegData|Path|str, tsout = None, level="pop_short"):
    config = FileConfig()
    if not type(dataname) == AegData:
        data = load_aegdata(config.refdir / dataname)
    filternodes = data.get_nodes(level, namelist)
    if not type(ts) == tskit.TreeSequence:
        ts=tskit.load(str(ts))
    ts_simp = ts.simplify(filternodes)
    if tsout is not None:
        ts_simp.dump(str(tsout))
    return ts_simp
    

def treemask_to_intervalmask(treemask: np.array, ts: tskit.TreeSequence, rev=False):
    """For conversion to interval masks in NSPOPE coaldecoder. 

        Note that their masks are by default for `excluding` info, ours are for `including` it. 
        This means that if rev=False (devfauls) sites from trees marked `0` will be included in the returned exclusion interval. 
        If we reverse, sites marked `1` will be included. 
        
    Return format is a numpy array (shape=(n, 2)) of [start stop) columns of sites to drop. These will correspond to tree bounds"""
    treemask = treemask.astype(bool) # this currently fits the reverse case. 
    if not rev: treemask = ~treemask # flip here as our semantics are opposite
    return(np.array([tuple(tree.interval) for (m, tree) in zip(treemask, tskit.trees()) if m], dtype=int))


def write_treemask(fhandle, treemask):
    np.savetxt(fhandle, treemask, newline='', fmt='%i')

def write_treebits(fhandle, treebits):
    np.savetxt(fhandle, treebits, fmt='%i') # for R... and graphing... 
    fhandle.write('\n')

def run_gnn_factors(ts: tskit.TreeSequence, config: FileConfig, data: AegData, k=2, thresh=0.9, nodes=5, exclude_individuals = False, 
                    exclude_poppairs = False, randomize=False, discrim=False, num_neighbours=1, makename=False):
    """ This one will pair with an iterator and run different kinds of gnn analyses all at once... """
    anclist = adx.adx_get_idxs_all(config.adx_dir / f"chr2a_apac_m10_{k}.qopt", thresh=thresh)
    focal_nodes = data.get_nodes_idv(None) # better here in future... 
    anclist = adx.adx_idx_downsample(anclist, nodes) if nodes is not None and not nodes=='a' else anclist
    if randomize: 
        anclist = adx.adx_idx_permute(anclist)
    anclist = data.get_group_nodes_idx(anclist)

    gnn = ts.genealogical_nearest_neighbours_advanced(focal_nodes, anclist, discrim = discrim, num_neighbours = num_neighbours)
    gnn_name = f"adx_k{k}_t{str(thresh)[2:]}_a{'a' if nodes is None else nodes}_{'x' if exclude_individuals else 'i'}i"
    gnn_name += f"_{'x' if exclude_poppairs else 'p'}p_{'r' if randomize else 'x'}r_{'d' if discrim else 'x'}d_n{num_neighbours}.gnn"
    print(gnn_name)
    if makename:
        return gnn, gnn_name
    else:
        return gnn
    
def run_gnn_factors_all(treefile, dataname: str|Path|AegData = "asiapac_ref.csv", kmin=2, kmax=20, thresh=0.9, nodenums=[2, 5, 10, 'a'], 
                        discrims = [True, False], neighbourlist = [1, 2, 4, 8], randomize=[False]):
    treefile = Path(treefile).resolve()
    workdir = treefile.parent
    config = FileConfig()
    if not type(dataname) == AegData:
        data = load_aegdata(config.refdir / dataname)
    else: data = dataname
    ts = tskit.load(str(treefile))

    for k in range(kmin, kmax+1):
        for discrim in discrims:
            for nodenum in nodenums:
                for rr in randomize:
                    for nn in neighbourlist:
                        ## skip step if nn >= nodenum
                        if type(nodenum) == int and nn >= nodenum:
                            continue
                        gnn, gnn_name = run_gnn_factors(ts, config, data, k=k, thresh=thresh, nodes=nodenum, 
                                                        discrim=discrim, num_neighbours=nn, randomize=rr, makename=True)
                        gnn_full = workdir / gnn_name
                        np.savetxt(gnn_full, gnn)


    # run_gnn_idv_ancestry... 

def get_gnn_alignment(gnn_file, k):
    gnn_file = Path(gnn_file).resolve()
    config = FileConfig()
    adx = pd.read_csv(config.adx_dir / f"chr2a_apac_m10_{k}.qopt", sep=" ", header=None)[[n for n in range(k)]]
    gnn = pd.read_csv(gnn_file, sep=" ", header=None)

    retval = []
    for d in range(adx.shape[0]): # number of rows
        al0 = 0
        al1 = 0
        for e in range(adx.shape[1]): # should be cols
            val0 = min([adx.iloc[d, e], gnn.iloc[d*2, e]]) # nodes... 
            al0 = al0 + val0
            val1 = min([adx.iloc[d, e], gnn.iloc[d*2+1, e]])
            al1 = al1 + val1
        align = np.mean([al0, al1]) # contains the alignment for an individual assuming all included... 
        retval.append(align)
    # at this point shouldhave list of alignmentsf or each individual
    return np.mean(retval)


def get_fst_paired(treefile, dataname: str|Path|AegData):
    treefile = Path(treefile).resolve()
    #workdir = treefile.parent
    config = FileConfig()
    if not type(dataname) == AegData:
        data = load_aegdata(config.refdir / dataname)
    else: data = dataname
    ts = tskit.load(str(treefile))
    pops = data.get_group_nodes("pop_short")
    #pops = data.get_group_nodes_pop() # this is the sample_sets list of lists... 
    print(pops)
    # get k tuples (all comparisons)
    indexes = make_ktuple_grid(len(pops))
    print(indexes)
    #fst = ts.Fst(pops, indexes=indexes)
    #np.save(treefile.with_suffix('fst'), fst)

def make_ktuple_grid(datalen):
    ktuple = []
    for n in range(datalen):
        for m in range(datalen):
            ktuple.append([n, m])
    return ktuple

def get_f3_stat(ts: tskit.TreeSequence, dataname: str|Path|AegData, p_a, p_b, p_c):
    """Implements F3 statistics on a population"""

    config = FileConfig()
    if not type(dataname) == AegData:
        data = load_aegdata(config.refdir / dataname)
    else: data = dataname
    sampsets = data.get_group_nodes("pop_short", [p_a, p_b, p_c])
    return ts.f3(sampsets)


# def run_gnn_idv_ancestry_level(ts: tskit.TreeSequence, data: AegData, focal: str|list=None, level: int=2, threads: int=1) -> np.ndarray:


### functions for matching to core individual... assume they are in some kind of central location... 



#### Now: need options for post-processing the GNN results in a meaninful way... preferably with the RefData that helped produce it
### Figures class ####

class Figure(object): # borrowed from Kelleher
    """
    Superclass of figures for the paper. Each figure is a concrete subclass.
    """
    name = None

    def __init__(self):
        datafile_name = "data/{}.csv".format(self.name)
        self.data = pd.read_csv(datafile_name)

    def save(self, figure_name=None, bbox_inches="tight"):
        if figure_name is None:
            figure_name = self.name
        print("Saving figure '{}'".format(figure_name))
        plt.savefig("figures/{}.pdf".format(figure_name), bbox_inches='tight', dpi=400)
        plt.savefig("figures/{}.png".format(figure_name), bbox_inches='tight', dpi=400)
        plt.close()

    def error_label(self, error, label_for_no_error = "No genotyping error"):
        """
        Make a nice label for an error parameter
        """
        try:
            error = float(error)
            return "Error rate = {}".format(error) if error else label_for_no_error
        except (ValueError, TypeError):
            try: # make a simplified label
                if "Empirical" in error:
                    error = "With genotyping"
            except:
                pass
            return "{} error".format(error) if error else label_for_no_error

class GnnFigure(Figure):

    name = "clustermap"
    def __init__(self):
        pass
        
    def plot_clustermap(self, gnn_file):
        df = pd.read_csv(gnn_file) #.set_index("c") # sets index to particular column that contains necessary data... 

        # Z score normalize
        for col in list(df): # note we are borrowing code here from Kelleher 2019 (Fig 4)
            df[col] = scipy.stats.zscore(df[col])

        row_linkage = scipy.cluster.hierarchy.linkage(df, method="average", optimal_ordering=True)
        # rotate linkage(row_linkate, -x)
        order = scipy.cluster.hierarchy.leaves_list(row_linkage)
        x_pop = df.index.values[order] # note that this might assume a specific index file in advance... 
        
        cg = sns.clustermap(df[x_pop], row_linkage=row_linkage, col_cluster=False, rasterized=True)
        cg.ax_heatmap.set_ylabel("")
        for tick in cg.ax_heatmap.get_xticklabels():
            tick.set_rotation(-45)
            tick.set_ha('lef')
            tick.set_rotation_mode("ancor")
        self.save("gnn_clustermap")

class GnnStructureFigure(Figure):
    name = "gnn_structure"

    def __init__(self):
        pass

    def plot_clustermap(self, df, pop_colours, region_colours, figsize=(10, 10)):

        dfg = df.groupby("pop").mean()
        #Zscore normalize
        for col in list(dfg):
            dfg[col] = scipy.stats.zscore(dfg[col])
        row_linkage = scipy.cluster.hierarchy.linkage(dfg, method="average")
        order = scipy.cluster.hierarchy.leaves_list(row_linkage)
        x_pop = dfg.index.values[order]

        colours = pd.Series(pop_colours)
        cg = sns.clustermap(
            dfg[x_pop], row_linkage=row_linkage, col_cluster=False, 
            row_colors=colours, figsize=figsize, rasterized=True)
        cg.ax_heatmap.set_ylabel("")

        for region, col in region_colours.items():
            cg.ax_col_dendrogram.bar(0, 0, color=col, label=region, linewidth=0)
        return cg


"""
Ingredients. 
    (1) I need to 'know' the order of individuals (idvs & nodes) in the tree sequences, and be able to access them. 
    (2) I need a hierarchy of population levels (pandas, etc.) that I can use to subselect elements
    (3) I need inputs that can take other population info and add it here. 
    (4) At core, I need an integrated data class? (built on pandas) that can extract core population subspecie etc. info and 
        apply them to tree sequences, etc. etc. We are relying too heavily on R. 

        Plan: R -> .csv -> pandas DataFrame -> Python class for dealing with things (assumes standard output)... 
        Need a central location (like 'genomes') that data can be drawn from (raw) for processign these things. 
        Need bespoke ways of classifying and splitting this data also. Combined from mine & other format. Interchange with 
        R. 
"""

####### DISTANCE METRICS  - Random Forests #######

"""
Plan is to take something like 1,000 random samples along the genome (points) & then compare trees at those points for 
similarity with each other. I can probably put all this together into some kind of distance matrix and graph them. 
"""

def kc_distance_samples(ts1: tskit.TreeSequence, ts2: tskit.TreeSequence, n: int=1000, left: int=None, right: int=None)-> float:
    outvals = np.zeros(n)
    print("simplifying")
    ts1 = ts1.simplify() # need to do this if things aren't already on track... provenance check? 
    print("and more")
    ts2 = ts2.simplify()
    if left is None:
        left = int(max(ts1.edges_left[0], ts2.edges_left[0]))
        print(f"LEFT BOUND: {left}")
    if right is None:
        right = int(min(ts1.edges_right[-1], ts2.edges_right[-1]))
        print(f"RIGHT BOUND: {right}")

    rng = Random() # instantiating class that doesn't share state

    # now assume left & right are appropriate distances... 
    for i in range(n):
        loc = rng.randrange(left, right) # both need to be inclusive
        loc_kc = ts1.at(loc, sample_lists=True).kc_distance(ts2.at(loc, sample_lists=True))
        print(f"POS: {loc}\tKC: {loc_kc}")
        outvals[i] = loc_kc # only uses topology... # assumes same underlying nodes... 
    return outvals.mean()



### STATISTISCS

# allele frequency spectrum
#  mean gene divergence (within + between nodes)
# diversity (pi)
# F2 (between sets of nodes)
# F3, F4, Fst (windowed)... 
# genetic relatedness... (!!!)... 
# general stat (windowed, weights.. (?))... 
# mean descendents (per node?)
# Tajima's D
# trait correlation (to regions... )
# Y2 statistic (pair of sets of nodes)
# Y3 statistic (triples of sets of nodes)

if __name__ == "__main__":
    from pathlib import Path
    # compare two ftsinfer forms... imp & not imp... 

    cdir = Path(__file__).parent
    treefile1 = cdir / "miss/nodate/tsinfer/inf/chrom2a_inf.trees"
    treefile2 = cdir / "impute/nodate/tsinfer/inf/chrom2a_apacest_i10_inf.trees"

    ts1 = tskit.load(treefile1)
    ts2 = tskit.load(treefile2)
    result = kc_distance_samples(ts1, ts2) # probably async-compatible... (?)... 
    print(f"---\nFINAL KC: {result}")

