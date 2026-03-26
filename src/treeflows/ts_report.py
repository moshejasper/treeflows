#!/bin/python

# probably run this as gendb?
import tskit
import numpy as np

from collections import Counter
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF
from pathlib import Path
import argparse
from functools import cached_property



def count_variables(variable_list, out=None):
    variable_counts = Counter(variable_list)
    if out is None:
        print('Roots:\t\tCount:\t')
    else:
        out.write('Roots:\tCount:\n')
    for variable, count in variable_counts.items():
        if out is None:
            print(f'{variable}\t\t{count}')
        else:
            out.write(f'{variable}\t{count:,}\n')

def get_root_distribs(ts, out=None):
    trees = [len(t.roots) for t in ts.trees()]
    count_variables(trees, out)

def assess_tree(ts): # rework as file attachment... 
    ntrees = ts.get_num_trees()
    sites = ts.get_num_sites()
    leng = ts.get_sequence_length()
    mutations = ts.get_num_mutations()
    print('----------------------')
    print(f'Number of Trees: {ntrees:g>20}')
    print(f'Number of Sites: {sites:g>20}')
    print(f'Sequence Length: {int(leng)}')
    print(f'Sites per Tree: {sites / ntrees:.2f}')
    print(f'Site density (/100bp): {(sites / leng) * 100:.2f}')
    print(f'bp per Tree: {leng / ntrees:.2f}')
    print('----------------------')
    get_root_distribs(ts)


class TreeStat:
    """Key class for putting information together about a tree sequence. Relies on very specific file structure. """

    def __init__(self, treefile):
        # we are going to assume we can get metadata from filepathway. 
        self.treepath = Path(treefile)
        tpath = self.treepath.resolve()
        components = tpath.parent.parts[tpath.parent.parts.index('tree_inference')+1:] # selects semantically relevant components
        self.ts = tskit.load(self.treepath)
        match components[0]:
            case 'est' | 'maj' | 'rule':
                self.ancestry = components[0]
            case _:
                self.ancestry = "undefined"
        match components[1]:
            case 'Aaa' | 'Aaa_in' | 'all' | 'apac' | 'apac_in':
                self.subset = components[1]
            case _: self.subset = "undefined" 
        match components[2]:
            case 'm10' | 'm25': self.missing = components[2]
            case _: self.missing = "undefined"
        match components[3]:
            case 'impute': self.imputed = True
            case 'miss': self.imputed = False
            case _: self.imputed = "undefined"
        match components[4]:
            case 'nodate' | 'date_rule' | 'rep_date': self.initdate = components[4]
            case _: self.initdate = "undefined"
        match components[5]:
            case 'relate' | 'singer' | 'tsinfer': self.program = components[5]
            case _: self.program = 'undefined'
        if len(components) == 6: self.params = "default"
        else: self.params = components[6]
    
    @property
    def num_trees(self):
        return self.ts.num_trees
    @property
    def num_edges(self):
        return self.ts.num_edges
    @property
    def num_nodes(self):
        return self.ts.num_nodes
    @property
    def num_sites(self):
        return self.ts.num_sites
    @cached_property
    def num_mutations_unique(self):
        return len(np.unique(self.ts.mutations_site))
    @property
    def num_mutations_total(self):
        return self.ts.num_mutations
    @property
    def num_individuals(self):
        return self.ts.num_individuals
    @property
    def num_samples(self):
        return self.ts.num_samples
    @property
    def sequence_length(self):
        return self.ts.sequence_length
    
    @property
    def frac_sites_present(self):
        return self.num_sites / self.sequence_length
    @property
    def mutation_redundancy(self):
        return self.num_mutations_total / self.num_mutations_unique
    
    @property
    def r_muts_tree(self):
        return self.num_mutations_total / self.num_trees # divide by mutation_reduncancy to get unique edges... 
    @property
    def r_muts_edge(self):
        return self.num_mutations_total / self.num_edges
    @property
    def r_edges_tree(self):
        return self.num_edges / self.num_trees
    
    @property
    def filesize_mb(self):
        return self.treepath.resolve().stat().st_size / (1024**2) # Megabytes
    
    @property
    def filesize_scaled(self):
        return self.treepath.resolve().stat().st_size / self.num_sites / self.num_individuals

    def __str__(self):
        ostr = f"Tree File: {self.treepath.stem:>26}\n" # string formatting - right-adjusted, 26... 
        ostr += f"Dataset: {self.subset:>28}\n"
        ostr += f"Ancestry: {self.ancestry:>22}\n"
        ostr += f"Missingness: {str(bool(self.missing)):>24}\n"
        ostr += f"Imputed: {self.imputed:28}\n"
        ostr += f"Initial Dating: {self.initdate:>21}\n"
        ostr += f"Inference Program: {self.program:>18}\n"
        ostr += f"Parameters: {self.params:>20}\n"
        return ostr
    
    def writelog(self, filename, append=False):
        mode = "a" if append else "w"
        with open(filename, mode) as logfile:
            if mode == "w":
                logline ="Treefile\tDataset\tAncestry\tMissingness\tImputed\tDate_Init\tProgram\tParams\tnum_trees"
                logline += "\tnum_edges\tnum_nodes\tnum_sites\tnum_mutations_total\tnum_mutations_unique\tnum_individuals"
                logline += "\tsum_samples\tsequence_length\tfrac_sites_present\tmutation_redundancy\tfilesize_mb"
                logline += "\tr_muts_tree\tr_muts_edge\tr_edges_tree\tfilesize_scaled"
                logfile.write(logline + "\n")
            logfile.write("\t".join(
                str(x) for x in [
                    self.treepath.stem, self.subset, self.ancestry, str(self.missing), str(bool(self.imputed)), 
                    self.initdate, self.program, self.params, self.num_trees, self.num_edges, self.num_nodes, self.num_sites, 
                    self.num_mutations_total, self.num_mutations_unique, self.num_individuals, self.num_samples, 
                    self.sequence_length, self.frac_sites_present, self.mutation_redundancy, self.filesize_mb, 
                    self.r_muts_tree, self.r_muts_edge, self.r_edges_tree, self.filesize_scaled]) + '\n')

def treestat_walker(root, logfile):
    root = Path(root).resolve()
    header = False
    with open(logfile, "w") as log:
        log.write("Treefile\tDataset\tAncestry\tMissingness\tImputed\tDate_Init\tProgram\tParams\tnum_trees\n")
    for t in root.rglob("*.trees"):
        print("--------")
        tree = TreeStat(t)
        if not header:
            tree.writelog(logfile, False)
            header = True
        else:
            tree.writelog(logfile, True)
        print(tree)

def assess_tree_report(ts, outfile, treefile=None): # rework as file attachment... 
    ntrees = ts.get_num_trees()
    sites = ts.get_num_sites()
    leng = ts.get_sequence_length()
    n_muts = len(np.unique(ts.mutations_site))
    with open(outfile, 'w') as ofi:
        ofi.write('----------------------\n')
        if treefile is not None:
            ofi.write(str(TreeStat(treefile)))
        ofi.write('----------------------\n')
        ofi.write(f'Number of Trees: {ntrees:>20,}\n')
        ofi.write(f'Number of Sites: {sites:>20,}\n')
        ofi.write(f'Number of Mut Sites: {n_muts:>16}\n')
        ofi.write(f'Number of Muts: {ts.num_mutations:>21,}\n')
        ofi.write(f'Number of Edges: {ts.num_edges:>20,}\n')
        ofi.write(f'Sequence Length: {int(leng):>20,}\n')
        ofi.write(f'Sites per Tree: {sites / ntrees:>21.2f}\n')
        ofi.write(f'Mutations per Tree: {ts.num_mutations / ntrees:>17.2f}\n')
        ofi.write(f'Mutations per Edge: {ts.num_mutations / ts.num_edges:17.2f}\n')
        ofi.write(f'Muts per Mut Site: {ts.num_mutations / n_muts:18.2f}\n')
        ofi.write(f'Site density (/100bp): {(sites / leng) * 100:>14.2f}\n')
        ofi.write(f'Mut Site density (/100bp): {(n_muts / leng) * 100:>10.2f}\n')
        ofi.write(f'bp per Tree: {leng / ntrees:>24.2f}\n')
        ofi.write('----------------------\n')
        get_root_distribs(ts, ofi)

def to_tree(ts, n):
    trees = [t for t in ts.trees()]
    return(trees[n])

def svg_to_pdf(svg, outfile, midfile = 'ueahfiaugh.svg'):
    midd = Path(midfile)
    with open(midd, "w") as mid:
        mid.write(svg)
    drawing = svg2rlg(midd)
    renderPDF.drawToFile(drawing, outfile + '.pdf')
    midd.unlink()

node_label_style = (
    '.node > .lab {"font-size": "10%"}'
    '.leaf > .lab {"text-anchor": "start"; "transform": "rotate(90deg)"}'
    '.sample > .lab {"transform": rotate(90deg)'
    '.site > .lab {"font-size": "10%"}'
)

def ts_to_pdf(ts, tname, pos=1, sized=(10000, 500)):
    tsweep = ts.at(pos)
    tsvg = tsweep.draw_svg(size=sized, y_gridlines = True, y_axis = True, style=node_label_style, x_axis=True, 
                           mutation_labels = {})
    svg_to_pdf(tsvg, tname)
    tsvg = tsweep.draw_svg(size=sized, y_gridlines = True, y_axis = True, style=node_label_style, x_axis=True, 
                           mutation_labels = {}, all_edge_mutations=True)
    svg_to_pdf(tsvg, tname + '_all')

def ts_make_report_file(tsfile, tname=None, pos=1):
    ts = tskit.load(tsfile + '.trees')
    if tname is None:
        tname = tsfile
    assess_tree(ts)
    ts_to_pdf(ts, tname, pos)

def ts_make_report_basic(tsfile):
    tsfile = Path(tsfile)
    ts = tskit.load(tsfile.with_suffix('.trees'))
    reportfile = tsfile.with_suffix('.ts_report')
    assess_tree_report(ts, reportfile, treefile=tsfile)


if __name__ == "__main__":

    treestat_walker(".", "test.txt")
    # intree = Path("./est/apac/m10/miss/nodate/tsinfer/inf/chrom2a_inf.trees")
    # mytree = TreeStat(intree)
    # print(mytree)
    # mytree.writelog("test.txt")
    # exit
    # parser = argparse.ArgumentParser("TS REPORT", "python ts_report ts_filename", "Generate basic tree sequence report")
    # parser.add_argument("filename", help="filename of tree sequence")
    # wd = Path(__file__).parent
    # args = parser.parse_args()
    # tstat = TreeStat(args.filename)
    # print(tstat)
    # ts_make_report_basic(args.filename)

    #ts_make_report_file('mj_init_may24_NC_035108.1:9451401-10451400_sub_d3m80c2_dummy_T_inf', 'report_test', pos=8000)