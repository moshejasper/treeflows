#!/bin/python3

import tskit
from collections import Counter
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF
from pathlib import Path

treeseq = 'loc1b3_filt75_ancestral_T_inf'


def count_variables(variable_list):
    """Print counts of items in `variable_list`."""
    variable_counts = Counter(variable_list)
    print('Roots:\t\t\tCount:\t')
    for variable, count in variable_counts.items():
        print(f'{variable}\t\t\t{count}')

def get_root_distrib(ts: tskit.TreeSequence)->None:
    """Count and print the distribution of number of roots per tree."""
    trees = [len(t.roots) for t in ts.trees()]
    count_variables(trees)


def assess_tree(ts: tskit.TreeSequence)->None:
    """Print a small summary of a tree sequence and its root distribution."""
    ntrees = ts.num_trees
    sites = ts.num_sites
    leng = ts.sequence_length
    print('----------------------')
    print(f'Number of Trees: {ntrees}')
    print(f'Number of Sites: {sites}')
    print(f'Sequence Length: {leng}')
    print(f'Sites per Tree: {sites / ntrees}')
    print(f'Site density (/100bp): {(sites / leng) * 100}')
    print(f'bp per Tree: {leng / ntrees}')
    print('----------------------')
    get_root_distrib(ts)
    
def to_tree(ts, n):
    """Return the nth tree from a tree sequence."""
    trees = [t for t in ts.trees()]
    return(trees[n])


def svg_to_pdf(svg, file, midfile = 'ueahfiauh.svg'):
    """Write an SVG string to a temporary file and convert it to a PDF."""
    midd = Path(midfile)
    with open(midd, 'w') as mid:
        mid.write(svg)
    drawing = svg2rlg(midd)
    renderPDF.drawToFile(drawing, file + '.pdf')
    midd.unlink()

def graph_tree_loc(tsname, loc=100000, timescale='rank'):
    """Render a tree at a genomic location to PDF(s) using `Tree.draw_svg`."""
    if not timescale in ['rank', 'time', 'log_time']:
        raise ValueError("timescale must be 'rank', 'time', or 'log_time'")
    ts = tskit.load(tsname + '.trees')
    assess_tree(ts)

    sized = (15000, 1000) # dimensions of final image (pixels)

    tree = ts.at(loc)
    node_label_style = (
        '.node > .lab {"font-size": "10%"}'
        '.leaf > .lab {"text-anchor": "start"; "transform": "rotate(90deg)"}'
        '.sample > .lab {"transform": rotate(90deg)'
        '.site > .lab {"font-size": "10%"}'
    )
    test_svg = tree.draw_svg(
        size = sized, 
        time_scale = timescale, 
        y_gridlines = True, 
        y_axis = True, 
        y_ticks = [0.001, 0.01, 0.1, 1], 
        x_axis = True,
        style=node_label_style, 
        mutation_labels={}
    )
    svg_to_pdf(test_svg, tsname)
    test_svg = tree.draw_svg(
        size = sized, 
        time_scale = timescale, 
        y_gridlines = True, 
        y_axis = True, 
        y_ticks = [0.001, 0.01, 0.1, 1], 
        x_axis = True,
        style=node_label_style, 
        #mutation_labels={}, 
        all_edge_mutations=True
    )
    svg_to_pdf(test_svg, tsname + '_all')
    return 0

if __name__ == '__main__':
    graph_tree_loc(treeseq)



