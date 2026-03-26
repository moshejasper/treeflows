# ts stats... 

import tskit
import argparse
import numpy as np

# assume treesequence... 

def ts_global_stats(ts):
    edges = ts.num_edges
    sites = ts.num_sites
    trees = ts.num_trees
    site_per_edge = sites / edges
    edgesum = np.sum(ts.edges_right - ts.edges_left)
    #edgesum = sum((n.right - n.left for n in ts.edges())) # probably can happen as array? 
    edgelength = edgesum / edges
    site_per_tree = sites / trees
    edges_per_tree = edges / trees
    # edge stats
    edgecoverage = edgesum / ts.sequence_length
    muts = ts.num_mutations

    print(f"Sites per edge: {site_per_edge:.2f}\tSites per tree: {site_per_tree:.2f}\tEdges per tree: {edges_per_tree:.2f}\t Edgelength: {edgelength:.2f}")
    print(f"Average edge coverage (site): {edgecoverage:.2f}\tMutation count: {muts}\tMutations per site: {muts/sites:.2f}")

# ts = tskit.load('chr2_basic_r7_simp.trees')
# ts_global_stats(ts)

def tree_local_stats(tree):
    # let's see
    edges = tree.num_edges
    # intv = tree.interval
    # span = tree.span
    muts = tree.num_mutations
    allmuts = tree_get_all_edge_mutations(tree)
    sites = tree.num_sites # probably same as mutations? (but we can compare)... (etc.)... 
    # rank = tree.rank()
    innernodes = tree_count_internal(tree)
    print(f"Tree length: {tree.span}\tInternal nodes: {innernodes}\tMutations: {muts}\tSites: {sites}\tM/S: {muts/sites:.2f}\tEdges: {edges}\tEdge Mutations: {allmuts}")

def tree_count_internal(tree):
    innercount = 0
    for node in tree.nodes():
        if tree.is_internal(node):
            innercount += 1
    return innercount

def tree_get_all_edge_mutations(tree):
    # adapted from draw function
    nodes = set(tree.nodes())
    tree_left = tree.interval.left
    tree_right = tree.interval.right
    ts = tree.tree_sequence # for reference
    edge_left = ts.tables.edges.left
    edge_right = ts.tables.edges.right
    node_edges = tree.edge_array
    #whittle mutations down
    mut_t = ts.tables.mutations
    focal_mutations = np.isin(mut_t.node, np.fromiter(nodes, mut_t.node.dtype)) # this is probably grabbing mutations at the same node but other edges!!!... 
    #return sum(focal_mutations)
    mutation_nodes = mut_t.node[focal_mutations]
    mutation_positions = ts.tables.sites.position[mut_t.site][focal_mutations]
    # print(mutation_positions)
    #print(edge_left[node_edges[mutation_nodes]])
    edges_with_muts = node_edges[mutation_nodes]
    mutleft = edge_left[edges_with_muts]
    mutright = edge_right[edges_with_muts]
    # print(len(mutleft))
    # print(len(mutation_positions))
    # print(len(mutright))
    #lb = mutleft <= mutation_positions
    #rb = mutation_positions <= mutright
    bb = (mutleft <= mutation_positions) & (mutation_positions <= mutright)
    num_muts = sum(bb)
    #tleftfilt = tree_left[edges_with_muts]
    #trightfilt = tree_right[edges_with_muts]
    #print(tleftfilt <= mutation_positions <= trightfilt)
    
    
    # mutation_ids = np.arange(ts.num_mutations, dtype=int)[focal_mutations]
    # there is more but this is probably good for now. 
    # num_muts = 0
    # for node, pos in zip(mutation_nodes, mutation_positions):
    #     curr_edge = node_edges[node]
    #     if curr_edge >= 0:
    #         if edge_left[curr_edge] <= pos <= edge_right[curr_edge]: #i.e. it is a legit mutation
    #             num_muts += 1
    return num_muts


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help="tree file for stats calculation")
    args = parser.parse_args()
    ts = tskit.load(args.file)
    print(f"Analysing file {args.file}")
    ts_global_stats(ts)
    tree = ts.first()
    for n in range(10):
        tree.next()
        tree_local_stats(tree)
    print('\n')