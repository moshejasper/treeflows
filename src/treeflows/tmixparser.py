#!/bin/python3

# this script will attempt to enable the parsing of treemix files... 


# file imports... 

import gzip
import numpy as np
import dendropy
from collections import defaultdict
import tsinfer
from tskit import MISSING_DATA
import pandas as pd
import cyvcf2
from pathlib import Path
from statistics import mean
import argparse


edgefile = 'mj_init_may24_total_bgl_tmx4R_mig3_boot20.edges.gz'
vertfile = 'mj_init_may24_total_bgl_tmx4R_mig3_boot20.vertices.gz'
treefile = 'mj_init_may24_total_bgl_tmx4R_mig3_boot20.treeout.gz'
filebase = 'mj_init_may24_total_bgl_oagKR_mig'


def get_elines(edgefile):
    print('starting')
    ef = pd.read_table(edgefile, names = ('pnode', 'cnode', 'length', 'weight', 'mig', 'mstart', 'blank'), sep=' ')
    print(ef[ef['mig'] == 'MIG'])


def get_vlines(vertfile):
    print('starting')
    vf = pd.read_table(vertfile, names = ('id', 'label', 'root', 'mig', 'tip', 'pid', 'c1id', 'nc1', 'c2id', 'nc2', 'subtree'), sep=' ')
    print(vf[vf['mig'] == 'MIG'])

def get_subtree(vf, id):
    if isinstance(id, pd.core.series.Series):
        id = int(id.iloc[0])
    if not isinstance(id, int):
        id = int(id)
    focus = vf[vf['id'] == id]
    subtree = focus['subtree']
    if len(subtree) > 1:
        exit("More than one vertex with same name!")
    arr_str = ', '.join(map(str, subtree.values.flatten()))
    arr_str = arr_str + ';'
    return(arr_str)


def get_treemaster(treefile):
    dtree = dendropy.Tree.get_from_path(treefile, schema='newick')
    dtree.infer_taxa()

def get_biparts(vertfile, edgefile, diff=False):
    vf = pd.read_table(vertfile, names = ('id', 'label', 'root', 'mig', 'tip', 'pid', 'c1id', 'nc1', 'c2id', 'nc2', 'subtree'), sep=' ')
    ef = pd.read_table(edgefile, names = ('pnode', 'cnode', 'length', 'weight', 'mig', 'mstart', 'blank'), sep=' ')

    mig_verts = vf[vf['mig'] == 'MIG']
    mig_edges = ef[ef['mig'] == 'MIG']
    migs = list()
    migweights = list()
    for i, row in mig_verts.iterrows():
        pid = int(row['pid'])
        c1id = int(row['c1id'])
        id = int(row['id'])
        migweight = float(mig_edges[mig_edges['pnode']==id]['weight'].iloc[0])
        c2id = mig_edges[mig_edges['pnode']==id]['cnode']
        ptree = get_subtree(vf, pid)
        c1tree = get_subtree(vf, c1id)
        c2tree = get_subtree(vf, c2id)
        trees = [ptree, c1tree, c2tree]
        taxas = [0,  0, 0]
        for n, tree in enumerate(trees):
            #print(str(tree))
            dtree = dendropy.Tree.get_from_string(tree, "newick")
            #dtree.print_plot()
            dtaxa = set()
            for node in dtree.leaf_node_iter():
                #print(type(node.taxon))
                #print(type(node.taxon.__str__()))
                dtaxa.add(node.taxon.__str__().strip("'"))
            #print(dtaxa)
            taxas[n] = dtaxa
        if diff:
            taxdiff = [taxa for taxa in taxas[0] if not taxa in taxas[1]]
            taxas[0] = taxdiff
        migs.append(taxas)
        migweights.append(migweight)

    return(migs, migweights)
    
def compare_biparts(filebase, start=0, stop=99, threshold = 0.5, mig=0, diff=False, summary=False):
    partlist = list()
    weightlist = list()
    if summary:
        diff = False

    treecount = stop - start + 1
    treeinit = treecount

    for n in range(start, stop+1):
        bootbase = filebase + str(mig) + "_boot" + str(n)
        efile = bootbase + '.edges.gz'
        vfile = bootbase + '.vertices.gz'
        if not Path(efile).exists():
            print(f'Missing bootstrap #{n}. skipping...')
            treecount -= 1
            continue
        biparts, migweights = get_biparts(vfile, efile, diff=diff)
        #print(migweights)
        partlist.extend(biparts)
        weightlist.extend(migweights)
    partdict = dict()
    weightdict = dict()
    #print(partlist)
    if len(partlist) == 0:
        print(f'No migration events detected across {treecount} successfully located bootstraps! ({treeinit - treecount} not located)')
        return(partdict)
    unpart = [partlist[0],] # first element of bipart list... 
    unweight = weightlist[0] # first element of weight list
    partdict[0] = 1
    weightdict[0] = [unweight,] ### finish here for the minute... 
    counter = 1
    for n, pp in enumerate(partlist[1:]):
        try:
            index = unpart.index(pp)
            partdict[index] += 1
            weightdict[index].append(weightlist[n])
        except ValueError:
            unpart.append(pp)
            partdict[counter] = 1
            weightdict[counter] = [weightlist[n],]
            counter += 1
    
    for k in weightdict.keys():
        weightdict[k] = mean(weightdict[k])

    # rename variables::: 
    if summary:
        cladenames, cladesets = get_clades()
        for m, parts in enumerate(unpart):
            for mm, part in enumerate(parts):
                for n, clade in enumerate(cladesets):
                    if all(v in part for v in clade):
                        part = [v for v in part if v not in clade]
                        #part = part[partidx]
                        part.append(cladenames[n])
                        #print(part)
                        #print(cladenames[n])
                parts[mm] = part
            #print(parts)
            unpart[m] = parts
    thresh = treecount * threshold
    print(f'Fractional threshold: {threshold}')
    print(f'Bootstraps located: {treecount}/{treeinit}')
    print(f'Count threshold: {thresh}')
    if diff:
        print('DIFF:')
    else:
        print('CLADES:')
    for n in range(len(partdict)):
        #print(partdict[n])
        if partdict[n] >= thresh:
            print('\n------')
            print(unpart[n][0])
            print(f'   |')
            print(f'   |--------------> {unpart[n][2]}')
            print(f'   V')
            print(unpart[n][1])
            print('------')
            print(f'Bootstraps: {partdict[n]}')
            print(f'Mean migration weight: {round(weightdict[n], 3)}\n')


def print_unpart(upart, boot, wdict):
    print('\n------')
    print(upart[0])
    print(f'   |')
    print(f'   |--------------> {upart[2]}')
    print(f'   V')
    print(upart[1])
    print('------')
    print(f'Bootstraps: {boot}')
    print(f'Mean migration weight: {round(wdict, 3)}\n')


def get_clades():
    cladenames = ['ROOT', 'All Aegypti (Africa samples)', 'Aegypti precursors (Ngoye)', 'Proto-aegypti (Paraguay)', 'Aegypti aegypti', 'Asia-Pacific', 'Asia-Pacific (xN)', 'Oceania', 
                  'Oceania (xN)', 'Papua New Guinea', 'Australia', 'Asia (xJ)', 'South-East Asia', 'South America', 'North-Central America (xSL)']
    cladesets = [['ANHM', 'ASNC', 'BALI', 'BNGK', 'CALI', 'CLMB', 'CPYR', 'CRNS', 'CRNV', 'DILI', 'DRNG', 'EFAT', 'ESNB', 'FRSN', 'JDDH', 'KHSN', 'KLLM', 'KRBT', 'KWAL', 'LAE', 'LNGR', 'MDNG', 'MLDV', 'MNTM', 'NADI', 'NGOY', 'NKLF', 'NOUM', 'OGDG', 'PHNX', 'PNNG', 'PRTM', 'RCLR', 'RDJN', 'SNGP', 'SNTR', 'STLC', 'THIS', 'TPCH', 'TRNI', 'TTLI', 'YGYK'],
                 ['ANHM', 'ASNC', 'BALI', 'BNGK', 'CALI', 'CLMB', 'CPYR', 'CRNS', 'CRNV', 'DILI', 'DRNG', 'EFAT', 'ESNB', 'FRSN', 'JDDH', 'KHSN', 'KLLM', 'KRBT', 'LAE', 'LNGR', 'MDNG', 'MLDV', 'MNTM', 'NADI', 'NGOY', 'NKLF', 'NOUM', 'OGDG', 'PHNX', 'PNNG', 'PRTM', 'RCLR', 'RDJN', 'SNGP', 'SNTR', 'STLC', 'THIS', 'TPCH', 'TRNI', 'TTLI', 'YGYK'],
                 ['ANHM', 'ASNC', 'BALI', 'BNGK', 'CALI', 'CLMB', 'CPYR', 'CRNS', 'CRNV', 'DILI', 'DRNG', 'EFAT', 'ESNB', 'FRSN', 'JDDH', 'KHSN', 'KLLM', 'KRBT', 'LAE', 'LNGR', 'MDNG', 'MLDV', 'MNTM', 'NADI', 'NGOY', 'NKLF', 'NOUM', 'PHNX', 'PNNG', 'PRTM', 'RCLR', 'RDJN', 'SNGP', 'SNTR', 'STLC', 'TPCH', 'TRNI', 'TTLI', 'YGYK'],
                 ['ANHM', 'ASNC', 'BALI', 'BNGK', 'CALI', 'CLMB', 'CPYR', 'CRNS', 'CRNV', 'DILI', 'DRNG', 'EFAT', 'ESNB', 'FRSN', 'JDDH', 'KHSN', 'KLLM', 'KRBT', 'LAE', 'LNGR', 'MDNG', 'MLDV', 'MNTM', 'NADI', 'NKLF', 'NOUM', 'PHNX', 'PNNG', 'PRTM', 'RCLR', 'RDJN', 'SNGP', 'SNTR', 'STLC', 'TPCH', 'TRNI', 'TTLI', 'YGYK'],
                 ['ANHM', 'BALI', 'BNGK', 'CALI', 'CLMB', 'CPYR', 'CRNS', 'CRNV', 'DILI', 'DRNG', 'EFAT', 'ESNB', 'FRSN', 'JDDH', 'KHSN', 'KLLM', 'KRBT', 'LAE', 'LNGR', 'MDNG', 'MLDV', 'MNTM', 'NADI', 'NKLF', 'NOUM', 'PHNX', 'PNNG', 'PRTM', 'RCLR', 'RDJN', 'SNGP', 'SNTR', 'STLC', 'TPCH', 'TRNI', 'TTLI', 'YGYK'],
                 ['BALI', 'BNGK', 'CLMB', 'CPYR', 'CRNS', 'DILI', 'DRNG', 'EFAT', 'ESNB', 'JDDH', 'KHSN', 'KLLM', 'KRBT', 'LAE', 'LNGR', 'MDNG', 'MLDV', 'MNTM', 'NADI', 'NKLF', 'NOUM', 'PNNG', 'PRTM', 'SNGP', 'TRNI', 'TTLI', 'YGYK'],
                 ['BALI', 'BNGK', 'CLMB', 'CPYR', 'CRNS', 'DILI', 'DRNG', 'EFAT', 'ESNB', 'JDDH', 'KHSN', 'KLLM', 'KRBT', 'LAE', 'LNGR', 'MDNG', 'MLDV', 'MNTM', 'NADI', 'NKLF', 'PNNG', 'PRTM', 'SNGP', 'TRNI', 'TTLI', 'YGYK'],
                 ['EFAT', 'KRBT', 'NADI', 'NKLF', 'NOUM', 'TTLI'],
                 ['EFAT', 'KRBT', 'NADI', 'NKLF', 'TTLI'],
                 ['ESNB', 'LAE', 'MDNG', 'PRTM'],
                 ['CPYR', 'CRNS', 'DRNG', 'LNGR', 'MNTM'],
                 ['BALI', 'BNGK', 'CLMB', 'KHSN', 'KLLM', 'MLDV', 'PNNG', 'SNGP', 'TRNI', 'YGYK', 'DILI'],
                 ['BALI', 'BNGK', 'KHSN', 'KLLM', 'PNNG', 'SNGP', 'TRNI', 'YGYK', 'DILI'],
                 ['CALI', 'RCLR', 'RDJN', 'SNTR'],
                 ['ANHM', 'CRNV', 'FRSN', 'PHNX', 'TPCH']]
    return(cladenames, cladesets)
#def make_

#get_biparts(vertfile, edgefile)
#compare_biparts(filebase, threshold = 0.05, mig=1   , diff=False, summary=True)
#get_treemaster(treefile)

# find parent, child, connector... --> load all three as tree (100 times) --> count equivalents --> return list... (hmm)... 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize bootstraps from Treemix-type Analyses")
    parser.add_argument("--fileprefix", type=str, help="Prefix of treemix bootstrap files")
    parser.add_argument("--threshold", type=float, help="Fraction of trees migration edge must be present in to be included")
    parser.add_argument("--mig", type=int, help="Migration edges subgroup to summarise")
    args = parser.parse_args()
    compare_biparts(args.fileprefix, threshold=args.threshold, mig=args.mig, diff=False, summary=True)