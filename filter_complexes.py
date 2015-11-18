#!/usr/bin/env python3

import ixntools as ixn
from itertools import combinations
from collections import namedtuple
import dbscan
from random import randint

def map_pfam():
    """Returns dictionary of PFAM domains mapped to uniprot ids."""
    pfam_dict = {}
    ribo = []
    with open('data/10090.tsv') as infile:
        for line in infile:
            line = line.strip().split('\t')
            if len(line) != 14:
                continue
            if line[0] not in pfam_dict:
                pfam_dict[line[0]] = [line[5]]
            else:
                pfam_dict[line[0]].append(line[5])
            if 'Ribosomal' in line[6]:
                ribo.append(line[5])
    pfam_dict = {p: {a for a in pfam_dict[p]} for p in pfam_dict}
    return pfam_dict, ribo

def load_table(n):
    """Select appropriate minimum chain length, and then return a list of
    complexes from the relevant table. i.e. table30.out for minimum chain
    length of 30 residues.
    """
    with open('data/ned_mapped_to_pdb.out') as infile:
        strucs = {l.split('_')[0] for l in infile}
    table = ixn.structure.Table(n)
    useable = table.filter(name=strucs, unq_chns=range(2, 500))
    return useable


def get_chains(pool, threshold):
    """For all structures in table, if at least x percent of that structure’s
    chains (x ~ 90) map to a mouse protein, then add it to a list of
    potentially useable structures.
    """
    chaindict = {}
    Chain = namedtuple('Chain', ['prot', 'chn', 'seqid', 'score', 'dclass'])
    with open('data/ned_mapped_to_pdb.out') as infile:
        for line in infile:
            struc, _, prot, seqid, score, dclass = line.split()
            struc, chn = struc.split('_')
            if struc not in chaindict:
                chaindict[struc] = [Chain(prot, chn, float(seqid),
                                    float(score), dclass)]
            else:
                chaindict[struc].append(Chain(prot, chn, float(seqid),
                                        float(score), dclass))
    useable = []
    # Allows for some missing chains, number indirectly defined by threshold.
    for struc in pool:
        chncount = len(chaindict[struc.name])
        trucount = struc.unq_chns
        if chncount >= threshold*trucount/100:
            struc.chains = chaindict[struc.name]
            useable.append(struc)
    return useable

def interfaces(pool):
    useable = []
    cmap = ixn.structure.chain_map()
    ints = ixn.structure.Interfaces()
    for struc in pool:
        try:
            pdb_chns = {c.chn for c in struc.chains}
            chns = {cmap[struc.name][c.chn] for c in struc.chains}
            print(chns)
        except:
            print()
            print(struc.name, pdb_chns, cmap[struc.name])
            print()

def filter_missing_classes(pool):
    useable = []
    for struc in pool:
        dclass = {c.dclass for c in struc.chains}
        if 'E' in dclass:
            useable.append(struc)
    return useable

def get_pfams(pool):
    """Adds PFAM domain annotations to proteins in cluster."""
    Pfam = namedtuple('Pfam', ['prot', 'pfam'])
    pfam_dict, ribo = map_pfam()
    for struc in pool:
        prots = [p.prot for p in struc.chains]
        pfams = [Pfam(p, pfam_dict.get(p, {randint(1, 10000)})) for p in prots]
        struc.pfams = pfams
    return pool

def filter_ribosomal(pool):
    useable = []
    pfam_dict, ribo = map_pfam()
    for struc in pool:
        pfams = [p.pfam for p in struc.pfams]
        pfams = {p for sublist in pfams for p in sublist}
        if len(set(ribo).intersection(pfams))/len(pfams) <= 0.4:
            useable.append(struc)
    return useable

def cluster_redundant(pool, threshold, key='Chains'):
    """Cluster groups of complexes together based on similarity of chains or
    pfam domains present, e.g. if less than 34 percent of chains are shared
    between two complexes, then they belong in different clusters.
    """
    points = []
    for struc in pool:
        chains = [p.prot for p in struc.chains]
        pfams = [p.pfam for p in struc.pfams]
        pfams = [p for sublist in pfams for p in sublist]
        if key == 'Domains':
            point = dbscan.Point(struc.name, pfams, struc)
        elif key == 'Chains':
            point = dbscan.Point(struc.name, chains, struc)
        points.append(point)
    scanr = dbscan.DBSCAN(points, threshold, 1)
    clusts = scanr.get_clusters()
    return clusts

def filter_redundant(clusters):
    """Sorts each cluster based on first the number of chains, and then the
    average sequence identity of mapped chains.
    """
    useable = []
    for k in clusters:
        # First sort by length of structure, then by seqid
        strucs = sorted(clusters[k].members,
                        key=lambda x: (len(x.value), x.original_attr.avgseqid),
                        reverse=True)
        useable.append(strucs[0].original_attr)
    return useable

def filter_ig(pool):
    """Removes structures containing mouse immunoglobin domains"""
    domains = {'ig','Ig_2', 'Ig_3' 'C1-set', 'C1-set_C', 'C2-set_2', 'C2-set',
               'V-set', 'V-set_2', 'V-set_CD46'}
    with open('data/10090.tsv') as infile:
        igs = {i.split()[5] for i in infile if i.split()[6] in domains}
    for struc in pool:
        for prot in struc.pfams:
            if len(prot.pfam.intersection(igs)) > 1:
                print(struc, prot.name)

def filter_pfam(pool, threshold):
    """Filters out complexes where proteins within complexes are paralogous.
    Suspect this is not currently a very good way of doing it.
    """
    useable = []
    for struc in pool:
        pfams = [p.pfam for p in struc.pfams]
        pfams = [p for sublist in pfams for p in sublist]
        count = 0
        for domain in set(pfams):
            c = pfams.count(domain) - 1
            count += c
        if count/len(pfams)*100 < threshold:
            useable.append(struc)
    return useable

def summary(pool):
    """Summarises distribution of different sized complexes in a pool."""
    all_lengths = [i.tot_chns for i in pool]
    lengths = {i.tot_chns for i in pool}
    for i in sorted(lengths):
        print(i, all_lengths.count(i))
    print('total:', len(pool))

def display_pairwise(pool):
    for struc in sorted(pool, key = lambda x: x.name):
        pairs = combinations(struc.chains, 2)
        for pair in pairs:
            info = '\t'.join([x.chn for x in pair] + [x.prot for x in pair] + [x.dclass for x in pair])
            print(struc, info, sep='\t')

if __name__ == '__main__':
    pool = load_table(1)
    pool =
    pool = get_chains(pool, 100)
    interfaces(pool)
    # pool = filter_missing_classes(pool)
    # pool = get_pfams(pool)
    # clusters = cluster_redundant(pool, 66, 'Chains')
    # pool = filter_redundant(clusters)
    # clusters = cluster_redundant(pool, 100, 'Domains')
    # pool = filter_redundant(clusters)
    # pool = filter_pfam(pool, 50)
    # pool = filter_ribosomal(pool)
    # display_pairwise(pool)
