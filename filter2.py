#!/usr/bin/env python3

import ixntools as ixn
from itertools import combinations
from collections import namedtuple
import dbscan

def map_pfam():
    pfam_dict = {}
    with open('data/10090.tsv') as infile:
        for line in infile:
            line = line.split()
            if line[0] not in pfam_dict:
                pfam_dict[line[0]] = [line[5]]
            else:
                pfam_dict[line[0]].append(line[5])
    # Possibly unhelpful?
    pfam_dict = {p: {a for a in pfam_dict[p]} for p in pfam_dict}
    return pfam_dict

def load_table(n):
    with open('data/ned_mapped_to_pdb.out') as infile:
        strucs = {l.split('_')[0] for l in infile}
    table = ixn.structure.Table(n)
    useable = table.filter(name=strucs, unq_chns=range(2, 200))
    return useable

def get_chains(pool, threshold):
    chaindict = {}
    Chain = namedtuple('Chain', ['prot', 'chn', 'seqid'])
    with open('data/ned_mapped_to_pdb.out') as infile:
        for line in infile:
            struc, _, prot, seqid, _, _ = line.split()
            struc, chn = struc.split('_')
            if struc not in chaindict:
                chaindict[struc] = [Chain(prot, chn, float(seqid))]
            else:
                chaindict[struc].append(Chain(prot, chn, float(seqid)))
    useable = []
    # Allows for some missing chains, number indirectly defined by threshold.
    for struc in pool:
        chncount = len(chaindict[struc.name])
        trucount = struc.unq_chns
        if chncount >= threshold*trucount/100 and chncount <= trucount:
            struc.chains = chaindict[struc.name]
            useable.append(struc)
    return useable

def get_pfams(pool):
    Pfam = namedtuple('Pfam', ['prot', 'pfam'])
    pfam_dict = map_pfam()
    for struc in pool:
        prots = [p.prot for p in struc.chains]
        pfams = [Pfam(p, pfam_dict.get(p, set())) for p in prots]
        struc.pfams = pfams
    return pool

def cluster_redundant(pool, threshold, key='Chains'):
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
    useable = []
    for k in clusters:
        # First sort by length of structure, then by seqid
        strucs = sorted(clusters[k].members,
                        key=lambda x: (len(x.value), x.original_attr.avgseqid),
                        reverse=True)
        useable.append(strucs[0].original_attr)
    return useable


def filter_pfam(pool, threshold):
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
    all_lengths = [i.unq_chns for i in pool]
    lengths = {i.unq_chns for i in pool}
    for i in sorted(lengths):
        print(i, all_lengths.count(i))
    print('total:', len(pool))


if __name__ == '__main__':
    pool = load_table(30)
    pool = get_chains(pool, 88)
    pool = get_pfams(pool)
    clusters = cluster_redundant(pool, 34, 'Domains')
    pool = filter_redundant(clusters)
    summary(pool)
    # clusters = cluster_redundant(pool, 34, 'Chains')
    # pool = filter_redundant(clusters)
    pool = filter_pfam(pool, 50)
    summary(pool)