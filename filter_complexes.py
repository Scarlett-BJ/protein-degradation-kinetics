#!/usr/bin/env python3

import ixntools as ixn
from itertools import combinations

def load_table(n):
    with open('data/ned_mapped_to_pdb.out') as infile:
        strucs = {l.split('_')[0] for l in infile}
    table = ixn.structure.Table(n)
    useable = table.filter(name=strucs, unq_chns=range(2, 200))
    return useable

def get_chains(structures, threshold):
    sdict = {}
    with open('data/ned_mapped_to_pdb.out') as infile:
        for line in infile:
            line = line.split()
            struc, chn = line[0].split('_')
            prot = line[2]
            seqid = line[3]
            if struc not in sdict:
                sdict[struc] = [(prot, chn, seqid)]
            else:
                sdict[struc].append((prot, chn, seqid))
    useable = []
    # Allows for some missing chains, number indirectly defined by threshold.
    for struc in structures:
        if len(sdict[struc.name]) >= threshold*struc.unq_chns:
            useable.append((struc, sdict[struc.name]))
    return useable

def compare_pair(pair, threshold):
    """Compare 2 complexes, if sufficiently similar return True."""
    pair = sorted(pair, key=lambda x: len(x[1]))
    p1 = {p[0] for p in pair[0][1]}
    p2 = {p[0] for p in pair[1][1]}
    if p1 == p2:
        return pair[0]
    intersect = len(p1.intersection(p2))
    if intersect >= threshold*len(p2):
        return pair[0]

def avg_seqid(struc):
    seqdids = [float(i[1][2]) for i in struc]
    return sum(seqids)/len(seqids)

def filter_redundant(structures, threshold):
    omit = []
    pairs = list(combinations(structures, 2))
    pairs = [sorted(p, key=lambda x: len(x[1])) for p in pairs]
    pairs = sorted(pairs, key=lambda x: len(pairs[0][1]))
    # print(pairs[0][0])
    for pair in pairs:
        red = compare_pair(pair, threshold)
        if red != None:
            omit.append(red[0].name)
    omit = set(omit)
    filtered = []
    for struc in structures:
        if struc[0].name not in omit:
            filtered.append(struc)
    for struc in filtered:
        struc[1].sort(key=lambda x: x[0])
    return sorted(filtered, key=lambda x: x[1][1])

def filter_red2(structures, threshold):
    omit = []
    used = list(combinations(structures, 2))
    used = {tuple(sorted((i[0][0].name, i[1][0].name))) for i in used}
    structures = sorted(structures, key=lambda x: len(x[1]))
    useable = []
    for s1 in  structures:
        name1 = s1[0].name
        redund = []
        for s2 in structures:
            if s1 == s2:
                continue
            pair = s1, s2
            npair = tuple(sorted([s1[0].name, s2[0].name]))
            if npair in used:
                red = compare_pair(pair, threshold)
                if red != None:
                    redund.append(red)
                else:
                    useable.append(pair)
                used.remove(npair)
        if len(redund) != 0:
            omit.append(redund)
    for reds in omit[:1]:
        print(reds)

def map_pfam():
    pfam_dict = {}
    with open('data/10090.tsv') as infile:
        for line in infile:
            line = line.split()
            if line[0] not in pfam_dict:
                pfam_dict[line[0]] = [line[5]]
            else:
                pfam_dict[line[0]].append(line[5])
    pfam_dict = {p: {a for a in pfam_dict[p]} for p in pfam_dict}
    return pfam_dict

# def map_upr():
#     with open('data/mouse_genes_ned.txt') as infile:
#         uprmap = {line.split(',')[3]: line.split(',')[0].split('-')[0]
#                   for line in infile}
#     return uprmap

def filter_pfam(strucs):
    pfam_dict = map_pfam()
    for struc in strucs:
        prots = [p[0] for p in struc[1]]
        try:
            pfams = [pfam_dict[p] for p in prots]
            if len(pfams[0].intersection(pfams[1])) > 0:
                print('ooops')
            # print(pfams)
        except:
            print(prots)

def print_table():
    strucs = load_table(30)
    strucs = get_chains(strucs, 0.8)
    strucs = filter_redundant(strucs, 0.66)
    for i in strucs:
        info = str(i[0])
        prots = ','.join([a[0] for a in i[1]])
        print(info, prots, sep='\t')

def summary():
    strucs = load_table(30)
    strucs = get_chains(strucs, 0.8)
    strucs = filter_redundant(strucs, 0.66)
    all_lengths = [i[0].unq_chns for i in strucs]
    lengths = {i[0].unq_chns for i in strucs}
    for i in sorted(lengths):
        print(i, all_lengths.count(i))
    print('total:', len(strucs))

if __name__ == '__main__':
    strucs = load_table(30)
    strucs = get_chains(strucs, 0.8)
    filter_red2(strucs, 0.66)


