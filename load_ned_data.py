#!/usr/bin/env python3

import ixntools as ix
from abundance import paxdb
from itertools import combinations


def load_mouse_genes_ned():
    data_dict = {}
    with open('data/mouse_genes_ned.txt') as infile:
        data = [line.strip().split(',') for line in infile.readlines()]
    for line in data:
        data_dict[line[3]] = {'ensg': line[2], 'gene': line[1],
                              'upr': line[0], 'score': float(line[4]),
                              'ab': float(line[5])}
    return data_dict

def load_ned_to_pdb():
    pdb_map = {}
    with open('data/ned_mapped_to_pdb.out') as infile:
        data = [line.strip().split() for line in infile.readlines()]
    for line in data:
        struc = line[0].split('_')
        if struc[0] not in pdb_map:
            pdb_map[struc[0]] = {line[1]: {'chain': struc[1],
                                           'upr': line[2],
                                           'seqid': float(line[3]),
                                           'score':float(line[4]),
                                           'cat': line[5]}}
        else:
            pdb_map[struc[0]][line[1]] = {'chain': struc[1],
                                          'upr': line[2],
                                          'seqid': float(line[3]),
                                          'score':float(line[4]),
                                          'cat': line[5]}
    return pdb_map


def load_mouse_orthologs():
    with open('data/mouse_human_orthologs.txt') as infile:
        data = [line.strip().split('\t') for line in infile.readlines()]
        orthomap = {line[0]: line[1] for line in data if len(line) == 2}
    return orthomap

def append_human_abundances():
    meta = paxdb.get_metadata()
    filename = meta['COLON']['filename']
    abunds = paxdb.Abundances(filename)
    omap = load_mouse_orthologs()
    mned = load_mouse_genes_ned()
    for m in list(mned.keys()):
        if m not in omap:
            mned.pop(m)
    for m in mned:
        if not abunds.ismember(omap[m]):
            print(m, omap[m], mned[m]['score'], mned[m]['ab'], 'NA',
              sep='\t')
            continue
        print(m, omap[m], mned[m]['score'], mned[m]['ab'], abunds[omap[m]],
              sep='\t')

def get_complexes(lbound, ubound):
    complexes = {}
    with open('data/ned_mapped_to_pdb.out') as infile:
        data = [line.split() for line in infile.readlines()]
    for line in data:
        comp = line[0].split('_')[0]
        if comp not in complexes:
            complexes[comp] = [line[1]]
        else:
            complexes[comp].append(line[1])
    # Filter by size
    for comp in list(complexes):
        if len(complexes[comp]) < lbound or len(complexes[comp]) > ubound:
            complexes.pop(comp)
            continue
        complexes[comp] = tuple(sorted(complexes[comp]))
    # Filter by duplicates
    cset = set(complexes.values())
    final = {}
    for comp in list(complexes):
        if complexes[comp] in cset:
            cset.remove(complexes[comp])
            final[comp] = complexes[comp]
    return final

def print_pairs(comps):
    ginfo = load_mouse_genes_ned()
    omap = load_mouse_orthologs()
    for comp in comps:
        all_present = True
        for ensp in comps[comp]:
            if ensp not in ginfo or ensp not in omap:
                all_present = False
                break
        if all_present == False:
            continue
        pairs = combinations(comps[comp], 2)
        for pair in pairs:
            score0 = str(round(ginfo[pair[0]]['score'], 3))
            score1 = str(round(ginfo[pair[1]]['score'], 3))
            print(comp, '\t'.join(pair), score0, score1, sep='\t')

def get_abundances():
    meta = paxdb.get_metadata()
    files = []
    for i in meta:
        if i == 'WHOLE_ORGANISM' or i == 'SALIVA':
            continue
        elif float(meta[i]['score']) >= 6.6:
            files.append(meta[i]['filename'])
    return files

def read_pairs():
    with open('data/pairs.tmp') as infile:
        data = [line.strip().split('\t') for line in infile.readlines()]
    abunds = [paxdb.Abundances(fname) for fname in get_abundances()]
    omap = load_mouse_orthologs()
    pdb_map = load_ned_to_pdb()
    ints = ix.Interface()
    for line in data:
        comp = line[0]
        chainA = pdb_map[comp][line[1]]['chain']
        chainB = pdb_map[comp][line[2]]['chain']
        try:
            interface = ints.get_interface(comp, chainA, chainB)
            line.append(str(interface))
        except:
            line.append('NA')
        print('\t'.join(line))


def extremely_ropey_test():
    """For each complex, are NED proteins expressed across a wider range of
    tissues that ED proteins? Seems to be unhealthily dependent on how many
    possible tissues are used. More tissues means more lo-quality data..."""
    comps = get_complexes(2, 10)
    ginfo = load_mouse_genes_ned()
    omap = load_mouse_orthologs()
    abunds = [paxdb.Abundances(fname) for fname in get_abundances()]
    nedcount_greater = 0
    trials = 0
    for comp in comps:
        all_present = True
        for ensp in comps[comp]:
            if ensp not in ginfo or ensp not in omap:
                all_present = False
                break
        if all_present == False:
            continue
        neds = [ensp for ensp in comps[comp] if ginfo[ensp]['score'] > 0.9]
        eds = [ensp for ensp in comps[comp] if ginfo[ensp]['score'] < 0.1]
        if len(eds) == 0 or len(neds) == 0:
            continue
        nedcounts = []
        edcounts = []
        for ned in neds:
            n = 0
            for tissue in abunds:
                if tissue.isexpressed(omap[ned]):
                    n += 1
            nedcounts.append(n)
        for ed in eds:
            n = 0
            for tissue in abunds:
                if tissue.isexpressed(omap[ed]):
                    n += 1
            edcounts.append(n)
        if sum(nedcounts)/len(nedcounts) > sum(edcounts)/len(edcounts):
            nedcount_greater += 1
        trials += 1
    print(nedcount_greater, trials)

if __name__ == '__main__':
    read_pairs()
