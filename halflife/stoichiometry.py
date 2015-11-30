#!/usr/bin/env python3

import numpy as np
from collections import namedtuple
from itertools import combinations
from scipy.stats import binom_test
from ixntools import structure

def load_mouse_prot_map():
    data_dict = {}
    with open('data/mouse_genes_ned.txt') as infile:
        data = [line.strip().split(',') for line in infile.readlines()]
    for line in data:
        data_dict[line[0]] = {'ensg': line[2], 'gene': line[1],
                              'ensp': line[3]}
    return data_dict

def load_data(filename):
    with open(filename) as infile:
        data = [line.strip().split('\t') for line in infile]
    header = data[0]
    data = data[1:]
    return header, data

def get_complexes(header, data):
    """Returns a dictionary of complexes. Each complex contains a list of
    namedtuples with information about the complex and proteins.
    """
    [str(i.strip('.')) for i in header[1:]]
    Info = namedtuple('Header', [str(i.strip('.')) for i in header[1:]])
    cdict = {}
    for line in data:
        if line[0] not in cdict:
            cdict[line[0]] = [Info(*line[1:])]
        else:
            cdict[line[0]].append(Info(*line[1:]))
    return cdict

def get_stoichiometry():
    header, data = load_data('data/protwise.tsv')
    cdict = get_complexes(header, data)
    chn_map = {}
    with open('data/ned_mapped_to_pdb.out') as infile:
        for line in infile:
            line = line.split()
            pdb, chn = line[0].split('_')
            if pdb not in chn_map:
                chn_map[pdb] = {line[2]: chn}
            else:
                chn_map[pdb][line[2]] = chn
    stoich_map = {}
    similar = structure.similar_chains()
    Info = namedtuple('Info', ['prot', 'chn', 'dcls', 'copy'])
    for pdb in cdict:
        new_info = []
        chns = {prot.prot: chn_map[pdb][prot.prot] for prot in cdict[pdb]}
        dcls = {prot.prot: prot.dcls for prot in cdict[pdb]}
        copies = {}
        for prot in chns:
            for unq in similar[pdb]:
                if chns[prot] in unq:
                    copies[prot] = len(unq)
            info = Info(prot, chn_map[pdb][prot], dcls[prot], copies[prot])
            new_info.append(info)
        cdict[pdb] = new_info
    return cdict

def analyse_stoichiometry_percomp():
    cdict = get_stoichiometry()
    success = 0
    trials = 0
    for pdb in cdict:
        comp = cdict[pdb]
        evals = []
        nvals = []
        uvals = []
        for prot in comp:
            if prot.dcls == 'E':
                evals.append(prot.copy)
            elif prot.dcls == 'N':
                nvals.append(prot.copy)
            elif prot.dcls == 'U':
                uvals.append(prot.copy)
        if len(nvals) > 0 and len(evals) > 0 and len(uvals) > 0:
            print(np.mean(nvals), np.mean(evals), np.mean(uvals))

def symmetry():
    header, data = load_data('data/protwise.tsv')
    cdict = get_complexes(header, data)
    symclasses = {line[3]: [] for line in data}
    for pdb in cdict:
        dclasses = {c.dcls for c in cdict[pdb]}
        sym = cdict[pdb][0].sym
        for d in dclasses:
            symclasses[sym].append(d)
    for c in ['C2', 'C2h', 'H', 'D']:
        symclasses['C'] += symclasses[c]
        symclasses.pop(c)
    for s in symclasses:
        print(s, symclasses[s].count('N') + symclasses[s].count('U'), symclasses[s].count('E'))

if __name__ == '__main__':
    symmetry()
