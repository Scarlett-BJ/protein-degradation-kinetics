#!/usr/bin/env python3

import numpy as np
from expression import paxdb, coexpressdb
from collections import namedtuple
from itertools import combinations
from scipy.stats import binom_test

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

def get_coexpression():
    gene_map = load_mouse_prot_map()
    coex = coexpressdb.Coexpression()
    data = load_pairwise()
    for line in data[1:]:
        if line[4] == 'NA' or line[7] == 'NA':
            print('\t'.join(line))
            continue
        ensp1 = gene_map[line[4]]['ensp']
        ensp2 = gene_map[line[7]]['ensp']
        if ensp1 not in coex.entrez or ensp2 not in coex.entrez:
            print('\t'.join(line))
            continue
        c = coex.get_coexpression(ensp1, ensp2)
        line.append(str(c))
        print('\t'.join(line))


def interfaces():
    """Uses measured interfaces, rather than all possible combinations of pairs
    in each complex.
    """
    header, data = load_data('data/pairwise_trimers.tsv')
    cdict = get_complexes(header, data)
    gene_map = load_mouse_prot_map()
    coex = coexpressdb.Coexpression()
    trials = 0
    success = 0
    for pdb in cdict:
        nvals = []
        evals = []
        comp = cdict[pdb]
        for pair in comp:
            ensp1 = gene_map[pair.p1]['ensp']
            ensp2 = gene_map[pair.p2]['ensp']
            # Remove homologous interactions, where cor == 1.0
            if ensp1 == ensp2:
                continue
            if float(pair.int) < 200:
                continue
            if ensp1 not in coex.entrez or ensp2 not in coex.entrez:
                continue
            c = coex.get_coexpression(ensp1, ensp2)
            if pair.p1cls != 'E' and pair.p2cls != 'E':
                nvals.append(c)
            elif pair.p1cls == 'E' or pair.p2cls == 'E':
                evals.append(c)
            if len(nvals) > 0 and len(evals) > 0:
                trials += 1
                if np.mean(nvals) > np.mean(evals):
                    success += 1
    print(success, trials, binom_test(success, trials))

def pairwise_combinations():
    header, data = load_data('data/protwise_trimers.tsv')
    cdict = get_complexes(header, data)
    gene_map = load_mouse_prot_map()
    coex = coexpressdb.Coexpression()
    trials = 0
    success = 0
    for pdb in cdict:
        nvals = []
        evals = []
        comp = list(combinations(cdict[pdb], 2))
        for pair in comp:
            ensp1 = gene_map[pair[0].prot]['ensp']
            ensp2 = gene_map[pair[1].prot]['ensp']
            if ensp1 not in coex.entrez or ensp2 not in coex.entrez:
                continue
            c = coex.get_coexpression(ensp1, ensp2)
            if pair[0].dcls != 'E' and pair[1].dcls != 'E':
                nvals.append(c)
            elif pair[0].dcls == 'E' or pair[1].dcls == 'E':
                evals.append(c)
            if len(nvals) > 0 and len(evals) > 0:
                trials += 1
                if np.mean(nvals) > np.mean(evals):
                    success += 1
    print(success, trials, binom_test(success, trials))

def avg_coexpression():
    header, data = load_data('data/protwise.tsv')
    cdict = get_complexes(header, data)
    gene_map = load_mouse_prot_map()
    coex = coexpressdb.Coexpression()
    success = 0
    trials = 0
    for pdb in cdict:
        nvals = []
        evals = []
        for prot1 in cdict[pdb]:
            avg_coex = []
            for prot2 in cdict[pdb]:
                if prot1 == prot2:
                    continue
                ensp1 = gene_map[prot1.prot]['ensp']
                ensp2 = gene_map[prot2.prot]['ensp']
                if ensp1 not in coex.entrez or ensp2 not in coex.entrez:
                    continue
                avg_coex.append(coex.get_coexpression(ensp1, ensp2))
            if len(avg_coex) == 0:
                continue
            if prot1.dcls == 'N':
                nvals.append(np.mean(avg_coex))
            elif prot1.dcls == 'E':
                evals.append(np.mean(avg_coex))
            if len(evals) > 0 and len(nvals) > 0:
                if np.mean(nvals) == np.mean(evals):
                    continue
                trials += 1
                if np.mean(nvals) > np.mean(evals):
                    success += 1
    print(success, trials, binom_test(success, trials))

def abundance_similarity():
    header, data = load_data('data/protwise_trimers.tsv')
    cdict = get_complexes(header, data)
    success = 0
    trials = 0
    for pdb in cdict:
        nvals = []
        evals = []
        comp = cdict[pdb]
        nprots = [prot for prot in comp if prot.dcls == 'N']
        eprots = [prot for prot in comp if prot.dcls == 'E' or prot.dcls == 'U']
        mean_n = []
        mean_e = []
        # if len(nprots) == 0 or len(eprots) == 0:
        #     continue
        for pn in nprots:
            for pn2 in nprots:
                if pn == pn2:
                    continue
                mean_n.append(((float(pn.ab) - float(pn2.ab))**2)**0.5)
            for pe in eprots:
                mean_e.append(((float(pn.ab) - float(pe.ab))**2)**0.5)
        if len(mean_n) == 0 or len(mean_e) == 0:
            continue
        trials += 1
        if np.mean(mean_n) > np.mean(mean_e):
            success += 1
    print(success, trials, binom_test(success, trials))

if __name__ == '__main__':
    abundance_similarity()