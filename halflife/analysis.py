#!/usr/bin/env python3

import numpy as np
from expression import paxdb, coexpressdb
from collections import namedtuple
from itertools import combinations

def load_pairwise():
    with open('data/pairwise.tsv') as infile:
        data = [line.split() for line in infile]
    return data

def load_protwise():
    with open('data/protwise.tsv') as infile:
        data = [line.split() for line in infile]
    return data

def load_mouse_prot_map():
    data_dict = {}
    with open('data/mouse_genes_ned.txt') as infile:
        data = [line.strip().split(',') for line in infile.readlines()]
    for line in data:
        data_dict[line[0]] = {'ensg': line[2], 'gene': line[1],
                              'ensp': line[3]}
    return data_dict

def interfaces():
    """Are interfaces between 'N' proteins larger than average in complex."""
    data = load_pairwise()
    cdict = {}
    for line in data[1:]:
        if line[1] == '2':
            continue
        if line[0] not in cdict:
            cdict[line[0]] = [(line[6], line[9], line[-1])]
        else:
            cdict[line[0]].append((line[6], line[9], line[-1]))
    success = 0
    trials = 0
    for c in cdict:
        classes = {i[0] for i in cdict[c]}.union({i[1] for i in cdict[c]})
        ints = cdict[c]
        nface = {float(i[2]) for i in ints if
                (i[1] == 'N' or i[0] == 'N') and float(i[2]) > 200}
        therest = {float(i[2]) for i in ints if float(i[2]) > 200} - nface
        if len(nface) > 0 and len(therest) > 0:
            trials += 1
            if np.mean(list(therest)) < np.mean(list(nface)):
                success += 1
    print(success, trials)
    print(len(cdict))


def interface_distribution():
    data = load_pairwise()
    efaces = []
    for line in data[1:]:
        if line[6] == 'U' and line[9] == 'U' and float(line[10]) > 200:
            print(line[10])

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

def analyse_coexpression():
    with open('data/coexpression.tsv') as infile:
        data = [line.strip().split('\t') for line in infile]
    nvals = []
    evals = []
    for line in data:
        if len(line) != 12:
            continue
        if line[6] != 'E' and line[9] != 'E':
            nvals.append(float(line[11]))
        elif line[6] == 'E' or line[9] == 'E':
            evals.append(float(line[11]))
    for i in nvals:
        print('n', i)
    for i in evals:
        print('e', i)

def analyse_coexpression_percomp():
    """This uses interfaces, probably a fairer way to do it is to use all
    possible pairwise combinations, as I suspect this is skewed by the fact
    that there are many homomeric interfaces present.
    """
    with open('data/coexpression.tsv') as infile:
        data = [line.strip().split('\t') for line in infile]
    cdict = {}
    coex = namedtuple('Int', ['c1', 'c2', 'coex'])
    for line in data:
        if len(line) != 12:
            continue
        if line[0] not in cdict:
            cdict[line[0]] = [coex(line[6], line[9], line[11])]
        else:
            cdict[line[0]].append(coex(line[6], line[9], line[11]))
    trials = 0
    success = 0
    for pdb in cdict:
        comp = cdict[pdb]
        nvals = []
        evals = []
        for i in comp:
            if i.c1 != 'E' and i.c2 != 'E':
                nvals.append(float(i.coex))
            elif i.c1 == 'E' or i.c2 == 'E':
                evals.append(float(i.coex))
            if len(nvals) > 0 and len(evals) > 0:
                trials += 1
                if np.mean(nvals) > np.mean(evals):
                    success += 1
    print(success, trials)

def get_pairwise_coexpression():
    gene_map = load_mouse_prot_map()
    coex = coexpressdb.Coexpression()
    data = load_protwise()
    cdict = {}
    for line in data[1:]:
        if line[0] not in cdict:
            cdict[line[0]] = [tuple(line)]
        else:
            cdict[line[0]].append(tuple(line))
    trials = 0
    success = 0
    # nvals = []
    # evals = []
    print(len(cdict))
    for pdb in cdict:
        # if len(cdict[pdb]) < 3:
        #     continue
        evals = []
        nvals = []
        pairs = list(combinations(cdict[pdb], 2))
        for p in pairs:
            c1 = p[0][7]
            c2 = p[1][7]
            ensp1 = gene_map[p[0][4]]['ensp']
            ensp2 = gene_map[p[1][4]]['ensp']
            if ensp1 not in coex.entrez or ensp2 not in coex.entrez:
                continue
            c = coex.get_coexpression(ensp1, ensp2)
            if c1 == 'N' or c2 == 'N':
                nvals.append(c)
            elif c1 == 'E' or c2 == 'E':
                if c1 == 'N' or c2 == 'N':
                    continue
                evals.append(c)
        if len(nvals) > 0 and len(evals) > 0:
            trials += 1
            if np.mean(nvals) > np.mean(evals):
                success += 1
    print(success, trials)
    # print(len(nvals), len(evals))
    # print(np.mean(nvals), np.mean(evals))


def analyse_coexpression_percomp_pairwise():
    """This uses interfaces, probably a fairer way to do it is to use all
    possible pairwise combinations, as I suspect this is skewed by the fact
    that there are many homomeric interfaces present.
    """
    with open('data/coexpression.tsv') as infile:
        data = [line.strip().split('\t') for line in infile]
    cdict = {}
    coex = namedtuple('Int', ['p1', 'p2', 'c1', 'c2', 'coex', 'int'])
    for line in data:
        if len(line) != 12:
            continue
        if line[0] not in cdict:
            cdict[line[0]] = [coex(line[4], line[7], line[6], line[9], line[11], line[10])]
        else:
            cdict[line[0]].append(coex(line[4], line[7], line[6], line[9], line[11], line[10]))
    trials = 0
    success = 0
    for pdb in cdict:
        comp = cdict[pdb]
        nvals = []
        evals = []
        for i in comp:
            if i.p1 == i.p2:
                continue
            if float(i.int) < 200:
                continue
            if i.c1 != 'E' and i.c2 != 'E':
                nvals.append(float(i.coex))
            elif i.c1 == 'E' or i.c2 == 'E':
                evals.append(float(i.coex))
            if len(nvals) > 0 and len(evals) > 0:
                trials += 1
                if np.mean(nvals) > np.mean(evals):
                    success += 1
    print(success, trials)

if __name__ == '__main__':
    # data = get_tissue_expression()
    # tissue_expression(data)
    # tissue_expression_percomp(data)
    # interfaces()
    # get_coexpression()
    analyse_coexpression_percomp_pairwise()