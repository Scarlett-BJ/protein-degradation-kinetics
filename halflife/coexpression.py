#!/usr/bin/env python3

import numpy as np
import ixntools as ix
from expression import paxdb, coexpressdb
from collections import namedtuple
from itertools import combinations
from scipy.stats import binom_test
import random
import re

# def load_mouse_prot_map():
#     data_dict = {}
#     with open('data/mouse_genes_ned.txt') as infile:
#         data = [line.strip().split(',') for line in infile.readlines()]
#     for line in data:
#         data_dict[line[0]] = {'ensg': line[2], 'gene': line[1],
#                               'ensp': line[3]}
#     return data_dict

def load_data(filename):
    with open(filename) as infile:
        data = [line.split() for line in infile]
    header = data[0]
    data = data[1:]
    return header, data

# def get_complexes(header, data):
#     """Returns a dictionary of complexes. Each complex contains a list of
#     namedtuples with information about the complex and proteins.

#     DEPRECATED
#     """
#     [str(i.strip('.')) for i in header[1:]]
#     Info = namedtuple('Header', [str(i.strip('.')) for i in header[1:]])
#     cdict = {}
#     for line in data:
#         if line[0] not in cdict:
#             cdict[line[0]] = [Info(*line[1:])]
#         else:
#             cdict[line[0]].append(Info(*line[1:]))
#     return cdict

# def get_coexpression():
#     """DEPRECATED"""
#     gene_map = load_mouse_prot_map()
#     coex = coexpressdb.Coexpression()
#     data = load_pairwise()
#     for line in data[1:]:
#         if line[4] == 'NA' or line[7] == 'NA':
#             print('\t'.join(line))
#             continue
#         ensp1 = gene_map[line[4]]['ensp']
#         ensp2 = gene_map[line[7]]['ensp']
#         if ensp1 not in coex.entrez or ensp2 not in coex.entrez:
#             print('\t'.join(line))
#             continue
#         c = coex.get_coexpression(ensp1, ensp2)
#         line.append(str(c))
#         print('\t'.join(line))

def calculate_avg_coexpression(coex, members):
    """Calculates the average coexpression for each protein in the complex.

    Correlation is taken from coexpressdb.jp data, which uses pearson
    correlations. Missing data is excluded from averages, but there is no limit
    to how few proteins will be used to calculate averages.

    Arguments:
        coex - the coexpressdb class instance containing coex info.
        members - a list of all proteins for which coexpression is required.

    Returns:
        dictionary of protein names, with values being average coexpression,
        or 'NA' if no avg. coexpression score could be calculated
    """
    avg_coex = {protein: [] for protein in members}
    for protein in avg_coex:
        for protein2 in avg_coex:
            if protein == protein2:
                continue
            cor = coex.get_coexpression(protein, protein2)
            if cor == None:
                continue
            avg_coex[protein].append(cor)
        if len(avg_coex[protein]) < 1:
            avg_coex[protein] = 'NA'
        else:
            avg_coex[protein] = np.mean(avg_coex[protein])
    return avg_coex

def match_entrez_to_uniprot(comp, entrezid):
    """Maps entrez ids to corresponding uniprot.

    Within CORUM files, subunits are given as a list of entrez and uniprot,
    with the index of each protein being equivalent. i.e. entrez[0] == upr[0]
    """
    for i in range(len(comp.entrez)):
        if entrezid in comp.entrez[i]:
            return comp.uniprot[i][comp.entrez[i].index(entrezid)]

def process_corum_data(species):
    """Calculates average coexpression per complex for each subunit."""
    core = ix.dbloader.LoadCorum(species.title(), 'core')
    coex = coexpressdb.Coexpression()
    if species == 'human':
        data = load_data('data/NED_human.txt')[1]
    elif species == 'mouse':
        data = load_data('data/NED_mouse_update.txt')[1]
    decay_defs = {line[0]: line[-1] for line in data}

    def get_best_from_multiple_subs(subunits):
        """In cases where 2 or more subunits cannot be distinguished from one
        another, CORUM encloses those subunits within brackets. This function
        returns the subunit for which most data is available, prioritising data
        on coexpression over data on decay classification.
        """
        possibilities = []
        for sub in subunits:
            upr = match_entrez_to_uniprot(comp, sub)
            coexpression = avg_coex[sub]
            uprdef = decay_defs.get(upr, 'NA')
            possibilities.append((struc, sub, upr, coexpression,
                                  uprdef, species))
        possibilities.sort(key=lambda x: (x.count('NA'), [x[3]].count('NA')))
        info = [str(item) for item in possibilities[0]]
        return info

    ofname = 'data/coexpressdb_corum_{0}.tsv'.format(species)
    header = '\t'.join(['comp', 'entrez', 'uniprot', 'avg.coex', 'def']) + '\n'
    with open(ofname, 'w') as outfile:
        outfile.write(header)
        for struc in core.strucs:
            comp = core[struc]
            uniprot = comp.uniprot
            entrez_list = [p for sublist in comp.entrez for p in sublist]
            avg_coex = calculate_avg_coexpression(coex, entrez_list)
            for subunit in comp.entrez:
                # Get coexpression data, then decay data. Else 'NA'
                if len(subunit) == 1:
                    upr = match_entrez_to_uniprot(comp, subunit[0])
                    coexpression = avg_coex[subunit[0]]
                    uprdef = decay_defs.get(upr, 'NA')
                    info = (struc, subunit[0], upr,
                            coexpression, uprdef, species)
                    info = [str(item) for item in info]
                else:
                    info = get_best_from_multiple_subs(subunit)
                outfile.write('\t'.join(info) + '\n')

def analyse_corum_data(filename, rand=False):
    """Predominantly bionmial test"""
    with open(filename) as infile:
        data = [line.strip().split('\t') for line in infile]
    strucs = {line[0]: [] for line in data}
    for line in data:
        strucs[line[0]].append(tuple(line[-3:-1]))
    success = 0
    trials = 0
    for struc in strucs:
        nvals = []
        evals = []
        comp = strucs[struc]
        # Shuffle if required for testing of null.
        if rand == True:
            def_dist = [subunit[1] for subunit in comp]
            random.shuffle(def_dist)
            for i in range(len(def_dist)):
                comp[i] = (comp[i][0], def_dist[i])
        if len(comp) <= 2:
            continue  # Because coex of subs in dimer are identical
        for subunit in comp:
            if subunit[0] == 'NA':
                continue
            if subunit[1] == 'UN':
                nvals.append(float(subunit[0]))
            elif subunit[1] == 'ED':
                evals.append(float(subunit[0]))
        if len(nvals) == 0 or len(evals) == 0 or nvals == evals:
            continue
        trials += 1
        if np.mean(nvals) > np.mean(evals):
            success += 1
    print(success, trials, binom_test(success, trials))
    all_neds = [float(i[-3]) for i in data if i[-2] == 'NED' if i[-3] != 'NA']
    all_eds = [float(i[-3]) for i in data if i[-2] == 'ED' if i[-3] != 'NA']
    all_uns = [float(i[-3]) for i in data if i[-2] == 'UN' if i[-3] != 'NA']
    all_nas = [float(i[-3]) for i in data if i[-2] == 'NA' if i[-3] != 'NA']
    print(np.var(all_neds))
    print(np.median(all_eds))
    print(np.median(all_uns))
    print(np.var(all_nas))
    print()


if __name__ == '__main__':
    # analyse_corum_data('data/coexpressdb_corum_human.tsv', rand=True)
    # analyse_corum_data('data/coexpressdb_corum_mouse.tsv', rand=True)
    analyse_corum_data('data/coexpressdb_corum_combined.tsv')
