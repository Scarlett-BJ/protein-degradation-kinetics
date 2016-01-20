#!/usr/bin/env python3

import numpy as np
import ixntools as ix
from expression import paxdb, coexpressdb
from collections import namedtuple
from itertools import combinations
from scipy.stats import binom_test
import random
import re

def load_data(filename):
    with open(filename) as infile:
        data = [line.split() for line in infile]
    header = data[0]
    data = data[1:]
    return header, data

def mouse_homologs(rev=False):
    with open('data/homologs.txt') as infile:
        homologs = {l.split()[0]: l.split()[1] for l in infile}
        homologs.pop('mouse')
    if rev:
        homologs = {val: key for key, val in homologs.items()}
    return homologs

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

def process_corum_data(species, homologs=False):
    """Calculates average coexpression per complex for each subunit."""
    if homologs:
        homologs = mouse_homologs()
        core = ix.dbloader.LoadCorum('Human', 'core')
        species = 'mouse'
    else:
        core = ix.dbloader.LoadCorum(species.title(), 'core')
    coex = coexpressdb.Coexpression()
    if species == 'human':
        data = load_data('data/NED_human_Abund.txt')[1]
    elif species == 'mouse':
        data = load_data('data/NED_mouse_Abund.txt')[1]
    decay_defs = {line[-2]: line[-3] for line in data}

    def get_best_from_multiple_subs(subunits):
        """In cases where 2 or more subunits cannot be distinguished from one
        another, CORUM encloses those subunits within brackets. This function
        returns the subunit for which most data is available, prioritising data
        on coexpression over data on decay classification.
        """
        possibilities = []
        for sub in subunits:
            upr = match_entrez_to_uniprot(comp, sub)
            if homologs:
                if upr in homologs:
                    upr = homologs[upr]
            coexpression = avg_coex[sub]
            uprdef = decay_defs.get(upr, 'NA')
            possibilities.append((struc, sub, upr, coexpression,
                                  uprdef, species))
        possibilities.sort(key=lambda x: (x.count('NA'), [x[3]].count('NA')))
        if len(possibilities) == 0:
            return
        info = [str(item) for item in possibilities[0]]
        return info

    if homologs:
        species += '_homologs'
    ofname = 'data/coexpressdb_corum_{0}.tsv'.format(species)
    header = '\t'.join(['comp', 'entrez', 'uniprot',
                        'avg.coex', 'def', 'species']) + '\n'

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
                if info:
                    outfile.write('\t'.join(info) + '\n')


def process_corum_data2(species, homologs=False):
    if homologs:
        homologs = mouse_homologs(rev = True)
        core = ix.dbloader.LoadCorum('Human', 'core')
    else:
        core = ix.dbloader.LoadCorum(species.title(), 'core')
    header, data = load_data('data/NED_{0}_Abund.txt'.format(species))
    decay = {line[-2]: line[-3] for line in data}
    coex = coexpressdb.Coexpression()

    def get_best_from_multiple_subs(subunits):
        """In cases where 2 or more subunits cannot be distinguished from one
        another, CORUM encloses those subunits within brackets. This function
        returns the subunit for which most data is available, prioritising data
        on coexpression over data on decay classification.
        """
        possibilities = []
        for entrez in subunits:
            upr = match_entrez_to_uniprot(core[struc], entrez)
            coexpression = avg_coex[entrez]
            if homologs and upr in homologs:
                upr = homologs[upr]  # swaps human id to mouse, or does nothing
            uprdecay = decay.get(upr, 'NA')
            possibilities.append((struc, entrez, upr, coexpression,
                                  uprdecay, species))
        possibilities.sort(key=lambda x: (x.count('NA'), [x[3]].count('NA')))
        if len(possibilities) == 0:  # possibly unused?
            return
        info = [str(item) for item in possibilities[0]]  # picks best option
        return info

    if homologs:
        species += '_homologs'
    ofname = 'data/coexpressdb_corum_{0}.tsv'.format(species)
    header = '\t'.join(['comp', 'entrez', 'uniprot',
                        'avg.coex', 'def', 'species']) + '\n'
    with open(ofname, 'w') as outfile:
        outfile.write(header)
        for struc in core.strucs:
            entrez = set(e for sublist in core[struc].entrez for e in sublist)
            avg_coex = calculate_avg_coexpression(coex, entrez)
            for subunit in core[struc].entrez:
                # Get coexpression data, then decay data. Else 'NA'
                info = get_best_from_multiple_subs(subunit)
                if info:
                    outfile.write('\t'.join(info) + '\n')

###############################################################################

def analyse_corum_data(filename, homologs=False):
    """Per complex binomial test for average subunit coexpression."""
    if homologs:
        homologs = mouse_homologs(rev = True)
        data = load_data('data/NED_mouse_Abund.txt')[-1]
        decay = {line[-2]: line[-3] for line in data[1:]}
    with open(filename) as infile:
        data = [line.strip().split('\t') for line in infile][1:]
    strucs = {line[0]: [] for line in data}
    for line in data:
        info = tuple(line[-3:-1])
        print(info)
        if homologs:
            if line[2] in homologs:
                info = (line[-3], decay.get(homologs[line[2]], 'NA'))
            else:
                continue
        strucs[line[0]].append(info)
    success, trials = 0, 0
    for struc in strucs:
        nvals, evals = [], []
        if len(strucs[struc]) <= 2:
            continue  # Because coex of subs in dimer are identical
        for subunit in strucs[struc]:
            if subunit[0] == 'NA':
                continue
            # Append avg coexpressions.
            if subunit[1] == 'NED':
                nvals.append(float(subunit[0]))
            elif subunit[1] == 'ED':
                evals.append(float(subunit[0]))
        if len(nvals) == 0 or len(evals) == 0 or nvals == evals:
            continue
        trials += 1
        if np.mean(nvals) > np.mean(evals):
            success += 1
    print(success, trials, binom_test(success, trials))

def main():
    # process_corum_data2('mouse', homologs=True)
    # analyse_corum_data('data/coexpressdb_corum_mouse_homologs.tsv')
    # process_corum_data2('human')
    analyse_corum_data('data/coexpressdb_corum_human.tsv', homologs=True)
    # process_corum_data2('mouse')
    # analyse_corum_data('data/coexpressdb_corum_mouse.tsv')

if __name__ == '__main__':
    main()


    # You need mouse entrez ids. Swap the subunits for mouse subunits and calculate MOUSE coexpressions. Sort your head out you silly twat!
