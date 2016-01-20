#!/usr/bin/env python3

import ixntools as ix
import numpy as np
from scipy.stats import binom_test
from collections import namedtuple

def load_NED_data(filename, sep='\t'):
    """Returns header and data from NED files."""
    with open(filename) as infile:
        data = [line.strip().split(sep) for line in infile]
    header = data[0]
    data = data[1:]
    return header, data

def map_mouse_homologs():
    with open('data/homologs.txt') as infile:
        homologs = {l.split()[0]: l.split()[1] for l in infile}
        homologs.pop('mouse')
    return homologs

def get_old_mouse_abundances():
    header, data = load_NED_data('data/NED_mouse_original.txt', ',')
    abunds = {line[0]: float(line[-1]) for line in data}
    return abunds

def get_coexpression():
    with open('data/coexpressdb_corum_combined.tsv') as infile:
        coexpress = {l.split()[2]: float(l.split()[-3]) for l in infile
                     if 'mouse' in l and l.split()[-3] != 'NA'}
    return coexpress

def tcount_binomial_test():
    fname = 'data/abundance/human_proteomicsdb_tissue_expression.txt'
    with open(fname) as infile:
        data = [line.strip().split('\t') for line in infile]
    data = {line[0]: [int(line[-2]),  line[-1]] for line in data[1:]}
    core = ix.dbloader.LoadCorum('Human', 'core')
    success, trials = 0, 0
    for struc in core.strucs:
        nvals, evals = [], []
        subunits = core[struc].uniprot
        for sub in subunits:
            if len(sub) == 1 and sub[0] in data:
                if data[sub[0]][1] == 'NED':
                    nvals.append(data[sub[0]][0])
                elif data[sub[0]][1] == 'ED':
                    evals.append(data[sub[0]][0])
            else:
                for p in sub:
                    if p in data and data[p][1] == 'NED':
                        nvals.append(data[p][0])
                        break
                    elif p in data and data[p][1] == 'ED':
                        evals.append(data[p][0])
                        break
        if len(evals) > 0 and len(nvals) > 0:
            trials += 1
            if np.mean(nvals) >= np.mean(evals):
                success += 1
    print(success, trials, binom_test(success, trials))

def abundance_binomial_test():
    core = ix.dbloader.LoadCorum('Human', 'core')
    header, data = load_NED_data('data/NED_mouse.txt')
    homologs = map_mouse_homologs()
    Prot = namedtuple('Prot', ['abundance', 'decay'])
    prots = {homologs[line[0]]: Prot(float(line[-1]), line[-2])
             for line in data if line[0] in homologs}
    success, trials = 0, 0
    for struc in core.strucs:
        nvals, evals = [], []
        subunits = core[struc].uniprot
        for sub in subunits:
            # Deals with ambiguous (bracketed) subunits
            if len(sub) == 1 and sub[0] in prots:
                if prots[sub[0]].decay == 'NED':
                    nvals.append(prots[sub[0]].abundance)
                elif prots[sub[0]].decay == 'ED':
                    evals.append(prots[sub[0]].abundance)
            else:
                for p in sub:
                    if p in prots and prots[p].decay == 'NED':
                        nvals.append(prots[p].abundance)
                        break
                    elif p in prots and prots[p].decay == 'ED':
                        evals.append(prots[p].abundance)
                        break
        if len(evals) > 0 and len(nvals) > 0:
            trials += 1
            if np.mean(nvals) >= np.mean(evals):
                success += 1
    print(success, trials, binom_test(success, trials))

def handle_corum_subs(subunits, NED_prots):
    final_subs = []
    for sublist in subunits:
        for prot in sublist:
                if prot in NED_prots:
                    final_subs.append(prot)
                    break
    return final_subs


def coexpression_corum_data(filename):
    """Per complex binomial test for average subunit coexpression. Different
    format to abundance and tissue count tests - takes its own dedicated file
    in which each complex is written in its entirety"""
    with open(filename) as infile:
        data = [line.strip().split('\t') for line in infile if 'mouse' in line]
    strucs = {line[0]: [] for line in data}
    for line in data:
        strucs[line[0]].append(tuple(line[-3:-1]))
    success = 0
    trials = 0
    strucs_for_plot = []
    for struc in strucs:
        nvals = []
        evals = []
        comp = strucs[struc]
        if len(comp) <= 2:
            continue  # Because coex of subs in dimer are identical
        for subunit in comp:
            if subunit[0] == 'NA':
                continue
            if subunit[1] == 'NED':
                nvals.append(float(subunit[0]))
            elif subunit[1] == 'ED':
                evals.append(float(subunit[0]))
        if len(nvals) == 0 or len(evals) == 0 or nvals == evals:
            continue
        strucs_for_plot.append(struc)
        trials += 1
        if np.mean(nvals) > np.mean(evals):
            success += 1
    print(success, trials, binom_test(success, trials))

def generic():
    core = ix.dbloader.LoadCorum('Mouse', 'core')
    header, data = load_NED_data('data/NED_mouse_update.txt')
    ned_prots = set(line[0] for line in data)
    homolog = map_mouse_homologs()
    abund = get_old_mouse_abundances()
    decay = {line[0]: line[-1] for line in data}
    coexpress = get_coexpression()
    success, trials = 0, 0
    for struc in core.strucs:
        nvals, evals = [], []
        subunits = handle_corum_subs(core[struc].uniprot, ned_prots)
        for sub in subunits:
            if decay.get(sub, None) == 'NED':
                nvals.append(coexpress.get(sub, None))
            elif decay.get(sub, None) == 'ED':
                evals.append(coexpress.get(sub, None))
        if None in nvals:
            nvals.remove(None)
        if None in evals:
            evals.remove(None)
        if len(nvals) == 0 or len(evals) == 0 or nvals == evals:
            continue
        trials += 1
        if np.mean(nvals) > np.mean(evals):
            success += 1
    print(success, trials, binom_test(success, trials))

def main():
    abundance_binomial_test()
    tcount_binomial_test()
    coexpression_corum_data('data/coexpressdb_corum_combined.tsv')

if __name__ == '__main__':
    # main()
    generic()
    # print(get_coexpression())
