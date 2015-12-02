#!/usr/bin/env python3

import numpy as np
from expression import paxdb, coexpressdb
from collections import namedtuple
from itertools import combinations
from scipy.stats import binom_test

def load_data(filename):
    with open(filename) as infile:
        data = [line.strip().split('\t') for line in infile]
    header = data[0]
    data = data[1:]
    return header, data

def load_mouse_prot_map():
    data_dict = {}
    with open('data/mouse_genes_ned.txt') as infile:
        data = [line.strip().split(',') for line in infile.readlines()]
    for line in data:
        data_dict[line[0]] = {'ensg': line[2], 'gene': line[1],
                              'ensp': line[3]}
    return data_dict

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

def relative_abundance():
    """Within complexes, are N proteins more abundant?"""
    header, data = load_data('data/protwise_trimers.tsv')
    cdict = get_complexes(header, data)
    trials = 0
    success = 0
    avg_distance = []
    for pdb in cdict:
        comp = cdict[pdb]
        nvals = [float(p.ab) for p in comp if p.dcls == 'N' or p.dcls == 'U']
        evals = [float(p.ab) for p in comp if p.dcls == 'E']
        if len(nvals) != 0 and len(evals) != 0:
            trials += 1
            emean = np.mean(evals)
            nmean = np.mean(nvals)
            avg_distance.append(nmean-emean)
            print(nmean, emean)
            if nmean > emean:
                success += 1
    print(success, trials, binom_test(success, trials), np.mean(avg_distance))

def get_tissue_expression():
    gene_map = load_mouse_prot_map()
    header, data = load_data('data/protwise.tsv')
    meta = paxdb.get_metadata('10090')
    abundances = []
    for i in meta:
        filename = meta[i]['filename']
        abundances.append(paxdb.Abundances(filename, 100))
    for line in data:
        ensp = gene_map[line[4]]['ensp']
        tcount = 0
        for tissue in abundances:
            if ensp in tissue.members:
                tcount += 1
        line.append(tcount)
    return data

def tissue_expression_percomp(data):
    """Are N proteins more widely expressed than E proteins"""
    cdict = {}
    avg_distance = []
    for line in data:
        if line[0] not in cdict:
            cdict[line[0]] = [(line[7], line[8])]
        else:
            cdict[line[0]].append((line[7], line[8]))
    trials = 0
    success = 0
    for pdb in cdict:
        comp = cdict[pdb]
        evals = []
        nvals = []
        for prot in comp:
            if prot[0] == 'E':
                evals.append(prot[1])
            elif prot[0] == 'N' or prot[0] == 'U':
                nvals.append(prot[1])
        if len(evals) > 0 and len(nvals) > 0:
            trials += 1
            emean = np.mean(evals)
            nmean = np.mean(nvals)
            avg_distance.append(((nmean-emean)**2)**0.5)
            if np.mean(nvals) >= np.mean(evals):
                success += 1
    print(success, trials, binom_test(success, trials))

def tissue_expression(data):
    """Are N proteins more widely expressed than E proteins"""
    nvals = []
    evals = []
    uvals = []
    for line in data:
        if line[7] == 'E':
            evals.append(line[8])
        elif line[7] == 'N':
            nvals.append(line[8])
        elif line[7] == 'U':
            uvals.append(line[8])
    print(np.mean(nvals), np.mean(uvals), np.mean(evals))

if __name__ == '__main__':
    # relative_abundance()
    data = get_tissue_expression()
    tissue_expression_percomp(data)
    tissue_expression(data)