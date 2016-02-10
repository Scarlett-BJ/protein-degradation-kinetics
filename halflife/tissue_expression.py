#!/usr/bin/env python3

import ixntools as ix
from expression import paxdb, proteomicsdb, proteomicsdb_requests as req
from halflife import utils
from halflife import abundance as ab

import numpy as np
from scipy.stats import binom_test

import logging
import sys

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
log = logging.getLogger('tissue')

def protein_map(species):
    """Returns dictionary of uniprot accession codes to ensp IDs."""
    if species == 'mouse':
        filename = 'data/mouse_uniprot_ensp.txt'
    elif species == 'human':
        filename = 'data/human_uniprot_ensp.txt'
    with open(filename) as infile:
        data = [line.split() for line in infile]
    prot_map = {}
    for line in data:
        if line[0] not in prot_map:
            prot_map[line[0]] = [line[1]]
        else:
            prot_map[line[0]].append(line[1])
    return prot_map

def load_paxdb_data(species):
    if species == 'mouse':
        meta = paxdb.get_metadata('10090')
    elif species == 'human':
        meta = paxdb.get_metadata('9606')
    data = utils.load_ned_data(species)[1]
    pmap = protein_map(species)
    skip = ['CELL_LINE']
    for i in skip:
        meta.pop(i)
    abunds = [paxdb.Abundances(meta[t]['filename'], 0) for t in sorted(meta)]
    header = '\t'.join(['prot', '\t'.join(sorted(meta)), 'tcount', 'def'])
    return header, data, abunds, pmap

def load_proteomicsdb_data():
    data = utils.load_ned_data('human')[1]
    isoab = proteomicsdb.Abundances('human_protdb_isoforms_expression.txt')
    abunds = proteomicsdb.Abundances('trembl_tissues.txt')
    tissues = abunds.tissues
    header = '{0}\t{1}\t{2}\t{3}\n'.format('prot', '\t'.join(tissues),
                                           'tcount', 'def')
    with open('data/human_proteomicsdb_tissue_expression.txt', 'w') as outfile:
        outfile.write(header)
        for line in data:
            prot = line[0]
            if prot in abunds.proteins:
                expression = [abunds.expression(prot, tis) for tis in ttargets]
            elif prot in isoab.proteins:
                expression = [isoab.expression(prot, tis) for tis in ttargets]
            else:
                continue
            tcount = str(len(expression) - expression.count('NA'))
            newline = '\t'.join([prot] + expression + [tcount, line[-1]])
            outfile.write(newline + '\n')

def tcount_binomial_test():
    fname = 'data/abundance/human_proteomicsdb_tissue_expression.txt'
    with open(fname) as infile:
        data = [line.strip().split('\t') for line in infile]
    data = {line[0]: [int(line[-2]),  line[-1]] for line in data[1:]}
    core = ix.dbloader.LoadCorum('human', 'core')
    success, trials = 0, 0
    for struc in core.strucs:
        comp = core[struc].uniprot
        nvals, evals = [], []
        for sub in comp:
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

def tcount_binomial_test2():
    fname = 'data/abundance/human_proteomicsdb_tissue_expression.txt'
    with open(fname) as infile:
        data = [line.strip().split('\t') for line in infile]
    tdict = {line[0]: [int(line[-2]),  line[-1]] for line in data[1:]}
    abdict = ab.abundance_dict('human')
    for prot in tdata:
        print(prot, tdata[prot][1], tdata[prot][0], abdict[prot][0])


def main():
    tcount_binomial_test()



if __name__ == '__main__':
    main()
