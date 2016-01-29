#!/usr/bin/env python3

from ixntools import dbloader
from halflife import utils
from numpy import mean
from scipy.stats import binom_test
import logging
import sys

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
log = logging.getLogger('abundance')

def abundance_dict(species):
    header, data = utils.load_ned_data(species)
    pi = header.index('proID')
    ai = header.index('Rel. abund')
    di = header.index('def')
    return {line[pi]: (float(line[ai]), line[di]) for line in data}

def abundance_binomial(species, homologs=False):
    adict = abundance_dict(species)
    if homologs:
        homologs = utils.get_uniprot_homologs()
    corum = dbloader.LoadCorum(version='core')
    success, trials = 0, 0
    for struc in corum.strucs:
        nvals, evals = [], []
        comp = corum[struc]
        for p in comp.uniprot:
            if homologs:
                if p[0] not in homologs or homologs[p[0]] not in adict:
                    continue
                abund, decay = adict[homologs[p[0]]]
            elif p[0] in adict:
                abund, decay = adict[p[0]]
            else:
                continue
            if decay == 'NED':
                nvals.append(abund)
            elif decay == 'ED':
                evals.append(abund)
        if len(nvals) == 0 or len(evals) == 0:
            continue
        trials += 1
        if mean(nvals) > mean(evals):
            success += 1
    pval = binom_test(success, trials)
    print(success, trials, pval)
    return success, trials, pval

def main():
    abundance_binomial('mouse', homologs=True)
    abundance_binomial('human')

if __name__ == '__main__':
    main()
