#!/usr/bin/env python3

from ixntools import dbloader
from halflife import utils
from numpy import mean, log
from scipy.stats import binom_test
from scipy.stats.mstats import gmean
import logging
import sys

# logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
# log = logging.getLogger('abundance')

def abundance_dict(species):
    """Returns dict - uniprot => (decay, abundance)"""
    header, data = utils.load_ned_data(species)
    pi = header.index('proID')
    ai = header.index('rel.abun')
    di = header.index('def')
    return {line[pi]: (float(line[ai]), line[di]) for line in data}

def maximise_subunit_data(comp, abdict, homologs):
    nvals, evals, uvals = [], [], []
    for sub in comp:
        decay = None
        # Pick one from bracketed subs, break as soon as data available
        for poss_sub in sub:
            if homologs:
                if poss_sub in homologs and homologs[poss_sub] in abdict:
                    abund, decay = abdict[homologs[poss_sub]]
                    break  # breaks out of possible subunits loop i.e. brackets
            else:
                if poss_sub in abdict:
                    abund, decay = abdict[poss_sub]
                    break
        if decay == None:
            continue
        if decay == 'NED':
            nvals.append(abund)
        elif decay == 'ED':
            evals.append(abund)
        elif decay == 'UN':
            uvals.append(abund)
    return nvals, evals, uvals

def abundance_binomial(species, homologs=False):
    abdict = abundance_dict(species)
    if homologs:
        homologs = utils.get_uniprot_homologs()
        corum = dbloader.LoadCorum(version='core')
    else:
        # Important to specify species to avoid homologs
        corum = dbloader.LoadCorum(species, 'core')
    success, trials = 0, 0
    for struc in corum.strucs:
        comp = corum[struc].uniprot
        nvals, evals, uvals = maximise_subunit_data(comp, abdict, homologs)
        # At least 1 NED protein and 1 ED protein must be present
        if len(nvals) == 0 or len(evals) == 0:
            continue
        trials += 1
        if gmean(nvals) > gmean(evals):
            success += 1
    pval = binom_test(success, trials)
    print(success, trials, pval)
    return success, trials, pval

def main():
    abundance_binomial('mouse', homologs=True)
    abundance_binomial('human')
    abundance_binomial('mouse')

if __name__ == '__main__':
    main()
