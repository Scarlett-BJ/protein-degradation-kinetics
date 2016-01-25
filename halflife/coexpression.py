#!/usr/bin/env python3

import numpy as np
import ixntools as ix
from expression import paxdb, coexpressdb
from collections import namedtuple
from itertools import combinations
from scipy.stats import binom_test
import random
import re

def load_NED_data(filename):
    with open(filename) as infile:
        data = [line.split() for line in infile]
    header = data[0]
    data = data[1:]
    return header, data

def mouse_homologs(rev=False):
    upr_homologs, entrez_homologs = {}, {}
    with open('data/homologs.txt') as infile:
        data = [line.split() for line in infile]
    for line in data[1:]:
        # if line[0] not in upr_homologs and line[1] not in entrez_homologs:
            upr_homologs[line[0]] = line[2]
            entrez_homologs[line[1]] = line[3]
    if rev:
        upr_homologs = {val: key for key, val in upr_homologs.items()}
        entrez_homologs = {val: key for key, val in entrez_homologs.items()}
    return upr_homologs, entrez_homologs


class CoexpressionCalculator(object):

    """Methods for calculating and tabulating average coexpression."""

    def __init__(self, species, homologs=False):
        if homologs:
            self.upr_homs, self.entrez_homs = mouse_homologs(rev = True)
            # Homologous species complexes, i.e human if using mouse homologs
            self.core = ix.dbloader.LoadCorum('Human', 'core')
        else:
            self.core = ix.dbloader.LoadCorum(species.title(), 'core')
        header, data = load_NED_data('data/NED_{0}_Abund.txt'.format(species))
        # True species decay, i.e. mouse if using mouse homologs
        self.decay = {line[-2]: line[-3] for line in data}
        self.coex = coexpressdb.Coexpression()
        self.homologs = homologs
        self.species = species

    def tabulate_results(self):
        self.tabulated = []
        for struc in sorted(self.core.strucs):
            self.struc = struc
            # Gets all possible subunit entrez_ids
            entrez = [e for sublst in self.core[struc].entrez for e in sublst]
            avg_coex = self._calculate_avg_coexpression(entrez)
            for subunit in self.core[struc].entrez:
                # Get coexpression data, then decay data. Else 'NA'
                info = self._get_best_from_multiple_subs(subunit, avg_coex)
                self.tabulated.append(info)

    def _calculate_avg_coexpression(self, subunits):
        """Calculates the average coexpression for each protein in the complex.

        Correlation is taken from coexpressdb.jp data, which uses pearson
        correlations. Missing data is excluded from averages, but there is no
        limit to how few proteins will be used to calculate averages.

        Arguments:
            coex - the coexpressdb class instance containing coex info.
            subunits - a list of all proteins for which coexpression is
            required.

        Returns:
            dictionary of protein names, with values being average
            coexpression, or 'NA' if no avg. coexpression score could be
            calculated.
        """
        if self.homologs:
            # Convert human entrez ids to mouse ids
            subunits = [self.entrez_homs[s] for s in subunits
                        if s in self.entrez_homs]
        avg_coex = {protein: [] for protein in subunits}
        for protein in avg_coex:
            for protein2 in avg_coex:
                if protein == protein2:
                    continue
                cor = self.coex.get_coexpression(protein, protein2)
                if cor == None:
                    continue
                avg_coex[protein].append(cor)
            if len(avg_coex[protein]) < 1:
                avg_coex[protein] = 'NA'
            else:
                avg_coex[protein] = np.mean(avg_coex[protein])
        avg_coex['NA'] = 'NA'
        return avg_coex

    def _match_entrez_to_uniprot(self, struc, entrezid):
        """Maps entrez ids to corresponding uniprot.

        Within CORUM files, subunits are given as a list of entrez and uniprot,
        with the ind of each protein being equivalent. i.e. entrez[0] == upr[0]
        """
        comp = self.core[struc]
        for i in range(len(comp.entrez)):
            if entrezid in comp.entrez[i]:
                return comp.uniprot[i][comp.entrez[i].index(entrezid)]

    def _get_best_from_multiple_subs(self, subunits, avg_coex):
        """In cases where 2 or more subunits cannot be distinguished from one
        another, CORUM encloses those subunits within brackets. This function
        returns the subunit for which most data is available, prioritising data
        on coexpression over data on decay classification.
        """
        possibilities = []
        for entrezid in subunits:
            upr = self._match_entrez_to_uniprot(self.struc, entrezid)
            if self.homologs:
                # Swap entrez and uniprot ids
                entrezid = self.entrez_homs.get(entrezid, 'NA')
                upr = self.upr_homs.get(upr, 'NA')
            coexpression = avg_coex[entrezid]
            decay = self.decay.get(upr, 'NA')
            possibilities.append((self.struc, entrezid, upr, coexpression,
                                  decay, self.species))
        possibilities.sort(key=lambda x: (x.count('NA'), [x[3]].count('NA')))
        if len(possibilities) == 0:  # possibly unused?
            return
        info = [str(item) for item in possibilities[0]]  # picks best option
        return info

    def write_to_file(self):
        if self.homologs:
            self.species += '_homologs'
        ofname = 'data/coexpressdb_corum_{0}.tsv'.format(self.species)
        header = '\t'.join(['comp', 'entrez', 'uniprot',
                            'avg.coex', 'def', 'species']) + '\n'
        with open(ofname, 'w') as outfile:
            outfile.write(header)
            for line in self.tabulated:
                outfile.write('\t'.join(line) + '\n')


###############################################################################

def analyse_corum_data(filename):
    """Per complex binomial test for average subunit coexpression."""
    with open(filename) as infile:
        data = [line.strip().split('\t') for line in infile][1:]
    strucs = {line[0]: [] for line in data}
    for line in data:
        info = tuple(line[-3:-1])
        if info == ('NA', 'NA'):
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
    calc = CoexpressionCalculator('mouse', homologs=True)
    calc.tabulate_results()
    calc.write_to_file()
    analyse_corum_data('data/coexpressdb_corum_mouse_homologs.tsv')

if __name__ == '__main__':
    main()


    # You need mouse entrez ids. Swap the subunits for mouse subunits and calculate MOUSE coexpressions. Sort your head out you silly twat!
