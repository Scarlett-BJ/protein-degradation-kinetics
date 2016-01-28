#!/usr/bin/env python3

from ixntools import dbloader
from expression import coexpressdb
from numpy import mean, median
from scipy.stats import binom_test
from collections import namedtuple
import logging
import sys

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
log = logging.getLogger('coexpression')

###############################################################################

def load_ned_data(filename):
    """Loads processed decay data from Selbach group"""
    with open(filename) as infile:
        data = [line.strip().split('\t') for line in infile]
    header = data[0]
    data = data[1:]
    return header, data

def get_homologs():
    """Returns dictionary mapping mouse homologs entrez <-> uniprot."""
    homologs = {}
    with open('data/homology/corum_mouse_homologs.txt') as infile:
        data = [line.strip().split('\t') for line in infile]
    # data must be sorted in order of sequence identity (high first)
    for line in data:
        original = line[1].split('|')[1]
        uniprot = line[0]
        entrez = line[3]
        # proteins map 1 to 1
        if original not in homologs:
            homologs[original] = [entrez]
        # genes map 1 to multiple
        if entrez not in homologs:
            homologs[entrez] = [uniprot]
        elif uniprot not in homologs[entrez]:
            homologs[entrez].append(uniprot)
    return homologs

###############################################################################


class CoexpressTable(object):
    """Combines data from NED, CORUM and CoexpressDB."""

    def __init__(self, species, homologs = False):
        """Mouse specific if using homologs, else human or mouse."""
        if homologs == True:
            neds = load_ned_data('data/NED_mouse_Abund.txt')[1]
            self.homologs = get_homologs()
            self.corum = dbloader.LoadCorum(version='core')
            self.coex = coexpressdb.Coexpression('mouse')
            self.decay = {line[-2]: line[-3] for line in neds}
            self.species = 'mouse_homologs'
        else:
            neds = load_ned_data('data/NED_{0}_Abund.txt'.format(species))[1]
            self.homologs = False
            self.corum = dbloader.LoadCorum(species.title(), 'core')
            self.coex = coexpressdb.Coexpression(species)
            self.decay = {line[-2]: line[-3] for line in neds}
            self.species = species
        self.outdata = []

    def _homologise_complex(self):
        """Swaps corum uniprot id with mouse homolog entrez id."""
        assert self.homologs
        entrezids = []
        for sub in self.complex.uniprot:
            for u in sub:
                if u in self.homologs:
                    entrezids.append(self.homologs[u][0])
                    break
        return entrezids

    def _avg_coexpression(self, subunits):
        """Uses CoexpressDB to calculate average coexpression of subunits"""
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
                avg_coex[protein] = mean(avg_coex[protein])
        return avg_coex

    def _convert_avgcoex_keys(self, avg_coex):
        """Maps CORUM entrez ids to corresponding uniprot id."""
        conversions = {}
        entrez = [e for sublist in self.complex.entrez for e in sublist]
        uniprot = [u for sublist in self.complex.uniprot for u in sublist]
        if len(entrez) != len(uniprot):
            log.debug('number of entrez subunits != number of uniprot')
        for i in range(len(entrez)):
            conversions[entrez[i]] = uniprot[i]
        avg_coex = {conversions[key]: val for key, val in avg_coex.items()}
        return avg_coex

    def _convert_homolog_avgcoex_keys(self, avg_coex):
        """As above, but using homology mappings to get uniprot from entrez."""
        assert self.homologs
        for entrez in list(avg_coex):
            if set(self.homologs[entrez]).intersection(self.decay) == set():
                avg_coex.pop(entrez)
            else:
                for upr in self.homologs[entrez]:
                    if upr in self.decay:
                        avg_coex[upr] = avg_coex[entrez]
                        avg_coex.pop(entrez)
                        break
        return avg_coex

    def process_data(self):
        """Gets avg. coexpression and decay of each subunit in CORUM complexes.

        Carries out required processing needed to map subunits from corum
        complexes to homologs (if required), decay-type classifications and
        average coexpression within a complex. In doing so builds outdata, a
        list of lists, where each inner list contains a the following
        attributes for a single subunit:
            corumid, unique_subs, uniprot_id, avg_coexpression, decay, species
        """
        for struc in self.corum.strucs:
            self.complex = self.corum[struc]
            usubs = len(self.complex.uniprot)
            if self.homologs:
                entrezids = self._homologise_complex()
                avg_coex = self._avg_coexpression(entrezids)
                avg_coex = self._convert_homolog_avgcoex_keys(avg_coex)
            else:
                entrezids = [sub[0] for sub in self.complex.entrez]
                avg_coex = self._avg_coexpression(entrezids)
                avg_coex = self._convert_avgcoex_keys(avg_coex)
            for subunit in avg_coex:
                info = [struc, str(usubs), subunit, str(avg_coex[subunit]),
                        self.decay.get(subunit, 'NA'), self.species]
                self.outdata.append(info)

    def write_to_file(self, filename):
        """Writes outdata to specified filename."""
        if self.outdata == []:
            log.warning('No data to write!')
        with open(filename, 'w') as outfile:
            header = ['comp', 'usubs', 'uniprot', 'avgcoex', 'def', 'species']
            outfile.write('{0}\n'.format('\t'.join(header)))
            for line in self.outdata:
                outfile.write('{0}\n'.format('\t'.join(line)))


###############################################################################

def coexpression_binomial(filename):
    """Per complex binomial test for average subunit coexpression.

    Returns:
        trials - number of complexes tested
        success - times median coexpression score of NEDs was higher than EDs.
        pval - exact, two-sided probability of success given p=0.5 under H0.
    """
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
        if median(nvals) > median(evals):
            success += 1
    pval = binom_test(success, trials)
    print(success, trials, pval)
    return success, trials, pval


def main():
    # tab = CoexpressTable('mouse', homologs=True)
    # tab.process_data()
    # tab.write_to_file('data/coexpressdb_corum_mouse_homologs.tsv')
    analyse_corum_data('data/coexpressdb_corum_mouse_homologs.tsv')

    # tab = CoexpressTable('mouse')
    # tab.process_data()
    # tab.write_to_file('data/coexpressdb_corum_mouse.tsv')
    analyse_corum_data('data/coexpressdb_corum_human.tsv')

if __name__ == '__main__':
    main()