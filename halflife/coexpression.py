#!/usr/bin/env python3

from ixntools import dbloader
from expression import coexpressdb
from halflife import utils
from numpy import mean, median
from scipy.stats import binom_test
import logging
import sys

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
log = logging.getLogger('coexpression')


class CoexpressTable(object):
    """Combines data from NED, CORUM and CoexpressDB."""

    def __init__(self, species, homologs = False):
        """Mouse specific if using homologs, else human or mouse."""
        neds = utils.load_ned_data(species)[1]
        if homologs == True:
            if species != 'mouse':
                raise ValueError('homologs only implemented for mouse')
            self._homologs = utils.get_homologs()
            self._corum = dbloader.LoadCorum(version='core')
        else:
            self._homologs = False
            self._corum = dbloader.LoadCorum(species, 'core')
        self._coex = coexpressdb.Coexpression(species)
        self._decay = {line[-2]: line[-3] for line in neds}
        self._species = species
        self._outdata = []

    def _homologise_complex(self):
        """Swaps corum uniprot id with mouse homolog entrez id."""
        assert self._homologs
        entrezids = []
        for sub in self._complex.uniprot:
            for u in sub:
                if u in self._homologs:
                    entrezids.append(self._homologs[u][0])
                    break
        return entrezids

    def _avg_coexpression(self, subunits):
        """Uses CoexpressDB to calculate average coexpression of subunits"""
        avg_coex = {protein: [] for protein in subunits}
        for protein in avg_coex:
            for protein2 in avg_coex:
                if protein == protein2:
                    continue
                cor = self._coex.get_coexpression(protein, protein2)
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
        entrez = [e for sublist in self._complex.entrez for e in sublist]
        uniprot = [u for sublist in self._complex.uniprot for u in sublist]
        if len(entrez) != len(uniprot):
            log.debug('number of entrez subunits != number of uniprot')
        for i in range(len(entrez)):
            conversions[entrez[i]] = uniprot[i]
        avg_coex2 = {conversions[key]: val for key, val in avg_coex.items()}
        return avg_coex2

    def _convert_homolog_avgcoex_keys(self, avg_coex):
        """As above, but using homology mappings to get uniprot from entrez."""
        assert self._homologs
        avg_coex2 = {}
        for entrez in list(avg_coex):
            for upr in self._homologs[entrez]:
                if upr in self._decay:
                    avg_coex2[upr] = avg_coex[entrez]
                    break
        return avg_coex2

    def process_data(self):
        """Gets avg. coexpression and decay of each subunit in CORUM complexes.

        Carries out required processing needed to map subunits from corum
        complexes to homologs (if required), decay-type classifications and
        average coexpression within a complex. In doing so builds outdata, a
        list of lists, where each inner list contains a the following
        attributes for a single subunit:
            corumid, unique_subs, uniprot_id, avg_coexpression, decay, species
        """
        for struc in self._corum.strucs:
            self._complex = self._corum[struc]
            usubs = len(self._complex.uniprot)
            if self._homologs:
                entrezids = self._homologise_complex()
                avg_coex = self._avg_coexpression(entrezids)
                avg_coex = self._convert_homolog_avgcoex_keys(avg_coex)
            else:
                entrezids = [sub[0] for sub in self._complex.entrez]
                avg_coex = self._avg_coexpression(entrezids)
                avg_coex = self._convert_avgcoex_keys(avg_coex)
            for subunit in avg_coex:
                info = [struc, str(usubs), subunit, str(avg_coex[subunit]),
                        self._decay.get(subunit, 'NA'), self._species]
                self._outdata.append(info)


    def write_to_file(self, filename):
        """Writes outdata to specified filename."""
        if self._outdata == []:
            log.warning('No data to write!')
        with open(filename, 'w') as outfile:
            header = ['comp', 'usubs', 'uniprot', 'avgcoex', 'def', 'species']
            outfile.write('{0}\n'.format('\t'.join(header)))
            for line in self._outdata:
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

def filter_topseqid_struc_coexpression():
    with open('data/coexpression/struc_coex.txt') as infile:
        new_data = [infile.readline()]
        sdict = []
        for line in  infile:
            sline = line.strip().split()
            if sline[2] not in sdict:
                sdict.append(sline[2])
                new_data.append(line)
    for line in new_data:
        print(line.strip())

def main():
    # # Mouse Homologs
    # tab = CoexpressTable('mouse', homologs=True)
    # tab.process_data()
    # tab.write_to_file('data/coexpressdb_corum_mouse_homologs.tsv')
    coexpression_binomial('data/coexpressdb_corum_mouse_homologs.tsv')
    # # Human complexes
    # tab = CoexpressTable('human')
    # tab.process_data()
    # tab.write_to_file('data/coexpressdb_corum_human.tsv')
    coexpression_binomial('data/coexpressdb_corum_human.tsv')
    # # Mouse complexes
    # tab = CoexpressTable('mouse')
    # tab.process_data()
    # tab.write_to_file('data/coexpressdb_corum_mouse.tsv')
    coexpression_binomial('data/coexpressdb_corum_mouse.tsv')

if __name__ == '__main__':
    main()
