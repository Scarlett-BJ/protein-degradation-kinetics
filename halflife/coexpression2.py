#!/usr/bin/env python3

from ixntools import dbloader
from expression import coexpressdb
from numpy import mean, median
from scipy.stats import binom_test
from collections import namedtuple
import logging
import sys

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
log = logging.getLogger()

###############################################################################

def load_ned_data(filename):
    with open(filename) as infile:
        data = [line.split() for line in infile]
    header = data[0]
    data = data[1:]
    return header, data

def load_corum_data(species):
    pass

def load_coexpression_data():
    pass

def get_homologs():
    entrez_homologs = {}
    uniprot_homologs = {}
    with open('data/corum_homologs.txt') as infile:
        data = [line.strip().split('\t') for line in infile]
    for line in data[1:]:
        human_entrez = line[0]
        mouse_entrez = line[1]
        if len(line) == 3:
            mouse_uniprot = line[2]
            if human_entrez not in uniprot_homologs:
                uniprot_homologs[human_entrez] = [mouse_uniprot]
            else:
                uniprot_homologs[human_entrez].append(mouse_uniprot)
        else:
            if human_entrez not in entrez_homologs:
                entrez_homologs[human_entrez] = [mouse_entrez]
            else:
                entrez_homologs[human_entrez].append(mouse_entrez)
    for key, val in entrez_homologs.items():
        entrez_homologs[key] = tuple(set(val))
    for key, val in uniprot_homologs.items():
        uniprot_homologs[key] = tuple(set(val))
    return entrez_homologs, uniprot_homologs

###############################################################################

class CoexpressCalculator(object):

    def __init__(self, species, homologs = False):
        if homologs == True:
            assert species != 'human'
            neddata = load_ned_data('data/NED_mouse_Abund.txt')[1]
            self.species = 'mouse_homologs'
            self.homolog_status = True
            self.entrezhoms, self.unihoms = get_homologs() # human to mouse
            self.corum = dbloader.LoadCorum('Human', 'core')
            self.coex = coexpressdb.Coexpression('mouse')
            self.decay = {line[-2]: line[-3] for line in neddata}
        else:
            filename = 'data/NED_{0}_Abund.txt'.format(species)
            neddata = load_ned_data(filename)[1]
            self.homolog_status = False
            self.corum = dbloader.LoadCorum(species.title(), 'core')
            self.coex = coexpressdb.Coexpression(species)
            self.decay = {line[-2]: line[-3] for line in neddata}
            self.species = species

    def process_data(self):
        self.outdata = []
        for struc in self.corum.strucs:
            self.complex = self.corum[struc]
            entrez = [e for sublist in self.complex.entrez for e in sublist]
            if self.homolog_status:
                # Convert human entrez to mouse entrez
                mouse_entrez = []
                reverse = {}
                for e in entrez:
                    if e in self.entrezhoms:
                        for eh in self.entrezhoms[e]:
                            if eh in self.coex.avail:
                                mouse_entrez.append(eh)
                                reverse[eh] = e
                                break
                # Get average coexpression, then convert back to human ids
                self.avg_coex = self._avg_coexpression(mouse_entrez)
                self.avg_coex = {reverse[key]: value
                                 for key, value in self.avg_coex.items()}
                for subunit in self.complex.entrez:
                    info = self._pick_best_subunits(subunit)
                    if info == None:
                        continue
                    info = [struc] + [str(i) for i in info] + [self.species]
                    self.outdata.append(info)
            else:
                self.avg_coex = self._avg_coexpression(entrez)
                for subunit in self.complex.entrez:
                    info = self._pick_best_subunits(subunit)
                    if info == None:
                        continue
                    info = [struc] + [str(i) for i in info] + [self.species]
                    self.outdata.append(info)

    def write_to_file(self, filename):
        with open(filename, 'w') as outfile:
            header = ['comp', 'entrez', 'uniprot', 'avgcoex', 'def', 'species']
            outfile.write('{0}\n'.format('\t'.join(header)))
            for line in self.outdata:
                outfile.write('{0}\n'.format('\t'.join(line)))

    def _entrez_to_uniprot(self):
        subunit_map = {}
        entrez = [e for sublist in self.complex.entrez for e in sublist]
        uniprot = [u for sublist in self.complex.uniprot for u in sublist]
        if len(entrez) != len(uniprot):
            print('number of entrez subunits != number of uniprot')
        for i in range(len(entrez)):
            subunit_map[entrez[i]] = uniprot[i]
        return subunit_map

    def _avg_coexpression(self, subunits):
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

    def _pick_best_subunits(self, subunits):
        """In cases where 2 or more subunits cannot be distinguished from one
        another, CORUM encloses those subunits within brackets. This function
        returns the subunit for which most data is available, prioritising data
        on coexpression over data on decay classification.
        """
        possibilities = []
        subunit_map = self._entrez_to_uniprot()
        for entrezid in subunits:
            if entrezid not in self.avg_coex:
                continue
            if self.homolog_status == True:
                for mouse_upr in self.unihoms.get(entrezid, ('NA')):
                    if mouse_upr in self.decay:
                        decay = self.decay[mouse_upr]
                        break
                    else:
                        decay = 'NA'
                mouse_entrezid = self.entrezhoms[entrezid]
                coexpression = self.avg_coex[entrezid]
                decay = self.decay.get(mouse_upr, 'NA')
                possibilities.append((mouse_entrezid, mouse_upr,
                                      coexpression, decay))
            else:
                upr = subunit_map[entrezid]
                coexpression = self.avg_coex[entrezid]
                decay = self.decay.get(upr, 'NA')
                possibilities.append((entrezid, upr, coexpression, decay))
        possibilities.sort(key=lambda x: (x.count('NA'), [x[2]].count('NA')))
        if len(possibilities) == 0:  # possibly unused?
            return
        return possibilities[0]


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
        if median(nvals) > median(evals):
            success += 1
    print(success, trials, binom_test(success, trials))


def main():
    # calc = CoexpressCalculator('mouse')
    # calc.process_data()
    # calc.write_to_file('data/coexpressdb_corum_mouse.tsv')
    analyse_corum_data('test.tmp')
    # calc = CoexpressCalculator('human')
    # calc.process_data()
    # calc.write_to_file('data/coexpressdb_corum_mouse.tsv')
    # analyse_corum_data('data/coexpressdb_corum_mouse.tsv')
    # analyse_corum_data('data/coexpressdb_corum_human.tsv')

if __name__ == '__main__':
    main()