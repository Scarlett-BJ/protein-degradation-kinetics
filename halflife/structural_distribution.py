#!/usr/bin/env python3

from collections import namedtuple
import logging
import sys

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
log = logging.getLogger('structural')

def load_genes(species):
    """Return dictionary contating genes mapped to pdb and decay."""
    genes = {}
    ribosome_blacklist = []
    filename = 'data/structural/all_{0}.out'.format(species)
    with open(filename) as infile:
        infile.readline()
        for line in infile:
            line = line.split()
            for i in ['Rps', 'Rpl', 'Mrpl', 'Mrps', # Mouse
                      'RPS', 'RPL', 'MRPL', 'MRPS']: # Human
                if i in line[1]:
                    ribosome_blacklist.append(line[0].split('_')[0])
            if float(line[2]) < 70.0:
                continue
            if line[1] in genes:
                # Ensure only top hit is selected. Assumes infile is sorted...
                continue
            else:
                pdb, chn = line[0].split('_')
                genes[line[1]] = (pdb, chn, line[3])
    return genes, set(ribosome_blacklist)

def load_structure_data():
    """Returns dictionary of structural information."""
    with open('data/structural/table1.out') as infile:
        strucs = {line.split()[0]: line.split()[1:] for line in infile}
    Stoich = namedtuple('Stoich', ['total', 'unique', 'symmetry', 'qtype'])
    for struc in strucs:
        i = strucs[struc]
        stoich = Stoich(int(i[0]), int(i[1]), i[2], 'NA')
        if stoich.total == 1 and stoich.unique == 1:
            qtype = 'mon'
        elif stoich.unique != 1:
            qtype = 'het'
        elif stoich.total != 1 and stoich.unique == 1:
            qtype = 'hom'
        strucs[struc] = Stoich(i[0], i[1], i[2], qtype)
    return strucs

def write_gene_info(species, filter_ribosomes=True):
    """Write out combined decay and structural info."""
    genes, rblist = load_genes(species)
    strucs = load_structure_data()
    if filter_ribosomes:
        fname = 'data/structural/NED_quaternary_{0}.txt'.format(species)
    else:
        fname = 'data/structural/NED_quaternary_{0}_ribo.txt'.format(species)
    with open(fname, 'w') as outfile:
        header = ['gene', 'struc', 'chn', 'decay.class', 'qtype',
                  'unq', 'tot', 'species', '\n']
        outfile.write('\t'.join(header))
        for gene in genes:
            pdb = genes[gene][0]
            if pdb not in strucs:
                log.warning('{0} {1} not in table1.out'.format(species, pdb))
                continue
            elif filter_ribosomes and pdb in rblist:
                continue
            struc = strucs[pdb]
            line = [gene, '\t'.join(genes[gene]), struc.qtype,
                    struc.unique, struc.total, species]
            outfile.write('\t'.join(line)+'\n')

def main():
    write_gene_info('mouse')
    write_gene_info('human')
    write_gene_info('mouse', False)
    write_gene_info('human', False)

if __name__ == '__main__':
    main()
