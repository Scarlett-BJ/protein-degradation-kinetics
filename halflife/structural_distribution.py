#!/usr/bin/env python3

from collections import namedtuple

def load_structures(filename):
    genes = {}
    with open(filename) as infile:
        for line in infile:
            line = line.split()
            if float(line[2]) < 70.0:
                break
            if line[1] in genes:
                continue
            if line[3] == 'X':
                continue
            else:
                pdb, chn = line[0].split('_')
                genes[line[1]] = (pdb, chn, line[3])
    return genes

def load_structure_data():
    with open('data/ned_structural_files/table1.out') as infile:
        strucs = {line.split()[0]: line.split()[1:] for line in infile}
    Stoich = namedtuple('Stoich', ['total', 'unique', 'symmetry'])
    for struc in strucs:
        i = strucs[struc]
        stoich = Stoich(int(i[0]), int(i[1]), i[2])
        if stoich.total == 1 and stoich.unique == 1:
            qtype = 'mon'
        elif stoich.unique != 1:
            qtype = 'het'
        elif stoich.total != 1 and stoich.unique == 1:
            qtype = 'hom'
        strucs[struc] = qtype
    return strucs

def write_gene_info(genes, strucs, species):
    fname = 'data/ned_structural_files/NED_quaternary_{0}.txt'.format(species)
    with open(fname, 'w') as outfile:
        header = ['gene', 'struc', 'chn', 'decay.class', 'qtype', '\n']
        outfile.write('\t'.join(header))
        for gene in genes:
            line = [gene, '\t'.join(genes[gene]), strucs[genes[gene][0]]]
            outfile.write('\t'.join(line)+'\n')


if __name__ == '__main__':
    genes = load_structures('data/ned_structural_files/all_human.out')
    strucs = load_structure_data()
    write_gene_info(genes, strucs, 'human')
