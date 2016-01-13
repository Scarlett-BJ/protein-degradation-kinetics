#!/usr/bin/env python3

from collections import namedtuple
from filter_strucs import load_pfam

def load_genes(filename):
    genes = {}
    with open(filename) as infile:
        for line in infile:
            line = line.split()
            if float(line[2]) < 70.0:
                break
            if line[1] in genes:
                continue
            # if line[3] == 'X':
            #     continue
            else:
                pdb, chn = line[0].split('_')
                genes[line[1]] = (pdb, chn, line[3])
    return genes

def load_structure_data():
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

def write_gene_info(genes, strucs, species):
    fname = 'data/structural/NED_quaternary_{0}.txt'.format(species)
    with open(fname, 'w') as outfile:
        header = ['gene', 'struc', 'chn', 'decay.class', 'qtype',
                  'unq', 'tot', '\n']
        outfile.write('\t'.join(header))
        for gene in genes:
            if genes[gene][0] not in strucs:
                continue
            struc = strucs[genes[gene][0]]
            line = [gene, '\t'.join(genes[gene]), struc.qtype,
                    struc.unique, struc.total]
            outfile.write('\t'.join(line)+'\n')

def main(species):
    genes = load_genes('data/structural/all_{0}.out'.format(species))
    strucs = load_structure_data()
    fstrucs = load_pfam('ibosomal', 'itochondrial ribo', 'roteasome')
    fstrucs = {s.lower() for s in fstrucs}
    fstrucs = fstrucs.intersection(set(strucs))
    # If you don't want to filter out bigger structures remove this loop
    for f in fstrucs:
        strucs.pop(f)
    write_gene_info(genes, strucs, species)

if __name__ == '__main__':
    main('mouse')
