#!/usr/bin/env python3

from collections import namedtuple

def load_genes(species):
    """Return dictionary contating genes mapped to pdb and decay."""
    genes = {}
    ribosome_blacklist = []
    filename = 'data/structural/all_{0}.out'.format(species)
    with open(filename) as infile:
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

def write_gene_info(genes, strucs, species):
    """Write out combined decay and structural info."""
    fname = 'data/structural/NED_quaternary_{0}.txt'.format(species)
    with open(fname, 'w') as outfile:
        header = ['gene', 'struc', 'chn', 'decay.class', 'qtype',
                  'unq', 'tot', 'species', '\n']
        outfile.write('\t'.join(header))
        for gene in genes:
            if genes[gene][0] not in strucs:
                continue
            struc = strucs[genes[gene][0]]
            line = [gene, '\t'.join(genes[gene]), struc.qtype,
                    struc.unique, struc.total, species]
            outfile.write('\t'.join(line)+'\n')

def main(species):
    genes, rblist = load_genes(species)
    strucs = load_structure_data()
    # If you don't want to filter out ribosomes remove this loop
    for r in rblist:
        strucs.pop(r)
    for s in strucs:
        if int(strucs[s].unique) > 64:
            print(s, strucs[s])
    write_gene_info(genes, strucs, species)

if __name__ == '__main__':
    main('mouse')
    main('human')
