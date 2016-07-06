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
    filename = 'data/revised_data/{0}_map.txt'.format(species)
    with open(filename) as infile:
        header = infile.readline()
        for line in infile:
            line = line.split()
            for i in ['Rps', 'Rpl', 'Mrpl', 'Mrps', # Mouse
                      'RPS', 'RPL', 'MRPL', 'MRPS']: # Human
                if i in line[0]:
                    ribosome_blacklist.append(line[2].split('_')[0])
            if float(line[3]) < 70.0 or int(line[6]) != 0:
                # Ensure only sufficiently similar and best option is selected.
                continue
            else:
                pdb, chn = line[2].split('_')
                genes[line[0]] = (pdb, chn, line[1])
    return genes, set(ribosome_blacklist)

def load_structure_data(species):
    """Returns dictionary of structural information."""
    with open('data/revised_data/{0}_QS.txt'.format(species)) as infile:
        header = infile.readline()
        qs = {line.split()[1]: line.strip().split()[2] for line in infile}
    with open('data/revised_data/{0}_map.txt'.format(species)) as infile:
        strucs = {line.split()[2]: line.split()[4] for line in infile}
    Stoich = namedtuple('Stoich', ['unique', 'qtype'])
    for struc in list(strucs):
        if struc in qs:
            stoich = Stoich(strucs[struc], qs[struc])
            strucs[struc.split('_')[0]] = stoich
            strucs.pop(struc)
        else:
            strucs.pop(struc)
    return strucs

def write_gene_info(species, filter_ribosomes=False):
    """Write out combined decay and structural info."""
    genes, rblist = load_genes(species)
    strucs = load_structure_data(species)
    if filter_ribosomes == True:
        fname = 'data/revised_data/NED_quaternary_{0}_worib.txt'.format(species)
    else:
        fname = 'data/revised_data/NED_quaternary_{0}.txt'.format(species)
    with open(fname, 'w') as outfile:
        header = ['gene', 'struc', 'chn', 'decay.class', 'qtype',
                  'unq', 'species', '\n']
        outfile.write('\t'.join(header))
        for gene in genes:
            pdb = genes[gene][0]
            if pdb not in strucs:
                log.warning('{0} {1} not available'.format(species, pdb))
                continue
            elif filter_ribosomes and pdb in rblist:
                log.warning('Ribosomal chain in {0} removed'.format(pdb))
                continue
            struc = strucs[pdb]
            line = [gene, '\t'.join(genes[gene]), struc.qtype,
                    struc.unique, species]
            outfile.write('\t'.join(line)+'\n')


# def filter_top_seqid_isize():
#     with open('data/old_data/isize.txt') as infile:
#         infile.readline()
#         print(len([i.split()[2] for in infile]))
#     with open('data/old_data/isize.txt') as infile:
#         new_data = [infile.readline()]
#         sdict = []
#         for line in  infile:
#             sline = line.strip().split('\t')
#             if sline[2] not in sdict:
#                 sdict.append(sline[2])
#                 new_data.append(line)
    # for line in new_data:
    #     print(line.strip())
    # print(len(new_data))

def filter_top_seqid_assembly():
    with open('data/structural/assembly.txt') as infile:
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
    write_gene_info('mouse')
    write_gene_info('human')
    write_gene_info('mouse', True)
    write_gene_info('human', True)
    # filter_top_seqid_isize()
    # a, b = load_genes('mouse')
    # print(b)

if __name__ == '__main__':
    main()
