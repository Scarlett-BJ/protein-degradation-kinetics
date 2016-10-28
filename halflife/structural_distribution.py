#!/usr/bin/env python3

from collections import namedtuple
import logging
import sys

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
log = logging.getLogger('structural')

# Figure 5, Panel A

def load_genes(species):
    """Return dictionary containing genes mapped to pdb and decay."""
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
                    ribosome_blacklist.append((line[0], line[2]))
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
    print(qs)
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
        for gene in qs:
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

################################################################################
## Above lies legacy code, below the good stuff

def load_strucs(species):
    strucs = {}
    Info = namedtuple('Info', ['decayclass', 'unique'])
    with open('data/revised_data/{0}_map.txt'.format(species)) as infile:
        infile.readline()
        for line in infile:
            line = line.split()
            # if line[6] == '0':
            strucs[(line[0], line[2])] = Info(line[1], line[4])
    return strucs

# Figure 5, Panel A
def qstype_data(species, remove_ribosomes=False):
    _, ribosome_blacklist = load_genes(species)
    strucs = load_strucs(species)
    with open ('data/revised_data/{0}_QS.txt'.format(species)) as infile:
        infile.readline()
        data = [line.strip().split('\t') for line in infile]
    header = 'gene\tpdb\tdecay.class\tqtype\tunq\n'
    if remove_ribosomes:
        outfilename = 'data/figdata/panela_{0}_noribo.txt'.format(species)
    else:
        outfilename = 'data/figdata/panela_{0}.txt'.format(species)
    with open(outfilename, 'w') as outfile:
        outfile.write(header)
        for l in data:
            if remove_ribosomes and (l[0], l[1]) in ribosome_blacklist:
                log.warning('ribosomal {0}, {1} removed'.format(l[0], l[1]))
                continue
            l = [l[0], l[1], strucs[(l[0], l[1])].decayclass,
                 l[2], strucs[(l[0], l[1])].unique]
            outfile.write('\t'.join(l)+'\n')

# Figure 5, Panel B
def interface_data(species):
    strucs = load_strucs(species)
    with open('data/revised_data/{0}_interfaces.txt'.format(species)) as infile:
        infile.readline()
        data = [line.strip().split('\t') for line in infile]
    header = 'Gene\tClass\tPDB chain\tInterface size\tUnique subunits\n'
    outfilename = 'data/figdata/panelb_{0}.txt'.format(species)
    with open(outfilename, 'w') as outfile:
        outfile.write(header)
        for l in data:
            l = [l[0], strucs[(l[0], l[1])].decayclass, l[1],
                 l[2], strucs[(l[0], l[1])].unique]
            outfile.write('\t'.join(l)+'\n')

# Figure 5, Panel C
def assembly_data(species):
    strucs = load_strucs(species)
    with open('data/revised_data/{0}_assembly.txt'.format(species)) as infile:
        infile.readline()
        data = [line.strip().split('\t') for line in infile]
    header = 'Gene\tClass\tPDB chain\tNormalised assembly order\tusubs\n'
    outfilename = 'data/figdata/panelc_{0}.txt'.format(species)
    with open(outfilename, 'w') as outfile:
        outfile.write(header)
        for l in data:
            l = [l[0], strucs[(l[0], l[1])].decayclass, l[1], l[2],
                 strucs[(l[0], l[1])].unique]
            outfile.write('\t'.join(l)+'\n')

# Figure 5, Panel D
def coexpression_data(species):
    strucs = load_strucs(species)
    filename = 'data/revised_data/{0}_coexpression.txt'.format(species)
    with open(filename) as infile:
        infile.readline()
        data = [line.strip().split('\t') for line in infile]
    header = 'Gene\tClass\tPDB chain\tCoexpression\tUnique subunits\n'
    outfilename = 'data/figdata/paneld_{0}.txt'.format(species)
    with open(outfilename, 'w') as outfile:
        outfile.write(header)
        for l in data:
            if (l[0], l[1]) not in strucs:
                log.warning('{0} missing'.format(l[1]))
                continue
            l = [l[0], strucs[(l[0], l[1])].decayclass, l[1],
                 l[2], strucs[(l[0], l[1])].unique]
            outfile.write('\t'.join(l)+'\n')

def main():
    for species in ('mouse', 'human'):
        qstype_data(species)
        qstype_data(species, remove_ribosomes=True)
        interface_data(species)
        assembly_data(species)
        coexpression_data(species)

if __name__ == '__main__':
    main()
