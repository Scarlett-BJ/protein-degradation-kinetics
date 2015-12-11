#!/usr/bin/env python3

from expression import paxdb, proteomicsdb, proteomicsdb_requests as req
import numpy as np

def load_data(filename):
    with open(filename) as infile:
        data = [line.split() for line in infile]
    header = data[0]
    data = data[1:]
    return header, data

def protein_list(species):
    if species == 'mouse':
        prots = [line[0] for line in load_data('data/NED_mouse_update.txt')[1]]
    elif species == 'human':
        prots = [line[0] for line in load_data('data/NED_human.txt')[1]]
    return prots

def protein_map(species):
    if species == 'mouse':
        filename = 'data/mouse_uniprot_ensp.txt'
    elif species == 'human':
        filename = 'data/human_uniprot_ensp.txt'
    with open(filename) as infile:
        data = [line.split() for line in infile]
    prot_map = {}
    for line in data:
        if line[0] not in prot_map:
            prot_map[line[0]] = [line[1]]
        else:
            prot_map[line[0]].append(line[1])
    return prot_map

def load_expression_data(species):
    if species == 'mouse':
        data = load_data('data/NED_mouse_update.txt')[1]
        meta = paxdb.get_metadata('10090')
        pmap = protein_map('mouse')
    elif species == 'human':
        data = load_data('data/NED_human.txt')[1]
        meta = paxdb.get_metadata('9606')
        pmap = protein_map('human')
    skip = ['CELL_LINE']
    for i in skip:
        meta.pop(i)
    abunds = [paxdb.Abundances(meta[t]['filename'], 0) for t in sorted(meta)]
    header = '\t'.join(['prot', '\t'.join(sorted(meta)), 'tcount', 'def'])
    return header, data, abunds, pmap

def print_expression_data(species):
    header, data, abunds, pmap = load_expression_data(species)
    ofname = 'data/{0}_paxdb_tissue_expression.txt'.format(species)
    with open(ofname, 'w') as outfile:
        outfile.write(header + '\n')
        for line in data:
            if line[0] not in pmap:
                continue
            # list of lists containing tissue expressions for each ensp isoform
            ensps = [[t.data.get(p, 'NA') for t in abunds]
                     for p in pmap[line[0]]]
            ensps = [i for i in ensps if i.count('NA') != len(i)]
            if len(ensps) == 0:
                continue
            # Sort by max tcount, then max avg abundance
            skey = lambda x: (x.count('NA'),
                              -np.mean([p for p in x if p != 'NA']))
            ensps.sort(key=skey)
            expression = '\t'.join([str(e) for e in ensps[0]])
            tcount = str(len(ensps[0]) - ensps[0].count('NA'))
            newline = '\t'.join([line[0], expression, tcount, line[-1]])
            outfile.write(newline + '\n')

def load_proteomicsdb_data():
    data = load_data('data/NED_human.txt')[1]
    isoab = proteomicsdb.Abundances('human_protdb_isoforms_expression.txt')
    abunds = proteomicsdb.Abundances('trembl_tissues.txt')
    tissues = abunds.tissues
    ttargets = ['B-lymphocyte', 'adipocyte', 'adrenal gland', 'brain',
                'colon', 'colonic epithelial cell', 'colorectal cancer cell',
                'cytotoxic T-lymphocyte', 'esophagus', 'gall bladder', 'gut',
                'heart', 'helper T-lymphocyte', 'ileum epithelial cell',
                'liver', 'lung', 'lymph node', 'monocyte', 'myometrium',
                'natural killer cell', 'ovary', 'pancreas',
                'pancreatic islet', 'placenta', 'prostate gland', 'rectum',
                'retina', 'stomach', 'testis', 'thyroid gland',
                'salivary gland', 'spinal cord', 'spleen', 'tonsil',
                'urinary bladder', 'uterine cervix', 'uterus']
    header = '{0}\t{1}\t{2}\t{3}\n'.format('prot', '\t'.join(ttargets),
                                           'tcount', 'def')
    with open('data/human_protdb_limited_tissues.txt', 'w') as outfile:
        outfile.write(header)
        for line in data:
            prot = line[0]
            if prot in abunds.proteins:
                expression = [abunds.expression(prot, tis) for tis in ttargets]
            elif prot in isoab.proteins:
                expression = [isoab.expression(prot, tis) for tis in ttargets]
            else:
                continue
            tcount = str(len(expression) - expression.count('NA'))
            newline = '\t'.join([prot] + expression + [tcount, line[-1]])
            outfile.write(newline + '\n')




def get_missing():
    abunds = proteomicsdb.Abundances('trembl_tissues.txt')
    abunds_iso = proteomicsdb.Abundances('human_protdb_isoforms_expression.txt')
    data = load_data('data/NED_human.txt')[1]
    prots = {line[0] for line in data}
    missing = prots.intersection(abunds_iso.proteins.union(abunds.proteins))
    print(len(prots) - len(missing))



if __name__ == '__main__':
    # print_expression_data('mouse')
    # print_expression_data('human')
    load_proteomicsdb_data()
    # help(req)
    # get_missing()
