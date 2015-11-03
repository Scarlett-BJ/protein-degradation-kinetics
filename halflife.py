#!/usr/bin/env python3

import ixntools as ix
from expression import paxdb
from itertools import combinations
from collections import namedtuple

"""To Do:
Think of a good way to format this so it's easily mutable and accessible.
Need to be able to filter quickly by complex attributes, e.g. size
Be able to add new information, therefore need class.

For all complexes:
    are NED proteins more abundant?

For all pairs:
    are NED pairs more likely to be coexpressed?

Develop further using new Hein data with stoichiometries...
"""

def get_gene_map():
    GeneID = namedtuple('GeneID', ['upr', 'hgnc', 'ensg'])
    with open('data/mouse_genes_ned.txt') as infile:
        data = [line.strip().split(',') for line in infile.readlines()]
    gene_map = {i[3]: GeneID(i[0], i[1], i[2]) for i in data}
    return gene_map


class MouseProtein(object):
    """Collection of info on mouse protein"""
    def __init__(self, name, score, relabund):
        self.name = name
        self.score = score
        self.relative_abund = relabund

    def set_ortholog(self, ortholog):
        self.otholog = ortholog

    def set_pdb(self, pdb, chain, seqid):
        self.pdb = pdb
        self.chain = chain
        self.seqid = seqid

    def get_chain(self, pdb):
        if pdb != self.pdb:
            raise ValueError('{0} not a member of {1}'.format(self.name, pdb))
        return self.chain

    def __repr__(self):
        return '\t'.join([str(self.score), str(self.relative_abund)])


class MouseTable(object):
    def __init__(self):
        fl1 = 'data/mouse_genes_ned.txt'
        fl2 = 'data/mouse_human_orthologs.txt'
        with open(fl1) as infile1, open(fl2) as infile2:
            data1 = [line.strip().split(',') for line in infile1.readlines()]
            data2 = [line.strip().split('\t') for line in infile2.readlines()]
        self._prots = {i[3]: MouseProtein(i[3], float(i[4]), float(i[5]))
                       for i in data1}
        for line in data2:
            if len(line) != 2:
                continue
            self._prots[line[0]].set_ortholog(line[1])
        self.members = sorted(self._prots)

    def __len__(self):
        return len(self.members)

    def __getitem__(self, protein):
        return self._prots[protein]

    def __repr__(self):
        return '\n'.join(self.members)

def get_pdb_info():
    pdb_info = {}
    Info = namedtuple('Info', ['chn', 'seqid'])
    with open('data/ned_mapped_to_pdb.out') as infile:
        data = [line.strip().split() for line in infile.readlines()]
    for line in data:
        pdb, chain = line[0].split('_')
        if pdb not in pdb_info:
            pdb_info[pdb] = {chain: Info(line[1], float(line[3]))}
        else:
            pdb_info[pdb][chain] = Info(line[1], float(line[3]))
    return pdb_info


if __name__ == '__main__':
    mouse = MouseTable()
    gmap = get_gene_map()
    for i in mouse.members:
        print(gmap[i].ensg)
