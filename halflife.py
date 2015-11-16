#!/usr/bin/env python3

import ixntools as ix
from expression import paxdb, coexpressdb
from itertools import combinations
from collections import namedtuple
import copy

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
    """Collection of info on individual mouse protein."""
    def __init__(self, name, score, relabund):
        self.name = name
        self.score = score
        self.relative_abund = relabund
        self.chain = None
        self.pdb = None
        self.seqid = None
        self.ortholog = None

    def set_pdb(self, pdb, chain, seqid):
        self.pdb = pdb
        self.chain = chain
        self.seqid = seqid

    @property
    def decay(self):
        if self.score >= 0.9:
            decay = 'N'
        elif self.score <= 0.1:
            decay = 'E'
        else:
            decay = 'U'
        return decay


class MouseTable(object):
    """Table of MouseProteins."""
    def __init__(self):
        with open('data/mouse_genes_ned.txt') as infile:
            data = [line.strip().split(',') for line in infile]
        self._prots = {i[3]: MouseProtein(i[3], float(i[4]), float(i[5]))
                       for i in data}
        self.members = set(self._prots)

    def set_orthologs(self):
        with open('data/mouse_human_orthologs.txt') as infile:
            for line in infile:
                if len(line) != 2 or line[0] not in self.members:
                    continue
            self._prots[line[0]].ortholog = line[1]

    def __len__(self):
        return len(self.members)

    def __getitem__(self, protein):
        return self._prots[protein]

    def __repr__(self):
        return '\n'.join(sorted(self.members))


class Complex(namedtuple('_Complex', ['name', 'members'])):
    """Simple representation of a protein complex"""
    def filter_members(self, clss):
        return [mem for mem in self.members if mem.get_class == clss]

    def __len__(self):
        return len(self.members)

    def __repr__(self):
        members = [mem.name for mem in self.members]
        return '{0}\n{1}'.format(self.name, '\n'.join(members))


class Complexes(MouseTable):
    """Table of complexes"""
    def __init__(self):
        super().__init__()
        self._init_raw_complexes()
        self._init_fleshed_complexes()
        self.members = set(self._complexes)

    def _init_raw_complexes(self):
        """Messy, see following 2 methods. _raw_complexes is especially bad
        because it more than doubles the amount of data being stored.
        """
        self._raw_complexes = {}
        Chain = namedtuple('Chain', ['prot', 'chain', 'seqid'])
        with open('data/ned_mapped_to_pdb.out') as infile:
            data = [line.strip().split() for line in infile]
        for line in data:
            pdb, chain = line[0].split('_')
            seqid = float(line[3])
            prot = Chain(line[1], chain, seqid)
            if pdb not in self._raw_complexes:
                self._raw_complexes[pdb] = [prot]
            else:
                self._raw_complexes[pdb].append(prot)

    def _filter_identical(self):
        """Deprecated, only useful in event of no prior reduncancy filtering
        of protein complexes.
        """
        cset = {tuple(sorted([p.prot for p in pdb]))
                for pdb in self._raw_complexes.values()}
        final = {}
        for comp in self._raw_complexes:
            val = tuple(sorted([p.prot for p in self._raw_complexes[comp]]))
            if val in cset:
                cset.remove(val)
                final[comp] = self._raw_complexes[comp]
        self._raw_complexes = final

    def _init_fleshed_complexes(self):
        """Will likely go the same way as _filter_identical."""
        self._complexes = {}
        for comp in self._raw_complexes:
            prots = {p.prot for p in self._raw_complexes[comp]}
            if len(prots.intersection(self.members)) != len(prots):
                continue
            fleshed = []
            for member in self._raw_complexes[comp]:
                prot = copy.copy(self._prots[member.prot])
                prot.set_pdb(comp, member.chain, member.seqid)
                fleshed.append(prot)
            self._complexes[comp] = Complex(comp, fleshed)

    def filter_by_list(self, filename):
        with open(filename) as infile:
            complist = {line.strip() for line in infile}
        for comp in list(self._complexes):
            if comp not in complist:
                self._complexes.pop(comp)
        self.members = set(self._complexes)

    def filter_by_size(self, lower, upper):
        for comp in list(self._complexes):
            size = len(self._complexes[comp])
            if size < lower or size > upper:
                self._complexes.pop(comp)
        self.members = set(self._complexes)

    def __getitem__(self, pdb):
        return self._complexes[pdb]

    def __len__(self):
        return len(self._complexes)


def hein():
    intome = ix.Interactome()
    comps = MouseProteins
    comps.filter_by_size(3, 3)


if __name__ == '__main__':
    gene_map = get_gene_map()
    complexes = Complexes()
    complexes.filter_by_list('complex_set.txt')
    for prot in complexes['5a1v'].members:
        print(prot.name, gene_map[prot.name].upr, prot.decay)
    # for comp in complexes.members:
    #     comp = complexes[comp]
    #     ncount = 0
    #     ucount = 0
    #     ecount = 0
    #     for prot in comp.members:
    #         dclass = prot.decay
    #         if dclass == 'N':
    #             ncount += 1
    #         elif dclass == 'U':
    #             ucount += 1
    #         elif dclass == 'E':
    #             ecount += 1
    #     ln = len(comp.members)
    #     print(comp.name, ln, ncount/ln, ucount/ln, ecount/ln)
