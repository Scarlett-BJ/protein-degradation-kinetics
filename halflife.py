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
        self.coex = None

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

    def get_class(self):
        if self.score >= 0.9:
            clss = 'N'
        elif self.score <= 0.1:
            clss = 'E'
        else:
            clss = 'U'
        return clss

    def set_avg_coexpression(self, coex):
        self.coex = coex

    def __repr__(self):
       return ' '.join([self.name, self.chain, str(self.score)])

class MouseTable(object):
    """Table of MouseProteins."""
    def __init__(self):
        coex = coexpressdb.Coexpression()
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
        self.members = set(self._prots)

    def __len__(self):
        return len(self.members)

    def __getitem__(self, protein):
        return self._prots[protein]

    def __repr__(self):
        return '\n'.join(sorted(self.members))


class Complex(object):
    def __init__(self, name, members):
        self.name = name
        self.members = members

    def filter_members(self, clss):
        subset = []
        for mem in self.members:
            if mem.get_class() == clss:
                subset.append(mem)
        return subset

    def __len__(self):
        return len(self.members)

    def __repr__(self):
        return self.name


class Complexes(MouseTable):
    def __init__(self):
        super().__init__()
        self._init_raw_complexes()
        self._filter_identical()
        self._init_fleshed_complexes()
        self.members = set(self._complexes)

    def _init_raw_complexes(self):
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


def coexpression_test():
    test = Complexes()
    test.filter_by_size(2, 50)
    coex = coexpressdb.Coexpression()
    mean = lambda x: sum(x)/len(x)
    successes = 0
    trials = 0
    for pdb in test.members:
        comp = test[pdb]
        avail = [p for p in comp.members if p.name in coex.entrez]
        navail = [p for p in avail if p.get_class() != 'E']
        eavail = [p for p in avail if p.get_class() != 'N']
        if len(avail) == 2 or len(navail) < 2 or len(eavail) < 2:
            continue
        # pairs = list(combinations(avail, 2))
        compcoex = mean([coex.get_coexpression(p[0].name, p[1].name)
                         for p in combinations(avail, 2)])
        ncoex = mean([coex.get_coexpression(p[0].name, p[1].name)
                      for p in combinations(navail, 2)])
        ecoex = mean([coex.get_coexpression(p[0].name, p[1].name)
                      for p in combinations(eavail, 2)])
        if ncoex > ecoex:
            successes += 1
        trials += 1
    print(successes, trials)
    #         avg_coex = []
    #         for p2 in avail:
    #             if p1.name == p2.name:
    #                 continue
    #             cor = round(coex.get_coexpression(p1.name, p2.name), 3)
    #             avg_coex.append(cor)
    #         avg_coex = sum(avg_coex)/len(avg_coex)
    #         p1.set_avg_coexpression(avg_coex)
    # wincount = 0
    # trials = 0
    # for pdb in test.members:
    #     comp = test[pdb]
    #     ncx = [p.coex for p in comp.filter_members('N') if p.coex != None]
    #     ecx = [p.coex for p in comp.filter_members('E') if p.coex != None]
    #     ncx += [p.coex for p in comp.filter_members('U') if p.coex != None]
    #     mean = lambda x: sum(x)/len(x)
    #     if len(ncx) == 0 or len(ecx) == 0:
    #         continue
    #     if mean(ecx) > mean(ncx):
    #         wincount += 1
    #     trials += 1
    # print(wincount, trials)

def interface_test():
    ints = ix.Interface()
    test = Complexes()
    test.filter_by_size(2, 20)
    mean = lambda x: sum(x)/len(x)
    wincount = 0
    trials = 0
    for pdb in test.members:
        comp = test[pdb]
        if len(comp) < 2:
            continue
        pairs = list(combinations(comp.members, 2))
        ue = []
        no_e = []
        # ee = []
        for p in pairs:
            p1c = p[0].get_class()
            p2c = p[1].get_class()
            chainA = p[0].chain
            chainB = p[1].chain
            if not ints.is_present(pdb, chainA, chainB):
                continue
            if p1c == 'E' or p2c == 'E':
                ue.append(ints.get_interface(pdb, chainA, chainB))
            elif p1c != 'E' and p2c != 'E':
                no_e.append(ints.get_interface(pdb, chainA, chainB))
            # elif p1c == 'E' and p2c == 'E':
            #     ee.append(ints.get_interface(pdb, chainA, chainB))
        if len(ue) == 0 or len(no_e) == 0:
            continue
        if mean(ue) < mean(no_e):
            wincount += 1
        trials += 1
    print(wincount, trials)

def ned_size_distribution():
    test = Complexes()
    test.filter_by_size(1, 10)
    # print('\t'.join(['pdb', 'len', 'len.N', 'len.E', 'len.U']))
    for pdb in test.members:
        comp = test[pdb]
        memscore = [mem.score for mem in comp.members]
        memabund = [mem.relative_abund for mem in comp.members]
        print(comp, len(comp), sum(memscore)/len(memscore),
              sum(memabund)/len(memabund))

def abundance_per_complex():
    test = Complexes()
    test.filter_by_size(2, 2)
    wincount = 0
    trials = 0
    for pdb in test.members:
        comp = test[pdb]
        n_mems = [m.relative_abund for m in comp.filter_members('N')]
        e_mems = [m.relative_abund for m in comp.filter_members('E')]
        u_mems = [m.relative_abund for m in comp.filter_members('U')]
        mean = lambda x: sum(x)/len(x)
        if len(n_mems) == 0 or len(e_mems) == 0:
            continue
        if mean(n_mems) > mean(e_mems):
            wincount += 1
        trials += 1
    print(wincount, trials)

def tcount(prot):
    meta = paxdb.get_metadata('10090')
    tcount = 0
    for tissue in meta:
        if tissue == 'WHOLE_ORGANISM' or tissue == 'CELL_LINE':
            continue
        tabund = paxdb.Abundances(meta[tissue]['filename'], 0)
        if tabund.isexpressed(prot):
            tcount += 1
    return tcount

def tissue_expression_percomp():
    test = Complexes()
    test.filter_by_size(2, 10)
    meta = paxdb.get_metadata('10090')
    wincount = 0
    trials = 0
    for pdb in test.members:
        comp = test[pdb]
        n_mems = [tcount(m.name) for m in comp.filter_members('N')]
        e_mems = [tcount(m.name) for m in comp.filter_members('E')]
        u_mems = [tcount(m.name) for m in comp.filter_members('U')]
        mems = e_mems + u_mems
        mean = lambda x: sum(x)/len(x)
        if len(n_mems) == 0 or len(mems) == 0:
            continue
        if mean(n_mems) > mean(mems):
            wincount += 1
        trials += 1
    print(wincount, trials)

def extremely_ropey_test():
    """For each complex, are NED proteins expressed across a wider range of
    tissues that ED proteins? Seems to be unhealthily dependent on how many
    possible tissues are used. More tissues means more lo-quality data..."""
    test = Complexes()
    test.filter_by_size(2, 10)
    meta = paxdb.get_metadata('10090')
    abunds = [paxdb.Abundances(meta[t]['filename'], 0) for t in meta
               if t != 'WHOLE_ORGANISM' and t != 'CELL_LINE']
    wincount = 0
    trials = 0
    def tcount(prot):
        t = 0
        for tissue in abunds:
            if tissue.isexpressed(prot):
                t += 1
        return t

    for pdb in test.members:
        comp = test[pdb]
        neds = [m.name for m in comp.filter_members('N')]
        eds = [m.name for m in comp.filter_members('E')]
        us = [m.name for m in comp.filter_members('U')]
        # eds += us
        if len(eds) == 0 or len(neds) == 0:
            continue
        nedcounts = [tcount(prot) for prot in neds if tcount != 0]
        edcounts = [tcount(prot) for prot in eds if tcount != 0]
        # ucounts = [tcount(prot) for prot in us]
        mean = lambda x: sum(x)/len(x)
        if mean(nedcounts) >= mean(edcounts):
            wincount += 1
        trials += 1
    print(wincount, trials)

def pairwise_coexpression():
    test = Complexes()
    test.filter_by_size(2, 20)
    coex = coexpressdb.Coexpression()
    mean = lambda x: sum(x)/len(x)
    wincount = 0
    trials = 0
    for pdb in test.members:
        comp = test[pdb]
        avail = [p for p in comp.members if p.name in coex.entrez]
        if len(avail) < 2:
            continue
        pairs = list(combinations(avail, 2))
        elist =  []
        no_elist = []
        for p in pairs:
            p1c = p[0].get_class()
            p2c = p[1].get_class()
            cor = coex.get_coexpression(p[0].name, p[1].name)
            if p1c != 'E' and p2c != 'E':
                no_elist.append(cor)
            if p1c == 'E' or p2c == 'E':
                elist.append(cor)
        if len(no_elist) == 0 or len(elist) == 0:
            continue
        if mean(no_elist) > mean(elist):
            wincount += 1
        trials += 1
    print(wincount, trials)



if __name__ == '__main__':
    interface_test()

