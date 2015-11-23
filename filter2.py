#!/usr/bin/env python3

from ixntools import structure
from collections import namedtuple
from itertools import combinations
from random import randint
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.neighbors import DistanceMetric
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.metrics import jaccard_similarity_score

def load_table(n):
    """Select appropriate minimum chain length, and then return a list of
    complexes from the relevant table. i.e. table30.out for minimum chain
    length of 30 residues.
    """
    with open('data/ned_mapped_to_pdb.out') as infile:
        strucs = {l.split('_')[0] for l in infile}
    table = structure.Table(n)
    useable = table.filter(name=strucs, unq_chns=range(2, 500))
    return useable

def get_chains(pool, threshold):
    """For all structures in table, if at least x percent of that structure’s
    chains (x ~ 90) map to a mouse protein, then add it to a list of
    potentially useable structures.
    """
    useable = []
    chaindict = {}
    Chain = namedtuple('Chain', ['prot', 'chn', 'seqid', 'score', 'dclass'])
    with open('data/ned_mapped_to_pdb.out') as infile:
        for line in infile:
            struc, _, prot, seqid, score, dclass = line.split()
            struc, chn = struc.split('_')
            chain = Chain(prot, chn, float(seqid), float(score), dclass)
            if struc not in chaindict:
                chaindict[struc] = [chain]
            else:
                chaindict[struc].append(chain)
    # Allows for some missing chains, number indirectly defined by threshold.
    for struc in pool:
        chncount = len(chaindict[struc.name])
        trucount = struc.unq_chns
        if chncount >= (threshold * struc.unq_chns)/100:
            struc.chains = chaindict[struc.name]
            useable.append(struc)
    # Important sort, maintains order for later clustering
    useable.sort(key=lambda x: x.name)
    return useable

def get_interfaces(pool):
    useable = []
    cmap = structure.chain_map()
    ints = structure.Interfaces()
    for struc in pool:
        if struc.name not in ints.members:
            continue
        unq_chns = {c.chn for c in struc.chains}
        smap = cmap[struc.name]
        struc_ints = {}
        for key in ints[struc.name]:
            if key[0] not in smap or key[1] not in smap:
                continue
            unq = tuple(sorted((smap[key[0]], smap[key[1]])))
            if unq[0] in unq_chns and unq[1] in unq_chns:
                interface = ints.get_interface(struc.name, key[0], key[1])
                if unq not in struc_ints:
                    struc_ints[unq] = interface
                else:
                    struc_ints[unq] = max((interface, struc_ints[unq]))
        struc.interfaces = struc_ints
        useable.append(struc)
    useable.sort(key=lambda x: x.name)
    return useable


class PFAMFilters(object):
    """Could probably stop using useable, adds unneccessary memory usage."""
    def __init__(self, pool):
        pfam_file = 'data/10090.tsv'
        self._pool = pool
        self._pfam, self._ribo, self._ig = structure.pfam_map(pfam_file)

    def assign_pfams(self):
        """Adds PFAM domain annotations to proteins in cluster."""
        Pfam = namedtuple('Pfam', ['prot', 'pfam'])
        i = 0
        for struc in self._pool:
            prots = [p.prot for p in struc.chains]
            pfams = [Pfam(p, self._pfam.get(p, {'u_'+str(i)})) for p in prots]
            struc.pfams = pfams
            i += 1

    def filter_immunoglobins(self):
        for struc in self._pool:
            pfams = [p.pfam for p in struc.pfams]
            pfams = {p for sublist in pfams for p in sublist}
            if len(pfams.intersection(self._ig)) > 0:
                self._pool.remove(struc)

    def filter_ribosomes(self, threshold=40):
        for struc in self._pool:
            pfams = [p.pfam for p in struc.pfams]
            pfams = {p for sublist in pfams for p in sublist}
            if len(pfams.intersection(self._ribo))/len(pfams)*100 > threshold:
                self._pool.remove(struc)

    def filter_paralogous(self, threshold=40):
        for struc in self._pool:
            jaccards = []
            pfams = [p.pfam for p in struc.pfams]
            for i in combinations(pfams, 2):
                jindex = len(i[0].intersection(i[1]))/len(i[0].union(i[1]))
                jaccards.append(jindex)
            median = np.median(jaccards)*100
            if median > threshold:
                self._pool.remove(struc)

    @property
    def pool(self):
        return sorted(self._pool, key=lambda x: x.name)

    @pool.setter
    def pool(self, pool):
        self._pool = sorted(pool, key=lambda x: x.name)


class Clusters(object):
    def __init__(self, pool):
        self._pool = pool
        self._clusters = []

    def cluster(self, eps, feature):
        # Nailed it... Might put this in cluster factory
        pool_dict = {c.name: c for c  in self._pool}
        if feature == 'prot':
            features = [sorted([p.prot for p in i.chains])
                        for i in pool_dict.values()]
        elif feature == 'pfam':
            features = [sorted([','.join(sorted(p.pfam)) for p in i.pfams])
                        for i in pool_dict.values()]
        features = MultiLabelBinarizer().fit_transform(features)
        # dist = DistanceMetric.get_metric('jaccard')
        # distance_matrix = dist.pairwise(features)
        # Not 100% sure this is correct yet.
        dbs = DBSCAN(eps=eps, min_samples=1, metric='jaccard')
        clust_labels = dbs.fit_predict(features)
        i = 0
        self._clusters = {c: [] for c in set(clust_labels)}
        for pdb in pool_dict:
            self._clusters[clust_labels[i]].append(pool_dict[pdb])
            i += 1
        return self._clusters

    def _chain(self, struc):
        return len(struc.chains)

    def _interface(self, struc):
        return len(struc.interfaces)

    def _seqid(self, struc):
        return struc.avgseqid

    def _name(self, struc):
        return struc.name

    def sort(self):
        sortfunc = lambda x: (self._interface(x), self._chain(x),
                              self._seqid(x), self._name(x))
        for label in self._clusters:
            if len(self._clusters[label]) == 1:
                continue
            self._clusters[label].sort(key=sortfunc, reverse=True)

    def get_complexes(self):
        complexes = []
        for label in self._clusters:
            complexes.append(self._clusters[label][0])
        return complexes

    @property
    def pool(self):
        return sorted(self._pool, key=lambda x: x.name)

    @pool.setter
    def pool(self, pool):
        self._pool = sorted(pool, key=lambda x: x.name)

    def __repr__(self):
        clist = [(str(k), ', '.join([c.name for c in v]))
                 for k, v in self._clusters.items()]
        clist = ['\t'.join(c) for c in clist]
        return '\n'.join(clist)

    def __getitem__(self, label):
        return self._clusters[label]

    def __len__(self):
        return len(self._clusters)


def summary(pool):
    """Summarises distribution of different sized complexes in a pool."""
    all_lengths = [i.tot_chns for i in pool]
    lengths = {i.unq_chns for i in pool}
    for i in sorted(lengths):
        print(i, all_lengths.count(i))
    print('total:', len(pool))

def pipeline():
    pool = load_table(30)
    pool = get_chains(pool, 80)
    pfam = PFAMFilters(pool)
    pfam.assign_pfams()
    pfam.filter_paralogous(50)
    pool = pfam.pool
    pool = get_interfaces(pool)
    clusters = Clusters(pool)
    clusters.cluster(0.34, 'prot')
    clusters.sort()
    pool = clusters.get_complexes()
    clusters.pool = pool
    clusters.cluster(0.1, 'pfam')
    clusters.sort()
    pool = clusters.get_complexes()
    pfam = PFAMFilters(pool)
    pfam.filter_immunoglobins()
    pfam.filter_ribosomes()
    pool = pfam.pool
    pool.sort(key=lambda x: x.name)
    return pool

def pairwise():
    pool = pipeline()
    header = ['struc', 'unq', 'tot', 'sym', 'p1', 'p1.score', 'p1.cls',
              'p2', 'p2.score', 'p2.cls', 'int',]
    print('\t'.join(header))
    for struc in pool:
        cdict = {c.chn: (c.prot, c.score, c.dclass) for c in struc.chains}
        for p in struc.interfaces:
            iface = struc.interfaces[p]
            info = [str(struc),
                    '\t'.join([str(i) for i in cdict[p[0]]]),
                    '\t'.join([str(i) for i in cdict[p[1]]]),
                    str(iface)]
            print('\t'.join(info))

def protwise():
    with open('data/mouse_genes_ned.txt') as infile:
        pdict = {i.split(',')[0]: i.split(',')[-1].strip() for i in infile}
    pool = pipeline()
    print('\t'.join(['struc', 'unq', 'tot', 'sym', 'prot', 's', 'ab', 'dcls']))
    for struc in pool:
        for prot in struc.chains:
            abund = pdict[prot.prot]
            print(struc, prot.prot, prot.score, abund, prot.dclass, sep='\t')

if __name__ == '__main__':
    # summary(pipeline())
    # pairwise()
    protwise()

