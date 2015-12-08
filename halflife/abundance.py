#!/usr/bin/env python3

import numpy as np
from expression import paxdb, coexpressdb
from collections import namedtuple
from itertools import combinations
from scipy.stats import binom_test

def load_data(filename):
    with open(filename) as infile:
        data = [line.strip().split('\t') for line in infile]
    header = data[0]
    data = data[1:]
    return header, data

def load_mouse_prot_map():
    data_dict = {}
    with open('data/mouse_genes_ned.txt') as infile:
        data = [line.strip().split(',') for line in infile.readlines()]
    for line in data:
        data_dict[line[0]] = {'ensg': line[2], 'gene': line[1],
                              'ensp': line[3]}
    return data_dict

def get_complexes(header, data):
    """Returns a dictionary of complexes. Each complex contains a list of
    namedtuples with information about the complex and proteins.
    """
    pmap = load_mouse_prot_map()
    Info = namedtuple('Info', [str(i.strip('.')) for i in header[1:]])
    cdict = {}
    for line in data:
        line[4] = pmap[line[4]]['ensp']
        if line[0] not in cdict:
            cdict[line[0]] = [Info(*line[1:])]
        else:
            cdict[line[0]].append(Info(*line[1:]))
    return cdict

def relative_abundance():
    """Within complexes, are N proteins more abundant?"""
    header, data = load_data('data/protwise_trimers.tsv')
    cdict = get_complexes(header, data)
    trials = 0
    success = 0
    avg_distance = []
    for pdb in cdict:
        comp = cdict[pdb]
        nvals = [float(p.ab) for p in comp if p.dcls == 'N' or p.dcls == 'U']
        evals = [float(p.ab) for p in comp if p.dcls == 'E']
        if len(nvals) != 0 and len(evals) != 0:
            trials += 1
            emean = np.mean(evals)
            nmean = np.mean(nvals)
            avg_distance.append(nmean-emean)
            print(nmean, emean)
            if nmean > emean:
                success += 1
    print(success, trials, binom_test(success, trials), np.mean(avg_distance))

# def get_tissue_expression(threshold1, threshold2):
#     gene_map = load_mouse_prot_map()
#     header, data = load_data('data/protwise.tsv')
#     meta = paxdb.get_metadata('10090')
#     for i in ['WHOLE_ORGANISM', 'CELL_LINE']:
#         meta.pop(i)
#     abundances = []
#     for i in meta:
#         filename = meta[i]['filename']
#         abundances.append(paxdb.Abundances(filename, 0, 10000000))
#     for line in data:
#         ensp = gene_map[line[4]]['ensp']
#         tcount = 0
#         for tissue in abundances:
#             if ensp not in tissue.members:
#                 continue
#             if tissue[ensp] < threshold1 or tissue[ensp] > threshold2:
#                 continue
#                 tcount += 1
#         line.append(tcount)
#     return data

def tissue_expression_percomp(t1, t2):
    header, data = load_data('data/protwise.tsv')
    cdict = get_complexes(header, data)
    meta = paxdb.get_metadata('10090')
    for i in ['WHOLE_ORGANISM', 'CELL_LINE']:
        if i == 'WHOLE_ORGANISM':
            worg = paxdb.Abundances(meta[i]['filename'], 0)
            meta.pop(i)
        else:
            meta.pop(i)
    tlist = sorted(meta)
    all_tissues = [paxdb.Abundances(meta[t]['filename'], 0) for t in tlist]
    success, trials = 0, 0
    for pdb in cdict:
        # nabund = [t[p.prot] for p in cdict[pdb] if p.dcls == 'N'
        #           for t in all_tissues if p.prot in t.members]
        # eabund = [t[p.prot] for p in cdict[pdb] if p.dcls != 'N'
        #           for t in all_tissues if p.prot in t.members]
        # if len(nabund) == 0 or len(eabund) == 0:
        #     continue
        # elif np.mean(nabund) > np.mean(eabund):
        #     continue
        evals = []
        nvals = []
        for prot in cdict[pdb]:
            tcount = 0
            for tis in all_tissues:
                if prot.prot in tis.members:
                    tcount += 1
            if tcount == 8 or tcount == 0:
                continue
            if prot.dcls == 'N':
                nvals.append(tcount)
            elif prot.dcls != 'N':
                evals.append(tcount)
        if len(evals) == 0 or len(nvals) == 0:
            continue
        # if np.mean(nvals) == np.mean(evals):
        #     continue
        trials += 1
        if np.mean(nvals) >= np.mean(evals):
            success += 1
    print(success, trials, binom_test(success, trials))



# def tissue_expression_percomp(data):
#     """Are N proteins more widely expressed than E proteins"""
#     cdict = {}
#     avg_distance = []
#     for line in data:
#         if line[0] not in cdict:
#             cdict[line[0]] = [(line[7], int(line[8]))]
#         else:
#             cdict[line[0]].append((line[7], int(line[8])))
#     trials = 0
#     success = 0
#     for pdb in cdict:
#         comp = cdict[pdb]
#         evals = []
#         nvals = []
#         for prot in comp:
#             # if prot[1] == 0:
#             #     continue
#             if prot[0] == 'E' or prot[0] == 'U':
#                 evals.append(prot[1])
#             elif prot[0] == 'N':
#                 nvals.append(prot[1])
#         if len(evals) > 0 and len(nvals) > 0:
#             # if np.mean(nvals) == np.mean(evals):
#             #     continue
#             trials += 1
#             emean = np.mean(evals)
#             nmean = np.mean(nvals)
#             avg_distance.append(((nmean-emean)**2)**0.5)
#             if np.mean(nvals) >= np.mean(evals):
#                 success += 1
#     print(success, trials, binom_test(success, trials))

# def tissue_expression(data):
#     """Are N proteins more widely expressed than E proteins"""
#     nvals = []
#     evals = []
#     uvals = []
#     for line in data:
#         if line[7] == 'E':
#             evals.append(line[8])
#         elif line[7] == 'N':
#             nvals.append(line[8])
#         elif line[7] == 'U':
#             uvals.append(line[8])
#     print(np.mean(nvals), np.mean(uvals), np.mean(evals))

def abundance_vs_tcount_paxdb():
    with open('data/mouse_genes_ned.txt') as infile:
        prots = {l.split(',')[3]: l.strip().split(',')[-2:] for l in infile}
    meta = paxdb.get_metadata('10090')
    for i in ['WHOLE_ORGANISM', 'CELL_LINE']:
        meta.pop(i)
    tlist = sorted(meta)
    all_tissues = [paxdb.Abundances(meta[t]['filename'], 0) for t in tlist]
    header = ['prot'] + tlist + ['score', 'tcount']
    print('\t'.join(header))
    for prot in prots:
        if 'ENSMUSP' not in prot:
            continue
        tcount = 0
        tis_abunds = []
        for abunds in all_tissues:
            if prot in abunds.members:
                tcount += 1
                tis_abunds.append(str(abunds[prot]))
            else:
                tis_abunds.append('NA')
        print(prot, '\t'.join(tis_abunds), prots[prot][0], tcount, sep='\t')

def abundance_vs_tcount_relabund():
    with open('data/mouse_genes_ned.txt') as infile:
        prots = {l.split(',')[3]: l.strip().split(',')[-2:] for l in infile}
    meta = paxdb.get_metadata('10090')
    for i in ['WHOLE_ORGANISM', 'CELL_LINE']:
        meta.pop(i)
    all_tissues = [paxdb.Abundances(meta[t]['filename'], 0) for t in meta]
    header = ['prot', 'abund', 'score', 'tcount']
    print('\t'.join(header))
    for prot in prots:
        if 'ENSMUSP' not in prot:
            continue
        tcount = 0
        for abunds in all_tissues:
            if prot in abunds.members:
                tcount += 1
        print(prot, prots[prot][1], prots[prot][0], tcount, sep='\t')
        # print(sprots)

if __name__ == '__main__':
    # relative_abundance()
    # data = get_tissue_expression()
    # tissue_expression_percomp(0, 10)


    # tissue_expression(data)
    # abundance_vs_tcount_paxdb()
    # abundance_vs_tcount_relabund()
    meta = paxdb.get_metadata('10090')
    for i in sorted(meta):
        print(i, meta[i]['coverage'], meta[i]['score'])

