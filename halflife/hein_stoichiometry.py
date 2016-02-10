#!/usr/bin/env python3

import ixntools as ix
from numpy import mean
from halflife import utils

def load_hein():
    """Loads dataset 3? from 3d interactome paper. Each line represents a
    single interaction between bait and prey.
    """
    filename = ix.set_data_dir('human/hein/hein.csv')
    with open(filename) as infile:
        data = [line.strip().split(',') for line in infile]
    header = data[0]
    data = data[1:]
    # Turn list of redundant proteins into set for searching.
    for line in data:
        line[2] = set(line[2].split(';'))  # bait
        line[3] = set(line[3].split(';'))  # prey
    return header, data

def load_corum_proteins(species):
    core = ix.dbloader.LoadCorum(species.title(), 'core')
    uniprots = []
    # Filter only proteins from CORUM strucs with x number of subunits.
    for s in list(core._strucs.values()):
        if len(s.uniprot) > 10 and len(s.uniprot) < 15:
            uniprots += [u for sublist in s.uniprot for u in sublist]
    return set(uniprots)

def load_structural():
    pass

def NED_core_interactor_test(protein_subset=None):
    header, data = utils.load_ned_data('human')
    # data = [line for line in data if line[0] in protein_subset]
    neds = set(line[-2] for line in data if line[-3] == 'NED')
    eds = set(line[-2] for line in data if line[-3] == 'ED')
    header, data = load_hein()
    baits = set(protein for line in data for protein in line[2])
    # Select set of proteins common to both NED data and stoich data.
    # baits = baits.intersection(neds.union(eds))
    # neds = neds.intersection(baits)
    # eds = eds.intersection(baits)
    # Test likelihood of finding NED proteins in Core
    print(len(neds), len(eds))
    ned_core, ed_core, ned_noncore, ed_noncore = 0, 0, 0, 0
    for line in data:
        # Ignore self-interactors
        if line[14] == '+' or line[15] == '+':
            continue
        if line[2].intersection(neds) != set():
            if line[13] == '+':
                ned_core += 1
            else:
                ned_noncore += 1
        elif line[2].intersection(eds) != set():
            if line[13] == '+':
                ed_core += 1
            else:
                ed_noncore += 1
    print(ned_core, ed_core, ned_noncore, ed_noncore)

def per_complex():
    core = ix.dbloader.LoadCorum('human'.title(), 'core')
    corum_uniprots = load_corum_proteins('human')
    # NED data
    header, data = load_ned()
    data = [line for line in data if line[0] in corum_uniprots]
    neds = set(line[0] for line in data if line[-1] == 'NED')
    eds = set(line[0] for line in data if line[-1] == 'ED')
    # Hein data
    header, data = load_hein()
    core_prots = []
    for line in data:
        # Ignore self-interactors
        if line[14] == '+' or line[15] == '+':
            continue
        if line[13] == '+':
            prots = list(line[2].intersection(neds))
            if prots != [] and (prots[0] in neds or prots[0] in eds):
                core_prots.append(prots[0])
    core_prots = set(core_prots)
    for s in core.strucs:
        for p in core[s].uniprot:
            if p[0] in neds and p[0] in core_prots:
                print('yes')
            elif p[0] in neds and p[0] in core_prots:
                print('no')
    print(len(neds.intersection(core_prots)))
    print(len(eds.intersection(core_prots)))

def main():
    NED_core_interactor_test()

if __name__ == '__main__':
    main()
