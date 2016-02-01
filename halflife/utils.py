#!/usr/bin/env python3

def load_ned_data(species):
    """Loads processed decay data from Selbach group"""
    if species == 'mouse':
        filename = 'data/NED_mouse_Abund.txt'
    elif species == 'human':
        filename = 'data/NED_human_Abund.txt'
    with open(filename) as infile:
        data = [line.strip().split('\t') for line in infile]
    header = data[0]
    data = data[1:]
    return header, data

def get_homologs():
    """Returns dictionary mapping mouse homologs entrez <-> uniprot."""
    homologs = {}
    with open('data/homology/corum_mouse_homologs.txt') as infile:
        data = [line.strip().split('\t') for line in infile]
    # data must be sorted in order of sequence identity (high first)
    for line in data:
        original = line[1].split('|')[1]
        uniprot = line[0]
        entrez = line[3]
        # proteins map 1 to 1
        if original not in homologs:
            homologs[original] = [entrez]
        # genes map 1 to multiple
        if entrez not in homologs:
            homologs[entrez] = [uniprot]
        elif uniprot not in homologs[entrez]:
            homologs[entrez].append(uniprot)
    return homologs

def get_uniprot_homologs(rev=False):
    """As above, but exclusively uniprot => mouse uniprot"""
    homologs = {}
    with open('data/homology/corum_mouse_homologs.txt') as infile:
        data = [line.strip().split('\t') for line in infile]
    for line in data:
        original = line[1].split('|')[1]
        uniprot = line[0]
        # Picks first, and subsequently best. Seqid must be in desc order!
        if original not in homologs:
            homologs[original] = uniprot
    if rev:
        homologs = {value: key for key, value in homologs.items()}
    return homologs


###############################################################################
## Following functions for processing blast homology data. Don't mess.

def get_homologs_from_blast():
    """Maps results of BLAST of CORUM subunits to mouse uniprot ids."""
    with open('data/homology/blasts.out') as infile:
        data = [line.strip().split('\t') for line in infile]
    homologs = []
    for line in data:
        line[0] = line[0].split('*')[0]
        if '|' in line[0]:
            line[0] = line[0].split('|')[1]
        # Can modify if neccesary
        if float(line[8]) >= 70.0:
            homologs.append(line[:2] + [line[8]])
    homologs.sort(key=lambda line: (line[1], -float(line[2])))
    return homologs

def map_entrez_to_homologs(homologs):
    """Writes out sorted homology mapping (high seqid first)"""
    with open('data/homology/uniprot_entrez.out') as infile:
        idmap = {l.split()[0]: l.split()[1] for l in infile}
    for line in homologs:
        uniprot = line[0].split('-')[0]
        if uniprot in idmap:
            line.append(idmap[uniprot])
    homologs = [line for line in homologs if len(line) == 4]
    with open('data/homology/corum_mouse_homologs.txt', 'w') as outfile:
        for line in homologs:
            outfile.write('\t'.join(line)+'\n')

