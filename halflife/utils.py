#!/usr/bin/env python3

def load_ned_data(filename):
    """Loads processed decay data from Selbach group"""
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
