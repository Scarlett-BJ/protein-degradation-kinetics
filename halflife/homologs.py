#!/usr/bin/env python3

def get_homologs():
    with open('data/homology/blasts.out') as infile:
        data = [line.strip().split('\t') for line in infile]
    homologs = []
    for line in data:
        line[0] = line[0].split('*')[0]
        if '|' in line[0]:
            line[0] = line[0].split('|')[1]
        if float(line[8]) >= 50.0:
            homologs.append(line[:2] + [line[8]])
    homologs.sort(key=lambda line: (line[1], -float(line[2])))
    return homologs

def map_entrez_to_homologs(homologs):
    with open('data/homology/uniprot_entrez.out') as infile:
        idmap = {l.split()[0]: l.split()[1] for l in infile}
    for line in homologs:
        uniprot = line[0].split('-')[0]
        if uniprot in idmap:
            line.append(idmap[uniprot])
    homologs = [line for line in homologs if len(line) == 4]
    for line in homologs:
        print('\t'.join(line))



if __name__ == '__main__':
    homologs = get_homologs()
    map_entrez_to_homologs(homologs)
