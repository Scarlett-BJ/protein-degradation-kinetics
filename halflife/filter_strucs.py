#!/usr/bin/env python3

def load_pfam(*args):
    strucs = {}
    with open('data/structural/pdb_pfam_mapping.txt') as infile:
        for line in infile:
            for query in args:
                if query in line:
                    line = line.split()
                    if line[0] not in strucs:
                        strucs[line[0]] = [line[1]]
                    else:
                        strucs[line[0]].append(line[1])
                    break
    strucs = {struc: set(strucs[struc]) for struc in strucs}
    return strucs
