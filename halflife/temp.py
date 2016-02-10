#!/usr/bin/env python3

with open('data/coexpression/struc_coex.txt') as infile:
    new_data = [infile.readline()]
    sdict = []
    for line in  infile:
        sline = line.strip().split()
        if sline[2] not in sdict:
            sdict.append(sline[2])
            new_data.append(line)
for line in new_data:
    print(line.strip())