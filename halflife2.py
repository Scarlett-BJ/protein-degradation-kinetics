#!/usr/bin/env python3

from ixntools import structure

def get_interfaces():
    cmap = structure.chain_map()
    ints = structure.Interfaces()
    with open('pairwise.tsv') as infile:
        data = [line.split() for line in infile.readlines()]
    for line in data:
        if line[0] not in ints.interfaces:
            continue
        try:
            pair = cmap[line[0]][line[4]], cmap[line[0]][line[5]]
            pair = tuple(sorted(pair))
            print(ints[line[0]][pair])
        except:
            print(line, ints[line[0]], cmap[line[0]])
            break


if __name__ == '__main__':
    get_interfaces()