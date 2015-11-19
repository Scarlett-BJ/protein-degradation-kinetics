#!/usr/bin/env python3

from ixntools import structure
from numpy import mean

def get_interfaces():
    cmap = structure.chain_map()
    ints = structure.Interfaces()
    with open('pairwise.tsv') as infile:
        data = [line.split() for line in infile.readlines()]
    for line in data:
        pair = cmap[line[0]][line[4]][0], cmap[line[0]][line[5]][0]
        pair = tuple(sorted(pair))
        line.append(ints[line[0]][pair])
        line = [str(i) for i in line]
    return data

def analyse(filename):
    with open(filename) as infile:
        data = [line.split() for line in infile]
    strucs = {line[0]: {'N': [], 'E': []} for line in data}
    for line in data:
        if line[8] != 'E' and line[9] != 'E':
            strucs[line[0]]['N'].append(float(line[10]))
        elif line[8] == 'E' or line[9] == 'E':
            strucs[line[0]]['E'].append(float(line[10]))
    trials = 0
    successes = 0
    for s in strucs:
        if strucs[s]['N'] == [] or strucs[s]['E'] == []:
            continue
        if mean(strucs[s]['N']) > mean(strucs[s]['E']):
            successes += 1
        trials += 1
    print(successes, trials)



if __name__ == '__main__':
    data = get_interfaces()
    with open('test.tmp', 'w') as outfile:
        for line in data:
            outfile.write('\t'.join([str(i) for i in line])+'\n')
    analyse('test.tmp')
