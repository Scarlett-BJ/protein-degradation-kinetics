#!/usr/bin/env python3

from ixntools import structure
from numpy import mean

def get_interfaces():
    cmap = structure.chain_map()
    ints = structure.Interfaces()
    with open('pairwise.tsv') as infile:
        data = [line.split() for line in infile.readlines()]
    for line in data:
        sints = ints[line[0]]
        for pair in sints:
            smap = cmap[line[0]]
            if pair[0] not in smap or pair[1] not in smap:
                continue
            npair = tuple(sorted((smap[pair[0]], smap[pair[1]])))
            if npair == (line[4], line[5]):
                line.append(sints[pair])
                break
        print('\t'.join([str(i) for i in line]))
    return data


def analyse(filename):
    with open(filename) as infile:
        data = [line.split() for line in infile]
    strucs = {line[0]: {'N': [], 'E': []} for line in data}
    for line in data:
        if len(line) != 11:
            continue
        if line[8] != 'E' and line[9] != 'E' and (line[8] != 'U' and line[9] != 'U'):
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
    # data = get_interfaces()
    # print(data)
    # with open('test.tmp', 'w') as outfile:
    #     for line in data:
    #         outfile.write('\t'.join([str(i) for i in line])+'\n')
    analyse('test.tmp')
