#!/usr/bin/env python3

import numpy as np

def load_pairwise():
    with open('data/pairwise.tsv') as infile:
        data = [line.split() for line in infile]
    cdict = {}
    for line in data[1:]:
        if line[0] not in cdict:
            cdict[line[0]] = [(line[6], line[9], line[-1])]
        else:
            cdict[line[0]].append((line[6], line[9], line[-1]))
    success = 0
    trials = 0
    for c in cdict:
        ints = cdict[c]
        nface = {float(i[2]) for i in ints if
                (i[1] == 'N' or i[0] == 'N') and float(i[2]) > 200}
        therest = {float(i[2]) for i in ints if float(i[2]) > 200}
        if len(therest) > 0 and len(nface) > 0:
            trials += 1
            if np.mean(list(therest)) < np.mean(list(nface)):
                success += 1
    print(success, trials)
    print(len(cdict))

load_pairwise()