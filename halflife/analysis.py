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
        classes = {i[0] for i in cdict[c]}.union({i[1] for i in cdict[c]})
        if 'E' in classes:
            ints = cdict[c]
            nface = {float(i[2]) for i in ints if
                    (i[1] == 'E' or i[0] == 'E') and float(i[2]) > 200}
            therest = {float(i[2]) for i in ints if float(i[2]) > 200} - nface
            if len(nface) > 0 and len(therest) > 0:
                trials += 1
                if np.mean(list(therest)) > np.mean(list(nface)):
                    success += 1
    print(success, trials)
    print(len(cdict))

def load_protwise():
    with open('data/protwise.tsv') as infile:
        data = [line.split() for line in infile]
    cdict = {}
    for line in data[1:]:
        if line[0] not in cdict:
            cdict[line[0]] = [line[1:]]
        else:
            cdict[line[0]].append(line[1:])
    s = 0
    t = 0
    for i in cdict:
        es = []
        ns = []
        for p in cdict[i]:
            if p[-1] == 'N':
                ns.append(float(p[-2]))
            elif p[-1] == 'E':
                es.append(float(p[-2]))
        if es != [] and ns != []:
            t += 1
            if np.mean(ns) > np.mean(es):
                s += 1
    print(s, t)

load_pairwise()