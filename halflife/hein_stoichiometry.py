#!/usr/bin/env python3

import ixntools as ix
from numpy import mean

def load_hein():
    filename = ix.set_data_dir('human/hein/hein.csv')
    with open(filename) as infile:
        data = [line.strip().split(',') for line in infile]
    header = data[0]
    data = data[1:]
    for line in data:
        line[3] = set(line[3].split(';'))
        line[2] = set(line[2].split(';'))
    return header, data

def load_ned():
    with open('data/NED_human.txt') as infile:
        data = [line.split() for line in infile]
    header = data[0]
    data = data[1:]
    return header, data

def compare_proteins():
    header, data = load_ned()
    neds = set(line[0] for line in data if line[-1] == 'NED')
    eds = set(line[0] for line in data if line[-1] == 'ED')
    uns = set(line[0] for line in data if line[-1] == 'UN')
    header, data = load_hein()
    ned_stoichs = []
    ed_stoichs = []
    un_stoichs = []
    ncount = 0
    ecount = 0
    ucount = 0
    for line in data:
        if line[14] == '+':
            continue
        if line[2].intersection(neds) != set():
            ned_stoichs.append(float(line[11]))
            if line[13] == '+':
                ncount += 1
        elif line[2].intersection(eds) != set():
            ed_stoichs.append(float(line[11]))
            if line[13] == '+':
                ecount += 1
        elif line[2].intersection(uns) != set():
            un_stoichs.append(float(line[11]))
            if line[13] == '+':
                ucount += 1
    # print(len(ned_stoichs), len(ed_stoichs), len(un_stoichs))
    # print(10 ** mean(ned_stoichs), 10 ** mean(ed_stoichs), 10 ** mean(un_stoichs))
        # print(line[6])
    print(len(ed_stoichs) - ecount, len(ned_stoichs) - ncount, len(un_stoichs) - ucount, ecount, ncount, ucount)
    # print(header)

def count_interactions():
    header, data = load_ned()
    neds = {line[0]: [] for line in data if line[-1] == 'NED'}
    eds = {line[0]: [] for line in data if line[-1] == 'ED'}
    uns = {line[0]: [] for line in data if line[-1] == 'UN'}
    header, data = load_hein()

    def ribo(x): return bool('RPS' in x or 'RPL' in x)

    def get_prot(*args): return [list(line[2].intersection(i)) for i in args]

    def get_interactors(c, d):
        if c != []:
            d[c[0]].append(float(line[11]))

    def pop_empty(d):
        for i in list(d):
            if d[i] == []:
                d.pop(i)

    for line in data:
        if ribo(line[4]) or ribo(line[5]):
            continue
        if line[14] == '+':
            continue
        if line[11] == 'NaN':
            continue
        if line[13] != '+':
            continue
        nedi, edi, uni = get_prot(neds, eds, uns)
        for c, d in [(nedi, neds), (edi, eds), (uni, uns)]:
            get_interactors(c, d)
    for d in [neds, eds, uns]:
        pop_empty(d)
    nedlens = [10 ** mean(i) for i in neds.values()]
    edlens = [10 ** mean(i) for i in eds.values()]
    unlens = [10 ** mean(i) for i in uns.values()]
    print(mean(nedlens), mean(edlens), mean(unlens))
    print(len(neds), len(eds), len(uns))
if __name__ == '__main__':
    count_interactions()
