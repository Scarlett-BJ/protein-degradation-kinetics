#!/usr/bin/env python3

with open('test.tmp') as f1, open('test1.tmp') as f2:
    d1 = sorted([set(l.strip().split('\t')[1].split(', ')) for l in f1])
    d2 = sorted([set(l.strip().split('\t')[1].split(', ')) for l in f2])
for i in d1:
    if i not in d2:
        print(d1[d1.index(i)], d2[d2.index(i)])