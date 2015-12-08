#!/usr/bin/env python3

from expression import paxdb

def load_strucs():
    with open('data/mouse_genes_ned.txt') as infile:
        prots = {l.split(',')[3]: l.strip().split(',')[-1] for l in infile}
    meta = paxdb.get_metadata('10090')
    for i in ['WHOLE_ORGANISM', 'CELL_LINE']:
        meta.pop(i)
    header = ['prot'] + sorted(meta) + ['tcount']
    print('\t'.join(header))
    for prot in prots:
        if 'ENSMUSP' not in prot:
            continue
        tcount = 0
        tis_abunds = []
        for tissue in sorted(meta):
            abunds = paxdb.Abundances(meta[tissue]['filename'], 0)
            if prot in abunds.members:
                tcount += 1
                tis_abunds.append(str(abunds[prot]))
            else:
                tis_abunds.append('NA')
        print(prot, '\t'.join(tis_abunds), tcount, sep='\t')


load_strucs()
# help(paxdb.Abundances)