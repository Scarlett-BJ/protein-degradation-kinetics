import pytest
import random
from halflife import structural_distribution as sd
import os

@pytest.fixture(params=['mouse', 'human'])
def species(request):
    return request.param

def test_load_genes(species):
    genes, rblist = sd.load_genes(species)
    strucs = set(i[0] for i in genes.values())
    # Bad test, but simple function so not sure how to test well
    assert len(set(strucs).intersection(rblist)) != set()

def test_load_structure_data():
    strucs = sd.load_structure_data()
    for pdb in strucs:
        assert int(strucs[pdb].total) >= int(strucs[pdb].unique)

def test_write_gene_info(species):
    genes, rblist = sd.load_genes(species)
    strucs = sd.load_structure_data()
    for r in rblist:
        strucs.pop(r)
    tmpfilename = '../tests/tmp/test_structural.tmp'
    with open(tmpfilename, 'w') as outfile:
        header = ['gene', 'struc', 'chn', 'decay.class', 'qtype',
                  'unq', 'tot', '\n']
        outfile.write('\t'.join(header))
        for gene in genes:
            if genes[gene][0] not in strucs:
                continue
            struc = strucs[genes[gene][0]]
            line = [gene, '\t'.join(genes[gene]), struc.qtype,
                    struc.unique, struc.total]
            outfile.write('\t'.join(line)+'\n')
    with open(tmpfilename) as infile:
        header = infile.readline().strip().split('\t')
        data = [line.strip().split('\t') for line in infile]
    for line in data:
        assert len(line) == len(header)
        assert line[1] not in rblist
    os.remove(tmpfilename)
