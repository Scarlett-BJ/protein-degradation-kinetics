import pytest
import halflife.utils as utils

@pytest.fixture(params=['mouse', 'human'])
def species(request):
    return request.param

def test_load_ned_data(species):
    header, data = utils.load_ned_data(species)
    for colname in ['proID', 'Rel. abund', 'def']:
        assert colname in header
    assert len(data[0]) == len(header)

def test_homolog_dict_builders():
    with open('data/homology/corum_mouse_homologs.txt') as infile:
        data = [line.strip().split('\t') for line in infile]
    # data must be sorted in order of sequence identity (high first)
    previous_original = data[0][1].split('|')[1]
    previous_seqid = float(data[0][2])
    for line in data:
        original = line[1].split('|')[1]
        seqid = float(line[2])
        assert seqid >= 70.0
        if original == previous_original:
            assert previous_seqid >= seqid
