import pytest
import halflife.abundance as ab
import halflife.utils as utils
from ixntools import dbloader

@pytest.fixture(params=['mouse', 'human', 'mouse_homologs'])
def species(request):
    return request.param


def test_abundance_dict(species):
    if species == 'mouse_homologs':
        return
    _, data = utils.load_ned_data(species)
    abdict = ab.abundance_dict(species)
    assert len(data) == len(abdict)

def test_maximise_subunit_data(species, homologs=False):
    if species == 'mouse_homologs':
        species = 'mouse'
        homologs = utils.get_uniprot_homologs()
    abdict = ab.abundance_dict(species)
    if homologs:
        homologs = utils.get_uniprot_homologs()
        corum = dbloader.LoadCorum(version='core')
    else:
        corum = dbloader.LoadCorum(species, 'core')
    success, trials = 0, 0
    for struc in corum.strucs:
        comp = corum[struc].uniprot
        nvals, evals, uvals = ab.maximise_subunit_data(comp, abdict, homologs)
        assert len(nvals + evals + uvals) <= len(comp)
