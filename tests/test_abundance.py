import pytest
import halflife.abundance as ab
import halflife.utils as utils

@pytest.fixture(params=['mouse', 'human'])
def species(request):
    return request.param

def test_abundance_dict(species):
    _, data = utils.load_ned_data(species)
    abdict = ab.abundance_dict(species)
    assert len(data) == len(abdict)
