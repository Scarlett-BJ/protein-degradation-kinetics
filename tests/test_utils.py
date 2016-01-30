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
