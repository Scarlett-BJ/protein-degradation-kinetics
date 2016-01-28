#!/usr/bin/env python3

import pytest
from halflife import coexpression as cx
import random

def test_load_ned_data():
    filenames = ['data/NED_mouse_Abund.txt',
                 'data/NED_human_Abund.txt']
    for filename in filenames:
        header, data = cx.load_ned_data(filename)
        assert len(data) != 0 and len(header) == len(data[0])

def test_homologs():
    homologs = cx.get_homologs()
    assert len(homologs) != 0


class TestCoexpressTable(object):

    @pytest.fixture
    def homologous(self):
        return cx.CoexpressTable(species='mouse', homologs=True)

    @pytest.fixture(params=['mouse', 'human'])
    def non_homologous(self, request):
        return cx.CoexpressTable(species=request.param)

    @pytest.fixture(params=[1])
    def set_complex(self, request, homologous, non_homologous):
        struc1 = homologous._corum.strucs[request.param]
        struc2 = non_homologous._corum.strucs[request.param]
        homologous._complex = homologous._corum[struc1]
        non_homologous._complex = non_homologous._corum[struc2]

    def test_init(self, homologous, non_homologous):
        attrs = {'_coex', '_homologs', '_corum', '_decay'}
        assert set(dir(homologous)).intersection(attrs) == attrs
        assert set(dir(non_homologous)).intersection(attrs) == attrs

    def test_homologise_complex(self, homologous, non_homologous, set_complex):
        entrezids = homologous._homologise_complex()
        assert len(entrezids) <= len(homologous._complex.entrez)
        with pytest.raises(AssertionError):
            non_homologous._homologise_complex()

    def test_avg_coexpression(self, homologous, non_homologous, set_complex):
        entrezids = homologous._homologise_complex()
        avg_coex = homologous._avg_coexpression(entrezids)
        assert len(avg_coex) == len(entrezids)
        entrezids = [e[0] for e in non_homologous._complex.entrez]
        avg_coex = non_homologous._avg_coexpression(entrezids)
        assert len(avg_coex) == len(entrezids)

    def test_convert_avgcoex_keys(self, non_homologous, set_complex):
        entrezids = [e[0] for e in non_homologous._complex.entrez]
        avg_coex1 = non_homologous._avg_coexpression(entrezids)
        avg_coex2 = non_homologous._convert_avgcoex_keys(avg_coex1)
        assert len(avg_coex1) == len(avg_coex2)

    def test_convert_homolog_avgcoex_keys(self, homologous, set_complex):
        entrezids = homologous._homologise_complex()
        avg_coex1 = homologous._avg_coexpression(entrezids)
        avg_coex2 = homologous._convert_homolog_avgcoex_keys(avg_coex1)
        assert avg_coex1 != avg_coex2
        assert len(avg_coex1) == len(avg_coex2)
