#!/usr/bin/env python3

import pytest
import random

import halflife.utils as ut
import halflife.coexpression as cx
import halflife.abundance as ab
import halflife.tissue_expression as te


class TestUtils(object):
    pass


class TestCoexpression(object):

    @pytest.fixture
    def homolog(self):
        return cx.CoexpressTable(species='mouse', homologs=True)

    @pytest.fixture(params=['mouse', 'human'])
    def non_homolog(self, request):
        return cx.CoexpressTable(species=request.param)

    @pytest.fixture(params=[random.randint(0, 100)])
    def set_complex(self, request, homolog, non_homolog):
        struc1 = homolog._corum.strucs[request.param]
        struc2 = non_homolog._corum.strucs[request.param]
        homolog._complex = homolog._corum[struc1]
        non_homolog._complex = non_homolog._corum[struc2]

    def test_init(self, homolog, non_homolog):
        attrs = {'_coex', '_homologs', '_corum', '_decay'}
        assert set(dir(homolog)).intersection(attrs) == attrs
        assert set(dir(non_homolog)).intersection(attrs) == attrs

    def test_homologise_complex(self, homolog, non_homolog, set_complex):
        """Most data is lost at this stage. Lots of homologs missing"""
        entrezids = homolog._homologise_complex()
        assert len(entrezids) <= len(homolog._complex.entrez)
        with pytest.raises(AssertionError):
            non_homolog._homologise_complex()

    def test_avg_coexpression(self, homolog, non_homolog, set_complex):
        entrezids = homolog._homologise_complex()
        avg_coex = homolog._avg_coexpression(entrezids)
        assert len(avg_coex) == len(entrezids)
        entrezids = [e[0] for e in non_homolog._complex.entrez]
        avg_coex = non_homolog._avg_coexpression(entrezids)
        assert len(avg_coex) == len(entrezids)

    def test_convert_avgcoex_keys(self, non_homolog, set_complex):
        entrezids = [e[0] for e in non_homolog._complex.entrez]
        avg_coex1 = non_homolog._avg_coexpression(entrezids)
        avg_coex2 = non_homolog._convert_avgcoex_keys(avg_coex1)
        assert len(avg_coex1) == len(avg_coex2)

    def test_convert_homolog_avgcoex_keys(self, homolog, set_complex):
        """It surprises me that this passes! be careful here, check it out"""
        entrezids = homolog._homologise_complex()
        avg_coex1 = homolog._avg_coexpression(entrezids)
        avg_coex2 = homolog._convert_homolog_avgcoex_keys(avg_coex1)
        # Should == pass instead of <= ??? implies all ids are present in decay
        assert len(avg_coex2) == len(avg_coex1)

    def test_process_data(self, homolog, non_homolog):
        homolog._corum.strucs = homolog._corum.strucs[:20]
        non_homolog._corum.strucs = non_homolog._corum.strucs[:20]
        homolog.process_data()
        non_homolog.process_data()
        assert len(homolog._outdata) > 0 and len(non_homolog._outdata) > 0

    def test_write_to_file(self):
        # Needs doing! Hard to figure out how temp directories work!
        pass


class TestAbundance(object):
    pass


class TestStructuralDistribution(object):
    pass


class TestTissueExpression(object):
    pass
