import pytest
import random
import halflife.coexpression as cx
import os

class TestCoexpressTable(object):

    @pytest.fixture(params=['mouse', 'human', 'mouse_homolog'])
    def coextable(self, request):
        if request.param == 'mouse':
            return cx.CoexpressTable(species='mouse')
        elif request.param == 'human':
            return cx.CoexpressTable(species='human')
        elif request.param == 'mouse_homolog':
            return cx.CoexpressTable(species='mouse', homologs=True)

    @pytest.fixture(params=[random.randint(0, 100)])
    def set_complex(self, request, coextable):
        struc = coextable._corum.strucs[request.param]
        coextable._complex = coextable._corum[struc]

    def test_init(self, coextable):
        attrs = {'_coex', '_homologs', '_corum', '_decay'}
        assert set(dir(coextable)).intersection(attrs) == attrs

    def test_homologise_complex(self, coextable, set_complex):
        """Most data is lost at this stage. Lots of homologs missing"""
        if coextable._homologs:
            entrezids = coextable._homologise_complex()
            assert len(entrezids) <= len(coextable._complex.entrez)
        else:
            with pytest.raises(AssertionError):
                coextable._homologise_complex()

    def test_avg_coexpression(self, coextable, set_complex):
        if coextable._homologs:
            entrezids = coextable._homologise_complex()
            avg_coex = coextable._avg_coexpression(entrezids)
            assert len(avg_coex) == len(entrezids)
        else:
            entrezids = [e[0] for e in coextable._complex.entrez]
            avg_coex = coextable._avg_coexpression(entrezids)
            assert len(avg_coex) == len(entrezids)

    def test_convert_avgcoex_keys(self, coextable, set_complex):
        if not coextable._homologs:
            entrezids = [e[0] for e in coextable._complex.entrez]
            avg_coex1 = coextable._avg_coexpression(entrezids)
            avg_coex2 = coextable._convert_avgcoex_keys(avg_coex1)
            assert len(avg_coex1) == len(avg_coex2)

    def test_convert_homolog_avgcoex_keys(self, coextable, set_complex):
        """It surprises me that this passes! be careful here, check it out"""
        if coextable._homologs:
            entrezids = coextable._homologise_complex()
            avg_coex1 = coextable._avg_coexpression(entrezids)
            avg_coex2 = coextable._convert_homolog_avgcoex_keys(avg_coex1)
            # Should == pass instead of <= ??? implies all ids are present in
            # decay dict?
            assert len(avg_coex2) == len(avg_coex1)

    def test_process_data(self, coextable):
        coextable._corum.strucs = coextable._corum.strucs[:20]
        coextable.process_data()
        assert len(coextable._outdata) > 0

    def test_write_to_file(self, tmpdir, coextable):
        coextable._corum.strucs = coextable._corum.strucs[:20]
        coextable.process_data()
        coextable.write_to_file('../tests/tmp/test_coex.tmp')
        with open('../tests/tmp/test_coex.tmp') as infile:
            data = infile.readlines()
            assert len(data) == len(coextable._outdata) + 1
        os.remove('../tests/tmp/test_coex.tmp')

