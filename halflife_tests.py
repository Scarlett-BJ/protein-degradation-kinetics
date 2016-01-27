#!/usr/bin/env python3

import os
import pytest
import halflife as hl

def test_load_ned_data():
    filenames = ['data/NED_mouse_Abund.txt',
                 'data/NED_human_Abund.txt']
    for filename in filenames:
        header, data = hl.coexpression.load_ned_data(filename)
        assert len(data) != 0 and len(header) == len(data[0])

def test_homologs():
    entrez_homologs, uniprot_homologs = hl.coexpression.get_homologs()
    assert len(entrez_homologs) != 0 and len(uniprot_homologs) != 0


if __name__ == '__main__':

    os.chdir('halflife')
    pytest.main(['-v', '../halflife_tests.py'])