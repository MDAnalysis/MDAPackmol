#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `mdapackmol` package."""
import os
import pytest
import MDAnalysis as mda

import mdapackmol


@pytest.fixture()
def in_tmpdir(tmpdir):
    os.chdir(str(tmpdir))

    yield str(tmpdir)


HERE = os.path.abspath(os.path.dirname(__file__))
WATER_PDB = os.path.join(HERE, 'water.pdb')
UREA_PDB = os.path.join(HERE, 'urea.pdb')

    
def test_mixture(in_tmpdir):
    water = mda.Universe(WATER_PDB)
    urea = mda.Universe(UREA_PDB)
    
    # PS(ag, number, instructions)
    mixture = mdapackmol.packmol(
        [mdapackmol.PackmolStructure(
            water,
            number=1000,
            instructions=['inside box 0. 0. 0. 40. 40. 40.'],
        ),
         mdapackmol.PackmolStructure(
             urea,
             number=400,
             instructions=['inside box 0. 0. 0. 40. 40. 40.'],
         ),
        ]
    )

    assert len(mixture.atoms) == 1000 * len(water.atoms) + 400 * len(urea.atoms)
