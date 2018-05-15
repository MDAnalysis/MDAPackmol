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

@pytest.fixture()
def water():
    return mda.Universe(WATER_PDB)


@pytest.fixture()
def urea():
    return mda.Universe(UREA_PDB)


class TestMixture(object):
    @staticmethod
    @pytest.fixture()
    def mixture(in_tmpdir, water, urea):
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

        return mixture

    def test_mixture_size(self, mixture, water, urea):
        assert len(mixture.atoms) == 1000 * len(water.atoms) + 400 * len(urea.atoms)

    def test_residue_size(self, mixture, water, urea):
        for res in mixture.residues[:3]:
            assert len(res.atoms) == len(water.atoms)

        for res in mixture.residues[1000:1003]:
            assert len(res.atoms) == len(urea.atoms)
    
    def test__resnames(self, mixture, water, urea):
        for res in mixture.residues[:3]:
            for atom_a, atom_b in zip(res.atoms, water.atoms):
                assert atom_a.resname == atom_b.resname

        for res in mixture.residues[1000:1003]:
            for atom_a, atom_b in zip(res.atoms, urea.atoms):
                assert atom_a.resname == atom_b.resname

def test_bonds(urea, water):
    urea.atoms.guess_bonds()
    water.atoms.guess_bonds()

    mixture = mdapackmol.packmol(
        [mdapackmol.PackmolStructure(
            water,
            number=50,
            instructions=['inside box 0. 0. 0. 40. 40. 40.'],
        ),
         mdapackmol.PackmolStructure(
             urea,
             number=50,
             instructions=['inside box 0. 0. 0. 40. 40. 40.'],
         ),
        ]
    )

    assert hasattr(mixture, 'bonds')
    assert len(mixture.bonds) == len(water.bonds) * 50 + len(urea.bonds) * 50
