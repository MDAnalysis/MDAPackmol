===========
MDA Packmol
===========


An MDAnalysis wrapper around Packmol


* Free software: GNU General Public License v3

Features
--------

* Allows combining MDAnalysis and Packmol

* Preserves the topology information (bonds etc) of your system after Packmol


Usage Example
-------------

.. code-block:: python

   import MDAnalysis as mda
   import mdapackmol
   
   # load individual molecule files
   water = mda.Universe('water.pdb')
   urea = mda.Universe('urea.pdb')
   
   # call Packmol with MDAnalysis objects as arguments
   # the 'instructions' allow for any valid Packmol commands
   system = mdapackmol.packmol(
       [mdapackmol.PackmolStructure(
           water, number=1000,
           instructions=['inside box 0. 0. 0. 40. 40. 40.']),
        mdapackmol.PackmolStructure(
           urea, number=400,
           instructions=['inside box 0. 0. 0. 40. 40. 40.'])]
   )
   
   # the returned system is a MDAnalysis Universe
   # with all topology information from building blocks retained
   # which can then be saved into any format
   # eg to Lammps data file:
   system.atoms.write('topology.data')
