==========
MDAPackmol
==========

An MDAnalysis_ wrapper around Packmol_

* Allows combining MDAnalysis and Packmol

* Preserves the topology information (bonds etc) of your system after Packmol

* Free software: GNU General Public License v3

.. _MDAnalysis: https://www.mdanalysis.org
.. _Packmol: http://m3g.iqm.unicamp.br/packmol/home.shtml

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


Citing
------

If you find mdapackmol useful for you, please cite the following sources:

 * L Martinez, R Andrade, E G Birgin, J M Martinez, "Packmol: A package for building initial configurations for molecular dynamics simulations". Journal of Computational Chemistry, 30, 2157-2164, 2009. 
 
 * R J Gowers, M Linke, J Barnoud, T J E Reddy, M N Melo, S L Seyler, D L Dotson, J Domanski, S Buchoux, I M Kenney, and O Beckstein. "MDAnalysis: A Python package for the rapid analysis of molecular dynamics simulations." In S. Benthall and S. Rostrup, editors, Proceedings of the 15th Python in Science Conference, pages 102-109, Austin, TX, 2016.
