"""Cool stuff which is maybe useful"""

from scipy.constants import N_A

def molecules_for_target_density(
    existing_molecules,
    solvent_molecule,
    target_density,
    boxsize):
    """Calculate how many solvent molecules to add to reach target density

    Parameters
    ----------
    existing_molecules : dict
       mapping of AtomGroups to number in box
    solvent_molecule : AtomGroup
       solvent molecule you want to add
    target_density : float
       target system density in kg/m3
    boxsize : 3 floats
       boxsize in each dimension in Angstrom

    Returns
    -------
    nsolvent, density
      number of solvent molecules to add, resulting density

    Example
    -------
    To find how many water molecules to solvate our protein to a density of
    985 kg/m^3.  We load AtomGroups "protein" and "water" (making sure that
    the mass is correct for these).  We specify that there will be 1 protein
    molecule in the box, and the solvent is the water AtomGroup.  We then
    pass the density and size of our box (20x20x20 is this example)::

       >>> molecules_for_target_density({protein: 1}, water,
                                        985.0, [20.0, 20.0, 20.0])
    """
    # system volume
    vol = boxsize[0] * boxsize[1] * boxsize[2] * 10 ** -30

    target_mass = target_density * vol  # kg

    existing_mass = sum(mol.total_mass() * quantity
                       for mol, quantity in existing_molecules.items())
    # from g/mol to kg
    existing_mass /= 1000 * N_A

    required_mass = target_mass - existing_mass

    solvent_mass = solvent_molecule.total_mass() / (1000 * N_A)

    nreq = int(required_mass / solvent_mass)
    # calculate resulting density as a check
    actual_density = (existing_mass + nreq * solvent_mass) / vol

    return nreq, actual_density
