def identify_functional_groups(molecule):

    groups = []

    bond_orders = molecule.get_bond_orders()

    # Alkene
    if 2 in bond_orders:
        groups.append("Alkene")

    # Alkyne
    if 3 in bond_orders:
        groups.append("Alkyne")

    # Alcohol
    oxygen_atoms = [
        atom_id
        for atom_id, atom in molecule.atoms.items()
        if atom.element == "O"
        and atom.hydrogens > 0
    ]

    if oxygen_atoms:
        groups.append("Hydroxyl group / possible alcohol")

    return groups