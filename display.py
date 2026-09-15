def display_molecule(molecule):

    print("\nMOLECULAR STRUCTURE")
    print("-" * 45)

    print(f"Formula: {molecule.formula()}")

    print("\nAtoms:")

    for atom in molecule.atoms.values():

        print(
            f"  {atom.id}: "
            f"{atom.element}"
            f"(H={atom.hydrogens}, "
            f"charge={atom.charge}, "
            f"radical={atom.radical_electrons}, "
            f"lone pairs={atom.lone_pairs})"
        )

    print("\nBonds:")

    for bond in molecule.bonds:

        print(
            f"  {bond.atom1} -- "
            f"{bond.order} -- {bond.atom2}"
            f"    "
            f"(σ={bond.sigma_bonds}, "
            f"π={bond.pi_bonds})"
        )