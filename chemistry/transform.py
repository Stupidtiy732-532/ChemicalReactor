from chemistry.molecule import Molecule
from chemistry.bonds import BondType


def find_bonds(molecule: Molecule, order: float):

    return [
        bond
        for bond in molecule.bonds
        if bond.order == order
    ]


def reduce_double_bond(molecule: Molecule):

    for bond in molecule.bonds:

        if bond.order == 2.0:

            bond.order = 1.0
            bond.bond_type = BondType.SINGLE

            return molecule

    raise ValueError("No double bond found")


def reduce_triple_bond(molecule: Molecule):

    for bond in molecule.bonds:

        if bond.order == 3.0:

            bond.order = 2.0
            bond.bond_type = BondType.DOUBLE

            return molecule

    raise ValueError("No triple bond found")


def add_atom_to_bond(
    molecule: Molecule,
    bond,
    element: str,
):

    from chemistry.atoms import Atom
    from chemistry.bonds import Bond

    atom_a = bond.atom_a
    atom_b = bond.atom_b

    new_atom = Atom(element=element)
    molecule.add_atom(new_atom)

    # Remove original bond
    molecule.bonds.remove(bond)

    # Split original bond
    first_bond = Bond(
        atom_a=atom_a,
        atom_b=new_atom,
        order=1.0,
        bond_type=BondType.SINGLE,
    )

    second_bond = Bond(
        atom_a=new_atom,
        atom_b=atom_b,
        order=1.0,
        bond_type=BondType.SINGLE,
    )

    molecule.add_bond(first_bond)
    molecule.add_bond(second_bond)

    return molecule