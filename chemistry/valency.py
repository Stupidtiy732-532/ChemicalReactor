from chemistry.atoms import Atom
from chemistry.molecule import Molecule


DEFAULT_VALENCY = {
    "H": 1,
    "D": 1,
    "T": 1,

    "C": 4,
    "N": 3,
    "O": 2,

    "F": 1,
    "Cl": 1,
    "Br": 1,
    "I": 1,

    "S": 2,
    "P": 3,
}


def bond_order_sum(molecule: Molecule, atom: Atom) -> float:

    total = 0.0

    for bond in molecule.get_bonds():

        if bond.atom_a is atom:
            total += bond.order

        elif bond.atom_b is atom:
            total += bond.order

    return total


def complete_hydrogens(molecule: Molecule) -> Molecule:

    atoms_to_add = []

    for atom in molecule.atoms:

        # Do not automatically complete hydrogen atoms.
        if atom.element in {"H", "D", "T"}:
            continue

        # Unknown elements are left untouched.
        if atom.element not in DEFAULT_VALENCY:
            continue

        # Charged/radical/lone-pair atoms need special rules later.
        if atom.charge != 0:
            continue

        if atom.lone_pairs != 0:
            continue

        if atom.unpaired_electrons != 0:
            continue

        valency = DEFAULT_VALENCY[atom.element]
        used_valency = bond_order_sum(molecule, atom)

        missing_hydrogens = int(valency - used_valency)

        if missing_hydrogens < 0:
            raise ValueError(
                f"Invalid valency for {atom}: "
                f"used bond order = {used_valency}, "
                f"allowed valency = {valency}"
            )

        for _ in range(missing_hydrogens):
            atoms_to_add.append(atom)

    for parent_atom in atoms_to_add:

        hydrogen = Atom(element="H")
        molecule.add_atom(hydrogen)

        from chemistry.bonds import Bond, BondType

        bond = Bond(
            atom_a=parent_atom,
            atom_b=hydrogen,
            order=1.0,
            bond_type=BondType.SINGLE,
        )

        molecule.add_bond(bond)

    return molecule