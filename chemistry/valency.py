from chemistry.atoms import Atom
from chemistry.bonds import Bond, BondType


VALENCY = {
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


def bond_order_sum(molecule, atom) -> float:

    total = 0.0

    for bond in molecule.bonds:

        if bond.atom_a is atom or bond.atom_b is atom:
            total += float(bond.order)

    return total


def complete_hydrogens(molecule):

    additions = []

    for atom in list(molecule.atoms):

        element = atom.element

        if element not in VALENCY:
            continue

        if element in {"H", "D", "T"}:
            continue

        # Only neutral atoms are handled by this basic rule.
        if getattr(atom, "charge", 0) != 0:
            continue

        used = bond_order_sum(molecule, atom)
        required = VALENCY[element]

        missing = required - used

        if missing < 0:
            raise ValueError(
                f"Invalid valency for {element}: "
                f"bond order {used} exceeds {required}"
            )

        if not missing.is_integer():
            raise ValueError(
                f"Non-integer hydrogen count for {element}: {missing}"
            )

        for _ in range(int(missing)):
            additions.append(atom)

    for parent in additions:

        hydrogen = Atom(element="H")

        molecule.add_atom(hydrogen)

        molecule.add_bond(
            Bond(
                atom_a=parent,
                atom_b=hydrogen,
                order=1.0,
                bond_type=BondType.SINGLE,
            )
        )

    return molecule