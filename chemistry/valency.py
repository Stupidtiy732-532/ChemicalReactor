"""
OrgReact — Valency and implicit hydrogen inference
"""

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


def bond_order_sum(molecule: Molecule, atom_index: int) -> float:
    return sum(
        bond.order
        for bond in molecule.get_bonds(atom_index)
    )


def infer_implicit_hydrogens(molecule: Molecule):
    """
    Infer missing hydrogen count from bond order.

    Explicit hydrogens already present in the molecule
    are not added again.

    Atoms with explicitly specified electron information
    are left unchanged for now.
    """

    for index, atom in enumerate(molecule.atoms):

        if atom.element not in DEFAULT_VALENCY:
            continue

        if atom.element in {"H", "D", "T"}:
            continue

        if atom.lone_pairs is not None:
            continue

        if atom.unpaired_electrons > 0:
            continue

        if atom.charge != 0:
            continue

        valency = DEFAULT_VALENCY[atom.element]

        used_valency = bond_order_sum(
            molecule,
            index,
        )

        missing = valency - used_valency

        if missing < 0:
            continue

        if missing != int(missing):
            continue

        for _ in range(int(missing)):
            h_index = molecule.add_atom(
                __import__(
                    "chemistry.atoms",
                    fromlist=["Atom"]
                ).Atom("H")
            )

            molecule.add_bond(
                __import__(
                    "chemistry.bonds",
                    fromlist=["Bond"]
                ).Bond(index, h_index)
            )