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


def atom_bond_order(molecule: Molecule, atom) -> float:

    total = 0.0

    for bond in molecule.bonds:

        if bond.atom_a is atom or bond.atom_b is atom:
            total += bond.order

    return total


def validate_molecule(molecule: Molecule) -> list[str]:

    errors = []

    for index, atom in enumerate(molecule.atoms):

        if atom.element not in DEFAULT_VALENCY:
            continue

        allowed = DEFAULT_VALENCY[atom.element]
        used = atom_bond_order(molecule, atom)

        # Explicit electrons modify the normal valency rules.
        if atom.charge != 0:
            continue

        if atom.lone_pairs != 0:
            continue

        if atom.unpaired_electrons != 0:
            continue

        if used > allowed:
            errors.append(
                f"Atom {index} ({atom}) exceeds valency: "
                f"{used} > {allowed}"
            )

    return errors


def assert_valid_molecule(molecule: Molecule):

    errors = validate_molecule(molecule)

    if errors:
        message = "Invalid molecule:\n"
        message += "\n".join(errors)
        raise ValueError(message)