# analyzer.py

import re
from molecule import Molecule


GROUP_TOKENS = [
    "COOH",
    "CHO",
    "CH3",
    "CH2",
    "NH2",
    "OH",
    "CO",
    "Cl",
    "Br",
    "CH",
    "C",
    "O",
    "N",
    "F",
    "I",
    "S",
    "P",
]


def tokenize(text):
    text = text.replace(" ", "")
    tokens = []
    i = 0

    while i < len(text):
        found = False

        for token in GROUP_TOKENS:
            if text.startswith(token, i):
                tokens.append(token)
                i += len(token)
                found = True
                break

        if not found:
            raise ValueError(
                f"Cannot understand structure near: {text[i:]}"
            )

    return tokens


def parse_group(molecule, token):
    """
    Returns:
        attachment_atom,
        created_atoms
    """

    if token == "CH3":
        atom = molecule.add_atom("C", 3)
        return atom, [atom]

    if token == "CH2":
        atom = molecule.add_atom("C", 2)
        return atom, [atom]

    if token == "CH":
        atom = molecule.add_atom("C", 1)
        return atom, [atom]

    if token == "C":
        atom = molecule.add_atom("C", 0)
        return atom, [atom]

    if token == "OH":
        atom = molecule.add_atom("O", 1)
        return atom, [atom]

    if token == "NH2":
        atom = molecule.add_atom("N", 2)
        return atom, [atom]

    if token in {"F", "Cl", "Br", "I"}:
        atom = molecule.add_atom(token, 0)
        return atom, [atom]

    if token == "O":
        atom = molecule.add_atom("O", 0)
        return atom, [atom]

    if token == "N":
        atom = molecule.add_atom("N", 0)
        return atom, [atom]

    if token == "S":
        atom = molecule.add_atom("S", 0)
        return atom, [atom]

    if token == "P":
        atom = molecule.add_atom("P", 0)
        return atom, [atom]

    if token == "CHO":
        carbon = molecule.add_atom("C", 1)
        oxygen = molecule.add_atom("O", 0)
        molecule.add_bond(carbon, oxygen, 2)
        return carbon, [carbon, oxygen]

    if token == "CO":
        carbon = molecule.add_atom("C", 0)
        oxygen = molecule.add_atom("O", 0)
        molecule.add_bond(carbon, oxygen, 2)
        return carbon, [carbon, oxygen]

    if token == "COOH":
        carbon = molecule.add_atom("C", 0)
        carbonyl_oxygen = molecule.add_atom("O", 0)
        hydroxyl_oxygen = molecule.add_atom("O", 1)

        molecule.add_bond(carbon, carbonyl_oxygen, 2)
        molecule.add_bond(carbon, hydroxyl_oxygen, 1)

        return carbon, [
            carbon,
            carbonyl_oxygen,
            hydroxyl_oxygen,
        ]

    raise ValueError(f"Unsupported group: {token}")


def parse_structure(text):
    tokens = tokenize(text)

    molecule = Molecule()
    previous_attachment = None

    for token in tokens:
        attachment, _ = parse_group(molecule, token)

        if previous_attachment is not None:
            molecule.add_bond(
                previous_attachment,
                attachment,
                1,
            )

        previous_attachment = attachment

    return molecule


def atom_is_carbon(molecule, atom_id):
    return molecule.atoms[atom_id].element == "C"


def detect_functional_groups(molecule):
    groups = []

    for atom in molecule.atoms.values():
        if atom.element != "C":
            continue

        neighbors = molecule.neighbors(atom.id)

        carbon_neighbors = [
            n for n, _ in neighbors
            if molecule.atoms[n].element == "C"
        ]

        oxygen_neighbors = [
            (n, order)
            for n, order in neighbors
            if molecule.atoms[n].element == "O"
        ]

        halogen_neighbors = [
            molecule.atoms[n].element
            for n, _ in neighbors
            if molecule.atoms[n].element in {"F", "Cl", "Br", "I"}
        ]

        nitrogen_neighbors = [
            n for n, _ in neighbors
            if molecule.atoms[n].element == "N"
        ]

        double_oxygen = any(order == 2 for _, order in oxygen_neighbors)

        if double_oxygen:
            single_oxygen = [
                n for n, order in oxygen_neighbors
                if order == 1
            ]

            if single_oxygen:
                oxygen = molecule.atoms[single_oxygen[0]]

                if oxygen.hydrogens:
                    groups.append("carboxylic acid")

            elif atom.hydrogens:
                groups.append("aldehyde")

            elif len(carbon_neighbors) == 2:
                groups.append("ketone")

            continue

        if oxygen_neighbors:
            oxygen_id, order = oxygen_neighbors[0]
            oxygen = molecule.atoms[oxygen_id]

            if oxygen.hydrogens:
                carbon_degree = len(carbon_neighbors)

                if carbon_degree == 1:
                    groups.append("primary alcohol")
                elif carbon_degree == 2:
                    groups.append("secondary alcohol")
                elif carbon_degree == 3:
                    groups.append("tertiary alcohol")
                else:
                    groups.append("alcohol")

        if halogen_neighbors:
            for halogen in halogen_neighbors:
                groups.append(f"haloalkane ({halogen})")

        if nitrogen_neighbors:
            groups.append("amine")

    for bond in molecule.bonds:
        a = molecule.atoms[bond.a]
        b = molecule.atoms[bond.b]

        if a.element == "C" and b.element == "C":
            if bond.order == 2:
                groups.append("alkene")
            elif bond.order == 3:
                groups.append("alkyne")

    if not groups:
        groups.append("hydrocarbon or unsupported structure")

    return sorted(set(groups))


def analyze(text):
    molecule = parse_structure(text)

    return {
        "molecule": molecule,
        "formula": molecule.formula(),
        "molar_mass": molecule.molar_mass(),
        "structure": molecule.condensed(),
        "functional_groups": detect_functional_groups(molecule),
    }