import re

from molecule import (
    Molecule,
    FunctionalGroup,
    Substituent,
    Atom,
    Bond,
)


TOKEN_PATTERN = re.compile(
    r"CH\d*|"
    r"NH\d*|"
    r"COOH|"
    r"CHO|"
    r"CO|"
    r"OH|"
    r"Cl|"
    r"Br|"
    r"[A-Z]|"
    r"\(|\)|=|#|-"
)


def tokenize_structure(text):
    text = text.replace(" ", "")

    tokens = TOKEN_PATTERN.findall(text)
    reconstructed = "".join(tokens)

    if reconstructed != text:
        raise ValueError(
            f"Unsupported notation: {text}"
        )

    return tokens


class StructureParser:
    def __init__(self, text):
        self.original_input = text
        self.tokens = tokenize_structure(text)
        self.position = 0

        self.atoms = []
        self.bonds = []

        self.previous_atom = None
        self.pending_bond_order = 1

    def new_atom(self, element, hydrogens=0):
        atom_id = len(self.atoms) + 1

        self.atoms.append(
            Atom(
                atom_id,
                element,
                hydrogens,
            )
        )

        return atom_id

    def add_bond(self, atom1, atom2, order=1):
        self.bonds.append(
            Bond(
                atom1,
                atom2,
                order,
            )
        )

    def attach_atom(self, element, hydrogens=0):
        atom_id = self.new_atom(
            element,
            hydrogens,
        )

        if self.previous_atom is not None:
            self.add_bond(
                self.previous_atom,
                atom_id,
                self.pending_bond_order,
            )

        self.previous_atom = atom_id
        self.pending_bond_order = 1

        return atom_id

    @staticmethod
    def hydrogen_count(token):
        match = re.search(r"H(\d*)", token)

        if not match:
            return 0

        number = match.group(1)

        if number == "":
            return 1

        return int(number)

    def parse_carbon(self, token):
        hydrogens = 0

        if token != "C":
            hydrogens = self.hydrogen_count(token)

        self.attach_atom(
            "C",
            hydrogens,
        )

    def parse_special_group(self, token):
        if token == "OH":
            self.attach_atom(
                "O",
                1,
            )
            return

        if token == "COOH":
            carbon_id = self.attach_atom(
                "C",
                0,
            )

            carbonyl_oxygen = self.new_atom(
                "O",
                0,
            )

            hydroxyl_oxygen = self.new_atom(
                "O",
                1,
            )

            self.add_bond(
                carbon_id,
                carbonyl_oxygen,
                2,
            )

            self.add_bond(
                carbon_id,
                hydroxyl_oxygen,
                1,
            )

            return

        if token == "CHO":
            carbon_id = self.attach_atom(
                "C",
                1,
            )

            oxygen_id = self.new_atom(
                "O",
                0,
            )

            self.add_bond(
                carbon_id,
                oxygen_id,
                2,
            )

            return

        if token == "CO":
            carbon_id = self.attach_atom(
                "C",
                0,
            )

            oxygen_id = self.new_atom(
                "O",
                0,
            )

            self.add_bond(
                carbon_id,
                oxygen_id,
                2,
            )

            return

        if token.startswith("NH"):
            hydrogens = self.hydrogen_count(token)

            self.attach_atom(
                "N",
                hydrogens,
            )

            return

        if token in {"Cl", "Br", "F", "I"}:
            self.attach_atom(
                token,
                0,
            )
            return

        raise ValueError(
            f"Unsupported token: {token}"
        )

    def parse(self):
        branch_stack = []

        while self.position < len(self.tokens):
            token = self.tokens[self.position]

            if token == "-":
                self.pending_bond_order = 1

            elif token == "=":
                self.pending_bond_order = 2

            elif token == "#":
                self.pending_bond_order = 3

            elif token == "(":
                if self.previous_atom is None:
                    raise ValueError(
                        "A branch cannot begin the molecule."
                    )

                branch_stack.append(
                    self.previous_atom
                )

            elif token == ")":
                if not branch_stack:
                    raise ValueError(
                        "Unmatched closing parenthesis."
                    )

                self.previous_atom = branch_stack.pop()

            elif token.startswith("CH") or token == "C":
                self.parse_carbon(token)

            else:
                self.parse_special_group(token)

            self.position += 1

        if branch_stack:
            raise ValueError(
                "Unmatched opening parenthesis."
            )

        molecule = Molecule(
            self.original_input
        )

        molecule.atoms = self.atoms
        molecule.bonds = self.bonds

        return molecule


class MoleculeAnalyzer:
    def analyze(self, molecule):
        self.find_parent_chain(molecule)
        self.find_functional_groups(molecule)
        self.find_substituents(molecule)

    def carbon_atoms(self, molecule):
        return [
            atom
            for atom in molecule.atoms
            if atom.element == "C"
        ]

    def carbon_neighbours(self, molecule, carbon_id):
        result = []

        for neighbour_id in molecule.neighbours(carbon_id):
            neighbour = molecule.atom_by_id(neighbour_id)

            if neighbour.element == "C":
                result.append(neighbour_id)

        return result

    def find_parent_chain(self, molecule):
        carbon_ids = {
            atom.id
            for atom in self.carbon_atoms(molecule)
        }

        best_path = []

        def dfs(current_id, visited, path):
            nonlocal best_path

            if len(path) > len(best_path):
                best_path = path.copy()

            for neighbour_id in self.carbon_neighbours(
                molecule,
                current_id,
            ):
                if neighbour_id not in visited:
                    visited.add(neighbour_id)
                    path.append(neighbour_id)

                    dfs(
                        neighbour_id,
                        visited,
                        path,
                    )

                    path.pop()
                    visited.remove(neighbour_id)

        for carbon_id in carbon_ids:
            dfs(
                carbon_id,
                {carbon_id},
                [carbon_id],
            )

        molecule.parent_chain = best_path

    def find_functional_groups(self, molecule):
        molecule.functional_groups = []

        for atom in molecule.atoms:
            if atom.element == "O":
                neighbours = molecule.neighbours(atom.id)

                if atom.hydrogens == 1:
                    carbon_neighbours = []

                    for neighbour_id in neighbours:
                        neighbour = molecule.atom_by_id(neighbour_id)

                        if neighbour.element == "C":
                            carbon_neighbours.append(neighbour_id)

                    if carbon_neighbours:
                        carbon_id = carbon_neighbours[0]

                        carbon_count = len(
                            self.carbon_neighbours(
                                molecule,
                                carbon_id,
                            )
                        )

                        if carbon_count == 1:
                            name = "primary alcohol"

                        elif carbon_count == 2:
                            name = "secondary alcohol"

                        elif carbon_count == 3:
                            name = "tertiary alcohol"

                        else:
                            name = "alcohol"

                        molecule.functional_groups.append(
                            FunctionalGroup(
                                name,
                                [carbon_id, atom.id],
                            )
                        )

                for neighbour_id in neighbours:
                    bond = molecule.bond_between(
                        atom.id,
                        neighbour_id,
                    )

                    if bond is None or bond.order != 2:
                        continue

                    neighbour = molecule.atom_by_id(neighbour_id)

                    if neighbour.element == "C":
                        carbon_id = neighbour_id

                        carbon_neighbours = molecule.neighbours(
                            carbon_id
                        )

                        has_oh = False

                        for other_id in carbon_neighbours:
                            other = molecule.atom_by_id(other_id)

                            if (
                                other.element == "O"
                                and other.id != atom.id
                                and other.hydrogens == 1
                            ):
                                has_oh = True

                        if has_oh:
                            name = "carboxylic acid"

                        else:
                            carbon_carbon_neighbours = [
                                other_id
                                for other_id in carbon_neighbours
                                if molecule.atom_by_id(
                                    other_id
                                ).element == "C"
                            ]

                            if len(carbon_carbon_neighbours) == 1:
                                name = "aldehyde"

                            else:
                                name = "ketone"

                        molecule.functional_groups.append(
                            FunctionalGroup(
                                name,
                                [carbon_id, atom.id],
                            )
                        )

            elif atom.element in {"Cl", "Br", "F", "I"}:
                neighbours = molecule.neighbours(atom.id)

                if neighbours:
                    molecule.functional_groups.append(
                        FunctionalGroup(
                            "halo group",
                            [neighbours[0], atom.id],
                        )
                    )

    def find_substituents(self, molecule):
        molecule.substituents = []

        parent_atoms = set(molecule.parent_chain)

        for atom in molecule.atoms:
            if atom.element != "C":
                continue

            if atom.id in parent_atoms:
                continue

            carbon_neighbours = self.carbon_neighbours(
                molecule,
                atom.id,
            )

            parent_neighbours = [
                neighbour_id
                for neighbour_id in carbon_neighbours
                if neighbour_id in parent_atoms
            ]

            if not parent_neighbours:
                continue

            parent_atom_id = parent_neighbours[0]

            parent_locant = (
                molecule.parent_chain.index(
                    parent_atom_id
                ) + 1
            )

            molecule.substituents.append(
                Substituent(
                    "alkyl substituent",
                    [atom.id],
                    parent_locant,
                )
            )