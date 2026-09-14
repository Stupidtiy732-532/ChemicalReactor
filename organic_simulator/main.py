from collections import Counter
import re


# ============================================================
# DATA STRUCTURES
# ============================================================

class Atom:
    def __init__(self, atom_id, element, hydrogens=0):
        self.id = atom_id
        self.element = element
        self.hydrogens = hydrogens


class Bond:
    def __init__(self, atom1, atom2, order=1):
        self.atom1 = atom1
        self.atom2 = atom2
        self.order = order


class FunctionalGroup:
    def __init__(
        self,
        name,
        atom_ids,
        locant=None,
        details="",
    ):
        self.name = name
        self.atom_ids = atom_ids
        self.locant = locant
        self.details = details


class Substituent:
    def __init__(
        self,
        name,
        atom_ids,
        parent_locant=None,
    ):
        self.name = name
        self.atom_ids = atom_ids
        self.parent_locant = parent_locant


class Molecule:
    def __init__(self, original_input):
        self.original_input = original_input

        self.atoms = []
        self.bonds = []

        self.parent_chain = []
        self.functional_groups = []
        self.substituents = []

    def atom_by_id(self, atom_id):
        for atom in self.atoms:
            if atom.id == atom_id:
                return atom

        raise ValueError(f"Atom {atom_id} does not exist.")

    def neighbours(self, atom_id):
        result = []

        for bond in self.bonds:
            if bond.atom1 == atom_id:
                result.append(bond.atom2)

            elif bond.atom2 == atom_id:
                result.append(bond.atom1)

        return result

    def bond_between(self, atom1, atom2):
        for bond in self.bonds:
            if {
                bond.atom1,
                bond.atom2,
            } == {atom1, atom2}:
                return bond

        return None

    def total_element_count(self, element):
        """
        Counts explicit atoms and implicit hydrogens.
        """

        total = 0

        for atom in self.atoms:
            if atom.element == element:
                total += 1

            if element == "H":
                total += atom.hydrogens

        return total

    @property
    def formula(self):
        counts = Counter()

        for atom in self.atoms:
            counts[atom.element] += 1
            counts["H"] += atom.hydrogens

        result = []

        # Hill-system order: C, H, then alphabetical elements.
        if counts["C"] > 0:
            result.append("C")

            if counts["C"] > 1:
                result.append(str(counts["C"]))

        if counts["H"] > 0:
            result.append("H")

            if counts["H"] > 1:
                result.append(str(counts["H"]))

        for element in sorted(counts):
            if element in {"C", "H"}:
                continue

            if counts[element] <= 0:
                continue

            result.append(element)

            if counts[element] > 1:
                result.append(str(counts[element]))

        return "".join(result)

    @property
    def molar_mass(self):
        atomic_masses = {
            "H": 1.008,
            "C": 12.011,
            "N": 14.007,
            "O": 15.999,
            "F": 18.998,
            "Cl": 35.45,
            "Br": 79.904,
            "I": 126.904,
        }

        mass = 0.0

        for atom in self.atoms:
            if atom.element not in atomic_masses:
                raise ValueError(
                    f"No atomic mass available for {atom.element}"
                )

            mass += atomic_masses[atom.element]

            mass += atom.hydrogens * atomic_masses["H"]

        return mass


# ============================================================
# TOKENIZATION
# ============================================================

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


# ============================================================
# STRUCTURE PARSER
# ============================================================

class StructureParser:
    """
    Supports basic open-chain aliphatic notation.

    Examples:

        CH3-CH2-OH
        CH3-CH(OH)-CH3
        CH3-CH(CH3)-CH2-OH
        CH2=CH-CH3
        CH3-COOH
        CH3-CHO
        CH3-CO-CH3
        CH3-CH2-Cl
    """

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

        atom = Atom(
            atom_id=atom_id,
            element=element,
            hydrogens=hydrogens,
        )

        self.atoms.append(atom)

        return atom_id

    def add_bond(self, atom1, atom2, order=1):
        self.bonds.append(
            Bond(
                atom1=atom1,
                atom2=atom2,
                order=order,
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
        if token == "C":
            hydrogens = 0
        else:
            hydrogens = self.hydrogen_count(token)

        self.attach_atom(
            element="C",
            hydrogens=hydrogens,
        )

    def parse_special_group(self, token):
        # Hydroxyl group: -OH
        if token == "OH":
            self.attach_atom(
                element="O",
                hydrogens=1,
            )
            return

        # Carboxylic acid group: -COOH
        if token == "COOH":
            carbon_id = self.attach_atom(
                element="C",
                hydrogens=0,
            )

            carbonyl_oxygen = self.new_atom(
                element="O",
                hydrogens=0,
            )

            hydroxyl_oxygen = self.new_atom(
                element="O",
                hydrogens=1,
            )

            self.add_bond(
                carbon_id,
                carbonyl_oxygen,
                order=2,
            )

            self.add_bond(
                carbon_id,
                hydroxyl_oxygen,
                order=1,
            )

            return

        # Aldehyde group: -CHO
        if token == "CHO":
            carbon_id = self.attach_atom(
                element="C",
                hydrogens=1,
            )

            oxygen_id = self.new_atom(
                element="O",
                hydrogens=0,
            )

            self.add_bond(
                carbon_id,
                oxygen_id,
                order=2,
            )

            return

        # Ketone carbonyl fragment: -CO-
        if token == "CO":
            carbon_id = self.attach_atom(
                element="C",
                hydrogens=0,
            )

            oxygen_id = self.new_atom(
                element="O",
                hydrogens=0,
            )

            self.add_bond(
                carbon_id,
                oxygen_id,
                order=2,
            )

            return

        # Amino group
        if token.startswith("NH"):
            hydrogens = self.hydrogen_count(token)

            self.attach_atom(
                element="N",
                hydrogens=hydrogens,
            )

            return

        # Halogens
        if token in {"Cl", "Br", "F", "I"}:
            self.attach_atom(
                element=token,
                hydrogens=0,
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
            original_input=self.original_input
        )

        molecule.atoms = self.atoms
        molecule.bonds = self.bonds

        return molecule


# ============================================================
# MOLECULE ANALYZER
# ============================================================

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

        for neighbour_id in molecule.neighbours(
            carbon_id
        ):
            neighbour = molecule.atom_by_id(
                neighbour_id
            )

            if neighbour.element == "C":
                result.append(neighbour_id)

        return result

    def find_parent_chain(self, molecule):
        """
        Finds the longest carbon path.

        This is currently the basic longest-chain algorithm.
        Functional-group priority will be added afterward.
        """

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
                current_id=carbon_id,
                visited={carbon_id},
                path=[carbon_id],
            )

        molecule.parent_chain = best_path

    def find_functional_groups(self, molecule):
        molecule.functional_groups = []

        for atom in molecule.atoms:

            # Alcohol and carbonyl oxygen
            if atom.element == "O":
                neighbours = molecule.neighbours(
                    atom.id
                )

                # Alcohol O-H
                if atom.hydrogens == 1:
                    carbon_neighbours = []

                    for neighbour_id in neighbours:
                        neighbour = molecule.atom_by_id(
                            neighbour_id
                        )

                        if neighbour.element == "C":
                            carbon_neighbours.append(
                                neighbour_id
                            )

                    if carbon_neighbours:
                        carbon_id = carbon_neighbours[0]

                        carbon_count = len(
                            self.carbon_neighbours(
                                molecule,
                                carbon_id,
                            )
                        )

                        if carbon_count == 1:
                            group_name = "primary alcohol"

                        elif carbon_count == 2:
                            group_name = "secondary alcohol"

                        elif carbon_count == 3:
                            group_name = "tertiary alcohol"

                        else:
                            group_name = "alcohol"

                        molecule.functional_groups.append(
                            FunctionalGroup(
                                name=group_name,
                                atom_ids=[
                                    carbon_id,
                                    atom.id,
                                ],
                            )
                        )

                # Carbonyl oxygen
                for neighbour_id in neighbours:
                    bond = molecule.bond_between(
                        atom.id,
                        neighbour_id,
                    )

                    if bond is None:
                        continue

                    if bond.order != 2:
                        continue

                    neighbour = molecule.atom_by_id(
                        neighbour_id
                    )

                    if neighbour.element == "C":
                        molecule.functional_groups.append(
                            FunctionalGroup(
                                name="carbonyl",
                                atom_ids=[
                                    neighbour_id,
                                    atom.id,
                                ],
                            )
                        )

            # Halogen group
            elif atom.element in {
                "Cl",
                "Br",
                "F",
                "I",
            }:
                neighbours = molecule.neighbours(
                    atom.id
                )

                if neighbours:
                    molecule.functional_groups.append(
                        FunctionalGroup(
                            name="halo group",
                            atom_ids=[
                                neighbours[0],
                                atom.id,
                            ],
                        )
                    )

    def find_substituents(self, molecule):
        molecule.substituents = []

        parent_atoms = set(
            molecule.parent_chain
        )

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
                    name="alkyl substituent",
                    atom_ids=[atom.id],
                    parent_locant=parent_locant,
                )
            )


# ============================================================
# DISPLAY
# ============================================================

def display_molecule(molecule):
    print()
    print("=" * 65)
    print("MOLECULE ANALYSIS")
    print("=" * 65)

    print(f"Input:      {molecule.original_input}")
    print(f"Formula:    {molecule.formula}")
    print(
        f"Molar mass: {molecule.molar_mass:.3f} g/mol"
    )

    print()
    print("Atoms:")

    for atom in molecule.atoms:
        print(
            f"  Atom {atom.id}: "
            f"{atom.element}, "
            f"implicit H = {atom.hydrogens}"
        )

    print()
    print("Bonds:")

    for bond in molecule.bonds:
        print(
            f"  {bond.atom1} - "
            f"{bond.atom2}, "
            f"order = {bond.order}"
        )

    print()
    print("Parent carbon chain:")

    if molecule.parent_chain:
        chain_text = " - ".join(
            f"C{index + 1}"
            for index in range(
                len(molecule.parent_chain)
            )
        )

        print(f"  {chain_text}")
        print(
            f"  Length: "
            f"{len(molecule.parent_chain)} carbon atoms"
        )

    else:
        print("  No carbon parent chain found.")

    print()
    print("Functional groups:")

    if molecule.functional_groups:
        for group in molecule.functional_groups:
            print(
                f"  {group.name}: "
                f"atoms {group.atom_ids}"
            )
    else:
        print("  None detected.")

    print()
    print("Substituents:")

    if molecule.substituents:
        for substituent in molecule.substituents:
            print(
                f"  {substituent.name}: "
                f"atoms {substituent.atom_ids}, "
                f"parent locant = "
                f"{substituent.parent_locant}"
            )
    else:
        print("  None detected.")

    print("=" * 65)


# ============================================================
# CLI
# ============================================================

def main():
    print("Organic Chemistry Simulator — Version 1")
    print("Open-chain aliphatic molecule analyzer")
    print("Type 'exit' to quit.")

    while True:
        print()

        text = input("Molecule > ").strip()

        if text.lower() == "exit":
            print("Exiting.")
            break

        if not text:
            continue

        try:
            parser = StructureParser(text)
            molecule = parser.parse()

            analyzer = MoleculeAnalyzer()
            analyzer.analyze(molecule)

            display_molecule(molecule)

        except Exception as error:
            print()
            print(f"Error: {error}")


if __name__ == "__main__":
    main()