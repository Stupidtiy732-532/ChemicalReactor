from chemistry.atoms import Atom
from chemistry.bonds import Bond, BondType
from chemistry.molecule import Molecule


class ChemicalParser:

    ELEMENTS = {
        "H", "D", "T",
        "B", "C", "N", "O", "F", "P", "S",
        "Cl", "Br", "I",
    }

    def __init__(self, formula: str):
        self.formula = formula
        self.position = 0
        self.length = len(formula)

        self.molecule = Molecule()

        self.current_atom = None
        self.pending_bond = BondType.SINGLE

        self.branch_stack = []

    def parse(self) -> Molecule:

        while self.position < self.length:

            char = self.formula[self.position]

            if char.isspace():
                self.position += 1
                continue

            # Bond symbols
            if char == "-":
                self.pending_bond = BondType.SINGLE
                self.position += 1
                continue

            if char == "=":
                self.pending_bond = BondType.DOUBLE
                self.position += 1
                continue

            if char == "#":
                self.pending_bond = BondType.TRIPLE
                self.position += 1
                continue

            # Branch opening
            if char == "(":
                if self.current_atom is None:
                    raise ValueError("Branch cannot start without an atom")

                self.branch_stack.append(self.current_atom)
                self.position += 1
                continue

            # Branch closing
            if char == ")":
                if not self.branch_stack:
                    raise ValueError("Unmatched ')'")

                self.current_atom = self.branch_stack.pop()
                self.position += 1
                continue

            # Isotope atom, e.g. [14C]
            if char == "[":
                atom = self.parse_isotope()
                self.connect_atom(atom)
                continue

            # Element atom
            if char.isupper():
                self.parse_element()
                continue

            raise ValueError(
                f"Unexpected character '{char}' at position {self.position}"
            )

        if self.branch_stack:
            raise ValueError("Unclosed branch '('")

        return self.molecule

    def parse_element(self):

        element = self.formula[self.position]
        self.position += 1

        # Two-letter elements: Cl and Br
        if self.position < self.length:

            possible_second = self.formula[self.position]

            candidate = element + possible_second

            if candidate in self.ELEMENTS:
                element = candidate
                self.position += 1

        # Number after element means attached hydrogen/deuterium count:
        #
        # C       -> carbon atom
        # CH3     -> carbon atom + 3 hydrogen atoms
        # CH2     -> carbon atom + 2 hydrogen atoms
        # OH      -> oxygen atom + hydrogen atom
        #
        # A number directly after an element is not atom repetition.

        self.connect_atom(Atom(element=element))

        if element in {"H", "D", "T"}:
            return

        if self.position >= self.length:
            return

        next_char = self.formula[self.position]

        # Explicit hydrogen count: CH3, CH2, CD3, etc.
        if next_char in {"H", "D", "T"}:

            isotope_hydrogen = next_char
            self.position += 1

            count = self.read_number()

            if count is None:
                count = 1

            for _ in range(count):
                hydrogen = Atom(element=isotope_hydrogen)
                self.connect_atom(hydrogen)

    def parse_isotope(self) -> Atom:

        # Skip '['
        self.position += 1

        isotope_start = self.position

        while (
            self.position < self.length
            and self.formula[self.position].isdigit()
        ):
            self.position += 1

        isotope_text = self.formula[isotope_start:self.position]

        if not isotope_text:
            raise ValueError("Missing isotope number")

        isotope = int(isotope_text)

        if self.position >= self.length:
            raise ValueError("Incomplete isotope notation")

        element = self.formula[self.position]
        self.position += 1

        if self.position < self.length:

            candidate = element + self.formula[self.position]

            if candidate in self.ELEMENTS:
                element = candidate
                self.position += 1

        if self.position >= self.length or self.formula[self.position] != "]":
            raise ValueError("Missing closing ']'")

        self.position += 1

        return Atom(
            element=element,
            isotope=isotope,
        )

    def connect_atom(self, atom: Atom):

        self.molecule.add_atom(atom)

        if self.current_atom is not None:

            bond = Bond(
                atom_a=self.current_atom,
                atom_b=atom,
                order=self.get_bond_order(self.pending_bond),
                bond_type=self.pending_bond,
            )

            self.molecule.add_bond(bond)

        self.current_atom = atom
        self.pending_bond = BondType.SINGLE

    def get_bond_order(self, bond_type: BondType) -> float:

        if bond_type == BondType.SINGLE:
            return 1.0

        if bond_type == BondType.DOUBLE:
            return 2.0

        if bond_type == BondType.TRIPLE:
            return 3.0

        return 1.0

    def read_number(self):

        start = self.position

        while (
            self.position < self.length
            and self.formula[self.position].isdigit()
        ):
            self.position += 1

        if start == self.position:
            return None

        return int(self.formula[start:self.position])