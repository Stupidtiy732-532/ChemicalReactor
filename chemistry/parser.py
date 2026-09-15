from chemistry.atoms import Atom
from chemistry.bonds import Bond, BondType
from chemistry.molecule import Molecule


class ChemicalParser:

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

            # Branch start
            if char == "(":
                self.branch_stack.append(self.current_atom)
                self.position += 1
                continue

            # Branch end
            if char == ")":
                self.current_atom = self.branch_stack.pop()
                self.position += 1

                # Handle repetition after branch: (...)2
                repeat_count = self.read_number()

                if repeat_count is not None:
                    self.repeat_last_branch(repeat_count)

                continue

            # Isotope notation: [14C]
            if char == "[":
                atom = self.parse_isotope_atom()
                self.connect_atom(atom)
                continue

            # Element symbol
            if char.isupper():

                element = char
                self.position += 1

                # Two-letter element symbols
                if self.position < self.length:
                    next_char = self.formula[self.position]

                    if next_char.islower():
                        element += next_char
                        self.position += 1

                atom = Atom(element=element)

                # Explicit hydrogen/count notation: CH3
                count = self.read_number()

                if count is None:
                    count = 1

                self.connect_atom(atom)

                # Add repeated atoms, e.g. H3
                if element not in {"H", "D", "T"} and count > 1:
                    for _ in range(count - 1):
                        repeated_atom = Atom(element=element)
                        self.connect_atom(repeated_atom)

                continue

            # Ignore whitespace
            if char.isspace():
                self.position += 1
                continue

            raise ValueError(
                f"Unexpected character '{char}' at position {self.position}"
            )

        return self.molecule

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

    def parse_isotope_atom(self):

        self.position += 1

        isotope_start = self.position

        while (
            self.position < self.length
            and self.formula[self.position].isdigit()
        ):
            self.position += 1

        isotope_text = self.formula[isotope_start:self.position]

        if not isotope_text:
            raise ValueError("Isotope number missing")

        isotope = int(isotope_text)

        if self.position >= self.length:
            raise ValueError("Incomplete isotope atom")

        element = self.formula[self.position]
        self.position += 1

        if self.position < self.length:
            next_char = self.formula[self.position]

            if next_char.islower():
                element += next_char
                self.position += 1

        if self.position >= self.length or self.formula[self.position] != "]":
            raise ValueError("Missing closing ']' in isotope notation")

        self.position += 1

        return Atom(
            element=element,
            isotope=isotope,
        )

    def repeat_last_branch(self, repeat_count: int):

        if not self.branch_stack:
            raise ValueError("Branch repetition has no valid branch context")

        branch_root = self.current_atom

        for _ in range(repeat_count - 1):
            clone = Atom(
                element=branch_root.element,
                isotope=branch_root.isotope,
                charge=branch_root.charge,
                lone_pairs=branch_root.lone_pairs,
                unpaired_electrons=branch_root.unpaired_electrons,
                label=branch_root.label,
            )

            self.molecule.add_atom(clone)

            bond = Bond(
                atom_a=self.current_atom,
                atom_b=clone,
                order=1.0,
                bond_type=BondType.SINGLE,
            )

            self.molecule.add_bond(bond)

            self.current_atom = clone