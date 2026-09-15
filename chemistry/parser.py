from chemistry.atoms import Atom
from chemistry.bonds import Bond, BondType
from chemistry.molecule import Molecule


class ChemicalParser:

    ELEMENTS = {
        "H", "D", "T",
        "B", "C", "N", "O", "F", "P", "S",
        "Cl", "Br", "I",
    }

    BOND_TYPES = {
        "-": BondType.SINGLE,
        "=": BondType.DOUBLE,
        "#": BondType.TRIPLE,
    }

    BOND_ORDERS = {
        BondType.SINGLE: 1.0,
        BondType.DOUBLE: 2.0,
        BondType.TRIPLE: 3.0,
    }

    HYDROGEN_ELEMENTS = {"H", "D", "T"}

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

            if char in self.BOND_TYPES:
                self.pending_bond = self.BOND_TYPES[char]
                self.position += 1
                continue

            if char == "(":
                self.open_branch()
                continue

            if char == ")":
                self.close_branch()
                continue

            if char == "[":
                self.connect_atom(self.parse_isotope())
                continue

            if char.isupper():
                self.parse_element()
                continue

            raise ValueError(
                f"Unexpected character '{char}' "
                f"at position {self.position}"
            )

        if self.branch_stack:
            raise ValueError("Unclosed branch '('")

        return self.molecule

    # ============================================================
    # BRANCHES
    # ============================================================

    def open_branch(self):
        if self.current_atom is None:
            raise ValueError("Branch cannot start without an atom")

        self.branch_stack.append(self.current_atom)
        self.position += 1

    def close_branch(self):
        if not self.branch_stack:
            raise ValueError("Unmatched ')'")

        self.current_atom = self.branch_stack.pop()
        self.position += 1

    # ============================================================
    # ELEMENTS
    # ============================================================

    def parse_element(self):
        element = self.read_element()

        self.connect_atom(Atom(element=element))

        if element in self.HYDROGEN_ELEMENTS:
            return

        if self.position >= self.length:
            return

        if self.formula[self.position] not in self.HYDROGEN_ELEMENTS:
            return

        hydrogen_element = self.formula[self.position]
        self.position += 1

        count = self.read_number()
        count = 1 if count is None else count

        for _ in range(count):
            self.connect_atom(Atom(element=hydrogen_element))

    def read_element(self) -> str:
        element = self.formula[self.position]
        self.position += 1

        if self.position < self.length:
            candidate = element + self.formula[self.position]

            if candidate in self.ELEMENTS:
                self.position += 1
                return candidate

        if element not in self.ELEMENTS:
            raise ValueError(
                f"Unknown element '{element}' "
                f"at position {self.position - 1}"
            )

        return element

    # ============================================================
    # ISOTOPES
    # ============================================================

    def parse_isotope(self) -> Atom:
        self.position += 1  # '['

        isotope = self.read_number()

        if isotope is None:
            raise ValueError("Missing isotope number")

        if self.position >= self.length:
            raise ValueError("Incomplete isotope notation")

        element = self.read_element()

        if self.position >= self.length or self.formula[self.position] != "]":
            raise ValueError("Missing closing ']'")

        self.position += 1

        return Atom(
            element=element,
            isotope=isotope,
        )

    # ============================================================
    # ATOMS / BONDS
    # ============================================================

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

    @classmethod
    def get_bond_order(cls, bond_type: BondType) -> float:
        return cls.BOND_ORDERS[bond_type]

    # ============================================================
    # NUMBERS
    # ============================================================

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