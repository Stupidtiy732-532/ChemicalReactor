from chemistry.atoms import Atom
from chemistry.bonds import Bond


class Molecule:

    def __init__(self):
        self.atoms: list[Atom] = []
        self.bonds: list[Bond] = []

    def add_atom(self, atom: Atom) -> int:
        """Add an atom and return its index."""
        index = len(self.atoms)
        self.atoms.append(atom)
        return index

    def add_bond(self, bond: Bond) -> None:
        """Add a bond to the molecule."""
        self.bonds.append(bond)

    def get_atom(self, index: int) -> Atom:
        """Return the atom at the given index."""
        return self.atoms[index]

    def get_bonds(self, atom_index: int) -> list[Bond]:
        """Return all bonds connected to an atom."""
        return [
            bond
            for bond in self.bonds
            if bond.atom_a == atom_index
            or bond.atom_b == atom_index
        ]

    def neighbors(self, atom_index: int) -> list[int]:
        """Return the indices of atoms directly bonded to an atom."""
        return [
            bond.other_atom(atom_index)
            for bond in self.get_bonds(atom_index)
        ]

    def formula(self) -> dict[str, int]:
        """Return element counts for the molecule."""
        result: dict[str, int] = {}

        for atom in self.atoms:
            result[atom.element] = result.get(atom.element, 0) + 1

        return result

    def __repr__(self) -> str:
        return f"Molecule(atoms={self.atoms!r}, bonds={self.bonds!r})"

    def __str__(self) -> str:
        atoms = ", ".join(str(atom) for atom in self.atoms)
        return f"Molecule({atoms})"