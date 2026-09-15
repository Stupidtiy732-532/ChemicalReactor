"""
OrgReact — Molecule representation

Atoms are vertices.
Bonds are edges.
"""

from chemistry.atoms import Atom
from chemistry.bonds import Bond


class Molecule:
    def __init__(self):
        self.atoms: list[Atom] = []
        self.bonds: list[Bond] = []

    def add_atom(self, atom: Atom) -> int:
        index = len(self.atoms)
        self.atoms.append(atom)
        return index

    def add_bond(self, bond: Bond):
        self.bonds.append(bond)

    def get_atom(self, index: int) -> Atom:
        return self.atoms[index]

    def get_bonds(self, atom_index: int) -> list[Bond]:
        return [
            bond
            for bond in self.bonds
            if bond.atom_a == atom_index
            or bond.atom_b == atom_index
        ]

    def neighbors(self, atom_index: int) -> list[int]:
        return [
            bond.other_atom(atom_index)
            for bond in self.get_bonds(atom_index)
        ]

    def formula(self) -> dict[str, int]:
        result = {}

        for atom in self.atoms:
            result[atom.element] = result.get(atom.element, 0) + 1

        return result

    def __repr__(self):
        return (
            f"Molecule("
            f"atoms={self.atoms!r}, "
            f"bonds={self.bonds!r}"
            f")"
        )

    def __str__(self):
        return "Molecule(" + ", ".join(
            str(atom) for atom in self.atoms
        ) + ")"