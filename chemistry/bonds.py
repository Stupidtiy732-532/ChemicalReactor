"""
OrgReact — Bond representation

No SMILES.
"""

from enum import Enum
from typing import Any


class BondType(Enum):
    SINGLE = "-"
    DOUBLE = "="
    TRIPLE = "#"
    PARTIAL = "partial"
    COORDINATE = "coordinate"


class Bond:
    """
    Represents a bond between two atoms.

    order:
        Mathematical bond order.
        1   = single
        2   = double
        3   = triple
        1.5 = partial

    bond_type:
        SINGLE, DOUBLE, TRIPLE, PARTIAL, COORDINATE

    expression:
        Optional mathematical expression for bond order.
    """

    def __init__(
        self,
        atom_a: int,
        atom_b: int,
        order: float = 1.0,
        bond_type: BondType = BondType.SINGLE,
        expression: Any = None,
        coordinate: bool = False,
    ):
        self.atom_a = atom_a
        self.atom_b = atom_b
        self.order = order
        self.bond_type = bond_type
        self.expression = expression
        self.coordinate = coordinate

    def __repr__(self):
        return (
            f"Bond("
            f"atom_a={self.atom_a}, "
            f"atom_b={self.atom_b}, "
            f"order={self.order}, "
            f"bond_type={self.bond_type}, "
            f"coordinate={self.coordinate}"
            f")"
        )

    def copy(self):
        return Bond(
            atom_a=self.atom_a,
            atom_b=self.atom_b,
            order=self.order,
            bond_type=self.bond_type,
            expression=self.expression,
            coordinate=self.coordinate,
        )

    def other_atom(self, atom_index: int) -> int:
        if atom_index == self.atom_a:
            return self.atom_b

        if atom_index == self.atom_b:
            return self.atom_a

        raise ValueError("Atom is not part of this bond")

    def __str__(self):
        if self.coordinate:
            return "coordinate"

        if self.bond_type == BondType.PARTIAL:
            return f"-_{{{self.expression}}}"

        return self.bond_type.value