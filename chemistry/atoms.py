"""Examples:
    C
    :C
    C.
    [14C]
    [2H]
    C^{+}
    C^{2}+
"""


class Atom:
    """
    Represents one atom in a chemical structure.

    Parameters
    ----------
    element : str
        Chemical element symbol, e.g. C, H, O, N.

    isotope : int | None
        Mass number, e.g. 14 for carbon-14.

    charge : int
        Net charge of this atom.

    lone_pairs : int | None
        Number of lone pairs.

    unpaired_electrons : int
        Number of unpaired electrons.

    label : str | None
        Optional identifier for the atom.
    """

    def __init__(
        self,
        element: str,
        isotope: int | None = None,
        charge: int = 0,
        lone_pairs: int | None = None,
        unpaired_electrons: int = 0,
        label: str | None = None,
    ):
        self.element = element
        self.isotope = isotope
        self.charge = charge
        self.lone_pairs = lone_pairs
        self.unpaired_electrons = unpaired_electrons
        self.label = label

    def __repr__(self):
        return (
            f"Atom("
            f"element={self.element!r}, "
            f"isotope={self.isotope!r}, "
            f"charge={self.charge!r}, "
            f"lone_pairs={self.lone_pairs!r}, "
            f"unpaired_electrons={self.unpaired_electrons!r}"
            f")"
        )

    def copy(self):
        return Atom(
            element=self.element,
            isotope=self.isotope,
            charge=self.charge,
            lone_pairs=self.lone_pairs,
            unpaired_electrons=self.unpaired_electrons,
            label=self.label,
        )

    def formal_charge(self):
        return self.charge

    def has_radical(self):
        return self.unpaired_electrons > 0

    def has_lone_pairs(self):
        return self.lone_pairs is not None and self.lone_pairs > 0

    def __str__(self):
        result = ""

        if self.isotope is not None:
            result += f"[{self.isotope}"

        result += self.element

        if self.isotope is not None:
            result += "]"

        if self.lone_pairs is not None:
            result += ":" * self.lone_pairs

        if self.unpaired_electrons > 0:
            result += "." * self.unpaired_electrons

        if self.charge != 0:
            sign = "+" if self.charge > 0 else "-"
            magnitude = abs(self.charge)

            if magnitude == 1:
                result += f"^{sign}"
            else:
                result += f"^{magnitude}{sign}"

        return result