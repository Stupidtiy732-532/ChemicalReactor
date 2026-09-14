"""
Structure parser and molecular analysis.
"""

import re
from collections import Counter


ELEMENTS = {
    "H", "He",
    "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
    "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se",
    "Br", "Kr", "Rb", "Sr", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "Pt", "Au",
    "Hg", "Pb", "Bi"
}


ATOMIC_MASSES = {
    "H": 1.008,
    "B": 10.81,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.998,
    "Na": 22.990,
    "Mg": 24.305,
    "Al": 26.982,
    "Si": 28.085,
    "P": 30.974,
    "S": 32.06,
    "Cl": 35.45,
    "K": 39.098,
    "Ca": 40.078,
    "Cr": 51.996,
    "Mn": 54.938,
    "Fe": 55.845,
    "Co": 58.933,
    "Ni": 58.693,
    "Cu": 63.546,
    "Zn": 65.38,
    "Br": 79.904,
    "Ag": 107.868,
    "Sn": 118.710,
    "I": 126.904,
    "Pt": 195.084,
    "Au": 196.967,
    "Hg": 200.592,
    "Pb": 207.2,
}


# Longest tokens must appear first.
TOKEN_PATTERN = re.compile(
    r"Cl|Br|Na|Mg|Al|Si|Ca|Cr|Mn|Fe|Co|Ni|Cu|Zn|"
    r"Ag|Sn|Pt|Au|Hg|Pb|"
    r"CH3|CH2|CH|NH2|NH|COOH|CHO|CO|OH|"
    r"[A-Z][a-z]?|"
    r"\d+|"
    r"\(|\)|=|#|-"
)


class StructureParser:

    def __init__(self, notation):
        self.notation = notation.replace(" ", "")
        self.tokens = TOKEN_PATTERN.findall(self.notation)
        self.position = 0

    def parse(self):
        if not self.tokens:
            raise ValueError("Empty molecular notation")

        atoms = []
        bonds = []

        current_atom = None
        pending_bond = 1

        while self.position < len(self.tokens):
            token = self.tokens[self.position]

            if token in {"=", "#", "-"}:
                pending_bond = {
                    "-": 1,
                    "=": 2,
                    "#": 3
                }[token]

                self.position += 1
                continue

            if token == "(":
                self.position += 1
                continue

            if token == ")":
                self.position += 1
                continue

            if token.isdigit():
                self.position += 1
                continue

            element, hydrogens = self.expand_token(token)

            atom_id = len(atoms)
            atoms.append({
                "id": atom_id,
                "element": element,
                "hydrogens": hydrogens
            })

            if current_atom is not None:
                bonds.append({
                    "atom1": current_atom,
                    "atom2": atom_id,
                    "order": pending_bond
                })

            current_atom = atom_id
            pending_bond = 1
            self.position += 1

        return atoms, bonds

    def expand_token(self, token):
        special = {
            "CH3": ("C", 3),
            "CH2": ("C", 2),
            "CH": ("C", 1),
            "C": ("C", 0),
            "NH2": ("N", 2),
            "NH": ("N", 1),
            "N": ("N", 0),
            "OH": ("O", 1),
            "O": ("O", 0),
            "COOH": ("C", 0),
            "CHO": ("C", 1),
            "CO": ("C", 0),
        }

        if token in special:
            return special[token]

        if token in ELEMENTS:
            return token, 0

        raise ValueError(f"Unsupported notation token: {token}")


class StructureAnalyzer:

    def __init__(self, atoms, bonds):
        self.atoms = atoms
        self.bonds = bonds

    def formula_counts(self):
        counts = Counter()

        for atom in self.atoms:
            counts[atom["element"]] += 1
            counts["H"] += atom.get("hydrogens", 0)

        return counts

    def formula(self):
        counts = self.formula_counts()

        order = []

        if "C" in counts:
            order.append("C")

        if "H" in counts:
            order.append("H")

        order.extend(
            element for element in sorted(counts)
            if element not in {"C", "H"}
        )

        result = []

        for element in order:
            amount = counts[element]

            if amount == 1:
                result.append(element)
            elif amount > 1:
                result.append(f"{element}{amount}")

        return "".join(result)

    def molar_mass(self):
        counts = self.formula_counts()
        mass = 0.0

        for element, amount in counts.items():
            mass += ATOMIC_MASSES.get(element, 0.0) * amount

        return mass

    def dbe(self):
        """
        Approximate degree of unsaturation:

        DBE = (2C + 2 + N - H - X) / 2

        X = F + Cl + Br + I

        This is most reliable for ordinary neutral organic molecules.
        """

        counts = self.formula_counts()

        carbon = counts.get("C", 0)
        nitrogen = counts.get("N", 0)
        hydrogen = counts.get("H", 0)

        halogens = sum(
            counts.get(element, 0)
            for element in ["F", "Cl", "Br", "I"]
        )

        return (2 * carbon + 2 + nitrogen - hydrogen - halogens) / 2

    def bond_summary(self):
        single = 0
        double = 0
        triple = 0

        for bond in self.bonds:
            if bond["order"] == 1:
                single += 1
            elif bond["order"] == 2:
                double += 1
            elif bond["order"] == 3:
                triple += 1

        return {
            "single": single,
            "double": double,
            "triple": triple
        }

    def functional_groups(self):
        groups = []

        elements = [atom["element"] for atom in self.atoms]

        if "OH" in self._raw_tokens():
            groups.append("Alcohol")

        if "COOH" in self._raw_tokens():
            groups.append("Carboxylic acid")

        if "CHO" in self._raw_tokens():
            groups.append("Aldehyde")

        if "NH2" in self._raw_tokens():
            groups.append("Primary amine")

        if any(element in elements for element in ["Cl", "Br", "I", "F"]):
            groups.append("Halo compound")

        if any(bond["order"] == 3 for bond in self.bonds):
            groups.append("Alkyne")

        if any(bond["order"] == 2 for bond in self.bonds):
            groups.append("Alkene or carbonyl")

        if not groups:
            groups.append("No recognised functional group")

        return sorted(set(groups))

    def _raw_tokens(self):
        return [
            atom.get("token", atom["element"])
            for atom in self.atoms
        ]

    def report(self):
        print("\nSTRUCTURE ANALYSIS")
        print("=" * 60)

        print(f"Formula       : {self.formula()}")
        print(f"Molar mass    : {self.molar_mass():.3f} g mol^-1")
        print(f"DBE           : {self.dbe():.2f}")

        print("\nFunctional groups:")
        for group in self.functional_groups():
            print(f"  - {group}")

        print("\nBond summary:")
        for name, amount in self.bond_summary().items():
            print(f"  {name.capitalize()} bonds: {amount}")

        print("\nAtoms:")
        for atom in self.atoms:
            print(
                f"  Atom {atom['id']}: "
                f"{atom['element']} "
                f"(implicit H = {atom.get('hydrogens', 0)})"
            )

        print("\nBonds:")
        for bond in self.bonds:
            print(
                f"  {bond['atom1']} "
                f"--{bond['order']}-- "
                f"{bond['atom2']}"
            )