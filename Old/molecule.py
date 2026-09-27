from dataclasses import dataclass, field
import json


@dataclass
class Vector3D:
    x: float = 0.0
    y: float = 0.0
    z: float = 0.0

    def as_tuple(self):
        return self.x, self.y, self.z

    def distance_to(self, other):
        dx = self.x - other.x
        dy = self.y - other.y
        dz = self.z - other.z

        return (dx * dx + dy * dy + dz * dz) ** 0.5

    def to_dict(self):
        return {
            "x": self.x,
            "y": self.y,
            "z": self.z,
        }


@dataclass
class Atom:
    id: int
    element: str
    hydrogens: int = 0
    charge: int = 0
    radical_electrons: int = 0
    lone_pairs: int = 0
    mass_number: int | None = None
    mark: str | None = None
    position: Vector3D = field(default_factory=Vector3D)

    def isotope_label(self):
        if self.mass_number is None:
            return self.element

        return f"{self.element}_{self.mass_number}"

    def display_label(self):
        label = self.isotope_label()

        if self.hydrogens:
            label += f"H{self.hydrogens}"

        if self.charge:
            sign = "+" if self.charge > 0 else "-"
            magnitude = abs(self.charge)
            label += f"^{'' if magnitude == 1 else magnitude}{sign}"

        if self.radical_electrons:
            label += "." * self.radical_electrons

        if self.lone_pairs:
            label += ":" * self.lone_pairs

        if self.mark is not None:
            label += f"X<{self.mark}>"

        return label

    def to_dict(self):
        return {
            "id": self.id,
            "element": self.element,
            "hydrogens": self.hydrogens,
            "charge": self.charge,
            "radical_electrons": self.radical_electrons,
            "lone_pairs": self.lone_pairs,
            "mass_number": self.mass_number,
            "mark": self.mark,
            "position": self.position.to_dict(),
        }

    def __repr__(self):
        return (
            f"Atom(id={self.id}, element={self.element}, "
            f"H={self.hydrogens}, charge={self.charge}, "
            f"radical={self.radical_electrons}, "
            f"lone_pairs={self.lone_pairs}, "
            f"isotope={self.mass_number}, mark={self.mark})"
        )


@dataclass
class Bond:
    atom1: int
    atom2: int
    order: int = 1

    @property
    def sigma_bonds(self):
        return 1

    @property
    def pi_bonds(self):
        return max(0, self.order - 1)

    def contains(self, atom_id):
        return atom_id in (self.atom1, self.atom2)

    def other_atom(self, atom_id):
        if atom_id == self.atom1:
            return self.atom2

        if atom_id == self.atom2:
            return self.atom1

        raise ValueError(f"Atom {atom_id} is not part of this bond.")

    def to_dict(self):
        return {
            "atom1": self.atom1,
            "atom2": self.atom2,
            "order": self.order,
        }

    def __repr__(self):
        return f"Bond({self.atom1}-{self.atom2}, order={self.order})"


class Molecule:
    ELEMENTS = {
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
        "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
        "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se",
        "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo",
        "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
        "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce",
        "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
        "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
        "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb",
        "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
        "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
        "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg",
        "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl",
        "Mc", "Lv", "Ts", "Og",
    }

    ATOMIC_MASSES = {
        "H": 1.008,
        "C": 12.011,
        "N": 14.007,
        "O": 15.999,
        "F": 18.998,
        "P": 30.974,
        "S": 32.06,
        "Cl": 35.45,
        "Br": 79.904,
        "I": 126.904,
        "Na": 22.990,
        "Mg": 24.305,
        "Al": 26.982,
        "Si": 28.085,
        "K": 39.098,
        "Ca": 40.078,
        "Fe": 55.845,
        "Cu": 63.546,
        "Zn": 65.38,
        "Ag": 107.868,
        "Au": 196.967,
    }

    def __init__(self):
        self.atoms = {}
        self.bonds = []

    def add_atom(
        self,
        element,
        hydrogens=0,
        charge=0,
        radical_electrons=0,
        lone_pairs=0,
        mass_number=None,
        mark=None,
        position=None,
    ):
        if element not in self.ELEMENTS:
            raise ValueError(f"Unknown element: {element}")

        if hydrogens < 0:
            raise ValueError("Hydrogen count cannot be negative.")

        if mass_number is not None and mass_number <= 0:
            raise ValueError("Mass number must be positive.")

        atom_id = len(self.atoms) + 1

        atom = Atom(
            id=atom_id,
            element=element,
            hydrogens=hydrogens,
            charge=charge,
            radical_electrons=radical_electrons,
            lone_pairs=lone_pairs,
            mass_number=mass_number,
            mark=mark,
            position=position or Vector3D(),
        )

        self.atoms[atom_id] = atom
        return atom_id

    def add_bond(self, atom1, atom2, order=1):
        if atom1 not in self.atoms or atom2 not in self.atoms:
            raise ValueError("Bond refers to an unknown atom.")

        if atom1 == atom2:
            raise ValueError("An atom cannot bond to itself.")

        if order not in (1, 2, 3):
            raise ValueError("Bond order must be 1, 2, or 3.")

        if self.bond_between(atom1, atom2) is not None:
            raise ValueError(
                f"Bond already exists between {atom1} and {atom2}."
            )

        self.bonds.append(Bond(atom1, atom2, order))

    def get_atom(self, atom_id):
        return self.atoms[atom_id]

    def bonds_of(self, atom_id):
        return [
            bond
            for bond in self.bonds
            if bond.contains(atom_id)
        ]

    def neighbors(self, atom_id):
        return [
            bond.other_atom(atom_id)
            for bond in self.bonds_of(atom_id)
        ]

    def bond_between(self, atom1, atom2):
        for bond in self.bonds:
            if (
                    (bond.atom1 == atom1 and bond.atom2 == atom2)
                    or
                    (bond.atom1 == atom2 and bond.atom2 == atom1)
            ):
                return bond

        return None

    def bond_order_sum(self, atom_id):
        return sum(
            bond.order
            for bond in self.bonds_of(atom_id)
        )

    def total_charge(self):
        return sum(atom.charge for atom in self.atoms.values())

    def molecular_mass(self):
        hydrogen_mass = self.ATOMIC_MASSES["H"]
        mass = 0.0

        for atom in self.atoms.values():
            if atom.mass_number is not None:
                mass += atom.mass_number
            else:
                mass += self.ATOMIC_MASSES.get(atom.element, 0.0)

            mass += atom.hydrogens * hydrogen_mass

        return mass

    def formula(self):
        counts = {}

        for atom in self.atoms.values():
            element = (
                atom.element
                if atom.mass_number is None
                else f"[{atom.mass_number}{atom.element}]"
            )

            counts[element] = counts.get(element, 0) + 1

            if atom.hydrogens:
                counts["H"] = counts.get("H", 0) + atom.hydrogens

        def hill_key(element):
            if element == "C" or (
                element.startswith("[")
                and element.endswith("C]")
            ):
                return 0, element

            if element == "H":
                return 1, element

            return 2, element

        return "".join(
            element + (str(count) if count != 1 else "")
            for element, count in sorted(
                counts.items(),
                key=lambda item: hill_key(item[0]),
            )
        )

    def to_dict(self):
        return {
            "formula": self.formula(),
            "total_charge": self.total_charge(),
            "molecular_mass_u": self.molecular_mass(),
            "atoms": [
                atom.to_dict()
                for atom in self.atoms.values()
            ],
            "bonds": [
                bond.to_dict()
                for bond in self.bonds
            ],
        }

    def to_json(self, indent=2):
        return json.dumps(self.to_dict(), indent=indent)

    def __repr__(self):
        return (
            f"Molecule(formula={self.formula()}, "
            f"atoms={len(self.atoms)}, "
            f"bonds={len(self.bonds)}, "
            f"charge={self.total_charge()})"
        )