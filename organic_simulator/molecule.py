# molecule.py

ATOMIC_MASS = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.998,
    "Cl": 35.45,
    "Br": 79.904,
    "I": 126.904,
    "Na": 22.990,
    "K": 39.098,
    "S": 32.06,
    "P": 30.974,
}

VALENCE = {
    "H": 1,
    "C": 4,
    "N": 3,
    "O": 2,
    "F": 1,
    "Cl": 1,
    "Br": 1,
    "I": 1,
    "Na": 1,
    "K": 1,
    "S": 2,
    "P": 3,
}


class Atom:
    def __init__(self, atom_id, element, hydrogens=0):
        self.id = atom_id
        self.element = element
        self.hydrogens = hydrogens

    def copy(self):
        return Atom(self.id, self.element, self.hydrogens)


class Bond:
    def __init__(self, a, b, order=1):
        self.a = a
        self.b = b
        self.order = order

    def connects(self, x, y):
        return (self.a == x and self.b == y) or (
            self.a == y and self.b == x
        )


class Molecule:
    def __init__(self):
        self.atoms = {}
        self.bonds = []
        self.next_id = 0

    def add_atom(self, element, hydrogens=0):
        atom = Atom(self.next_id, element, hydrogens)
        self.atoms[self.next_id] = atom
        self.next_id += 1
        return atom.id

    def add_bond(self, a, b, order=1):
        self.bonds.append(Bond(a, b, order))

    def get_bond(self, a, b):
        for bond in self.bonds:
            if bond.connects(a, b):
                return bond
        return None

    def remove_bond(self, a, b):
        self.bonds = [
            bond for bond in self.bonds
            if not bond.connects(a, b)
        ]

    def neighbors(self, atom_id):
        result = []

        for bond in self.bonds:
            if bond.a == atom_id:
                result.append((bond.b, bond.order))
            elif bond.b == atom_id:
                result.append((bond.a, bond.order))

        return result

    def clone(self):
        new = Molecule()

        for atom in self.atoms.values():
            new.atoms[atom.id] = atom.copy()

        new.bonds = [
            Bond(bond.a, bond.b, bond.order)
            for bond in self.bonds
        ]

        new.next_id = self.next_id
        return new

    def formula(self):
        counts = {}

        for atom in self.atoms.values():
            counts[atom.element] = counts.get(atom.element, 0) + 1

            if atom.hydrogens:
                counts["H"] = counts.get("H", 0) + atom.hydrogens

        order = ["C", "H", "N", "O", "F", "Cl", "Br", "I", "S", "P", "Na", "K"]

        result = ""

        for element in order:
            if element in counts:
                number = counts[element]
                result += element
                if number != 1:
                    result += str(number)

        for element, number in counts.items():
            if element not in order:
                result += element
                if number != 1:
                    result += str(number)

        return result

    def molar_mass(self):
        mass = 0.0

        for atom in self.atoms.values():
            mass += ATOMIC_MASS.get(atom.element, 0.0)
            mass += atom.hydrogens * ATOMIC_MASS["H"]

        return mass

    def condensed(self):
        """
        Produces a useful structural representation.
        This is not a full IUPAC canonical SMILES generator.
        """

        carbon_atoms = [
            atom.id
            for atom in self.atoms.values()
            if atom.element == "C"
        ]

        if not carbon_atoms:
            return self.formula()

        start = None

        for carbon in carbon_atoms:
            carbon_neighbors = [
                self.atoms[n]
                for n, _ in self.neighbors(carbon)
                if self.atoms[n].element == "C"
            ]

            if len(carbon_neighbors) <= 1:
                start = carbon
                break

        if start is None:
            start = carbon_atoms[0]

        chain = []
        visited = set()
        current = start

        while current is not None and current not in visited:
            visited.add(current)
            atom = self.atoms[current]

            if atom.element != "C":
                break

            carbon_neighbors = [
                (n, order)
                for n, order in self.neighbors(current)
                if self.atoms[n].element == "C"
                and n not in visited
            ]

            oxygen_neighbors = [
                (self.atoms[n], order)
                for n, order in self.neighbors(current)
                if self.atoms[n].element == "O"
            ]

            nitrogen_neighbors = [
                (self.atoms[n], order)
                for n, order in self.neighbors(current)
                if self.atoms[n].element == "N"
            ]

            halogen_neighbors = [
                self.atoms[n].element
                for n, order in self.neighbors(current)
                if self.atoms[n].element in {"F", "Cl", "Br", "I"}
            ]

            carbonyl = any(order == 2 for _, order in oxygen_neighbors)

            if carbonyl:
                if any(
                    self.atoms[n].element == "O"
                    and order == 1
                    for n, order in self.neighbors(current)
                ):
                    group = "COOH"
                elif atom.hydrogens >= 1:
                    group = "CHO"
                else:
                    group = "CO"
            elif oxygen_neighbors:
                oxygen, order = oxygen_neighbors[0]

                if oxygen.hydrogens:
                    group = "CH" + str(atom.hydrogens) + "OH"
                else:
                    group = "CH" + str(atom.hydrogens) + "O"
            elif nitrogen_neighbors:
                group = "CH" + str(atom.hydrogens) + "NH2"
            elif halogen_neighbors:
                group = "CH" + str(atom.hydrogens) + halogen_neighbors[0]
            else:
                if atom.hydrogens == 3:
                    group = "CH3"
                elif atom.hydrogens == 2:
                    group = "CH2"
                elif atom.hydrogens == 1:
                    group = "CH"
                else:
                    group = "C"

            chain.append(group)

            if carbon_neighbors:
                current = carbon_neighbors[0][0]
            else:
                current = None

        if chain:
            return "-".join(chain)

        return self.formula()

    def describe(self):
        return {
            "formula": self.formula(),
            "molar_mass": round(self.molar_mass(), 3),
            "structure": self.condensed(),
        }