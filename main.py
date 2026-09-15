from chemistry.atoms import Atom
from chemistry.bonds import Bond, BondType
from chemistry.molecule import Molecule


def main():
    molecule = Molecule()

    carbon_1 = molecule.add_atom(Atom("C"))
    carbon_2 = molecule.add_atom(Atom("C"))

    molecule.add_bond(
        Bond(
            atom_a=carbon_1,
            atom_b=carbon_2,
            order=1.0,
            bond_type=BondType.SINGLE,
        )
    )

    print(molecule)
    print("Formula:", molecule.formula())
    print("Neighbors of carbon 1:", molecule.neighbors(0))


if __name__ == "__main__":
    main()