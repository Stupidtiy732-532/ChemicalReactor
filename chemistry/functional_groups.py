from dataclasses import dataclass

from chemistry.molecule import Molecule
from chemistry.bonds import BondType


@dataclass(frozen=True)
class FunctionalGroup:
    name: str
    atoms: tuple
    description: str


def bond_between(molecule: Molecule, atom_a, atom_b):

    for bond in molecule.bonds:

        if (
            bond.atom_a is atom_a and bond.atom_b is atom_b
        ) or (
            bond.atom_a is atom_b and bond.atom_b is atom_a
        ):
            return bond

    return None


def neighbors_of(molecule: Molecule, atom):

    return molecule.neighbors(atom)


def find_functional_groups(molecule: Molecule) -> list[FunctionalGroup]:

    groups = []

    for atom in molecule.atoms:

        # -------------------------------------------------
        # ALCOHOL: C-O-H
        # -------------------------------------------------

        if atom.element == "O":

            oxygen_neighbors = neighbors_of(molecule, atom)

            if len(oxygen_neighbors) == 1:

                carbon = oxygen_neighbors[0]

                if carbon.element == "C":

                    groups.append(
                        FunctionalGroup(
                            name="alcohol",
                            atoms=(carbon, atom),
                            description="Carbon bonded to hydroxyl oxygen",
                        )
                    )

            # -------------------------------------------------
            # ETHER: C-O-C
            # -------------------------------------------------

            elif len(oxygen_neighbors) == 2:

                first, second = oxygen_neighbors

                if (
                    first.element == "C"
                    and second.element == "C"
                ):

                    groups.append(
                        FunctionalGroup(
                            name="ether",
                            atoms=(first, atom, second),
                            description="Oxygen bonded to two carbon atoms",
                        )
                    )

        # -------------------------------------------------
        # CARBON-BASED GROUPS
        # -------------------------------------------------

        if atom.element != "C":
            continue

        carbon_neighbors = neighbors_of(molecule, atom)

        double_bonded_oxygen = None
        single_bonded_oxygen = None
        single_bonded_nitrogen = None
        single_bonded_halogen = None

        for neighbor in carbon_neighbors:

            bond = bond_between(molecule, atom, neighbor)

            if bond is None:
                continue

            if neighbor.element == "O":

                if bond.order == 2:
                    double_bonded_oxygen = neighbor

                elif bond.order == 1:
                    single_bonded_oxygen = neighbor

            elif neighbor.element == "N" and bond.order == 1:

                single_bonded_nitrogen = neighbor

            elif neighbor.element in {"F", "Cl", "Br", "I"}:

                if bond.order == 1:
                    single_bonded_halogen = neighbor

        # -------------------------------------------------
        # CARBONYL: C=O
        # -------------------------------------------------

        if double_bonded_oxygen is not None:

            groups.append(
                FunctionalGroup(
                    name="carbonyl",
                    atoms=(atom, double_bonded_oxygen),
                    description="Carbon double-bonded to oxygen",
                )
            )

            # -------------------------------------------------
            # CARBOXYLIC ACID: C(=O)-O
            # -------------------------------------------------

            if single_bonded_oxygen is not None:

                groups.append(
                    FunctionalGroup(
                        name="carboxylic_acid",
                        atoms=(
                            atom,
                            double_bonded_oxygen,
                            single_bonded_oxygen,
                        ),
                        description="Carbon bonded to C=O and C-O",
                    )
                )

            # -------------------------------------------------
            # AMIDE: C(=O)-N
            # -------------------------------------------------

            if single_bonded_nitrogen is not None:

                groups.append(
                    FunctionalGroup(
                        name="amide",
                        atoms=(
                            atom,
                            double_bonded_oxygen,
                            single_bonded_nitrogen,
                        ),
                        description="Carbonyl carbon bonded to nitrogen",
                    )
                )

        # -------------------------------------------------
        # HALO COMPOUND: C-X
        # -------------------------------------------------

        if single_bonded_halogen is not None:

            groups.append(
                FunctionalGroup(
                    name="halo_compound",
                    atoms=(atom, single_bonded_halogen),
                    description="Carbon bonded to halogen",
                )
            )

    return groups