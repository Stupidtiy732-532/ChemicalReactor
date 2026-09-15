"""
OrgReact — Molecular Classification

Classifies:
- Primary / secondary / tertiary / quaternary carbon
- Primary / secondary / tertiary alcohol
- Aldehydes / ketones
- Carboxylic acids
- Ethers
- Haloalkanes
- Alkenes / alkynes
"""

from dataclasses import dataclass, field
from collections import defaultdict


# ============================================================
# DATA STRUCTURES
# ============================================================

@dataclass
class FunctionalGroup:
    name: str
    atom_ids: tuple
    details: str = ""


@dataclass
class MoleculeClassification:
    carbon_types: dict = field(default_factory=dict)
    alcohol_sites: list = field(default_factory=list)
    functional_groups: list = field(default_factory=list)

    def all_labels(self):
        labels = []

        labels.extend(self.carbon_types.values())

        for alcohol in self.alcohol_sites:
            labels.append(alcohol["classification"])

        for group in self.functional_groups:
            labels.append(group.name)

        return labels


# ============================================================
# HELPERS
# ============================================================

def element(atom):
    return str(atom.element).capitalize()


def atom_by_id(molecule, atom_id):
    return molecule.atoms[atom_id]


def bond_other_atom(bond, atom_id):
    if bond.atom1 == atom_id:
        return bond.atom2

    if bond.atom2 == atom_id:
        return bond.atom1

    return None


def connected_atoms(molecule, atom_id):
    """
    Returns:
        list[(neighbor_id, bond_order)]
    """

    result = []

    for bond in molecule.bonds:
        other = bond_other_atom(bond, atom_id)

        if other is not None:
            result.append((other, bond.order))

    return result


def carbon_neighbors(molecule, atom_id):
    """
    Returns directly bonded carbon atoms.

    Important:
    Carbon classification depends on the number of carbon
    neighbors, not the total bond order.
    """

    result = []

    for neighbor_id, bond_order in connected_atoms(molecule, atom_id):
        neighbor = atom_by_id(molecule, neighbor_id)

        if element(neighbor) == "C":
            result.append(neighbor_id)

    return result


def carbon_class(molecule, carbon_id):
    """
    Classifies one carbon atom.

    1 carbon neighbor -> primary carbon
    2 carbon neighbors -> secondary carbon
    3 carbon neighbors -> tertiary carbon
    4 carbon neighbors -> quaternary carbon
    """

    count = len(carbon_neighbors(molecule, carbon_id))

    classes = {
        0: "methyl/carbon with no carbon neighbors",
        1: "primary carbon",
        2: "secondary carbon",
        3: "tertiary carbon",
        4: "quaternary carbon",
    }

    return classes.get(count, f"hyper-substituted carbon ({count} carbon neighbors)")


def atom_has_hydrogen(atom):
    return getattr(atom, "hydrogens", 0) > 0


# ============================================================
# ALCOHOL CLASSIFICATION
# ============================================================

def classify_alcohol_site(molecule, oxygen_id):
    """
    Detects and classifies an alcohol oxygen.

    Alcohol:
        O-H
        O attached to carbon

    The alcohol type depends on the carbon attached to OH.
    """

    oxygen = atom_by_id(molecule, oxygen_id)

    if element(oxygen) != "O":
        return None

    if not atom_has_hydrogen(oxygen):
        return None

    neighbors = connected_atoms(molecule, oxygen_id)

    carbon_neighbors_of_oxygen = []

    for neighbor_id, bond_order in neighbors:
        neighbor = atom_by_id(molecule, neighbor_id)

        if element(neighbor) == "C":
            carbon_neighbors_of_oxygen.append(neighbor_id)

    # Ordinary alcohol must have one carbon attached to oxygen.
    if len(carbon_neighbors_of_oxygen) != 1:
        return None

    carbon_id = carbon_neighbors_of_oxygen[0]
    carbon_type = carbon_class(molecule, carbon_id)

    alcohol_type_map = {
        "primary carbon": "primary alcohol",
        "secondary carbon": "secondary alcohol",
        "tertiary carbon": "tertiary alcohol",
    }

    classification = alcohol_type_map.get(
        carbon_type,
        f"alcohol attached to {carbon_type}"
    )

    return {
        "oxygen_id": oxygen_id,
        "carbon_id": carbon_id,
        "classification": classification,
    }


# ============================================================
# CARBONYL CLASSIFICATION
# ============================================================

def carbonyl_oxygen(molecule, carbon_id):
    """
    Returns oxygen atoms double-bonded to a carbon atom.
    """

    result = []

    for neighbor_id, bond_order in connected_atoms(molecule, carbon_id):
        neighbor = atom_by_id(molecule, neighbor_id)

        if element(neighbor) == "O" and bond_order == 2:
            result.append(neighbor_id)

    return result


def classify_carbonyl(molecule, carbon_id):
    """
    Distinguishes:
    - Aldehyde
    - Ketone
    - Generic carbonyl
    """

    carbon = atom_by_id(molecule, carbon_id)

    if element(carbon) != "C":
        return None

    double_bonded_oxygen = carbonyl_oxygen(molecule, carbon_id)

    if not double_bonded_oxygen:
        return None

    carbon_neighbors = carbon_neighbors_of_carbonyl(molecule, carbon_id)

    # Carbonyl carbon attached to two carbons = ketone.
    if len(carbon_neighbors) == 2:
        return "ketone"

    # Carbonyl carbon attached to one carbon and hydrogen = aldehyde.
    if len(carbon_neighbors) == 1 and atom_has_hydrogen(carbon):
        return "aldehyde"

    # Formaldehyde: H-C(=O)-H
    if len(carbon_neighbors) == 0 and atom_has_hydrogen(carbon):
        return "aldehyde"

    return "carbonyl compound"


def carbon_neighbors_of_carbonyl(molecule, carbon_id):
    return carbon_neighbors(molecule, carbon_id)


# ============================================================
# OTHER FUNCTIONAL GROUPS
# ============================================================

def detect_ether(molecule, oxygen_id):
    """
    Ether:
        C-O-C

    Oxygen must:
        - Have no hydrogen
        - Be attached to two carbons
    """

    oxygen = atom_by_id(molecule, oxygen_id)

    if element(oxygen) != "O":
        return False

    if atom_has_hydrogen(oxygen):
        return False

    neighbors = connected_atoms(molecule, oxygen_id)

    carbon_count = 0

    for neighbor_id, bond_order in neighbors:
        neighbor = atom_by_id(molecule, neighbor_id)

        if element(neighbor) == "C" and bond_order == 1:
            carbon_count += 1

    return carbon_count == 2


def detect_carboxylic_acid(molecule, carbon_id):
    """
    Carboxylic acid:
        C(=O)-O-H
    """

    carbon = atom_by_id(molecule, carbon_id)

    if element(carbon) != "C":
        return False

    has_carbonyl = False
    has_hydroxyl_oxygen = False

    for neighbor_id, bond_order in connected_atoms(molecule, carbon_id):
        neighbor = atom_by_id(molecule, neighbor_id)

        if element(neighbor) != "O":
            continue

        if bond_order == 2:
            has_carbonyl = True

        elif bond_order == 1 and atom_has_hydrogen(neighbor):
            has_hydroxyl_oxygen = True

    return has_carbonyl and has_hydroxyl_oxygen


def detect_haloalkane(molecule, carbon_id):
    """
    Detects carbon-halogen bonds:
        C-F
        C-Cl
        C-Br
        C-I
    """

    carbon = atom_by_id(molecule, carbon_id)

    if element(carbon) != "C":
        return False

    halogens = {"F", "Cl", "Br", "I"}

    for neighbor_id, bond_order in connected_atoms(molecule, carbon_id):
        neighbor = atom_by_id(molecule, neighbor_id)

        if element(neighbor) in halogens and bond_order == 1:
            return True

    return False


def detect_unsaturation(molecule):
    """
    Detects C=C and C≡C bonds.
    """

    groups = []

    for bond in molecule.bonds:
        atom1 = atom_by_id(molecule, bond.atom1)
        atom2 = atom_by_id(molecule, bond.atom2)

        if element(atom1) != "C" or element(atom2) != "C":
            continue

        if bond.order == 2:
            groups.append(
                FunctionalGroup(
                    name="alkene",
                    atom_ids=(bond.atom1, bond.atom2),
                    details="carbon-carbon double bond"
                )
            )

        elif bond.order == 3:
            groups.append(
                FunctionalGroup(
                    name="alkyne",
                    atom_ids=(bond.atom1, bond.atom2),
                    details="carbon-carbon triple bond"
                )
            )

    return groups


# ============================================================
# MAIN CLASSIFIER
# ============================================================

def classify_molecule(molecule):
    """
    Analyze the entire molecule graph.
    """

    result = MoleculeClassification()

    # --------------------------------------------------------
    # Carbon classification
    # --------------------------------------------------------

    for atom_id, atom in molecule.atoms.items():
        if element(atom) == "C":
            result.carbon_types[atom_id] = carbon_class(
                molecule,
                atom_id
            )

    # --------------------------------------------------------
    # Oxygen-based groups
    # --------------------------------------------------------

    for atom_id, atom in molecule.atoms.items():
        atom_element = element(atom)

        if atom_element == "O":
            alcohol = classify_alcohol_site(molecule, atom_id)

            if alcohol is not None:
                result.alcohol_sites.append(alcohol)

            if detect_ether(molecule, atom_id):
                result.functional_groups.append(
                    FunctionalGroup(
                        name="ether",
                        atom_ids=(atom_id,),
                        details="C-O-C"
                    )
                )

    # --------------------------------------------------------
    # Carbon-based groups
    # --------------------------------------------------------

    for atom_id, atom in molecule.atoms.items():
        if element(atom) != "C":
            continue

        carbonyl_type = classify_carbonyl(molecule, atom_id)

        if carbonyl_type is not None:
            result.functional_groups.append(
                FunctionalGroup(
                    name=carbonyl_type,
                    atom_ids=(atom_id,),
                    details="carbon-oxygen double bond"
                )
            )

        if detect_carboxylic_acid(molecule, atom_id):
            result.functional_groups.append(
                FunctionalGroup(
                    name="carboxylic acid",
                    atom_ids=(atom_id,),
                    details="C(=O)-O-H"
                )
            )

        if detect_haloalkane(molecule, atom_id):
            result.functional_groups.append(
                FunctionalGroup(
                    name="haloalkane",
                    atom_ids=(atom_id,),
                    details="carbon-halogen bond"
                )
            )

    # --------------------------------------------------------
    # C=C and C≡C
    # --------------------------------------------------------

    result.functional_groups.extend(
        detect_unsaturation(molecule)
    )

    return result


# ============================================================
# OUTPUT
# ============================================================

def classification_summary(classification):
    lines = []

    lines.append("MOLECULAR CLASSIFICATION")
    lines.append("=" * 40)

    if classification.carbon_types:
        lines.append("")
        lines.append("CARBON TYPES")

        for atom_id, carbon_type in classification.carbon_types.items():
            lines.append(
                f"  Atom {atom_id}: {carbon_type}"
            )

    if classification.alcohol_sites:
        lines.append("")
        lines.append("ALCOHOL SITES")

        for alcohol in classification.alcohol_sites:
            lines.append(
                f"  O atom {alcohol['oxygen_id']} "
                f"attached to C atom {alcohol['carbon_id']}: "
                f"{alcohol['classification']}"
            )

    if classification.functional_groups:
        lines.append("")
        lines.append("FUNCTIONAL GROUPS")

        for group in classification.functional_groups:
            atom_text = ", ".join(map(str, group.atom_ids))

            lines.append(
                f"  {group.name} "
                f"(atoms: {atom_text})"
            )

            if group.details:
                lines.append(
                    f"    {group.details}"
                )

    if not classification.all_labels():
        lines.append("")
        lines.append("No recognized functional groups.")

    return "\n".join(lines)