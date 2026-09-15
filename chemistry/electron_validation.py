from dataclasses import dataclass

from chemistry.molecule import Molecule


@dataclass(frozen=True)
class AtomAudit:
    atom_index: int
    element: str
    bond_order_sum: float
    lone_pairs: int
    unpaired_electrons: int
    charge: int
    valid: bool
    message: str


# Maximum ordinary bond-order occupancy.
# This is only a structural limit, not a complete quantum-chemical model.
MAX_VALENCE = {
    "H": 1,
    "D": 1,
    "T": 1,
    "C": 4,
    "N": 4,
    "O": 2,
    "F": 1,
    "Cl": 1,
    "Br": 1,
    "I": 1,
    "S": 6,
    "P": 5,
}


def audit_atom(molecule: Molecule, atom_index: int) -> AtomAudit:
    atom = molecule.atoms[atom_index]

    bond_order_sum = sum(
        float(bond.order)
        for bond in molecule.bonds
        if bond.atom_a == atom_index or bond.atom_b == atom_index
    )

    lone_pairs = atom.lone_pairs
    unpaired_electrons = atom.unpaired_electrons
    charge = atom.charge

    valid = True
    message = "OK"

    if lone_pairs < 0:
        valid = False
        message = "Lone-pair count cannot be negative."

    elif unpaired_electrons < 0:
        valid = False
        message = "Unpaired-electron count cannot be negative."

    elif atom.element not in MAX_VALENCE:
        valid = False
        message = f"Unknown valence rules for element {atom.element!r}."

    elif bond_order_sum > MAX_VALENCE[atom.element]:
        valid = False
        message = (
            f"Bond-order sum {bond_order_sum} exceeds "
            f"maximum valence {MAX_VALENCE[atom.element]}."
        )

    return AtomAudit(
        atom_index=atom_index,
        element=atom.element,
        bond_order_sum=bond_order_sum,
        lone_pairs=lone_pairs,
        unpaired_electrons=unpaired_electrons,
        charge=charge,
        valid=valid,
        message=message,
    )


def audit_molecule(molecule: Molecule) -> list[AtomAudit]:
    return [
        audit_atom(molecule, index)
        for index in range(len(molecule.atoms))
    ]


def validate_molecule(
    molecule: Molecule,
    *,
    strict: bool = False,
) -> tuple[bool, list[AtomAudit]]:
    audits = audit_molecule(molecule)

    valid = all(audit.valid for audit in audits)

    if strict:
        for audit in audits:
            expected = MAX_VALENCE.get(audit.element)

            if expected is None:
                continue

            # Strict mode rejects incomplete ordinary valence.
            # Explicit lone pairs and radicals are retained as supplied.
            if (
                audit.bond_order_sum < expected
                and audit.lone_pairs == 0
                and audit.unpaired_electrons == 0
                and audit.charge == 0
            ):
                valid = False

    return valid, audits