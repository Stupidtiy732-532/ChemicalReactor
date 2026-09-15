from chemistry.parser import parse
from chemistry.electron_validation import validate_molecule
from chemistry.reactions import (
    alkene_hydrogenation,
    alkyne_partial_hydrogenation,
)


def show_audit(label, molecule):
    valid, audits = validate_molecule(molecule)

    print(f"\n{label}: {molecule}")
    print("Formula:", molecule.formula())
    print("Valid:", valid)

    for audit in audits:
        print(
            f"  atom {audit.atom_index}: "
            f"{audit.element}, "
            f"bond-order={audit.bond_order_sum}, "
            f"lone-pairs={audit.lone_pairs}, "
            f"unpaired={audit.unpaired_electrons}, "
            f"charge={audit.charge} "
            f"-> {audit.message}"
        )


alkene = parse("C=C")
show_audit("Alkene before", alkene)

alkene_product = alkene_hydrogenation(alkene).products[0]
show_audit("Alkene after", alkene_product)


alkyne = parse("C#C")
show_audit("Alkyne before", alkyne)

alkyne_product = alkyne_partial_hydrogenation(alkyne).products[0]
show_audit("Alkyne after", alkyne_product)