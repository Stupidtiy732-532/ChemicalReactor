from analyzer import StructureParser, MoleculeAnalyzer
from reactions import ReactionEngine


def display_molecule(molecule):
    print()
    print("=" * 65)
    print("MOLECULE ANALYSIS")
    print("=" * 65)

    print(f"Input:      {molecule.original_input}")

    print(f"Formula:    {molecule.formula}")

    print(
        f"Molar mass: "
        f"{molecule.molar_mass:.3f} g/mol"
    )

    print()
    print("Atoms:")

    for atom in molecule.atoms:
        print(
            f"  Atom {atom.id}: "
            f"{atom.element}, "
            f"implicit H = {atom.hydrogens}"
        )

    print()
    print("Bonds:")

    for bond in molecule.bonds:
        print(
            f"  {bond.atom1} - "
            f"{bond.atom2}, "
            f"order = {bond.order}"
        )

    print()
    print("Parent carbon chain:")

    if molecule.parent_chain:
        chain_text = " - ".join(
            f"C{index + 1}"
            for index in range(
                len(molecule.parent_chain)
            )
        )

        print(f"  {chain_text}")

        print(
            f"  Length: "
            f"{len(molecule.parent_chain)} "
            f"carbon atoms"
        )

    else:
        print("  No carbon parent chain found.")

    print()
    print("Functional groups:")

    if molecule.functional_groups:
        for group in molecule.functional_groups:
            print(
                f"  {group.name}: "
                f"atoms {group.atom_ids}"
            )

    else:
        print("  None detected.")

    print()
    print("Substituents:")

    if molecule.substituents:
        for substituent in molecule.substituents:
            print(
                f"  {substituent.name}: "
                f"atoms {substituent.atom_ids}, "
                f"parent locant = "
                f"{substituent.parent_locant}"
            )
    else:
        print("  None detected.")

    print("=" * 65)


def run_reaction_loop(molecule):
    reaction_engine = ReactionEngine()

    while True:
        reagent_input = input(
            "\nEnter reagent, or press Enter to skip reaction: "
        ).strip()

        if reagent_input == "":
            break

        reagent = reaction_engine.database.find(
            reagent_input
        )

        if reagent is None:
            print("Unknown reagent.")
            continue

        print(
            f"Recognized reagent: "
            f"{reagent.name}"
        )

        if reagent.conditions:
            print(
                f"Conditions: "
                f"{reagent.conditions}"
            )

        excess_input = input(
            "Use reagent in excess? (y/n): "
        ).strip().lower()

        excess = excess_input == "y"

        reagent_moles = None

        if not excess:
            while True:
                try:
                    reagent_moles = float(input("Enter reagent amount in mol: "))
                    if reagent_moles < 0:
                        print("Amount cannot be negative.")
                        continue
                    break
                except ValueError:
                    print("Enter a valid numerical amount.")

        try:
            result = reaction_engine.react(
                molecule,
                reagent_input,
                reagent_moles,
                excess,
            )
            result.display()

            if result.product_molecule is not None:
                use_product = input("\nUse product as next molecule? (y/n): ").strip().lower()

                if use_product == "y":
                    molecule = result.product_molecule
                    analyzer = MoleculeAnalyzer()
                    analyzer.analyze(molecule)

                    print("\nUpdated molecule:")

                    display_molecule(molecule)
        except Exception as error:
            print()
            print(f"Reaction error: {error}")


def main():
    print("Organic Chemistry Simulator — Version 2")
    print("Open-chain aliphatic molecule analyzer")
    print("Type 'exit' to quit.")

    while True:
        print()
        text = input("Molecule > ").strip()

        if text.lower() == "exit":
            print("Exiting.")
            break

        if not text:
            continue

        try:
            parser = StructureParser(text)

            molecule = parser.parse()
            analyzer = MoleculeAnalyzer()
            analyzer.analyze(molecule)

            display_molecule(molecule)
            run_reaction_loop(molecule)

        except Exception as error:
            print()
            print(
                f"Error: {error}"
            )


if __name__ == "__main__":
    main()