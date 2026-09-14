# main.py

import os

from analyzer import analyze, parse_structure
from reactions import ReactionEngine


def clear():
    os.system("cls" if os.name == "nt" else "clear")


def pause():
    input("\nPress Enter to continue...")


def structure_mode():
    clear()

    print("STRUCTURE ANALYSIS")
    print("=" * 60)
    print("Examples:")
    print("  CH3CH2OH")
    print("  CH3CHO")
    print("  CH3COCH3")
    print("  CH3COOH")
    print("  CH3CH2Br")
    print()

    text = input("Enter structure: ").strip()

    try:
        result = analyze(text)

        print("\nRESULT")
        print("-" * 60)
        print("Formula:", result["formula"])
        print("Molar mass:", round(result["molar_mass"], 3), "g/mol")
        print("Structure:", result["structure"])
        print("Functional groups:")

        for group in result["functional_groups"]:
            print("  -", group)

    except Exception as error:
        print("\nERROR:", error)

    pause()


def reaction_mode():
    engine = ReactionEngine()

    while True:
        clear()

        print("REACTION MODE")
        print("=" * 60)
        print("Enter 'back' to return to the main menu.")
        print()

        reactant = input("Reactant structure: ").strip()

        if reactant.lower() == "back":
            return

        reagent = input("Reagent: ").strip()

        if reagent.lower() == "back":
            return

        conditions = input(
            "Conditions, if any "
            "(aq / alc / heat / acidic): "
        ).strip()

        try:
            molecule = parse_structure(reactant)
            result = engine.react(
                molecule,
                reagent,
                conditions,
            )

            clear()
            result.display()

        except Exception as error:
            print("\nREACTION ERROR")
            print("-" * 60)
            print(error)

        choice = input(
            "\nPress Enter for another reaction, "
            "or type 'back': "
        ).strip()

        if choice.lower() == "back":
            return


def quick_examples():
    clear()

    print("QUICK EXAMPLES")
    print("=" * 60)

    engine = ReactionEngine()

    examples = [
        ("CH3CH2OH", "PCC", ""),
        ("CH3CH2OH", "KMnO4", ""),
        ("CH3CHO", "NaBH4", ""),
        ("CH3COCH3", "LiAlH4", ""),
        ("CH3CH2Br", "KOH", "aq"),
        ("CH3CH2Br", "KOH", "alc heat"),
        ("CH3CH=CH2", "Br2", ""),
        ("CH3CH=CH2", "H2/Ni", ""),
    ]

    for reactant, reagent, conditions in examples:
        try:
            molecule = parse_structure(reactant)
            result = engine.react(
                molecule,
                reagent,
                conditions,
            )

            print(result.equation)

        except Exception as error:
            print(
                f"{reactant} + {reagent}: "
                f"ERROR: {error}"
            )

    pause()


def main():
    while True:
        clear()

        print("CHEMICALIST — ORGANIC CHEMISTRY SIMULATOR")
        print("=" * 60)
        print("1. Structure analysis")
        print("2. Reaction mode")
        print("3. Quick reaction examples")
        print("4. Exit")
        print()

        choice = input("Choose: ").strip()

        if choice == "1":
            structure_mode()

        elif choice == "2":
            reaction_mode()

        elif choice == "3":
            quick_examples()

        elif choice == "4":
            clear()
            print("Exiting Chemicalist.")
            break

        else:
            print("Invalid choice.")
            pause()


if __name__ == "__main__":
    main()