import os

from analyzer import StructureParser, StructureAnalyzer
from reagents import ReagentDatabase


def clear_screen():
    os.system("cls" if os.name == "nt" else "clear")


def pause():
    input("\nPress Enter to continue...")


def header(title):
    clear_screen()

    print("=" * 80)
    print(f"{title:^80}")
    print("=" * 80)


def read_number(prompt, minimum=None, maximum=None):
    while True:
        try:
            value = int(input(prompt))

            if minimum is not None and value < minimum:
                print(f"Enter a value >= {minimum}")
                continue

            if maximum is not None and value > maximum:
                print(f"Enter a value <= {maximum}")
                continue

            return value

        except ValueError:
            print("Enter a valid number.")


def read_structure():
    while True:
        notation = input(
            "\nEnter condensed structure "
            "(or type 'back'): "
        ).strip()

        if notation.lower() == "back":
            return None

        try:
            parser = StructureParser(notation)
            atoms, bonds = parser.parse()

            return atoms, bonds

        except ValueError as error:
            print(f"\nError: {error}")


def structure_analysis_mode():
    while True:
        header("STRUCTURE ANALYSIS MODE")

        print("1. Analyse a structure")
        print("2. Back")

        choice = input("\nSelect: ").strip()

        if choice == "2":
            return

        if choice != "1":
            continue

        structure = read_structure()

        if structure is None:
            continue

        atoms, bonds = structure

        analyzer = StructureAnalyzer(atoms, bonds)
        analyzer.report()

        pause()


def reagent_database_mode(database):
    while True:
        header("REAGENT DATABASE")

        print("1. List all reagents")
        print("2. List reagent categories")
        print("3. Search reagents")
        print("4. View reagent details")
        print("5. Back")

        choice = input("\nSelect: ").strip()

        if choice == "1":
            header("ALL AVAILABLE REAGENTS")
            database.display_all()
            pause()

        elif choice == "2":
            header("REAGENT CATEGORIES")
            database.display_categories()
            pause()

        elif choice == "3":
            query = input("\nSearch term: ")
            results = database.search(query)

            print("\nSEARCH RESULTS")
            print("=" * 80)

            if not results:
                print("No matching reagents.")

            for index, reagent in enumerate(results, start=1):
                print(
                    f"{index}. {reagent.name} "
                    f"({reagent.formula})"
                )

            pause()

        elif choice == "4":
            database.display_all()

            number = read_number(
                "\nEnter reagent number: ",
                1,
                len(database.all())
            )

            reagent = database.get_by_number(number)

            print("\nREAGENT DETAILS")
            print("=" * 60)
            print(f"Name       : {reagent.name}")
            print(f"Formula    : {reagent.formula}")
            print(f"Category   : {reagent.category}")
            print(f"Uses       : {reagent.uses}")
            print(f"Conditions : {reagent.conditions}")

            pause()

        elif choice == "5":
            return


def reaction_mode(database):
    while True:
        header("REACTION MODE")

        print("1. Enter molecule and reagent")
        print("2. View available reagents")
        print("3. Back")

        choice = input("\nSelect: ").strip()

        if choice == "3":
            return

        if choice == "2":
            database.display_all()
            pause()
            continue

        if choice != "1":
            continue

        structure = read_structure()

        if structure is None:
            continue

        atoms, bonds = structure

        analyzer = StructureAnalyzer(atoms, bonds)

        print("\nStarting molecule:")
        print(f"Formula: {analyzer.formula()}")
        print(f"Molar mass: {analyzer.molar_mass():.3f} g mol^-1")

        database.display_all()

        reagent_number = read_number(
            "\nSelect reagent number: ",
            1,
            len(database.all())
        )

        reagent = database.get_by_number(reagent_number)

        print("\nSELECTED REAGENT")
        print("=" * 60)
        print(f"Name     : {reagent.name}")
        print(f"Formula  : {reagent.formula}")
        print(f"Category : {reagent.category}")
        print(f"Uses     : {reagent.uses}")
        print(f"Conditions: {reagent.conditions}")

        print("\nREACTION ENGINE STATUS")
        print("=" * 60)

        # Current engine connection point.
        # Existing ReactionEngine can be connected here.
        #
        # Example:
        #
        # engine = ReactionEngine()
        # result = engine.react(molecule, reagent)
        #
        # At present this displays the selected reaction
        # and leaves actual transformation to reactions.py.

        print("\nReaction selected successfully.")

        print(
            "\nNote: The reaction engine must now map "
            "this reagent category to a transformation."
        )

        pause()


def functional_group_database_mode():
    functional_groups = {
        "Alkane": "C-C single bonds",
        "Alkene": "C=C",
        "Alkyne": "C≡C",
        "Haloalkane": "C-X, X = F, Cl, Br, I",
        "Alcohol": "R-OH",
        "Phenol": "Ar-OH",
        "Ether": "R-O-R",
        "Aldehyde": "R-CHO",
        "Ketone": "R-CO-R",
        "Carboxylic acid": "R-COOH",
        "Ester": "R-COOR",
        "Acid chloride": "R-COCl",
        "Amine": "R-NH2, R2NH, R3N",
        "Amide": "R-CONH2",
        "Nitrile": "R-CN",
        "Nitro compound": "R-NO2",
        "Azo compound": "R-N=N-R",
        "Diazonium salt": "Ar-N2+",
        "Thiol": "R-SH",
        "Sulphide": "R-S-R",
        "Disulphide": "R-S-S-R",
        "Grignard reagent": "R-Mg-X",
        "Organolithium": "R-Li",
        "Enol": "C=C-OH",
        "Enolate": "C=C-O−",
        "Acetal": "R-CH(OR)2",
        "Hemiacetal": "R-CH(OH)(OR)",
    }

    while True:
        header("FUNCTIONAL-GROUP DATABASE")

        names = list(functional_groups)

        for index, name in enumerate(names, start=1):
            print(f"{index:02d}. {name}")

        print(f"{len(names) + 1:02d}. Back")

        choice = read_number(
            "\nSelect functional group: ",
            1,
            len(names) + 1
        )

        if choice == len(names) + 1:
            return

        name = names[choice - 1]

        print("\nFUNCTIONAL GROUP")
        print("=" * 60)
        print(f"Name       : {name}")
        print(f"Structure   : {functional_groups[name]}")

        pause()


def main_menu():
    database = ReagentDatabase()

    while True:
        header("CHEMICALIST — ORGANIC CHEMISTRY SIMULATOR")

        print("1. Reaction mode")
        print("2. Structure analysis mode")
        print("3. Reagent database")
        print("4. Functional-group database")
        print("5. Exit")

        choice = input("\nSelect mode: ").strip()

        if choice == "1":
            reaction_mode(database)

        elif choice == "2":
            structure_analysis_mode()

        elif choice == "3":
            reagent_database_mode(database)

        elif choice == "4":
            functional_group_database_mode()

        elif choice == "5":
            clear_screen()
            print("Chemicalist closed.")
            break

        else:
            print("Invalid selection.")
            pause()


if __name__ == "__main__":
    main_menu()