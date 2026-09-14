import os

from analyzer import StructureParser, StructureAnalyzer
from reagents import ReagentDatabase
from reactions import ReactionEngine
from species import SpeciesParser


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
    engine = ReactionEngine()
    current_formula = None
    current_groups = []

    while True:
        header("REACTION MODE")

        print(f"Current molecule: {current_formula or 'None'}")
        print(
            "Current reagents: "
            + (engine.reagent_text() or "None")
        )

        print()
        print("1. Enter molecule")
        print("2. Enter any reagent manually")
        print("3. Choose reagent from database")
        print("4. Remove reagent")
        print("5. Clear reagents")
        print("6. Run reaction")
        print("7. Back")

        choice = input("\nSelect: ").strip()

        if choice == "1":
            formula = input(
                "\nEnter molecular formula/structure: "
            ).strip()

            if formula:
                current_formula = formula

                groups_text = input(
                    "\nEnter recognised functional groups "
                    "separated by commas:\n"
                    "Example: primary alcohol, alkene\n"
                    "Groups: "
                ).strip()

                current_groups = [
                    group.strip().lower()
                    for group in groups_text.split(",")
                    if group.strip()
                ]

                print("\nMolecule loaded.")

            pause()

        elif choice == "2":
            reagent = input(
                "\nEnter any reagent, acid, base or species: "
            ).strip()

            if reagent:
                engine.add_reagent(reagent)
                print(f"\nAdded reagent: {reagent}")

            pause()

        elif choice == "3":
            database.display_all()

            number = read_number(
                "\nSelect reagent number: ",
                1,
                len(database.all())
            )

            reagent = database.get_by_number(number)
            engine.add_reagent(reagent.formula)

            print(f"\nAdded reagent: {reagent.formula}")
            pause()

        elif choice == "4":
            if not engine.reagents:
                print("\nNo reagents loaded.")
                pause()
                continue

            for index, reagent in enumerate(
                engine.reagents,
                start=1
            ):
                print(f"{index}. {reagent}")

            number = read_number(
                "\nSelect reagent to remove: ",
                1,
                len(engine.reagents)
            )

            removed = engine.reagents[number - 1]
            engine.remove_reagent(removed)

            print(f"\nRemoved: {removed}")
            pause()

        elif choice == "5":
            engine.clear_reagents()
            print("\nAll reagents cleared.")
            pause()

        elif choice == "6":
            if current_formula is None:
                print("\nEnter a molecule first.")
                pause()
                continue

            if not engine.reagents:
                print("\nEnter at least one reagent.")
                pause()
                continue

            result = engine.react(
                current_formula,
                current_groups
            )

            print("\nREACTION RESULT")
            print("=" * 80)
            print(result)

            pause()

        elif choice == "7":
            return


def species_analysis_mode():
    parser = SpeciesParser()

    while True:
        header("CHEMICAL SPECIES MODE")

        print("1. Dissociate species")
        print("2. Classify species")
        print("3. Back")

        choice = input("\nSelect: ").strip()

        if choice == "3":
            return

        if choice not in {"1", "2"}:
            continue

        formula = input(
            "\nEnter any chemical species: "
        ).strip()

        if not formula:
            continue

        if choice == "1":
            parser.explain(formula)

        elif choice == "2":
            print(
                f"\nClassification: "
                f"{parser.classify(formula)}"
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
        print("5. Chemical species / acid-base mode")
        print("6. Exit")

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
            species_analysis_mode()

        elif choice == "6":
            clear_screen()
            print("Chemicalist closed.")
            break

        else:
            print("Invalid selection.")
            pause()


if __name__ == "__main__":
    main_menu()