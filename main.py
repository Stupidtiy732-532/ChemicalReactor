import os
import sys
import json

from parser import parse_structure
from geometry import GeometryEngine
from analyzer import MolecularAnalyzer


# ==========================================================
# TERMINAL UTILITIES
# ==========================================================

def clear_screen():
    """
    Clear the terminal screen.
    """

    os.system("cls" if os.name == "nt" else "clear")


def pause():
    """
    Wait for the user before continuing.
    """

    input("\nPress Enter to continue...")


def print_header(title):
    """
    Print a consistent section header.
    """

    print("=" * 70)
    print(title.center(70))
    print("=" * 70)


def print_banner():
    """
    Print the main application banner.
    """

    print("=" * 70)
    print("ORGREACT".center(70))
    print("Organic Chemistry Molecular Simulator".center(70))
    print("=" * 70)
    print()


def print_error(error):
    """
    Print errors in a readable format.
    """

    print()
    print("[ERROR]")
    print(str(error))
    print()


# ==========================================================
# MOLECULE INPUT
# ==========================================================

def request_structure():
    """
    Ask the user for a molecular structure.
    """

    print("Enter a molecular structure.")
    print("Examples:")
    print("  CH3-CH2-OH")
    print("  C=C")
    print("  C1=CC=CC=C1")
    print("  [NH4^+]")
    print()

    return input("Structure > ").strip()


def parse_user_structure():
    """
    Read and parse a structure from the terminal.
    """

    structure = request_structure()

    if not structure:
        print_error("No structure entered.")
        pause()
        return None

    try:
        molecule = parse_structure(structure)
        return molecule

    except Exception as error:
        print_error(
            f"Could not parse structure:\n{error}"
        )
        pause()
        return None


# ==========================================================
# MOLECULE DISPLAY
# ==========================================================

def show_basic_information(molecule):
    """
    Display basic molecular information.
    """

    print_header("MOLECULE INFORMATION")

    print(f"Formula:       {molecule.formula()}")
    print(f"Atoms:         {len(molecule.atoms)}")
    print(f"Bonds:         {len(molecule.bonds)}")
    print(f"Total charge:  {molecule.total_charge()}")
    print(f"Molar mass:    {molecule.molecular_mass():.6f} u")

    print()
    print("Atoms")
    print("-" * 70)

    for atom_id, atom in molecule.atoms.items():

        print(
            f"ID={atom_id} | "
            f"Element={atom.element} | "
            f"H={atom.hydrogens} | "
            f"Charge={atom.charge} | "
            f"Radicals={atom.radical_electrons} | "
            f"Lone pairs={atom.lone_pairs}"
        )

    print()
    print("Bonds")
    print("-" * 70)

    if not molecule.bonds:
        print("No bonds.")

    else:
        for bond in molecule.bonds:

            print(
                f"{bond.atom1} - {bond.atom2} | "
                f"order={bond.order} | "
                f"sigma={bond.sigma_bonds} | "
                f"pi={bond.pi_bonds}"
            )


def show_structural_analysis(molecule):
    """
    Display the complete structural analysis.
    """

    print_header("STRUCTURAL ANALYSIS")

    analyzer = MolecularAnalyzer(molecule)
    analyzer.print_summary()


def show_geometry(molecule):
    """
    Generate and display molecular coordinates.
    """

    print_header("MOLECULAR GEOMETRY")

    try:
        geometry = GeometryEngine(molecule)
        geometry.generate_layout()

        print("Coordinates in metres:")
        print("-" * 70)

        for atom_id, coordinates in geometry.coordinates().items():

            print(
                f"Atom {atom_id}: "
                f"x={coordinates[0]:.6e} m, "
                f"y={coordinates[1]:.6e} m, "
                f"z={coordinates[2]:.6e} m"
            )

        print()
        print("XYZ FORMAT")
        print("-" * 70)
        print(geometry.export_xyz())

    except Exception as error:
        print_error(
            f"Could not generate geometry:\n{error}"
        )


def show_json(molecule):
    """
    Display the molecule as JSON.
    """

    print_header("MOLECULE JSON")

    try:
        print(molecule.to_json())

    except Exception as error:
        print_error(
            f"Could not export JSON:\n{error}"
        )


def save_json(molecule):
    """
    Save the molecule as a JSON file.
    """

    print_header("SAVE MOLECULE")

    filename = input(
        "Filename [molecule.json] > "
    ).strip()

    if not filename:
        filename = "molecule.json"

    if not filename.lower().endswith(".json"):
        filename += ".json"

    try:
        with open(
            filename,
            "w",
            encoding="utf-8"
        ) as file:

            json.dump(
                molecule.to_dict(),
                file,
                indent=4
            )

        print()
        print(f"Saved molecule to: {filename}")

    except Exception as error:
        print_error(
            f"Could not save file:\n{error}"
        )


def save_xyz(molecule):
    """
    Generate and save an XYZ file.
    """

    print_header("SAVE XYZ")

    filename = input(
        "Filename [molecule.xyz] > "
    ).strip()

    if not filename:
        filename = "molecule.xyz"

    if not filename.lower().endswith(".xyz"):
        filename += ".xyz"

    try:
        geometry = GeometryEngine(molecule)
        geometry.generate_layout()

        xyz_text = geometry.export_xyz()

        with open(
            filename,
            "w",
            encoding="utf-8"
        ) as file:

            file.write(xyz_text)

        print()
        print(f"Saved XYZ file to: {filename}")

    except Exception as error:
        print_error(
            f"Could not save XYZ file:\n{error}"
        )


# ==========================================================
# MOLECULE MENU
# ==========================================================

def molecule_menu(molecule):
    """
    Menu for operations on the currently loaded molecule.
    """

    while True:

        clear_screen()
        print_banner()

        print(
            f"Current molecule: "
            f"{molecule.formula()}"
        )

        print(
            f"Atoms: {len(molecule.atoms)} | "
            f"Bonds: {len(molecule.bonds)} | "
            f"Charge: {molecule.total_charge()}"
        )

        print()
        print("1. Basic molecular information")
        print("2. Structural analysis")
        print("3. Generate geometry / XYZ")
        print("4. Show JSON")
        print("5. Save JSON")
        print("6. Save XYZ")
        print("7. Return to main menu")
        print()

        choice = input("Select option > ").strip()

        if choice == "1":

            clear_screen()
            show_basic_information(molecule)
            pause()

        elif choice == "2":

            clear_screen()
            show_structural_analysis(molecule)
            pause()

        elif choice == "3":

            clear_screen()
            show_geometry(molecule)
            pause()

        elif choice == "4":

            clear_screen()
            show_json(molecule)
            pause()

        elif choice == "5":

            clear_screen()
            save_json(molecule)
            pause()

        elif choice == "6":

            clear_screen()
            save_xyz(molecule)
            pause()

        elif choice == "7":
            return

        else:

            print()
            print("Invalid option.")
            pause()


# ==========================================================
# MAIN MENU
# ==========================================================

def main_menu():
    """
    Main OrgReact application loop.
    """

    while True:

        clear_screen()
        print_banner()

        print("1. Create / analyze molecule")
        print("2. Exit")
        print()

        choice = input("Select option > ").strip()

        if choice == "1":

            molecule = parse_user_structure()

            if molecule is not None:
                molecule_menu(molecule)

        elif choice == "2":

            clear_screen()
            print_banner()
            print("Goodbye.")
            sys.exit(0)

        else:

            print()
            print("Invalid option.")
            pause()


# ==========================================================
# PROGRAM ENTRY POINT
# ==========================================================

if __name__ == "__main__":

    try:
        main_menu()

    except KeyboardInterrupt:

        clear_screen()
        print()
        print("OrgReact closed by user.")
        sys.exit(0)