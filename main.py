import json
import os
import sys

from analyzer import MolecularAnalyzer
from geometry import GeometryEngine
from parser import parse_structure


WIDTH = 70


def clear_screen():
    os.system("cls" if os.name == "nt" else "clear")


def pause():
    input("\nPress Enter to continue...")


def header(title):
    print("=" * WIDTH)
    print(title.center(WIDTH))
    print("=" * WIDTH)


def banner():
    header("CHEMICAL REACTOR")
    print("Organic Chemistry Molecular Structure Simulator".center(WIDTH))
    print()


def show_error(message):
    print(f"\n[ERROR]\n{message}\n")


def request_structure():
    print("Enter a molecular structure.")
    print("Type 'help' for notation help.")
    print("Structure > ", end="")
    return input().strip()


def parse_user_structure():
    structure = request_structure()

    if not structure:
        show_error("No structure entered.")
        pause()
        return None

    try:
        return parse_structure(structure)
    except Exception as error:
        show_error(f"Could not parse structure:\n{error}")
        pause()
        return None


def show_basic_information(molecule):
    header("MOLECULE INFORMATION")

    print(f"Formula:        {molecule.formula()}")
    print(f"Atoms:          {len(molecule.atoms)}")
    print(f"Bonds:          {len(molecule.bonds)}")
    print(f"Total charge:   {molecule.total_charge()}")
    print(f"Molecular mass: {molecule.molecular_mass():.6f} u")

    print("\nAtoms")
    print("-" * WIDTH)

    for atom_id, atom in molecule.atoms.items():
        print(
            f"ID={atom_id} | "
            f"Element={atom.element} | "
            f"H={atom.hydrogens} | "
            f"Charge={atom.charge} | "
            f"Radicals={atom.radical_electrons} | "
            f"Lone pairs={atom.lone_pairs}"
        )

    print("\nBonds")
    print("-" * WIDTH)

    if not molecule.bonds:
        print("No bonds.")
        return

    for bond in molecule.bonds:
        print(
            f"{bond.atom1} - {bond.atom2} | "
            f"order={bond.order} | "
            f"sigma={bond.sigma_bonds} | "
            f"pi={bond.pi_bonds}"
        )


def show_structural_analysis(molecule):
    header("STRUCTURAL ANALYSIS")
    MolecularAnalyzer(molecule).print_summary()


def create_geometry(molecule):
    geometry = GeometryEngine(molecule)
    geometry.generate_layout()
    return geometry


def show_geometry(molecule):
    header("MOLECULAR GEOMETRY")

    try:
        geometry = GeometryEngine(molecule)
        geometry.generate_layout()

        print("Coordinates in metres")
        print("-" * WIDTH)

        for atom_id, (x, y, z) in geometry.coordinates().items():
            print(
                f"Atom {atom_id}: "
                f"x={x:.6e} m, "
                f"y={y:.6e} m, "
                f"z={z:.6e} m"
            )

        print("\nXYZ FORMAT")
        print("-" * WIDTH)
        print(geometry.export_xyz())

    except Exception as error:
        show_error(f"Could not generate geometry:\n{error}")


def show_json(molecule):
    header("MOLECULE JSON")

    try:
        print(molecule.to_json())
    except Exception as error:
        show_error(f"Could not export JSON:\n{error}")


def get_filename(prompt, extension):
    filename = input(prompt).strip()

    if not filename:
        filename = f"molecule.{extension}"

    if not filename.lower().endswith(f".{extension}"):
        filename += f".{extension}"

    return filename


def save_json(molecule):
    header("SAVE MOLECULE")

    filename = get_filename(
        "Filename [molecule.json] > ",
        "json",
    )

    try:
        with open(filename, "w", encoding="utf-8") as file:
            json.dump(
                molecule.to_dict(),
                file,
                indent=4,
            )

        print(f"\nSaved molecule to: {filename}")

    except OSError as error:
        show_error(f"Could not save file:\n{error}")


def save_xyz(molecule):
    header("SAVE XYZ")

    filename = get_filename(
        "Filename [molecule.xyz] > ",
        "xyz",
    )

    try:
        geometry = create_geometry(molecule)

        with open(filename, "w", encoding="utf-8") as file:
            file.write(geometry.export_xyz())

        print(f"\nSaved XYZ file to: {filename}")

    except OSError as error:
        show_error(f"Could not save XYZ file:\n{error}")
    except Exception as error:
        show_error(f"Could not generate XYZ file:\n{error}")


MOLECULE_ACTIONS = {
    "1": show_basic_information,
    "2": show_structural_analysis,
    "3": show_geometry,
    "4": show_json,
    "5": save_json,
    "6": save_xyz,
}


def molecule_menu(molecule):
    actions = {
        "1": show_basic_information,
        "2": show_structural_analysis,
        "3": show_geometry,
        "4": show_json,
        "5": save_json,
        "6": save_xyz,
    }

    while True:
        clear_screen()
        banner()

        print(f"Current molecule: {molecule.formula()}")
        print(
            f"Atoms: {len(molecule.atoms)} | "
            f"Bonds: {len(molecule.bonds)} | "
            f"Charge: {molecule.total_charge()}"
        )

        print(
            "\n"
            "1. Basic molecular information\n"
            "2. Structural analysis\n"
            "3. Generate geometry / XYZ\n"
            "4. Show JSON\n"
            "5. Save JSON\n"
            "6. Save XYZ\n"
            "7. Return to main menu\n"
        )

        choice = input("Select option > ").strip()

        if choice == "7":
            return

        action = actions.get(choice)

        if action is None:
            print("\nInvalid option.")
            pause()
            continue

        clear_screen()
        action(molecule)
        pause()


def show_help():
    clear_screen()
    banner()
    header("STRUCTURE NOTATION")

    print("""
Atoms
-----
C        Carbon
O        Oxygen
N        Nitrogen
Cl       Chlorine
Na       Sodium

Hydrogens
---------
C H      One hydrogen
CH3      Three hydrogens
[NH4^+]  Four hydrogens inside brackets

Bonds
-----
-        Single bond
=        Double bond
#        Triple bond
≡        Triple bond
-=       Triple bond
=-       Triple bond

Branches
--------
( ... )  Branch from the previous atom

Rings
-----
@number  Ring connection marker

Example:
C@1-C-C-C-C-C@1

Isotopes
--------
_element_mass

Example:
C_13
O_18

Charges
-------
^+       +1
^-       -1
^2+      +2
^+2      +2
^++      +2
^--      -2

Radicals
--------
.        One radical electron

Examples:
C.
..C

Lone pairs
----------
:        One lone pair

Examples:
O:
O::

Markers
-------
<name>       Attach a marker
X<name>      Attach a marker

Examples:
C<central>
X<central>
[O<X>]

Bracket atoms
-------------
[ ... ]  Explicit atom with modifiers

Examples:
[NH4^+]
[O^2-]
[C_13H4]

Important
---------
This notation uses explicit atoms.
It does NOT use SMILES-style implicit
carbon backbones.

Spaces are ignored.

Press Enter to return.
""")

    input()


def main_menu():
    while True:
        clear_screen()
        banner()

        print(
            "1. Create / analyze molecule\n"
            "2. Notation help\n"
            "3. Exit\n"
        )

        choice = input("Select option > ").strip().lower()

        if choice == "1":
            molecule = parse_user_structure()

            if molecule is not None:
                molecule_menu(molecule)

        elif choice in ("2", "help", "h"):
            show_help()

        elif choice in ("3", "exit", "q"):
            clear_screen()
            banner()
            print("Goodbye.")
            return

        else:
            print("\nInvalid option.")
            pause()


if __name__ == "__main__":
    try:
        main_menu()
    except KeyboardInterrupt:
        clear_screen()
        print("\nChemical Reactor closed by user.")
        sys.exit(0)