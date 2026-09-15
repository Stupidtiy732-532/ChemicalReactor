from chemistry.parser import ChemicalParser
from chemistry.valency import complete_hydrogens
from chemistry.validator import assert_valid_molecule


def main():

    examples = [
        "C-C",
        "C=C",
        "C#C",
        "C-O",
        "C-N",
        "CH3-CH2-OH",
        "CH3-CH=CH2",
        "CH3-C#CH",
    ]

    for notation in examples:

        print(f"Input: {notation}")

        try:
            molecule = ChemicalParser(notation).parse()

            assert_valid_molecule(molecule)

            complete_hydrogens(molecule)

            print(f"Formula: {molecule.formula()}")
            print(f"Atoms:   {len(molecule.atoms)}")
            print(f"Bonds:   {len(molecule.bonds)}")
            print("Status:  valid")

        except Exception as error:

            print(f"ERROR:   {error}")

        print()


if __name__ == "__main__":
    main()