from chemistry.parser import ChemicalParser
from chemistry.valency import complete_hydrogens


def main():

    examples = [
        "CH3-CH(CH3)-CH3",
        "CH3-C(CH3)2-CH3",
        "CH3-CH2-CH2-CH3",
        "CH3-CH=CH2",
        "CH3-C#CH",
        "CH3-CH2-OH",
    ]

    for formula in examples:

        try:
            molecule = ChemicalParser(formula).parse()

            explicit_formula = molecule.formula()

            complete_hydrogens(molecule)

            completed_formula = molecule.formula()

            print(f"Input:              {formula}")
            print(f"Before completion:  {explicit_formula}")
            print(f"After completion:   {completed_formula}")
            print(f"Atoms:              {len(molecule.atoms)}")
            print(f"Bonds:              {len(molecule.bonds)}")
            print()

        except Exception as error:

            print(f"Input: {formula}")
            print(f"ERROR: {error}")
            print()


if __name__ == "__main__":
    main()