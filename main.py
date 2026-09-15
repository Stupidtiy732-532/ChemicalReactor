from chemistry.parser import ChemicalParser
from chemistry.valency import infer_implicit_hydrogens


def main():

    examples = [
        "C-C",
        "C=C",
        "C#C",
        "CH4",
        "CH3-CH(CH3)-CH3",
        "CH3-C(CH3)2-CH3",
        "C(CH3)4",
    ]

    for formula in examples:

        molecule = ChemicalParser(formula).parse()

        infer_implicit_hydrogens(molecule)

        print(formula)
        print(molecule.formula())
        print()


if __name__ == "__main__":
    main()