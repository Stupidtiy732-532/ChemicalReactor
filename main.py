from chemistry.parser import ChemicalParser
from chemistry.functional_groups import find_functional_groups


def main():

    examples = [
        "C-O-H",
        "C-O-C",
        "C=O",
        "C(=O)-O",
        "C(=O)-N",
        "C-Cl",
    ]

    for notation in examples:

        print(f"Input: {notation}")

        try:
            molecule = ChemicalParser(notation).parse()
            groups = find_functional_groups(molecule)

            print(f"Formula: {molecule.formula()}")

            if not groups:
                print("Groups:  none")
            else:
                for group in groups:
                    print(
                        f"Group:   {group.name} — "
                        f"{group.description}"
                    )

        except Exception as error:
            print(f"ERROR: {error}")

        print()


if __name__ == "__main__":
    main()