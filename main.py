from chemistry.parser import ChemicalParser
from chemistry.engine import ReactionEngine


def test(notation: str, reaction_id: str):

    reactant = ChemicalParser(notation).parse()

    print(f"Reactant: {notation}")
    print(f"Before:  {reactant.formula()}")

    engine = ReactionEngine()
    result = engine.run(reaction_id, reactant)

    product = result.products[0]

    print(f"After:   {product.formula()}")

    for observation in result.observations:
        print(f"- {observation}")

    print()


def main():

    test("C=C", "alkene_hydrogenation")
    test("C#C", "alkyne_partial_hydrogenation")


if __name__ == "__main__":
    main()