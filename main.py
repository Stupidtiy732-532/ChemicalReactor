from chemistry.parser import ChemicalParser
from chemistry.engine import ReactionEngine


def main():

    engine = ReactionEngine()

    print("Available reactions:\n")
    engine.list_reactions()

    reactant = ChemicalParser("C=C").parse()

    result = engine.run(
        "alkene_hydrogenation",
        reactant,
    )

    print("Reaction:", result.reaction_id)
    print("Reactants:", len(result.reactants))
    print("Products:", len(result.products))

    print("\nObservations:")

    for observation in result.observations:
        print("-", observation)


if __name__ == "__main__":
    main()