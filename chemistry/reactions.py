from chemistry.parser import ChemicalParser
from chemistry.reaction import Reaction, ReactionResult


def parse(formula: str):
    return ChemicalParser(formula).parse()


def alkane_chlorination(reactant):

    product = parse("C-Cl")

    return ReactionResult(
        reaction_id="alkane_chlorination",
        reactants=(reactant,),
        products=(product,),
        observations=(
            "Free-radical substitution",
            "Requires UV light or heat",
            "Hydrogen is substituted by chlorine",
        ),
    )


def alkene_hydrogenation(reactant):

    product = parse("C-C")

    return ReactionResult(
        reaction_id="alkene_hydrogenation",
        reactants=(reactant,),
        products=(product,),
        observations=(
            "Addition reaction",
            "Requires H2 and Ni/Pt/Pd catalyst",
            "C=C becomes C-C",
        ),
    )


def alkene_bromination(reactant):

    product = parse("C-C")

    return ReactionResult(
        reaction_id="alkene_bromination",
        reactants=(reactant,),
        products=(product,),
        observations=(
            "Electrophilic addition",
            "Bromine adds across the double bond",
        ),
    )


REACTIONS = {
    "alkane_chlorination": Reaction(
        id="alkane_chlorination",
        name="Free-radical chlorination of alkane",
        reagent="Cl2",
        conditions=("UV light", "heat"),
        apply=alkane_chlorination,
    ),

    "alkene_hydrogenation": Reaction(
        id="alkene_hydrogenation",
        name="Catalytic hydrogenation of alkene",
        reagent="H2",
        conditions=("Ni/Pt/Pd", "heat if required"),
        apply=alkene_hydrogenation,
    ),

    "alkene_bromination": Reaction(
        id="alkene_bromination",
        name="Bromination of alkene",
        reagent="Br2",
        conditions=("inert solvent",),
        apply=alkene_bromination,
    ),
}