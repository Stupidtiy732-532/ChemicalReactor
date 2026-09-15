from chemistry.reaction import Reaction, ReactionResult
from chemistry.transform import (
    reduce_double_bond,
    reduce_triple_bond,
)
from chemistry.valency import complete_hydrogens


def alkene_hydrogenation(reactant):

    product = reduce_double_bond(reactant)
    complete_hydrogens(product)

    return ReactionResult(
        reaction_id="alkene_hydrogenation",
        reactants=(reactant,),
        products=(product,),
        observations=(
            "C=C converted to C-C",
            "H2 added across the double bond",
            "Catalyst: Ni, Pt, or Pd",
        ),
    )


def alkyne_partial_hydrogenation(reactant):

    product = reduce_triple_bond(reactant)
    complete_hydrogens(product)

    return ReactionResult(
        reaction_id="alkyne_partial_hydrogenation",
        reactants=(reactant,),
        products=(product,),
        observations=(
            "C#C converted to C=C",
            "Partial hydrogenation",
        ),
    )


REACTIONS = {
    "alkene_hydrogenation": Reaction(
        id="alkene_hydrogenation",
        name="Alkene hydrogenation",
        reagent="H2",
        conditions=("Ni/Pt/Pd catalyst",),
        apply=alkene_hydrogenation,
    ),

    "alkyne_partial_hydrogenation": Reaction(
        id="alkyne_partial_hydrogenation",
        name="Alkyne partial hydrogenation",
        reagent="H2",
        conditions=("catalyst",),
        apply=alkyne_partial_hydrogenation,
    ),
}