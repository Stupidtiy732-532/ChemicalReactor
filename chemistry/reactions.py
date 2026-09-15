from chemistry.reaction import Reaction, ReactionResult
from chemistry.transform import (
    reduce_double_bond,
    reduce_triple_bond,
)


def alkene_hydrogenation(molecule):
    product = reduce_double_bond(molecule)

    return ReactionResult(
        success=True,
        reactants=[molecule],
        products=[product],
        message="Reduced one explicit C=C bond to C-C.",
    )


def alkyne_partial_hydrogenation(molecule):
    product = reduce_triple_bond(molecule)

    return ReactionResult(
        success=True,
        reactants=[molecule],
        products=[product],
        message="Reduced one explicit C#C bond to C=C.",
    )


REACTIONS = [
    Reaction(
        id="alkene_hydrogenation",
        name="Alkene hydrogenation",
        reactants=("C=C",),
        products=("C-C",),
        handler=alkene_hydrogenation,
    ),
    Reaction(
        id="alkyne_partial_hydrogenation",
        name="Partial alkyne hydrogenation",
        reactants=("C#C",),
        products=("C=C",),
        handler=alkyne_partial_hydrogenation,
    ),
]