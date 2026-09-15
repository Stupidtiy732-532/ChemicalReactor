from chemistry.reactions import REACTIONS


class ReactionEngine:

    def __init__(self):
        self.reactions = REACTIONS

    def list_reactions(self):

        for reaction in self.reactions.values():

            print(f"{reaction.id}: {reaction.name}")
            print(f"  Reagent: {reaction.reagent}")
            print(f"  Conditions: {', '.join(reaction.conditions)}")
            print()

    def run(self, reaction_id: str, reactant):

        if reaction_id not in self.reactions:
            raise ValueError(
                f"Unknown reaction: {reaction_id}"
            )

        reaction = self.reactions[reaction_id]

        return reaction.apply(reactant)