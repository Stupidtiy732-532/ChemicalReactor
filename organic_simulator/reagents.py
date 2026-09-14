class Reagent:
    def __init__(
        self,
        name,
        aliases,
        reaction_type,
        conditions="",
    ):
        self.name = name
        self.aliases = aliases
        self.reaction_type = reaction_type
        self.conditions = conditions

    def matches(self, text):
        text = text.lower().strip()

        return (
            text == self.name.lower()
            or text in [
                alias.lower()
                for alias in self.aliases
            ]
        )


class ReagentDatabase:
    def __init__(self):
        self.reagents = [
            Reagent(
                "Na",
                [
                    "sodium",
                    "metallic sodium",
                ],
                "sodium",
            ),

            Reagent(
                "NaHCO3",
                [
                    "nahco3",
                    "sodium bicarbonate",
                    "baking soda",
                ],
                "bicarbonate",
            ),

            Reagent(
                "H2/Ni",
                [
                    "h2",
                    "hydrogen",
                    "h2 ni",
                    "hydrogen nickel",
                ],
                "hydrogenation",
                "Nickel catalyst",
            ),

            Reagent(
                "Br2",
                [
                    "br2",
                    "bromine",
                ],
                "bromination",
            ),

            Reagent(
                "KOH(aq)",
                [
                    "koh aq",
                    "aqueous koh",
                    "aq koh",
                ],
                "aqueous_koh",
            ),

            Reagent(
                "KOH(alc)",
                [
                    "koh alc",
                    "alcoholic koh",
                    "alc koh",
                ],
                "alcoholic_koh",
            ),

            Reagent(
                "KMnO4",
                [
                    "kmno4",
                    "potassium permanganate",
                ],
                "permanganate",
            ),

            Reagent(
                "K2Cr2O7/H+",
                [
                    "k2cr2o7",
                    "acidified dichromate",
                    "acidified potassium dichromate",
                ],
                "dichromate",
            ),

            Reagent(
                "PCC",
                [
                    "pcc",
                    "pyridinium chlorochromate",
                ],
                "pcc",
            ),

            Reagent(
                "Tollens",
                [
                    "tollens reagent",
                    "ammoniacal silver nitrate",
                ],
                "tollens",
            ),

            Reagent(
                "Fehling",
                [
                    "fehling solution",
                    "fehlings solution",
                ],
                "fehling",
            ),

            Reagent(
                "ROH/H+",
                [
                    "alcohol acid",
                    "esterification",
                    "roh h+",
                ],
                "esterification",
            ),
        ]

    def find(self, text):
        for reagent in self.reagents:
            if reagent.matches(text):
                return reagent

        return None