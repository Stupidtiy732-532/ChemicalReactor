"""
Reaction engine.

Known reactions are handled explicitly.
Unknown reagents are still accepted and analysed.
"""

from species import SpeciesParser


class ReactionResult:

    def __init__(
        self,
        success=False,
        reaction="",
        explanation="",
        products=None,
        product=None
    ):
        self.success = success
        self.reaction = reaction
        self.explanation = explanation
        self.products = products or []
        self.product = product

    def __str__(self):
        output = []

        output.append(
            "Reaction detected: "
            + ("YES" if self.success else "NO")
        )

        if self.reaction:
            output.append(f"\nReaction:\n{self.reaction}")

        if self.explanation:
            output.append(f"\nExplanation:\n{self.explanation}")

        if self.products:
            output.append(
                "\nProducts:\n"
                + "\n".join(f"  - {product}" for product in self.products)
            )

        return "\n".join(output)


class ReactionEngine:

    def __init__(self):
        self.reagents = []
        self.species_parser = SpeciesParser()

    def add_reagent(self, reagent):
        reagent = reagent.strip()

        if reagent and reagent not in self.reagents:
            self.reagents.append(reagent)

    def remove_reagent(self, reagent):
        if reagent in self.reagents:
            self.reagents.remove(reagent)

    def clear_reagents(self):
        self.reagents.clear()

    def reagent_text(self):
        return ", ".join(self.reagents)

    def react(self, formula, functional_groups=None):
        functional_groups = functional_groups or []

        if not self.reagents:
            return ReactionResult(
                success=False,
                explanation="No reagents have been entered."
            )

        all_results = []

        for reagent in self.reagents:
            result = self.react_single(
                formula,
                reagent,
                functional_groups
            )

            all_results.append(result)

        return self.combine_results(all_results)

    def react_single(
        self,
        formula,
        reagent,
        functional_groups
    ):
        formula = formula.strip()
        reagent = reagent.strip()

        reagent_lower = reagent.lower()

        # -------------------------------------------------
        # ACID-BASE DISSOCIATION
        # -------------------------------------------------

        if reagent in self.species_parser.ACIDS:
            ions = self.species_parser.dissociate(reagent)

            return ReactionResult(
                success=True,
                reaction=f"{reagent} → {' + '.join(ions)}",
                explanation=(
                    "The entered reagent is recognised as an acid. "
                    "Its acidic proton is represented as H+."
                ),
                products=ions
            )

        if reagent in self.species_parser.BASES:
            ions = self.species_parser.dissociate(reagent)

            return ReactionResult(
                success=True,
                reaction=f"{reagent} → {' + '.join(ions)}",
                explanation=(
                    "The entered reagent is recognised as a base. "
                    "It produces hydroxide or another basic species."
                ),
                products=ions
            )

        if reagent in self.species_parser.SALTS:
            ions = self.species_parser.dissociate(reagent)

            return ReactionResult(
                success=True,
                reaction=f"{reagent} → {' + '.join(ions)}",
                explanation="The entered salt was split into its ions.",
                products=ions
            )

        # -------------------------------------------------
        # SODIUM METAL + ALCOHOL
        # -------------------------------------------------

        if reagent in {"Na", "sodium"}:
            if "alcohol" in functional_groups:
                return ReactionResult(
                    success=True,
                    reaction=(
                        f"{formula} + Na → "
                        f"{formula.replace('OH', 'ONa')} + H2"
                    ),
                    explanation=(
                        "Sodium removes the acidic hydrogen of the "
                        "alcoholic OH group."
                    ),
                    products=[
                        formula.replace("OH", "ONa"),
                        "H2"
                    ]
                )

            if "carboxylic acid" in functional_groups:
                return ReactionResult(
                    success=True,
                    reaction=(
                        f"{formula} + Na → "
                        f"{formula.replace('COOH', 'COONa')} + H2"
                    ),
                    explanation=(
                        "Sodium reacts with the acidic proton of "
                        "the carboxylic acid."
                    ),
                    products=[
                        formula.replace("COOH", "COONa"),
                        "H2"
                    ]
                )

        # -------------------------------------------------
        # BICARBONATE TEST
        # -------------------------------------------------

        if reagent in {"NaHCO3", "sodium bicarbonate"}:
            if "carboxylic acid" in functional_groups:
                salt = formula.replace("COOH", "COONa")

                return ReactionResult(
                    success=True,
                    reaction=(
                        f"{formula} + NaHCO3 → "
                        f"{salt} + CO2 + H2O"
                    ),
                    explanation=(
                        "Carboxylic acids release carbon dioxide "
                        "with sodium bicarbonate."
                    ),
                    products=[salt, "CO2", "H2O"]
                )

        # -------------------------------------------------
        # OXIDATION
        # -------------------------------------------------

        oxidants = {
            "KMnO4",
            "K2Cr2O7",
            "K2Cr2O7/H+",
            "PCC",
            "PDC",
            "CrO3",
            "H2CrO4",
            "Jones reagent",
            "O3",
            "Tollens",
            "Fehling",
        }

        if reagent in oxidants or "oxid" in reagent_lower:
            if "primary alcohol" in functional_groups:
                return ReactionResult(
                    success=True,
                    reaction=f"{formula} → aldehyde / carboxylic acid",
                    explanation=(
                        "A primary alcohol can be oxidised first "
                        "to an aldehyde and further to a carboxylic acid."
                    ),
                    products=["aldehyde", "carboxylic acid"]
                )

            if "secondary alcohol" in functional_groups:
                return ReactionResult(
                    success=True,
                    reaction=f"{formula} → ketone",
                    explanation="Secondary alcohol oxidation gives a ketone.",
                    products=["ketone"]
                )

            if "aldehyde" in functional_groups:
                return ReactionResult(
                    success=True,
                    reaction=f"{formula} → carboxylic acid",
                    explanation="Aldehydes can be oxidised to carboxylic acids.",
                    products=["carboxylic acid"]
                )

        # -------------------------------------------------
        # REDUCTION
        # -------------------------------------------------

        reductants = {
            "H2/Ni",
            "H2/Pd",
            "H2/Pt",
            "NaBH4",
            "LiAlH4",
            "DIBAL-H",
            "Zn/HCl",
            "Zn-Hg/HCl",
            "NH2NH2/KOH",
        }

        if reagent in reductants or "reduc" in reagent_lower:
            if "alkene" in functional_groups:
                return ReactionResult(
                    success=True,
                    reaction=f"{formula} + H2 → alkane",
                    explanation="Catalytic hydrogenation reduces C=C.",
                    products=["alkane"]
                )

            if "alkyne" in functional_groups:
                return ReactionResult(
                    success=True,
                    reaction=f"{formula} + 2H2 → alkane",
                    explanation="Complete hydrogenation reduces C≡C.",
                    products=["alkane"]
                )

            if "aldehyde" in functional_groups:
                return ReactionResult(
                    success=True,
                    reaction=f"{formula} → primary alcohol",
                    explanation="Aldehydes reduce to primary alcohols.",
                    products=["primary alcohol"]
                )

            if "ketone" in functional_groups:
                return ReactionResult(
                    success=True,
                    reaction=f"{formula} → secondary alcohol",
                    explanation="Ketones reduce to secondary alcohols.",
                    products=["secondary alcohol"]
                )

        # -------------------------------------------------
        # HALOGENATION
        # -------------------------------------------------

        if reagent in {"Br2", "Cl2"}:
            if "alkene" in functional_groups:
                return ReactionResult(
                    success=True,
                    reaction=f"{formula} + {reagent} → vicinal dihalo compound",
                    explanation=(
                        "Halogen adds across the carbon-carbon double bond."
                    ),
                    products=["vicinal dihalo compound"]
                )

        # -------------------------------------------------
        # SUBSTITUTION
        # -------------------------------------------------

        if reagent in {
            "KOH(aq)",
            "NaOH(aq)",
            "KOH",
            "NaOH"
        }:
            if "halo compound" in functional_groups:
                return ReactionResult(
                    success=True,
                    reaction=f"{formula} → alcohol",
                    explanation=(
                        "Aqueous hydroxide promotes nucleophilic "
                        "substitution of a haloalkane."
                    ),
                    products=["alcohol"]
                )

        # -------------------------------------------------
        # ELIMINATION
        # -------------------------------------------------

        if reagent in {
            "KOH(alc)",
            "NaOH(alc)",
            "alcoholic KOH",
            "alcoholic NaOH"
        }:
            if "halo compound" in functional_groups:
                return ReactionResult(
                    success=True,
                    reaction=f"{formula} → alkene",
                    explanation=(
                        "Alcoholic hydroxide and heat promote "
                        "β-elimination."
                    ),
                    products=["alkene"]
                )

        # -------------------------------------------------
        # UNKNOWN REAGENT
        # -------------------------------------------------

        return ReactionResult(
            success=False,
            reaction="No implemented transformation",
            explanation=(
                f"'{reagent}' was accepted as an arbitrary reagent, "
                "but no reaction rule currently matches it."
            ),
            products=[]
        )

    def combine_results(self, results):
        successful = [
            result for result in results
            if result.success
        ]

        if not successful:
            return ReactionResult(
                success=False,
                explanation="\n\n".join(
                    result.explanation for result in results
                )
            )

        reaction_text = "\n".join(
            result.reaction for result in successful
        )

        explanations = "\n".join(
            result.explanation for result in successful
        )

        products = []

        for result in successful:
            products.extend(result.products)

        return ReactionResult(
            success=True,
            reaction=reaction_text,
            explanation=explanations,
            products=products
        )