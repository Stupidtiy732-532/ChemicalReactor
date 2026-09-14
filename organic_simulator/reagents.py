"""
Chemicalist reagent database.

Focused on NCERT/JEE organic chemistry.
"""

from dataclasses import dataclass


@dataclass
class Reagent:
    name: str
    formula: str
    category: str
    uses: str
    conditions: str = ""
    danger: str = ""


class ReagentDatabase:

    def __init__(self):
        self.reagents = [
            # -------------------------------------------------
            # OXIDISING AGENTS
            # -------------------------------------------------

            Reagent(
                "Acidified potassium dichromate",
                "K2Cr2O7/H+",
                "Oxidising agent",
                "Primary alcohol → aldehyde → carboxylic acid; "
                "secondary alcohol → ketone",
                "Usually aqueous acidic medium"
            ),

            Reagent(
                "Potassium permanganate",
                "KMnO4",
                "Oxidising agent",
                "Oxidation of alcohols, alkenes and alkynes; "
                "strong oxidative cleavage",
                "Acidic, neutral or alkaline medium"
            ),

            Reagent(
                "PCC",
                "C5H5NHCrO3Cl",
                "Oxidising agent",
                "Primary alcohol → aldehyde; "
                "secondary alcohol → ketone",
                "Anhydrous medium"
            ),

            Reagent(
                "PDC",
                "PDC",
                "Oxidising agent",
                "Alcohol oxidation",
                "Usually non-aqueous"
            ),

            Reagent(
                "Jones reagent",
                "CrO3/H2SO4",
                "Oxidising agent",
                "Strong oxidation of primary and secondary alcohols",
                "Aqueous acidic medium"
            ),

            Reagent(
                "Ozone",
                "O3",
                "Oxidising agent",
                "Ozonolysis of alkenes and alkynes",
                "Followed by reductive or oxidative workup"
            ),

            Reagent(
                "Hot alkaline KMnO4",
                "KMnO4/OH−/heat",
                "Oxidising agent",
                "Strong oxidation and cleavage of unsaturation",
                "Heat"
            ),

            Reagent(
                "Cold dilute alkaline KMnO4",
                "KMnO4/OH−",
                "Oxidising agent",
                "Alkene → vicinal diol",
                "Baeyer test"
            ),

            Reagent(
                "Chromic acid",
                "H2CrO4",
                "Oxidising agent",
                "Oxidation of alcohols",
                "Strong oxidising conditions"
            ),

            Reagent(
                "Tollens reagent",
                "[Ag(NH3)2]+",
                "Oxidation test",
                "Aldehyde → carboxylate; silver mirror test",
                "Ammoniacal silver nitrate"
            ),

            Reagent(
                "Fehling solution",
                "Cu2+/alkaline",
                "Oxidation test",
                "Aliphatic aldehyde → carboxylate; brick-red Cu2O",
                "Alkaline heating"
            ),

            Reagent(
                "Benedict solution",
                "Cu2+/alkaline",
                "Oxidation test",
                "Reducing aldehydes and reducing sugars",
                "Heating"
            ),

            # -------------------------------------------------
            # REDUCING AGENTS
            # -------------------------------------------------

            Reagent(
                "Catalytic hydrogenation",
                "H2/Ni",
                "Reducing agent",
                "Alkene → alkane; alkyne → alkane",
                "Hydrogen gas and nickel catalyst"
            ),

            Reagent(
                "Hydrogenation over palladium",
                "H2/Pd",
                "Reducing agent",
                "Reduction of unsaturated compounds",
                "Catalytic hydrogenation"
            ),

            Reagent(
                "Hydrogenation over platinum",
                "H2/Pt",
                "Reducing agent",
                "Reduction of unsaturated compounds",
                "Catalytic hydrogenation"
            ),

            Reagent(
                "Sodium borohydride",
                "NaBH4",
                "Reducing agent",
                "Aldehyde → primary alcohol; "
                "ketone → secondary alcohol",
                "Mild reducing agent"
            ),

            Reagent(
                "Lithium aluminium hydride",
                "LiAlH4",
                "Reducing agent",
                "Reduction of aldehydes, ketones, acids, esters, "
                "amides and nitriles",
                "Dry ether, followed by hydrolysis"
            ),

            Reagent(
                "DIBAL-H",
                "DIBAL-H",
                "Reducing agent",
                "Ester or nitrile → aldehyde under controlled conditions",
                "Low temperature"
            ),

            Reagent(
                "Sodium in alcohol",
                "Na/C2H5OH",
                "Reducing agent",
                "Reduction of some unsaturated and carbonyl systems",
                "Proton source present"
            ),

            Reagent(
                "Clemmensen reduction",
                "Zn-Hg/HCl",
                "Reducing agent",
                "Aldehyde or ketone → hydrocarbon",
                "Strongly acidic medium"
            ),

            Reagent(
                "Wolff-Kishner reduction",
                "NH2NH2/KOH",
                "Reducing agent",
                "Aldehyde or ketone → hydrocarbon",
                "Strongly basic, high temperature"
            ),

            Reagent(
                "Zn/HCl",
                "Zn/HCl",
                "Reducing agent",
                "Reduction of several functional groups",
                "Acidic medium"
            ),

            # -------------------------------------------------
            # HALOGENATION
            # -------------------------------------------------

            Reagent(
                "Bromine",
                "Br2",
                "Addition / halogenation",
                "Alkene or alkyne bromination; unsaturation test",
                "Usually inert solvent"
            ),

            Reagent(
                "Chlorine",
                "Cl2",
                "Addition / halogenation",
                "Halogenation of alkenes and alkanes",
                "Light or suitable solvent"
            ),

            Reagent(
                "NBS",
                "N-bromosuccinimide",
                "Halogenation",
                "Allylic and benzylic bromination",
                "Often light or radical initiator"
            ),

            Reagent(
                "Cl2/hν",
                "Cl2/hν",
                "Free-radical substitution",
                "Alkane chlorination",
                "Ultraviolet light"
            ),

            Reagent(
                "Br2/hν",
                "Br2/hν",
                "Free-radical substitution",
                "Alkane bromination",
                "Ultraviolet light"
            ),

            # -------------------------------------------------
            # SUBSTITUTION AND ELIMINATION
            # -------------------------------------------------

            Reagent(
                "Aqueous KOH",
                "KOH(aq)",
                "Nucleophilic substitution",
                "Haloalkane → alcohol",
                "Aqueous medium"
            ),

            Reagent(
                "Aqueous NaOH",
                "NaOH(aq)",
                "Nucleophilic substitution",
                "Haloalkane → alcohol",
                "Aqueous medium"
            ),

            Reagent(
                "Alcoholic KOH",
                "KOH(alc)",
                "Elimination",
                "Haloalkane → alkene",
                "Heat"
            ),

            Reagent(
                "Alcoholic NaOH",
                "NaOH(alc)",
                "Elimination",
                "Haloalkane → alkene",
                "Heat"
            ),

            Reagent(
                "Sodium ethoxide",
                "C2H5ONa",
                "Elimination / substitution",
                "E2 elimination or SN2 substitution",
                "Usually alcoholic medium"
            ),

            Reagent(
                "Sodium amide",
                "NaNH2",
                "Strong base",
                "Alkyne formation by double dehydrohalogenation",
                "Liquid ammonia or suitable solvent"
            ),

            # -------------------------------------------------
            # ADDITION
            # -------------------------------------------------

            Reagent(
                "Hydrogen halide",
                "HCl",
                "Electrophilic addition",
                "Alkene → haloalkane",
                "Markovnikov addition generally"
            ),

            Reagent(
                "Hydrogen bromide",
                "HBr",
                "Electrophilic addition",
                "Alkene → bromoalkane",
                "Peroxide effect possible"
            ),

            Reagent(
                "Hydrogen iodide",
                "HI",
                "Electrophilic addition",
                "Alkene → iodoalkane",
                ""
            ),

            Reagent(
                "Water / acid",
                "H2O/H+",
                "Hydration",
                "Alkene → alcohol",
                "Acid-catalysed"
            ),

            Reagent(
                "H2SO4 followed by water",
                "H2SO4/H2O",
                "Hydration",
                "Alkene → alcohol",
                "Acid-catalysed hydration"
            ),

            Reagent(
                "Hydroboration oxidation",
                "BH3; H2O2/OH−",
                "Hydration",
                "Alkene → anti-Markovnikov alcohol",
                "Syn addition"
            ),

            # -------------------------------------------------
            # ALCOHOL / ETHER / ESTER CHEMISTRY
            # -------------------------------------------------

            Reagent(
                "Concentrated sulphuric acid",
                "conc. H2SO4",
                "Dehydrating / acidic reagent",
                "Alcohol dehydration; esterification catalyst",
                "Heat where required"
            ),

            Reagent(
                "Sodium metal",
                "Na",
                "Acid-base / reducing agent",
                "Alcohol → sodium alkoxide + H2",
                "Dry conditions preferred"
            ),

            Reagent(
                "Sodium bicarbonate",
                "NaHCO3",
                "Acid-base test",
                "Carboxylic acid → CO2 effervescence",
                "Aqueous medium"
            ),

            Reagent(
                "Acetyl chloride",
                "CH3COCl",
                "Acylation",
                "Alcohol or amine acetylation",
                "Usually base present"
            ),

            Reagent(
                "Acetic anhydride",
                "(CH3CO)2O",
                "Acylation",
                "Alcohol or amine acetylation",
                ""
            ),

            Reagent(
                "Ethanol / acid",
                "C2H5OH/H+",
                "Esterification",
                "Carboxylic acid → ester",
                "Fischer esterification"
            ),

            # -------------------------------------------------
            # AROMATIC CHEMISTRY
            # -------------------------------------------------

            Reagent(
                "Conc. HNO3/H2SO4",
                "HNO3/H2SO4",
                "Electrophilic aromatic substitution",
                "Benzene nitration",
                "Nitrating mixture"
            ),

            Reagent(
                "Cl2/FeCl3",
                "Cl2/FeCl3",
                "Electrophilic aromatic substitution",
                "Benzene chlorination",
                ""
            ),

            Reagent(
                "Br2/FeBr3",
                "Br2/FeBr3",
                "Electrophilic aromatic substitution",
                "Benzene bromination",
                ""
            ),

            Reagent(
                "CH3Cl/AlCl3",
                "CH3Cl/AlCl3",
                "Friedel-Crafts alkylation",
                "Alkylation of benzene",
                ""
            ),

            Reagent(
                "RCOCl/AlCl3",
                "RCOCl/AlCl3",
                "Friedel-Crafts acylation",
                "Acylation of benzene",
                ""
            ),

            # -------------------------------------------------
            # TESTS
            # -------------------------------------------------

            Reagent(
                "Lucas reagent",
                "conc. HCl/ZnCl2",
                "Functional-group test",
                "Distinguishes primary, secondary and tertiary alcohols",
                ""
            ),

            Reagent(
                "Iodoform reagent",
                "I2/NaOH",
                "Functional-group test",
                "Test for CH3CO− or CH3CH(OH)− group",
                "Yellow CHI3 precipitate"
            ),

            Reagent(
                "Neutral FeCl3",
                "FeCl3",
                "Functional-group test",
                "Phenol test",
                "Coloured complex"
            ),

            Reagent(
                "Sodium nitroprusside",
                "Na2[Fe(CN)5NO]",
                "Functional-group test",
                "Certain ketone tests",
                ""
            ),
        ]

    def all(self):
        return self.reagents

    def by_category(self, category):
        return [
            reagent for reagent in self.reagents
            if reagent.category.lower() == category.lower()
        ]

    def categories(self):
        return sorted(set(reagent.category for reagent in self.reagents))

    def search(self, text):
        text = text.lower()

        return [
            reagent for reagent in self.reagents
            if text in reagent.name.lower()
            or text in reagent.formula.lower()
            or text in reagent.category.lower()
        ]

    def display_all(self):
        print("\nAVAILABLE REAGENTS")
        print("=" * 80)

        for index, reagent in enumerate(self.reagents, start=1):
            print(
                f"{index:02d}. {reagent.name:<32} "
                f"{reagent.formula:<24} "
                f"[{reagent.category}]"
            )

    def display_categories(self):
        print("\nREAGENT CATEGORIES")
        print("=" * 40)

        for index, category in enumerate(self.categories(), start=1):
            print(f"{index:02d}. {category}")

    def get_by_number(self, number):
        if 1 <= number <= len(self.reagents):
            return self.reagents[number - 1]

        return None