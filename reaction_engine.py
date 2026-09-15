import re
from dataclasses import dataclass
from parser import parse_structure


@dataclass
class ReactionResult:
    reaction: str
    products: list[str]
    conditions: str = ""
    notes: str = ""

    def display(self):
        print("\n" + "=" * 70)
        print("REACTION RESULT".center(70))
        print("=" * 70)

        print(f"\nReaction: {self.reaction}")

        if self.conditions:
            print(f"Conditions: {self.conditions}")

        print("\nProducts:")
        for index, product in enumerate(self.products, 1):
            try:
                molecule = parse_structure(product)
                print(
                    f"  {index}. {product}"
                    f"    [{molecule.formula()}]"
                )
            except Exception:
                print(f"  {index}. {product}")

        if self.notes:
            print(f"\nNote: {self.notes}")


class ReactionEngine:

    def react(self, structure, reagent, conditions=""):
        structure = structure.replace(" ", "")
        reagent = reagent.strip().lower()
        conditions = conditions.strip().lower()

        # =====================================================
        # ALKENES
        # =====================================================

        if self.is_alkene(structure):

            if reagent in {"h2/ni", "h2/pt", "h2/pd", "h2"}:
                return ReactionResult(
                    f"{structure} + H2",
                    [self.hydrogenate_alkene(structure)],
                    "Ni / Pt / Pd",
                    "Catalytic hydrogenation"
                )

            if reagent in {"br2", "br2/ccl4", "bromine"}:
                return ReactionResult(
                    f"{structure} + Br2",
                    [self.add_halogen(structure, "Br")],
                    "CCl4",
                    "Vicinal dibromide formation"
                )

            if reagent in {"cl2", "chlorine"}:
                return ReactionResult(
                    f"{structure} + Cl2",
                    [self.add_halogen(structure, "Cl")],
                    "",
                    "Vicinal dichloride formation"
                )

            if reagent in {"h2o", "h2o/h+", "steam"}:
                return ReactionResult(
                    f"{structure} + H2O",
                    [self.hydrate_alkene(structure)],
                    "dilute H2SO4 / H+",
                    "Markovnikov hydration"
                )

            if reagent in {"hbr", "hbr/h+"}:
                return ReactionResult(
                    f"{structure} + HBr",
                    [self.add_hx_markovnikov(structure, "Br")],
                    "",
                    "Markovnikov addition"
                )

            if reagent in {"hbr/peroxide", "hbr/roor", "hbr/roo"}:
                return ReactionResult(
                    f"{structure} + HBr",
                    [self.add_hx_anti_markovnikov(structure, "Br")],
                    "Peroxide",
                    "Anti-Markovnikov radical addition"
                )

            if reagent in {"hcl", "hcl/h+"}:
                return ReactionResult(
                    f"{structure} + HCl",
                    [self.add_hx_markovnikov(structure, "Cl")],
                    "",
                    "Markovnikov addition"
                )

            if reagent in {"hi", "hi/h+"}:
                return ReactionResult(
                    f"{structure} + HI",
                    [self.add_hx_markovnikov(structure, "I")],
                    "",
                    "Markovnikov addition"
                )

            if reagent in {"cold kmno4", "kmno4/cold", "alkaline kmno4"}:
                return ReactionResult(
                    f"{structure} + KMnO4",
                    [self.dihydroxylate_alkene(structure)],
                    "Cold dilute alkaline KMnO4",
                    "Vicinal glycol formation"
                )

            if reagent in {"o3", "ozone"}:
                return ReactionResult(
                    f"{structure} + O3",
                    self.ozonolysis(structure),
                    "O3 followed by Zn/H2O",
                    "Reductive ozonolysis"
                )

        # =====================================================
        # ALKYNES
        # =====================================================

        if self.is_alkyne(structure):

            if reagent in {"h2/ni", "h2/pt", "h2/pd", "h2 excess"}:
                return ReactionResult(
                    f"{structure} + 2H2",
                    [self.hydrogenate_alkyne(structure)],
                    "Excess H2 / Ni",
                    "Complete reduction to alkane"
                )

            if reagent in {"h2/lin", "h2/lindlar", "lindlar"}:
                return ReactionResult(
                    f"{structure} + H2",
                    [self.partial_alkyne_reduction(structure)],
                    "Lindlar catalyst",
                    "Cis-alkene formation"
                )

            if reagent in {"na/nh3", "li/nh3", "na/liquid nh3"}:
                return ReactionResult(
                    f"{structure} + H2",
                    [self.trans_alkene(structure)],
                    "Na / liquid NH3",
                    "Trans-alkene formation"
                )

            if reagent in {"h2o/hgso4", "hgso4/h2so4", "hydration"}:
                return ReactionResult(
                    f"{structure} + H2O",
                    [self.hydrate_alkyne(structure)],
                    "HgSO4 / H2SO4",
                    "Enol-keto tautomerism"
                )

        # =====================================================
        # HALOALKANES
        # =====================================================

        if self.is_haloalkane(structure):

            if reagent in {"koh alc", "koh(alc)", "alc koh", "alcoholic koh"}:
                return ReactionResult(
                    f"{structure} + KOH",
                    [self.eliminate_haloalkane(structure)],
                    "Alcoholic KOH / heat",
                    "β-elimination"
                )

            if reagent in {"koh aq", "koh(aq)", "aqueous koh"}:
                return ReactionResult(
                    f"{structure} + KOH",
                    [self.substitute_halide(structure, "OH")],
                    "Aqueous KOH",
                    "Nucleophilic substitution"
                )

            if reagent in {"kcn", "nacn", "cn-"}:
                return ReactionResult(
                    f"{structure} + KCN",
                    [self.substitute_halide(structure, "CN")],
                    "Alcoholic KCN",
                    "Cyanide substitution"
                )

            if reagent in {"agno2"}:
                return ReactionResult(
                    f"{structure} + AgNO2",
                    [self.substitute_halide(structure, "NO2")],
                    "",
                    "Nitroalkane formation"
                )

            if reagent in {"kno2", "nano2"}:
                return ReactionResult(
                    f"{structure} + KNO2",
                    [self.substitute_halide(structure, "ONO")],
                    "",
                    "Alkyl nitrite formation"
                )

            if reagent in {"nh3", "alc nh3"}:
                return ReactionResult(
                    f"{structure} + NH3",
                    [self.substitute_halide(structure, "NH2")],
                    "Alcoholic NH3",
                    "Amine formation"
                )

        # =====================================================
        # ALCOHOLS
        # =====================================================

        if self.is_alcohol(structure):

            if reagent in {"pcc"}:
                return ReactionResult(
                    f"{structure} + PCC",
                    [self.oxidize_alcohol(structure, "pcc")],
                    "PCC / CH2Cl2",
                    "Primary alcohol → aldehyde; secondary alcohol → ketone"
                )

            if reagent in {"kmno4", "k2cr2o7", "acidified k2cr2o7"}:
                return ReactionResult(
                    f"{structure} + oxidant",
                    [self.oxidize_alcohol(structure, "strong")],
                    "Acidified KMnO4 / K2Cr2O7",
                    "Strong oxidation"
                )

            if reagent in {"conc h2so4", "h2so4 heat", "dehydration"}:
                return ReactionResult(
                    f"{structure} → alkene",
                    [self.dehydrate_alcohol(structure)],
                    "Conc. H2SO4 / heat",
                    "Dehydration"
                )

            if reagent in {"soCl2", "soc l2", "thionyl chloride"}:
                return ReactionResult(
                    f"{structure} + SOCl2",
                    [self.replace_oh(structure, "Cl")],
                    "SOCl2",
                    "Alcohol → chloroalkane"
                )

            if reagent in {"pcl5"}:
                return ReactionResult(
                    f"{structure} + PCl5",
                    [self.replace_oh(structure, "Cl")],
                    "PCl5",
                    "Alcohol → chloroalkane"
                )

            if reagent in {"pbr3"}:
                return ReactionResult(
                    f"{structure} + PBr3",
                    [self.replace_oh(structure, "Br")],
                    "PBr3",
                    "Alcohol → bromoalkane"
                )

            if reagent in {"na", "sodium"}:
                return ReactionResult(
                    f"{structure} + Na",
                    [self.replace_oh(structure, "ONa"), "H2"],
                    "",
                    "Sodium alkoxide formation"
                )

        # =====================================================
        # CARBOXYLIC ACIDS
        # =====================================================

        if self.is_carboxylic_acid(structure):

            if reagent in {"socl2", "thionyl chloride"}:
                return ReactionResult(
                    f"{structure} + SOCl2",
                    [structure.replace("COOH", "COCl")],
                    "SOCl2",
                    "Acid chloride formation"
                )

            if reagent in {"nahco3", "nahco3 aq", "sodium bicarbonate"}:
                return ReactionResult(
                    f"{structure} + NaHCO3",
                    ["CO2", "H2O", "Na salt"],
                    "Aqueous",
                    "CO2 effervescence"
                )

            if reagent in {"li alh4", "lialh4", "alh4"}:
                return ReactionResult(
                    f"{structure} + LiAlH4",
                    [self.reduce_acid(structure)],
                    "LiAlH4 / ether",
                    "Carboxylic acid → primary alcohol"
                )

        # =====================================================
        # ALDEHYDES / KETONES
        # =====================================================

        if self.is_carbonyl(structure):

            if reagent in {"hcn", "kcn/hcn"}:
                return ReactionResult(
                    f"{structure} + HCN",
                    [self.add_hcn(structure)],
                    "KCN / HCN",
                    "Cyanohydrin formation"
                )

            if reagent in {"nh2oh", "hydroxylamine"}:
                return ReactionResult(
                    f"{structure} + NH2OH",
                    [structure + "=NOH"],
                    "",
                    "Oxime formation"
                )

            if reagent in {"2,4-dnp", "2,4-dnph"}:
                return ReactionResult(
                    f"{structure} + 2,4-DNP",
                    [structure + "=NNH-C6H3(NO2)2"],
                    "",
                    "2,4-DNP derivative"
                )

            if reagent in {"nabh4", "nab h4"}:
                return ReactionResult(
                    f"{structure} + NaBH4",
                    [self.reduce_carbonyl(structure)],
                    "NaBH4",
                    "Carbonyl → alcohol"
                )

        return ReactionResult(
            f"{structure} + {reagent}",
            [],
            conditions,
            "No rule implemented for this combination yet."
        )

    # =========================================================
    # DETECTION
    # =========================================================

    @staticmethod
    def is_alkene(s):
        return "=" in s and "C=C" in s and "#" not in s

    @staticmethod
    def is_alkyne(s):
        return "#" in s or "C≡C" in s

    @staticmethod
    def is_haloalkane(s):
        return any(x in s for x in ("Cl", "Br", "I", "F")) and "C" in s

    @staticmethod
    def is_alcohol(s):
        return "OH" in s and "COOH" not in s

    @staticmethod
    def is_carboxylic_acid(s):
        return "COOH" in s

    @staticmethod
    def is_carbonyl(s):
        return "CHO" in s or "CO" in s

    # =========================================================
    # SIMPLE TRANSFORMATIONS
    # =========================================================

    @staticmethod
    def hydrogenate_alkene(s):
        return s.replace("=", "-")

    @staticmethod
    def hydrogenate_alkyne(s):
        return s.replace("#", "-").replace("≡", "-")

    @staticmethod
    def partial_alkyne_reduction(s):
        return s.replace("#", "=").replace("≡", "=")

    @staticmethod
    def trans_alkene(s):
        return s.replace("#", "=").replace("≡", "=")

    @staticmethod
    def add_halogen(s, halogen):
        return s.replace("=", f"(-{halogen})-")

    @staticmethod
    def hydrate_alkene(s):
        return s.replace("=C", "-C(OH)")

    @staticmethod
    def add_hx_markovnikov(s, halogen):
        return s.replace("=C", f"-C(H)(-{halogen})")

    @staticmethod
    def add_hx_anti_markovnikov(s, halogen):
        return s.replace("C=", f"C(-{halogen})-")

    @staticmethod
    def dihydroxylate_alkene(s):
        return s.replace("=C", "-C(OH)")

    @staticmethod
    def ozonolysis(s):
        parts = s.split("=")
        if len(parts) == 2:
            left = parts[0].replace("C", "CHO")
            right = parts[1].replace("C", "CHO")
            return [left, right]
        return ["aldehydes / ketones"]

    @staticmethod
    def eliminate_haloalkane(s):
        return re.sub(r"(CH2|CH|C)(Cl|Br|I|F)", "C=", s, count=1)

    @staticmethod
    def substitute_halide(s, group):
        return re.sub(r"(Cl|Br|I|F)", group, s, count=1)

    @staticmethod
    def replace_oh(s, group):
        return s.replace("OH", group, 1)

    @staticmethod
    def oxidize_alcohol(s, mode):
        if "CH2OH" in s:
            return s.replace("CH2OH", "CHO", 1)

        if "CHOH" in s:
            return s.replace("CHOH", "CO", 1)

        if mode == "strong" and "CH2OH" in s:
            return s.replace("CH2OH", "COOH", 1)

        return s

    @staticmethod
    def dehydrate_alcohol(s):
        return s.replace("CH2OH", "CH=", 1)

    @staticmethod
    def reduce_acid(s):
        return s.replace("COOH", "CH2OH", 1)

    @staticmethod
    def add_hcn(s):
        return s.replace("CO", "C(OH)(CN)", 1)

    @staticmethod
    def reduce_carbonyl(s):
        return s.replace("CHO", "CH2OH").replace("CO", "CHOH", 1)