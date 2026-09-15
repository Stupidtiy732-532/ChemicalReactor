from dataclasses import dataclass


@dataclass(frozen=True)
class ReactionProfile:
    id: str
    reagent: str
    conditions: tuple[str, ...]
    substrate: tuple[str, ...]
    reaction_type: str
    product_rule: str
    notes: str = ""


REACTION_PROFILES = [

    # =========================================================
    # ALKANES
    # =========================================================

    ReactionProfile(
        id="alkane_free_radical_halogenation",
        reagent="Cl2/Br2",
        conditions=("sunlight", "UV light", "heat"),
        substrate=("alkane",),
        reaction_type="free-radical substitution",
        product_rule="C-H → C-X",
        notes="Produces a mixture depending on selectivity.",
    ),

    # =========================================================
    # ALKENES
    # =========================================================

    ReactionProfile(
        id="alkene_hydrogenation",
        reagent="H2/Ni",
        conditions=("heat", "catalyst"),
        substrate=("alkene",),
        reaction_type="addition",
        product_rule="C=C → C-C",
    ),

    ReactionProfile(
        id="alkene_bromination",
        reagent="Br2",
        conditions=("CCl4", "inert solvent"),
        substrate=("alkene",),
        reaction_type="electrophilic addition",
        product_rule="C=C → vicinal dibromide",
    ),

    ReactionProfile(
        id="alkene_chlorination",
        reagent="Cl2",
        conditions=("inert solvent",),
        substrate=("alkene",),
        reaction_type="electrophilic addition",
        product_rule="C=C → vicinal dichloride",
    ),

    ReactionProfile(
        id="alkene_hydrohalogenation",
        reagent="HCl/HBr/HI",
        conditions=("normal conditions",),
        substrate=("alkene",),
        reaction_type="electrophilic addition",
        product_rule="C=C → alkyl halide",
        notes="Usually Markovnikov addition.",
    ),

    ReactionProfile(
        id="alkene_hbr_peroxide",
        reagent="HBr",
        conditions=("peroxide", "ROOR", "organic peroxide"),
        substrate=("alkene",),
        reaction_type="free-radical addition",
        product_rule="anti-Markovnikov bromoalkane",
        notes="Peroxide effect is applicable mainly to HBr.",
    ),

    ReactionProfile(
        id="alkene_hydration",
        reagent="H2O/H+",
        conditions=("acidic", "dilute H2SO4"),
        substrate=("alkene",),
        reaction_type="electrophilic addition",
        product_rule="C=C → alcohol",
        notes="Usually Markovnikov hydration.",
    ),

    ReactionProfile(
        id="alkene_cold_kmno4",
        reagent="KMnO4",
        conditions=("cold", "dilute", "alkaline"),
        substrate=("alkene",),
        reaction_type="oxidation",
        product_rule="alkene → vicinal diol",
        notes="Baeyer oxidation.",
    ),

    ReactionProfile(
        id="alkene_ozonolysis",
        reagent="O3",
        conditions=("ozonolysis", "Zn/H2O", "reductive workup"),
        substrate=("alkene",),
        reaction_type="oxidative cleavage",
        product_rule="C=C → aldehydes/ketones",
    ),

    ReactionProfile(
        id="alkene_hot_kmno4",
        reagent="KMnO4",
        conditions=("hot", "concentrated", "acidic", "strong oxidation"),
        substrate=("alkene",),
        reaction_type="oxidative cleavage",
        product_rule="C=C → carboxylic acids/ketones/CO2",
    ),

    # =========================================================
    # ALKYNES
    # =========================================================

    ReactionProfile(
        id="alkyne_complete_hydrogenation",
        reagent="H2/Ni",
        conditions=("excess hydrogen", "catalyst"),
        substrate=("alkyne",),
        reaction_type="complete hydrogenation",
        product_rule="C≡C → C-C",
    ),

    ReactionProfile(
        id="alkyne_partial_hydrogenation",
        reagent="H2/Lindlar",
        conditions=("controlled",),
        substrate=("alkyne",),
        reaction_type="partial hydrogenation",
        product_rule="alkyne → cis-alkene",
    ),

    ReactionProfile(
        id="alkyne_dissolving_metal_reduction",
        reagent="Na/NH3",
        conditions=("liquid ammonia",),
        substrate=("alkyne",),
        reaction_type="partial reduction",
        product_rule="alkyne → trans-alkene",
    ),

    ReactionProfile(
        id="alkyne_hydration",
        reagent="HgSO4/H2SO4/H2O",
        conditions=("acidic",),
        substrate=("alkyne",),
        reaction_type="addition followed by tautomerism",
        product_rule="alkyne → ketone",
        notes="Terminal alkynes generally give methyl ketones.",
    ),

    ReactionProfile(
        id="alkyne_hydroboration",
        reagent="BH3·THF/H2O2/OH-",
        conditions=("basic workup",),
        substrate=("alkyne",),
        reaction_type="anti-Markovnikov hydration",
        product_rule="terminal alkyne → aldehyde",
    ),

    # =========================================================
    # HALOALKANES
    # =========================================================

    ReactionProfile(
        id="haloalkane_aqueous_koh",
        reagent="KOH",
        conditions=("aqueous", "water"),
        substrate=("haloalkane",),
        reaction_type="nucleophilic substitution",
        product_rule="R-X → R-OH",
    ),

    ReactionProfile(
        id="haloalkane_alcoholic_koh",
        reagent="KOH",
        conditions=("alcoholic", "ethanolic", "heat"),
        substrate=("haloalkane",),
        reaction_type="β-elimination",
        product_rule="R-X → alkene",
    ),

    ReactionProfile(
        id="haloalkane_kcn",
        reagent="KCN",
        conditions=("alcoholic medium",),
        substrate=("haloalkane",),
        reaction_type="nucleophilic substitution",
        product_rule="R-X → R-CN",
    ),

    ReactionProfile(
        id="haloalkane_agcn",
        reagent="AgCN",
        conditions=("normal conditions",),
        substrate=("haloalkane",),
        reaction_type="nucleophilic substitution",
        product_rule="R-X → R-NC",
    ),

    ReactionProfile(
        id="haloalkane_nh3",
        reagent="NH3",
        conditions=("alcoholic ammonia", "excess ammonia"),
        substrate=("haloalkane",),
        reaction_type="nucleophilic substitution",
        product_rule="R-X → R-NH2",
    ),

    ReactionProfile(
        id="haloalkane_kno2",
        reagent="KNO2",
        conditions=("normal conditions",),
        substrate=("haloalkane",),
        reaction_type="nucleophilic substitution",
        product_rule="R-X → R-ONO",
    ),

    ReactionProfile(
        id="haloalkane_agno2",
        reagent="AgNO2",
        conditions=("normal conditions",),
        substrate=("haloalkane",),
        reaction_type="nucleophilic substitution",
        product_rule="R-X → R-NO2",
    ),

    # =========================================================
    # ALCOHOLS
    # =========================================================

    ReactionProfile(
        id="primary_alcohol_mild_oxidation",
        reagent="PCC",
        conditions=("anhydrous", "controlled oxidation"),
        substrate=("primary alcohol",),
        reaction_type="oxidation",
        product_rule="primary alcohol → aldehyde",
    ),

    ReactionProfile(
        id="primary_alcohol_strong_oxidation",
        reagent="KMnO4/K2Cr2O7",
        conditions=("acidic", "hot", "strong oxidation"),
        substrate=("primary alcohol",),
        reaction_type="oxidation",
        product_rule="primary alcohol → carboxylic acid",
    ),

    ReactionProfile(
        id="secondary_alcohol_oxidation",
        reagent="PCC/KMnO4/K2Cr2O7",
        conditions=("oxidation conditions",),
        substrate=("secondary alcohol",),
        reaction_type="oxidation",
        product_rule="secondary alcohol → ketone",
    ),

    ReactionProfile(
        id="tertiary_alcohol_oxidation",
        reagent="oxidant",
        conditions=("ordinary conditions",),
        substrate=("tertiary alcohol",),
        reaction_type="no ordinary reaction",
        product_rule="no reaction",
        notes="No α-hydrogen on the carbon bearing OH.",
    ),

    ReactionProfile(
        id="alcohol_dehydration",
        reagent="conc. H2SO4",
        conditions=("heat", "acidic"),
        substrate=("alcohol",),
        reaction_type="β-elimination",
        product_rule="alcohol → alkene",
    ),

    ReactionProfile(
        id="alcohol_sodium",
        reagent="Na",
        conditions=("normal conditions",),
        substrate=("alcohol",),
        reaction_type="acid-base reaction",
        product_rule="alcohol → sodium alkoxide + H2",
    ),

    ReactionProfile(
        id="alcohol_haloalkane",
        reagent="SOCl2/PCl5/PBr3",
        conditions=("normal conditions",),
        substrate=("alcohol",),
        reaction_type="substitution",
        product_rule="alcohol → haloalkane",
    ),

    # =========================================================
    # PHENOLS
    # =========================================================

    ReactionProfile(
        id="phenol_bromination",
        reagent="Br2",
        conditions=("bromine water",),
        substrate=("phenol",),
        reaction_type="electrophilic aromatic substitution",
        product_rule="phenol → 2,4,6-tribromophenol",
    ),

    ReactionProfile(
        id="phenol_nitration",
        reagent="HNO3",
        conditions=("dilute",),
        substrate=("phenol",),
        reaction_type="electrophilic aromatic substitution",
        product_rule="phenol → nitrophenol",
    ),

    ReactionProfile(
        id="phenol_kolbe",
        reagent="CO2/NaOH",
        conditions=("pressure", "acidification"),
        substrate=("sodium phenoxide",),
        reaction_type="carboxylation",
        product_rule="phenoxide → salicylic acid",
    ),

    ReactionProfile(
        id="phenol_reimer_tiemann",
        reagent="CHCl3/NaOH",
        conditions=("heat", "basic"),
        substrate=("phenol",),
        reaction_type="formylation",
        product_rule="phenol → salicylaldehyde",
    ),

    # =========================================================
    # ALDEHYDES AND KETONES
    # =========================================================

    ReactionProfile(
        id="carbonyl_nabh4",
        reagent="NaBH4",
        conditions=("alcoholic medium",),
        substrate=("aldehyde", "ketone"),
        reaction_type="reduction",
        product_rule="carbonyl → alcohol",
    ),

    ReactionProfile(
        id="carbonyl_lialh4",
        reagent="LiAlH4",
        conditions=("dry ether", "hydrolysis"),
        substrate=("aldehyde", "ketone"),
        reaction_type="reduction",
        product_rule="carbonyl → alcohol",
    ),

    ReactionProfile(
        id="carbonyl_hcn",
        reagent="HCN",
        conditions=("trace CN-",),
        substrate=("aldehyde", "ketone"),
        reaction_type="nucleophilic addition",
        product_rule="carbonyl → cyanohydrin",
    ),

    ReactionProfile(
        id="aldehyde_tollens",
        reagent="Tollens reagent",
        conditions=("ammoniacal",),
        substrate=("aldehyde",),
        reaction_type="oxidation test",
        product_rule="aldehyde → carboxylate + silver mirror",
    ),

    ReactionProfile(
        id="aldehyde_fehling",
        reagent="Fehling solution",
        conditions=("alkaline", "heating"),
        substrate=("aldehyde",),
        reaction_type="oxidation test",
        product_rule="aldehyde → carboxylate + Cu2O",
    ),

    ReactionProfile(
        id="carbonyl_2_4_dnp",
        reagent="2,4-DNP",
        conditions=("acidic",),
        substrate=("aldehyde", "ketone"),
        reaction_type="condensation",
        product_rule="carbonyl → 2,4-DNP hydrazone",
    ),

    ReactionProfile(
        id="carbonyl_hydroxylamine",
        reagent="NH2OH",
        conditions=("acidic",),
        substrate=("aldehyde", "ketone"),
        reaction_type="condensation",
        product_rule="carbonyl → oxime",
    ),

    ReactionProfile(
        id="carbonyl_clemmensen",
        reagent="Zn-Hg/HCl",
        conditions=("acidic",),
        substrate=("aldehyde", "ketone"),
        reaction_type="reduction",
        product_rule="carbonyl → CH2",
    ),

    ReactionProfile(
        id="carbonyl_wolff_kishner",
        reagent="NH2NH2/KOH",
        conditions=("strong base", "heat"),
        substrate=("aldehyde", "ketone"),
        reaction_type="reduction",
        product_rule="carbonyl → CH2",
    ),

    # =========================================================
    # CARBOXYLIC ACIDS AND DERIVATIVES
    # =========================================================

    ReactionProfile(
        id="acid_nahco3",
        reagent="NaHCO3",
        conditions=("aqueous",),
        substrate=("carboxylic acid",),
        reaction_type="acid-base reaction",
        product_rule="acid → carboxylate + CO2 + H2O",
    ),

    ReactionProfile(
        id="acid_socl2",
        reagent="SOCl2",
        conditions=("anhydrous",),
        substrate=("carboxylic acid",),
        reaction_type="substitution",
        product_rule="carboxylic acid → acyl chloride",
    ),

    ReactionProfile(
        id="acid_lialh4",
        reagent="LiAlH4",
        conditions=("dry ether", "hydrolysis"),
        substrate=("carboxylic acid",),
        reaction_type="reduction",
        product_rule="carboxylic acid → primary alcohol",
    ),

    ReactionProfile(
        id="ester_hydrolysis",
        reagent="NaOH/H2O",
        conditions=("heating", "saponification"),
        substrate=("ester",),
        reaction_type="hydrolysis",
        product_rule="ester → carboxylate + alcohol",
    ),

    # =========================================================
    # AMINES
    # =========================================================

    ReactionProfile(
        id="amine_diazotisation",
        reagent="NaNO2/HCl",
        conditions=("0–5 °C", "cold acidic"),
        substrate=("primary aromatic amine",),
        reaction_type="diazotisation",
        product_rule="amine → diazonium salt",
    ),

    ReactionProfile(
        id="diazonium_sandmeyer_cl",
        reagent="CuCl/HCl",
        conditions=("aqueous",),
        substrate=("diazonium salt",),
        reaction_type="Sandmeyer reaction",
        product_rule="diazonium group → Cl",
    ),

    ReactionProfile(
        id="diazonium_sandmeyer_br",
        reagent="CuBr/HBr",
        conditions=("aqueous",),
        substrate=("diazonium salt",),
        reaction_type="Sandmeyer reaction",
        product_rule="diazonium group → Br",
    ),

    ReactionProfile(
        id="diazonium_iodide",
        reagent="KI",
        conditions=("normal conditions",),
        substrate=("diazonium salt",),
        reaction_type="substitution",
        product_rule="diazonium group → I",
    ),

    ReactionProfile(
        id="diazonium_reduction",
        reagent="H3PO2",
        conditions=("aqueous",),
        substrate=("diazonium salt",),
        reaction_type="reduction",
        product_rule="diazonium group → H",
    ),
]


def normalize(text: str) -> str:
    return (
        text.lower()
        .strip()
        .replace(" ", "")
        .replace("-", "")
        .replace("–", "")
        .replace(".", "")
    )


def find_profiles(
    reagent: str | None = None,
    substrate: str | None = None,
    reaction_type: str | None = None,
) -> list[ReactionProfile]:

    results = []

    for item in REACTION_PROFILES:

        if reagent is not None:
            reagent_key = normalize(reagent)
            profile_key = normalize(item.reagent)

            if reagent_key not in profile_key and profile_key not in reagent_key:
                continue

        if substrate is not None:
            if substrate.lower() not in item.substrate:
                continue

        if reaction_type is not None:
            if normalize(reaction_type) != normalize(item.reaction_type):
                continue

        results.append(item)

    return results


def profile_summary(profile: ReactionProfile) -> str:
    lines = [
        f"ID: {profile.id}",
        f"Reagent: {profile.reagent}",
        f"Conditions: {', '.join(profile.conditions)}",
        f"Substrate: {', '.join(profile.substrate)}",
        f"Reaction type: {profile.reaction_type}",
        f"Product rule: {profile.product_rule}",
    ]

    if profile.notes:
        lines.append(f"Notes: {profile.notes}")

    return "\n".join(lines)