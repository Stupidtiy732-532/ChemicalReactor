from dataclasses import dataclass, field


@dataclass(frozen=True)
class ReagentProfile:
    name: str
    formula: str
    aliases: tuple[str, ...]
    category: tuple[str, ...]
    strength: str | None = None
    medium: str | None = None
    conditions: tuple[str, ...] = ()
    acts_on: tuple[str, ...] = ()
    products: str | None = None
    notes: str = ""


def profile(
    name,
    formula,
    aliases,
    category,
    strength=None,
    medium=None,
    conditions=(),
    acts_on=(),
    products=None,
    notes="",
):
    return ReagentProfile(
        name=name,
        formula=formula,
        aliases=tuple(aliases),
        category=tuple(category),
        strength=strength,
        medium=medium,
        conditions=tuple(conditions),
        acts_on=tuple(acts_on),
        products=products,
        notes=notes,
    )


NCERT_REAGENTS = [

    # =========================================================
    # OXIDISING AGENTS
    # =========================================================

    profile(
        "Potassium permanganate",
        "KMnO4",
        ["kmno4", "potassium permanganate", "permanganate"],
        ["oxidant"],
        "condition-dependent",
        conditions=[
            "cold dilute alkaline",
            "hot concentrated",
            "acidic",
        ],
        acts_on=[
            "alkene",
            "alkyne",
            "primary alcohol",
            "secondary alcohol",
            "aldehyde",
        ],
        notes="Products depend strongly on temperature, concentration and medium.",
    ),

    profile(
        "Potassium dichromate",
        "K2Cr2O7",
        ["k2cr2o7", "potassium dichromate", "dichromate"],
        ["oxidant"],
        "strong",
        medium="acidic",
        acts_on=[
            "primary alcohol",
            "secondary alcohol",
            "aldehyde",
        ],
        notes="Usually used as acidified potassium dichromate.",
    ),

    profile(
        "Chromic acid",
        "H2CrO4",
        ["h2cro4", "chromic acid", "jones reagent", "jones"],
        ["oxidant"],
        "strong",
        medium="acidic",
        acts_on=[
            "primary alcohol",
            "secondary alcohol",
            "aldehyde",
        ],
    ),

    profile(
        "Pyridinium chlorochromate",
        "PCC",
        ["pcc", "pyridinium chlorochromate"],
        ["oxidant"],
        "mild",
        acts_on=[
            "primary alcohol",
            "secondary alcohol",
        ],
        products="primary alcohol → aldehyde; secondary alcohol → ketone",
    ),

    profile(
        "Pyridinium dichromate",
        "PDC",
        ["pdc", "pyridinium dichromate"],
        ["oxidant"],
        "mild",
        acts_on=[
            "primary alcohol",
            "secondary alcohol",
        ],
    ),

    profile(
        "Ozone",
        "O3",
        ["o3", "ozone"],
        ["oxidant", "cleaving agent"],
        "strong",
        acts_on=[
            "alkene",
            "alkyne",
        ],
        products="ozonolysis products",
    ),

    profile(
        "Manganese dioxide",
        "MnO2",
        ["mno2", "manganese dioxide"],
        ["oxidant"],
        "selective",
        acts_on=[
            "allylic alcohol",
            "benzylic alcohol",
        ],
    ),

    profile(
        "Copper(II) oxide",
        "CuO",
        ["cuo", "copper oxide"],
        ["oxidant"],
        "mild",
        acts_on=[
            "primary alcohol",
            "secondary alcohol",
        ],
    ),

    profile(
        "Nitric acid",
        "HNO3",
        ["hno3", "nitric acid"],
        ["oxidant", "nitrating agent", "acid"],
        acts_on=[
            "alcohol",
            "aromatic compound",
        ],
    ),

    # =========================================================
    # REDUCING AGENTS
    # =========================================================

    profile(
        "Hydrogen over nickel",
        "H2/Ni",
        ["h2/ni", "h2ni", "hydrogen nickel"],
        ["reductant", "hydrogenation catalyst"],
        "catalytic",
        acts_on=[
            "alkene",
            "alkyne",
            "carbonyl",
        ],
    ),

    profile(
        "Hydrogen over palladium",
        "H2/Pd",
        ["h2/pd", "h2pd", "hydrogen palladium"],
        ["reductant", "hydrogenation catalyst"],
        "catalytic",
        acts_on=[
            "alkene",
            "alkyne",
            "carbonyl",
        ],
    ),

    profile(
        "Hydrogen over platinum",
        "H2/Pt",
        ["h2/pt", "h2pt", "hydrogen platinum"],
        ["reductant", "hydrogenation catalyst"],
        "catalytic",
        acts_on=[
            "alkene",
            "alkyne",
        ],
    ),

    profile(
        "Sodium borohydride",
        "NaBH4",
        ["nabh4", "sodium borohydride"],
        ["reductant"],
        "mild",
        acts_on=[
            "aldehyde",
            "ketone",
        ],
        products="aldehyde/ketone → alcohol",
    ),

    profile(
        "Lithium aluminium hydride",
        "LiAlH4",
        ["lialh4", "lithium aluminium hydride", "lithium aluminum hydride"],
        ["reductant"],
        "strong",
        acts_on=[
            "aldehyde",
            "ketone",
            "carboxylic acid",
            "ester",
            "amide",
            "nitrile",
        ],
    ),

    profile(
        "Lindlar catalyst",
        "H2/Lindlar",
        ["lindlar", "lindlar catalyst", "h2/lindlar"],
        ["reductant", "selective hydrogenation"],
        "selective",
        acts_on=["alkyne"],
        products="alkyne → cis-alkene",
    ),

    profile(
        "Clemmensen reduction",
        "Zn-Hg/HCl",
        ["zn-hg", "znhg", "zn/hg", "clemmensen"],
        ["reductant"],
        "strong",
        medium="acidic",
        acts_on=["aldehyde", "ketone"],
        products="carbonyl → hydrocarbon",
    ),

    profile(
        "Wolff–Kishner reduction",
        "NH2NH2/KOH",
        ["nh2nh2", "hydrazine", "wolff kishner", "wolff-kishner"],
        ["reductant"],
        "strong",
        medium="basic",
        acts_on=["aldehyde", "ketone"],
        products="carbonyl → hydrocarbon",
    ),

    profile(
        "Tin and hydrochloric acid",
        "Sn/HCl",
        ["sn/hcl", "snhcl", "tin hcl"],
        ["reductant"],
        "strong",
        medium="acidic",
        acts_on=["nitro compound"],
        products="nitro compound → amine",
    ),

    profile(
        "Iron and hydrochloric acid",
        "Fe/HCl",
        ["fe/hcl", "fehcl", "iron hcl"],
        ["reductant"],
        "strong",
        medium="acidic",
        acts_on=["nitro compound"],
        products="nitro compound → amine",
    ),

    # =========================================================
    # HALOGENATING AGENTS
    # =========================================================

    profile(
        "Chlorine",
        "Cl2",
        ["cl2", "chlorine"],
        ["halogenating agent", "electrophile", "radical reagent"],
        acts_on=["alkene", "alkane", "aromatic compound"],
    ),

    profile(
        "Bromine",
        "Br2",
        ["br2", "bromine"],
        ["halogenating agent", "electrophile", "radical reagent"],
        acts_on=["alkene", "alkane", "aromatic compound"],
    ),

    profile(
        "Hydrogen chloride",
        "HCl",
        ["hcl", "hydrogen chloride"],
        ["acid", "halogenating agent"],
        acts_on=["alkene", "alcohol"],
    ),

    profile(
        "Hydrogen bromide",
        "HBr",
        ["hbr", "hydrogen bromide"],
        ["acid", "halogenating agent"],
        acts_on=["alkene", "alcohol"],
    ),

    profile(
        "Hydrogen iodide",
        "HI",
        ["hi", "hydrogen iodide"],
        ["acid", "halogenating agent", "reductant"],
        acts_on=["alkene", "alcohol"],
    ),

    profile(
        "Thionyl chloride",
        "SOCl2",
        ["socl2", "thionyl chloride"],
        ["halogenating agent"],
        acts_on=["alcohol", "carboxylic acid"],
        products="alcohol → chloroalkane; acid → acyl chloride",
    ),

    profile(
        "Phosphorus pentachloride",
        "PCl5",
        ["pcl5", "phosphorus pentachloride"],
        ["halogenating agent"],
        acts_on=["alcohol", "carboxylic acid"],
    ),

    profile(
        "Phosphorus tribromide",
        "PBr3",
        ["pbr3", "phosphorus tribromide"],
        ["halogenating agent"],
        acts_on=["alcohol"],
        products="alcohol → bromoalkane",
    ),

    profile(
        "N-Bromosuccinimide",
        "NBS",
        ["nbs", "n-bromosuccinimide"],
        ["halogenating agent", "radical reagent"],
        acts_on=["allylic compound", "benzylic compound"],
    ),

    # =========================================================
    # DEHYDRATING AGENTS
    # =========================================================

    profile(
        "Concentrated sulphuric acid",
        "conc. H2SO4",
        ["conc h2so4", "concentrated sulphuric acid", "h2so4"],
        ["acid", "dehydrating agent"],
        "strong",
        acts_on=["alcohol", "alkene", "carboxylic acid"],
    ),

    profile(
        "Phosphoric acid",
        "H3PO4",
        ["h3po4", "phosphoric acid"],
        ["acid", "dehydrating agent"],
        acts_on=["alcohol"],
    ),

    profile(
        "Aluminium oxide",
        "Al2O3",
        ["al2o3", "aluminium oxide", "alumina"],
        ["dehydrating agent"],
        acts_on=["alcohol"],
        conditions=["heat"],
        products="alcohol → alkene",
    ),

    # =========================================================
    # NUCLEOPHILES / SUBSTITUTION REAGENTS
    # =========================================================

    profile(
        "Potassium cyanide",
        "KCN",
        ["kcn", "potassium cyanide"],
        ["nucleophile", "substitution reagent"],
        acts_on=["haloalkane"],
        products="R-X → R-CN",
    ),

    profile(
        "Silver cyanide",
        "AgCN",
        ["agcn", "silver cyanide"],
        ["nucleophile", "substitution reagent"],
        acts_on=["haloalkane"],
        products="R-X → R-NC",
    ),

    profile(
        "Ammonia",
        "NH3",
        ["nh3", "ammonia"],
        ["nucleophile", "base"],
        acts_on=["haloalkane", "acyl derivative"],
    ),

    profile(
        "Potassium hydroxide aqueous",
        "KOH(aq)",
        ["koh aq", "aqueous koh", "koh water"],
        ["base", "nucleophile"],
        acts_on=["haloalkane"],
        products="haloalkane → alcohol",
    ),

    profile(
        "Potassium hydroxide alcoholic",
        "KOH(alc.)",
        ["koh alc", "alcoholic koh", "koh ethanol"],
        ["base", "elimination reagent"],
        acts_on=["haloalkane"],
        products="haloalkane → alkene",
    ),

    profile(
        "Sodium ethoxide",
        "C2H5ONa",
        ["c2h5ona", "sodium ethoxide"],
        ["base", "nucleophile"],
        acts_on=["haloalkane"],
    ),

    profile(
        "Sodium nitrite",
        "NaNO2",
        ["nano2", "sodium nitrite"],
        ["nucleophile", "nitrosating agent"],
        acts_on=["haloalkane", "amine"],
    ),

    profile(
        "Silver nitrite",
        "AgNO2",
        ["agno2", "silver nitrite"],
        ["nucleophile", "substitution reagent"],
        acts_on=["haloalkane"],
        products="R-X → R-NO2",
    ),

    # =========================================================
    # CARBONYL REAGENTS
    # =========================================================

    profile(
        "Hydrogen cyanide",
        "HCN",
        ["hcn", "hydrogen cyanide"],
        ["nucleophile", "addition reagent"],
        acts_on=["aldehyde", "ketone"],
        products="carbonyl → cyanohydrin",
    ),

    profile(
        "Sodium bisulphite",
        "NaHSO3",
        ["nahso3", "sodium bisulphite"],
        ["addition reagent"],
        acts_on=["aldehyde", "ketone"],
    ),

    profile(
        "Hydroxylamine",
        "NH2OH",
        ["nh2oh", "hydroxylamine"],
        ["condensation reagent"],
        acts_on=["aldehyde", "ketone"],
        products="carbonyl → oxime",
    ),

    profile(
        "2,4-Dinitrophenylhydrazine",
        "2,4-DNP",
        ["2,4-dnp", "dnph", "2,4-dinitrophenylhydrazine"],
        ["condensation reagent", "test reagent"],
        acts_on=["aldehyde", "ketone"],
        products="carbonyl → 2,4-DNP derivative",
    ),

    profile(
        "Tollens reagent",
        "[Ag(NH3)2]+",
        ["tollens", "tollens reagent", "ammoniacal silver nitrate"],
        ["oxidant", "test reagent"],
        acts_on=["aldehyde"],
        products="aldehyde → carboxylate; silver mirror",
    ),

    profile(
        "Fehling solution",
        "Cu2+ complex",
        ["fehling", "fehling solution"],
        ["oxidant", "test reagent"],
        acts_on=["aldehyde"],
        products="aldehyde → carboxylate; Cu2O precipitate",
    ),

    profile(
        "Iodine and sodium hydroxide",
        "I2/NaOH",
        ["i2/naoh", "iodoform reagent", "iodine naoh"],
        ["oxidant", "test reagent"],
        acts_on=["methyl ketone", "ethanol", "ethanal"],
        products="iodoform reaction",
    ),

    # =========================================================
    # AROMATIC REAGENTS
    # =========================================================

    profile(
        "Nitrating mixture",
        "conc. HNO3/conc. H2SO4",
        ["nitrating mixture", "hno3 h2so4", "conc hno3 h2so4"],
        ["electrophile generator", "aromatic substitution reagent"],
        acts_on=["benzene", "aromatic compound"],
        products="aromatic compound → nitro compound",
    ),

    profile(
        "Ferric chloride",
        "FeCl3",
        ["fecl3", "ferric chloride"],
        ["Lewis acid", "catalyst"],
        acts_on=["aromatic compound"],
    ),

    profile(
        "Ferric bromide",
        "FeBr3",
        ["febr3", "ferric bromide"],
        ["Lewis acid", "catalyst"],
        acts_on=["aromatic compound"],
    ),

    profile(
        "Aluminium chloride",
        "AlCl3",
        ["alcl3", "aluminium chloride", "aluminum chloride"],
        ["Lewis acid", "catalyst"],
        acts_on=["aromatic compound", "alkyl halide", "acyl chloride"],
    ),

    profile(
        "Sodium nitrite and hydrochloric acid",
        "NaNO2/HCl",
        ["nano2/hcl", "diazotisation mixture"],
        ["diazotising reagent"],
        medium="acidic",
        acts_on=["primary aromatic amine"],
        products="aromatic amine → diazonium salt",
    ),

    profile(
        "Hypophosphorous acid",
        "H3PO2",
        ["h3po2", "hypophosphorous acid"],
        ["reductant"],
        acts_on=["diazonium salt"],
        products="diazonium group → hydrogen",
    ),
]


def normalize(text: str) -> str:
    return (
        text.strip()
        .lower()
        .replace(" ", "")
        .replace("-", "")
        .replace("–", "")
        .replace("_", "")
        .replace(".", "")
    )


def find_reagent(text: str) -> ReagentProfile | None:
    key = normalize(text)

    for reagent in NCERT_REAGENTS:
        possible_names = (
            reagent.name,
            reagent.formula,
            *reagent.aliases,
        )

        for alias in possible_names:
            if normalize(alias) == key:
                return reagent

    return None


def search_reagents(category: str) -> list[ReagentProfile]:
    key = normalize(category)

    return [
        reagent
        for reagent in NCERT_REAGENTS
        if any(normalize(item) == key for item in reagent.category)
    ]


def reagent_summary(text: str) -> str:
    reagent = find_reagent(text)

    if reagent is None:
        return f"Unknown reagent: {text}"

    lines = [
        f"Name: {reagent.name}",
        f"Formula: {reagent.formula}",
        f"Category: {', '.join(reagent.category)}",
    ]

    if reagent.strength:
        lines.append(f"Strength: {reagent.strength}")

    if reagent.medium:
        lines.append(f"Medium: {reagent.medium}")

    if reagent.conditions:
        lines.append(f"Conditions: {', '.join(reagent.conditions)}")

    if reagent.acts_on:
        lines.append(f"Acts on: {', '.join(reagent.acts_on)}")

    if reagent.products:
        lines.append(f"Products: {reagent.products}")

    if reagent.notes:
        lines.append(f"Notes: {reagent.notes}")

    return "\n".join(lines)