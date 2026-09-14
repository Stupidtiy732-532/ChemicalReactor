"""
Chemical species handling.

This is a simplified acid/base dissociation system.
It does not calculate equilibrium constants yet.
"""


class Species:

    def __init__(
        self,
        formula,
        charge=0,
        phase="unknown"
    ):
        self.formula = formula.strip()
        self.charge = charge
        self.phase = phase

    def __str__(self):
        if self.charge == 0:
            return self.formula

        if self.charge == 1:
            return self.formula + "+"

        if self.charge == -1:
            return self.formula + "-"

        return f"{self.formula}{self.charge:+d}"


class SpeciesParser:

    # Strong acids and common weak acids.
    ACIDS = {
        "HCl": ["H+", "Cl-"],
        "HBr": ["H+", "Br-"],
        "HI": ["H+", "I-"],
        "HNO3": ["H+", "NO3-"],
        "HClO4": ["H+", "ClO4-"],
        "HClO3": ["H+", "ClO3-"],
        "H2SO4": ["H+", "HSO4-"],
        "HSO4-": ["H+", "SO4-2"],
        "H2CO3": ["H+", "HCO3-"],
        "HCO3-": ["H+", "CO3-2"],
        "CH3COOH": ["H+", "CH3COO-"],
        "HCOOH": ["H+", "HCOO-"],
        "C6H5COOH": ["H+", "C6H5COO-"],
        "NH4+": ["H+", "NH3"],
    }

    BASES = {
        "NaOH": ["Na+", "OH-"],
        "KOH": ["K+", "OH-"],
        "LiOH": ["Li+", "OH-"],
        "NH3": ["NH4+", "OH-"],
        "NH2-": ["NH3"],
    }

    SALTS = {
        "NaCl": ["Na+", "Cl-"],
        "KCl": ["K+", "Cl-"],
        "NaBr": ["Na+", "Br-"],
        "KBr": ["K+", "Br-"],
        "NaI": ["Na+", "I-"],
        "KI": ["K+", "I-"],
        "NaNO3": ["Na+", "NO3-"],
        "KNO3": ["K+", "NO3-"],
        "NaHCO3": ["Na+", "HCO3-"],
        "Na2CO3": ["Na+", "Na+", "CO3-2"],
        "Na2SO4": ["Na+", "Na+", "SO4-2"],
        "K2SO4": ["K+", "K+", "SO4-2"],
    }

    def dissociate(self, formula):
        formula = formula.strip()

        if formula in self.ACIDS:
            return self.ACIDS[formula]

        if formula in self.BASES:
            return self.BASES[formula]

        if formula in self.SALTS:
            return self.SALTS[formula]

        return [formula]

    def classify(self, formula):
        formula = formula.strip()

        if formula in self.ACIDS:
            return "acid"

        if formula in self.BASES:
            return "base"

        if formula in self.SALTS:
            return "salt"

        if formula in {"H+", "H3O+"}:
            return "acidic ion"

        if formula in {"OH-", "NH2-"}:
            return "basic ion"

        return "unknown"

    def explain(self, formula):
        species_type = self.classify(formula)
        products = self.dissociate(formula)

        print("\nSPECIES ANALYSIS")
        print("=" * 60)
        print(f"Input      : {formula}")
        print(f"Type       : {species_type}")
        print(f"Dissociates: {' + '.join(products)}")