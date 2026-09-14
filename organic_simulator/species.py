"""
Chemical species and simplified acid/base dissociation.
"""


class Species:

    def __init__(
        self,
        formula,
        charge=0,
        phase="unknown"
    ):
        self.formula = formula
        self.charge = charge
        self.phase = phase

    def __repr__(self):
        if self.charge == 0:
            charge = ""
        elif self.charge == 1:
            charge = "+"
        elif self.charge == -1:
            charge = "-"
        else:
            charge = f"{self.charge:+d}"

        return f"{self.formula}{charge}"


class SpeciesParser:

    ACIDS = {
        "HCl": ["H+", "Cl-"],
        "HBr": ["H+", "Br-"],
        "HI": ["H+", "I-"],
        "HNO3": ["H+", "NO3-"],
        "H2SO4": ["H+", "HSO4-"],
        "HSO4": ["H+", "SO4-2"],
        "HClO4": ["H+", "ClO4-"],
        "HClO3": ["H+", "ClO3-"],
        "CH3COOH": ["H+", "CH3COO-"],
        "HCOOH": ["H+", "HCOO-"],
        "H2CO3": ["H+", "HCO3-"],
        "NH4+": ["H+", "NH3"],
    }

    BASES = {
        "NaOH": ["Na+", "OH-"],
        "KOH": ["K+", "OH-"],
        "LiOH": ["Li+", "OH-"],
        "Ca(OH)2": ["Ca+2", "OH-", "OH-"],
        "Ba(OH)2": ["Ba+2", "OH-", "OH-"],
        "NH3": ["NH4+", "OH-"],
    }

    SALTS = {
        "NaCl": ["Na+", "Cl-"],
        "KCl": ["K+", "Cl-"],
        "NaBr": ["Na+", "Br-"],
        "KBr": ["K+", "Br-"],
        "NaNO3": ["Na+", "NO3-"],
        "KNO3": ["K+", "NO3-"],
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