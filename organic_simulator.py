
"""
NCERT + JEE Organic Chemistry Reaction Simulator
Custom notation — NOT SMILES
Python 3.10+
"""

import re
from dataclasses import dataclass, field


# ============================================================
# CHEMICAL DEFINITIONS
# ============================================================

VALENCIES = {
    "H": {1},
    "D": {1},
    "T": {1},
    "C": {4},
    "N": {3, 5},
    "O": {2},
    "F": {1},
    "Cl": {1},
    "Br": {1},
    "I": {1},
    "S": {2, 4, 6},
    "P": {3, 5},
}

BOND_ORDERS = {
    "-": 1.0,
    "=": 2.0,
    "#": 3.0,
}


# ============================================================
# STRUCTURE DATA
# ============================================================

@dataclass
class Atom:
    element: str
    isotope: int | None = None
    charge: int = 0
    lone_pairs: int | None = None
    radicals: int = 0
    index: int = 0
    bonds: list = field(default_factory=list)

    def __str__(self):
        isotope = (
            f"[{self.isotope}{self.element}]"
            if self.isotope is not None
            else self.element
        )

        charge = ""

        if self.charge > 0:
            charge = f"^{self.charge}+"
        elif self.charge < 0:
            charge = f"^{abs(self.charge)}-"

        radical = "." * self.radicals

        lone_pairs = (
            ":" * self.lone_pairs
            if self.lone_pairs is not None
            else ""
        )

        return isotope + charge + radical + lone_pairs


@dataclass
class Bond:
    atom1: int
    atom2: int
    order: float = 1.0

    def symbol(self):
        if self.order == 1:
            return "-"
        if self.order == 2:
            return "="
        if self.order == 3:
            return "#"

        return f"-_{{{self.order:g}}}"


@dataclass
class Molecule:
    atoms: list[Atom] = field(default_factory=list)
    bonds: list[Bond] = field(default_factory=list)

    def add_atom(self, atom):
        atom.index = len(self.atoms)
        self.atoms.append(atom)
        return atom.index

    def add_bond(self, a, b, order=1.0):
        bond = Bond(a, b, order)
        self.bonds.append(bond)

        self.atoms[a].bonds.append(bond)
        self.atoms[b].bonds.append(bond)

    def formula(self):
        counts = {}

        for atom in self.atoms:
            counts[atom.element] = (
                counts.get(atom.element, 0) + 1
            )

        ordered = []

        if "C" in counts:
            ordered.append("C")

        if "H" in counts:
            ordered.append("H")

        ordered += sorted(
            element
            for element in counts
            if element not in ("C", "H")
        )

        return "".join(
            element + (str(count) if count > 1 else "")
            for element in ordered
            for count in [counts[element]]
        )

    def show(self):
        print("\nAtoms:")

        for atom in self.atoms:
            print(f"  {atom.index}: {atom}")

        print("\nBonds:")

        for bond in self.bonds:
            print(
                f"  {bond.atom1} "
                f"{bond.symbol()} "
                f"{bond.atom2}"
            )

        print(f"\nFormula: {self.formula()}")

    def validate(self):
        errors = []

        for atom in self.atoms:

            if atom.element not in VALENCIES:
                continue

            bond_order = sum(
                bond.order
                for bond in atom.bonds
            )

            allowed = VALENCIES[atom.element]

            # Explicit charge/radical/lone-pair chemistry
            # is not fully inferred in this MVP.
            if atom.charge != 0:
                continue

            if atom.radicals != 0:
                continue

            if atom.lone_pairs is not None:
                continue

            if not any(
                abs(bond_order - val) < 1e-9
                for val in allowed
            ):
                errors.append(
                    f"Atom {atom.index} ({atom.element}) "
                    f"has bond order {bond_order:g}; "
                    f"allowed valencies: {allowed}"
                )

        return errors


# ============================================================
# CUSTOM NOTATION PARSER
# ============================================================

class Parser:

    def __init__(self, text):
        self.text = text.replace(" ", "")
        self.pos = 0
        self.molecule = Molecule()
        self.previous_atom = None
        self.pending_bond = 1.0

    def peek(self):
        if self.pos >= len(self.text):
            return ""

        return self.text[self.pos]

    def consume(self):
        char = self.peek()
        self.pos += 1
        return char

    def parse(self):
        while self.pos < len(self.text):

            char = self.peek()

            if char == "(":
                self.consume()
                self.parse_branch()

            elif char == ")":
                raise ValueError("Unexpected ')'")

            elif char in "-=#":
                self.parse_bond()

            elif char == "^":
                self.parse_charge()

            elif char == ".":
                self.parse_radical()

            elif char == ":":
                self.parse_lone_pair()

            elif char == "[":
                self.parse_isotope_atom()

            elif char.isupper():
                self.parse_atom()

            else:
                raise ValueError(
                    f"Unexpected character '{char}' "
                    f"at position {self.pos}"
                )

        return self.molecule

    def parse_atom(self):

        match = re.match(
            r"(Cl|Br|[A-Z][a-z]?)",
            self.text[self.pos:]
        )

        if not match:
            raise ValueError("Invalid atom")

        element = match.group(1)
        self.pos += len(element)

        atom = Atom(element=element)

        index = self.molecule.add_atom(atom)

        if self.previous_atom is not None:
            self.molecule.add_bond(
                self.previous_atom,
                index,
                self.pending_bond
            )

        self.previous_atom = index
        self.pending_bond = 1.0

    def parse_isotope_atom(self):

        end = self.text.find("]", self.pos)

        if end == -1:
            raise ValueError("Unclosed isotope bracket")

        content = self.text[self.pos + 1:end]

        match = re.fullmatch(
            r"(\d+)(Cl|Br|[A-Z][a-z]?)",
            content
        )

        if not match:
            raise ValueError(
                f"Invalid isotope notation: [{content}]"
            )

        isotope = int(match.group(1))
        element = match.group(2)

        self.pos = end + 1

        atom = Atom(
            element=element,
            isotope=isotope
        )

        index = self.molecule.add_atom(atom)

        if self.previous_atom is not None:
            self.molecule.add_bond(
                self.previous_atom,
                index,
                self.pending_bond
            )

        self.previous_atom = index
        self.pending_bond = 1.0

    def parse_bond(self):
        char = self.consume()
        self.pending_bond = BOND_ORDERS[char]

    def parse_branch(self):

        saved_atom = self.previous_atom

        self.parse_until_branch_end()

        self.previous_atom = saved_atom

    def parse_until_branch_end(self):

        while self.pos < len(self.text):

            char = self.peek()

            if char == ")":
                self.consume()
                return

            if char == "(":
                self.consume()
                self.parse_branch()

            elif char in "-=#":
                self.parse_bond()

            elif char == "^":
                self.parse_charge()

            elif char == ".":
                self.parse_radical()

            elif char == ":":
                self.parse_lone_pair()

            elif char == "[":
                self.parse_isotope_atom()

            elif char.isupper():
                self.parse_atom()

            else:
                raise ValueError(
                    f"Unexpected character '{char}' "
                    f"at position {self.pos}"
                )

    def parse_charge(self):

        self.consume()

        start = self.pos

        while self.peek().isdigit():
            self.consume()

        number = self.text[start:self.pos] or "1"

        sign = self.consume()

        if sign not in "+-":
            raise ValueError("Invalid charge")

        if self.previous_atom is None:
            raise ValueError("Charge without atom")

        charge = int(number)

        if sign == "-":
            charge = -charge

        self.molecule.atoms[
            self.previous_atom
        ].charge = charge

    def parse_radical(self):

        count = 0

        while self.peek() == ".":
            self.consume()
            count += 1

        if self.previous_atom is None:
            raise ValueError("Radical without atom")

        self.molecule.atoms[
            self.previous_atom
        ].radicals += count

    def parse_lone_pair(self):

        count = 0

        while self.peek() == ":":
            self.consume()
            count += 1

        if self.previous_atom is None:
            raise ValueError("Lone pair without atom")

        self.molecule.atoms[
            self.previous_atom
        ].lone_pairs = count


# ============================================================
# REACTION DATABASE
# ============================================================

REACTIONS = [

    {
        "name": "Alkene hydrogenation",
        "category": "Hydrocarbons",
        "reactants": ["ethene"],
        "reagents": "H2 / Ni, Pt or Pd",
        "equation": "CH2=CH2 + H2 → CH3-CH3",
        "products": "Ethane",
        "type": "Addition",
        "conditions": "Catalytic hydrogenation",
        "observation": "Double bond is reduced.",
    },

    {
        "name": "Markovnikov hydration",
        "category": "Alkenes",
        "reactants": ["propene"],
        "reagents": "H2O / H+",
        "equation": "CH3-CH=CH2 + H2O → CH3-CH(OH)-CH3",
        "products": "Propan-2-ol",
        "type": "Electrophilic addition",
        "conditions": "Acid-catalysed hydration",
        "observation": "Alcohol is formed.",
    },

    {
        "name": "Anti-Markovnikov addition",
        "category": "Alkenes",
        "reactants": ["propene"],
        "reagents": "HBr / peroxide",
        "equation": "CH3-CH=CH2 + HBr → CH3-CH2-CH2Br",
        "products": "1-Bromopropane",
        "type": "Free-radical addition",
        "conditions": "Presence of peroxide",
        "observation": "Anti-Markovnikov product.",
    },

    {
        "name": "Alcohol oxidation",
        "category": "Alcohols",
        "reactants": ["ethanol"],
        "reagents": "K2Cr2O7 / H+",
        "equation": "CH3CH2OH → CH3COOH",
        "products": "Ethanoic acid",
        "type": "Oxidation",
        "conditions": "Acidified potassium dichromate",
        "observation": "Primary alcohol oxidizes to acid.",
    },

    {
        "name": "Ethanol dehydration",
        "category": "Alcohols",
        "reactants": ["ethanol"],
        "reagents": "Conc. H2SO4, 443 K",
        "equation": "CH3CH2OH → CH2=CH2 + H2O",
        "products": "Ethene + water",
        "type": "Elimination",
        "conditions": "443 K, concentrated sulfuric acid",
        "observation": "Water is eliminated.",
    },

    {
        "name": "Haloalkane hydrolysis",
        "category": "Haloalkanes",
        "reactants": ["chloroethane"],
        "reagents": "Aqueous KOH",
        "equation": "CH3CH2Cl + KOH → CH3CH2OH + KCl",
        "products": "Ethanol",
        "type": "Nucleophilic substitution",
        "conditions": "Aqueous KOH",
        "observation": "Cl is replaced by OH.",
    },

    {
        "name": "Wurtz reaction",
        "category": "Haloalkanes",
        "reactants": ["chloroethane"],
        "reagents": "2Na / dry ether",
        "equation": "2CH3CH2Cl + 2Na → CH3CH2CH2CH3 + 2NaCl",
        "products": "Butane",
        "type": "Coupling",
        "conditions": "Sodium in dry ether",
        "observation": "Two alkyl groups couple.",
    },

    {
        "name": "Williamson ether synthesis",
        "category": "Alcohols and ethers",
        "reactants": ["sodium ethoxide", "chloroethane"],
        "reagents": "Dry ether",
        "equation": "C2H5ONa + C2H5Cl → C2H5OC2H5 + NaCl",
        "products": "Diethyl ether",
        "type": "SN2 nucleophilic substitution",
        "conditions": "Dry ether",
        "observation": "Ether formation.",
    },

    {
        "name": "Esterification",
        "category": "Carboxylic acids",
        "reactants": ["ethanoic acid", "ethanol"],
        "reagents": "Conc. H2SO4",
        "equation": "CH3COOH + C2H5OH ⇌ CH3COOC2H5 + H2O",
        "products": "Ethyl ethanoate",
        "type": "Esterification",
        "conditions": "Acid catalyst",
        "observation": "Fruity-smelling ester.",
    },

    {
        "name": "Aldol addition",
        "category": "Aldehydes and ketones",
        "reactants": ["ethanal"],
        "reagents": "Dilute NaOH",
        "equation": "2CH3CHO → CH3CH(OH)CH2CHO",
        "products": "3-Hydroxybutanal",
        "type": "Aldol addition",
        "conditions": "Dilute base",
        "observation": "Beta-hydroxy aldehyde forms.",
    },

    {
        "name": "Cannizzaro reaction",
        "category": "Aldehydes",
        "reactants": ["benzaldehyde"],
        "reagents": "Conc. NaOH",
        "equation": "2C6H5CHO → C6H5CH2OH + C6H5COONa",
        "products": "Benzyl alcohol + sodium benzoate",
        "type": "Disproportionation",
        "conditions": "Concentrated alkali",
        "observation": "One molecule oxidized, one reduced.",
    },

    {
        "name": "Hoffmann bromamide reaction",
        "category": "Amines",
        "reactants": ["acetamide"],
        "reagents": "Br2 / NaOH",
        "equation": "CH3CONH2 → CH3NH2",
        "products": "Methylamine",
        "type": "Rearrangement",
        "conditions": "Bromine and sodium hydroxide",
        "observation": "Amine has one less carbon.",
    },

    {
        "name": "Diazotisation",
        "category": "Amines",
        "reactants": ["aniline"],
        "reagents": "NaNO2 + HCl, 273–278 K",
        "equation": "C6H5NH2 → C6H5N2+Cl−",
        "products": "Benzenediazonium chloride",
        "type": "Diazotisation",
        "conditions": "Cold acidic solution",
        "observation": "Diazonium salt forms.",
    },

    {
        "name": "Sandmeyer reaction",
        "category": "Amines",
        "reactants": ["benzenediazonium chloride"],
        "reagents": "CuCl / HCl",
        "equation": "C6H5N2+Cl− → C6H5Cl + N2",
        "products": "Chlorobenzene",
        "type": "Substitution",
        "conditions": "CuCl",
        "observation": "Diazonium group replaced.",
    },

]


# ============================================================
# REACTION FUNCTIONS
# ============================================================

def show_reaction(reaction):

    print("\n" + "=" * 60)
    print(reaction["name"])
    print("=" * 60)

    print("Category:", reaction["category"])
    print("Reactants:", ", ".join(reaction["reactants"]))
    print("Reagents:", reaction["reagents"])
    print("Equation:", reaction["equation"])
    print("Products:", reaction["products"])
    print("Reaction type:", reaction["type"])
    print("Conditions:", reaction["conditions"])
    print("Observation:", reaction["observation"])


def search_reactions():

    query = input("\nSearch reaction: ").strip().lower()

    results = []

    for reaction in REACTIONS:

        text = " ".join(
            str(value)
            for value in reaction.values()
        ).lower()

        if query in text:
            results.append(reaction)

    if not results:
        print("No reactions found.")
        return

    for i, reaction in enumerate(results, 1):
        print(f"{i}. {reaction['name']}")

    choice = input(
        "\nEnter number to view, or press Enter to return: "
    ).strip()

    if choice.isdigit():

        index = int(choice) - 1

        if 0 <= index < len(results):
            show_reaction(results[index])


def browse_reactions():

    categories = sorted(
        set(reaction["category"] for reaction in REACTIONS)
    )

    print("\nCategories:")

    for i, category in enumerate(categories, 1):
        print(f"{i}. {category}")

    choice = input(
        "\nSelect category, or press Enter to return: "
    ).strip()

    if not choice.isdigit():
        return

    index = int(choice) - 1

    if not 0 <= index < len(categories):
        return

    category = categories[index]

    results = [
        reaction
        for reaction in REACTIONS
        if reaction["category"] == category
    ]

    print(f"\n{category}")

    for i, reaction in enumerate(results, 1):
        print(f"{i}. {reaction['name']}")

    choice = input("\nSelect reaction: ").strip()

    if choice.isdigit():

        index = int(choice) - 1

        if 0 <= index < len(results):
            show_reaction(results[index])


def parse_structure():

    text = input(
        "\nEnter structure in custom notation: "
    ).strip()

    if not text:
        return

    try:

        molecule = Parser(text).parse()

        molecule.show()

        errors = molecule.validate()

        print("\nValidation:")

        if errors:
            print("INVALID STRUCTURE")

            for error in errors:
                print(" -", error)

        else:
            print("Basic valency validation passed.")

    except Exception as error:
        print("\nPARSER ERROR:", error)

5

def show_help():

    print("\n")
    print("=" * 70)
    print("NCERT + JEE ORGANIC CHEMISTRY SIMULATOR — HELP")
    print("=" * 70)

    print("""
COMMANDS
--------

1  Parse and validate custom structure
   Enter a chemical structure using the custom notation.

2  Search reaction database
   Search by reaction name, reactant, reagent, or category.

3  Browse reactions by category
   View reactions grouped by organic chemistry topic.

4  List all reactions
   Display every reaction currently in the database.

5  Exit
   Close the simulator.

H  Help
   Display this help menu.

CUSTOM CHEMICAL NOTATION
------------------------

ATOMS
-----

C       Carbon
H       Hydrogen
N       Nitrogen
O       Oxygen
F       Fluorine
Cl      Chlorine
Br      Bromine
I       Iodine
S       Sulfur
P       Phosphorus

ISOTOPES
--------

[14C]       Carbon-14
[2H]        Deuterium
[3H]        Tritium

Examples:

[14C]H4
[14C]D4
[14C]D3H
[14C]([14C]T3H)2([14C]D3H)H

D and T may be used as hydrogen isotope symbols.

BONDS
-----

-           Single bond
=           Double bond
#           Triple bond

Examples:

C-C
C=C
C#C

Implicit single bonds:

CC
CH3CH3
H4C

ATOMS
-----

An atom can be written without explicitly specifying
all of its hydrogens.

Examples:

C
CH4
CH3CH3
CH3CH2OH

STRICT NOTATION:

Explicit atom valencies must be satisfied.
Incomplete structures may be rejected.

LONE PAIRS
----------

:C          Lone pair notation
C:          Lone pair notation

::C         Two lone-pair notation
C::         Two lone-pair notation

Examples:

:O
O:
:N
N:

RADICALS
--------

.C          One unpaired electron
C.          One unpaired electron

..C         Two unpaired electrons
C..         Two unpaired electrons

Examples:

CH3.
.C
C..

CHARGES
-------

^+          Charge +1
^-          Charge -1
^2+         Charge +2
^2-         Charge -2

Examples:

C^+
C^-
O^-
N^+
C^2+
O^2-

A charge belongs to the atom immediately before it.

PARENTHESES
-----------

Parentheses represent branches or groups.

Examples:

CH3(CH3)CH3
C(C)(C)C
C(OH)C

ISOTOPIC BRANCHES
-----------------

[14C]([14C]T3H)2([14C]D3H)H

This represents a structure containing isotopically
specified carbon and hydrogen atoms.

RINGS
-----

@n is a ring marker.

Planned notation:

@1
@2
@3

Ring closure and ring bonds are part of the extended
notation system.

PARTIAL BONDS
-------------

Planned notation:

-_{3/2}

Example:

C-_{3/2}C

The mathematical expression inside {} specifies
the desired bond order.

COORDINATE BONDS
----------------

Planned notation:

A->B
A<-B
A=>B
A<=B
A<#B
A#>B

The direction identifies the donor and acceptor.

RESONANCE
---------

Planned notation:

{structure1 | structure2}

Example:

{CH3-C(=O)-O^- | CH3-C(-O^-)=O}

Resonance contributors represent alternative
Lewis structures of the same species.

REACTION DATABASE
-----------------

Search examples:

ethanol
alkene
oxidation
substitution
benzaldehyde
amines
diazotisation
esterification

Reaction information includes:

- Reaction name
- Category
- Reactants
- Reagents
- Equation
- Products
- Reaction type
- Conditions
- Observation

IMPORTANT
---------

This is a chemical notation and reaction simulator.

It is NOT SMILES.

The parser validates basic valency and structure.
It does not yet perform complete quantum chemistry,
automatic mechanism generation, or full resonance analysis.
""")




# ============================================================
# MENU
# ============================================================

def menu():

    while True:

        print("\n")
        print("=" * 60)
        print("NCERT + JEE ORGANIC CHEMISTRY SIMULATOR")
        print("=" * 60)

        print("1. Parse and validate custom structure")
        print("2. Search reaction database")
        print("3. Browse reactions by category")
        print("4. List all reactions")
        print("5. Help")
        print("6. Exit")

        choice = input("\nChoose option: ").strip()

        if choice == "1":
            parse_structure()

        elif choice == "2":
            search_reactions()

        elif choice == "3":
            browse_reactions()

        elif choice == "4":

            for i, reaction in enumerate(REACTIONS, 1):
                print(f"{i}. {reaction['name']}")

        elif choice == "5":
            show_help()

        elif choice == "6":
            print("Exiting.")
            break

        elif choice.lower() == "h":
            show_help()

        else:
            print("Invalid choice.")


if __name__ == "__main__":
    menu()