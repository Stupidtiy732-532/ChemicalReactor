# reactions.py

from analyzer import detect_functional_groups
from molecule import Molecule


class ReactionResult:
    def __init__(self, reactant, products, equation, explanation):
        self.reactant = reactant
        self.products = products
        self.equation = equation
        self.explanation = explanation

    def display(self):
        print("\nREACTION")
        print("-" * 60)
        print(self.equation)
        print()
        print(self.explanation)
        print()
        print("PRODUCTS")

        for product in self.products:
            if isinstance(product, Molecule):
                print(
                    f"{product.condensed()}   "
                    f"[{product.formula()}]"
                )
            else:
                print(product)


class ReactionEngine:
    def react(self, molecule, reagent, conditions=""):
        reagent_clean = (
            reagent.lower()
            .replace(" ", "")
            .replace("₄", "4")
            .replace("₂", "2")
        )

        conditions_clean = conditions.lower().replace(" ", "")

        # Hydrogenation must be checked BEFORE generic reduction.
        if reagent_clean in {
            "h2/ni",
            "h2/pd",
            "h2/pt",
            "h2",
            "h₂/ni",
            "h₂",
        }:
            return self.hydrogenation(molecule)

        if self.is_oxidant(reagent_clean):
            return self.oxidation(molecule)

        if self.is_reducing_agent(reagent_clean):
            return self.reduction(molecule)

        if self.is_aqueous_koh(reagent_clean, conditions_clean):
            return self.aqueous_substitution(molecule)

        if self.is_alcoholic_koh(reagent_clean, conditions_clean):
            try:
                return self.elimination(molecule)
            except ValueError:
                raise ValueError(
                    "KOH(alc) eliminates HX from haloalkanes. "
                    "It does not normally dehydrate alcohols. "
                    "For alcohol dehydration use conc. H2SO4 or H3PO4 with heat."
                )

        if reagent_clean in {"br2", "bromine", "br₂"}:
            return self.bromination(molecule)

        if reagent_clean in {"na", "sodium"}:
            return self.sodium_reaction(molecule)

        if reagent_clean in {
            "nahco3",
            "nahco₃",
            "sodiumbicarbonate",
        }:
            return self.bicarbonate_reaction(molecule)

        if reagent_clean in {"hbr", "hcl", "hi"}:
            return self.hydrohalogenation(
                molecule,
                reagent_clean,
            )

        if reagent_clean in {"h2o/h+", "h2o", "dil.h2so4"}:
            return self.hydration(molecule)

        if reagent_clean in {
            "conc.h2so4",
            "h2so4",
            "h3po4",
        }:
            return self.dehydration(molecule)

        raise ValueError(
            "No reaction rule implemented for this reagent yet."
        )

    def is_oxidant(self, reagent):
        oxidants = {
            "kmno4",
            "kmno4/h+",
            "k2cr2o7",
            "k2cr2o7/h+",
            "pcc",
            "jones",
            "cro3",
            "oxidation",
        }

        return reagent in oxidants

    def is_reducing_agent(self, reagent):
        reducers = {
            "nab h4".replace(" ", ""),
            "nabh4",
            "lialh4",
            "h2/ni",
            "h2/pd",
            "h2/pt",
            "reduction",
        }

        return reagent in reducers

    def is_aqueous_koh(self, reagent, conditions):
        return (
            reagent in {"koh", "naoh"}
            and (
                "aq" in conditions
                or "aqueous" in conditions
                or conditions == ""
            )
        )

    def is_alcoholic_koh(self, reagent, conditions):
        return (
                reagent in {"koh", "naoh"}
                and (
                        "alc" in conditions
                        or "alcoholic" in conditions
                )
        )

    def oxidation(self, molecule):
        product = molecule.clone()
        groups = detect_functional_groups(product)

        if "primary alcohol" in groups:
            self.convert_primary_alcohol_to_aldehyde(product)

            return ReactionResult(
                molecule,
                [product],
                f"{molecule.condensed()} + [O] → "
                f"{product.condensed()}",
                "Primary alcohol oxidised to aldehyde.",
            )

        if "secondary alcohol" in groups:
            self.convert_secondary_alcohol_to_ketone(product)

            return ReactionResult(
                molecule,
                [product],
                f"{molecule.condensed()} + [O] → "
                f"{product.condensed()}",
                "Secondary alcohol oxidised to ketone.",
            )

        if "aldehyde" in groups:
            self.convert_aldehyde_to_acid(product)

            return ReactionResult(
                molecule,
                [product],
                f"{molecule.condensed()} + [O] → "
                f"{product.condensed()}",
                "Aldehyde oxidised to carboxylic acid.",
            )

        raise ValueError("No oxidisable functional group found.")

    def reduction(self, molecule):
        product = molecule.clone()
        groups = detect_functional_groups(product)

        if "aldehyde" in groups:
            self.convert_aldehyde_to_primary_alcohol(product)

            return ReactionResult(
                molecule,
                [product],
                f"{molecule.condensed()} + 2[H] → "
                f"{product.condensed()}",
                "Aldehyde reduced to primary alcohol.",
            )

        if "ketone" in groups:
            self.convert_ketone_to_secondary_alcohol(product)

            return ReactionResult(
                molecule,
                [product],
                f"{molecule.condensed()} + 2[H] → "
                f"{product.condensed()}",
                "Ketone reduced to secondary alcohol.",
            )

        raise ValueError("No reducible carbonyl group found.")

    def convert_primary_alcohol_to_aldehyde(self, molecule):
        for atom in molecule.atoms.values():
            if atom.element != "C":
                continue

            neighbors = molecule.neighbors(atom.id)

            oxygen = None

            for n, order in neighbors:
                if molecule.atoms[n].element == "O" and order == 1:
                    oxygen = molecule.atoms[n]
                    break

            if oxygen and oxygen.hydrogens == 1:
                oxygen.hydrogens = 0
                atom.hydrogens += 1

                carbon_neighbors = [
                    n for n, _ in neighbors
                    if molecule.atoms[n].element == "C"
                ]

                if len(carbon_neighbors) == 1:
                    molecule.remove_bond(atom.id, oxygen.id)
                    molecule.add_bond(atom.id, oxygen.id, 2)
                    return

    def convert_secondary_alcohol_to_ketone(self, molecule):
        for atom in molecule.atoms.values():
            if atom.element != "C":
                continue

            carbon_neighbors = [
                n for n, _ in molecule.neighbors(atom.id)
                if molecule.atoms[n].element == "C"
            ]

            oxygen = None

            for n, order in molecule.neighbors(atom.id):
                if molecule.atoms[n].element == "O" and order == 1:
                    oxygen = molecule.atoms[n]
                    break

            if oxygen and oxygen.hydrogens == 1:
                if len(carbon_neighbors) == 2:
                    oxygen.hydrogens = 0
                    molecule.remove_bond(atom.id, oxygen.id)
                    molecule.add_bond(atom.id, oxygen.id, 2)
                    atom.hydrogens = 0
                    return

    def convert_aldehyde_to_acid(self, molecule):
        for atom in molecule.atoms.values():
            if atom.element != "C" or atom.hydrogens < 1:
                continue

            oxygen = None

            for n, order in molecule.neighbors(atom.id):
                if molecule.atoms[n].element == "O" and order == 2:
                    oxygen = molecule.atoms[n]
                    break

            if oxygen:
                atom.hydrogens -= 1
                hydroxyl = molecule.add_atom("O", 1)
                molecule.add_bond(atom.id, hydroxyl, 1)
                return

    def convert_aldehyde_to_primary_alcohol(self, molecule):
        for atom in molecule.atoms.values():
            if atom.element != "C":
                continue

            oxygen = None

            for n, order in molecule.neighbors(atom.id):
                if molecule.atoms[n].element == "O" and order == 2:
                    oxygen = molecule.atoms[n]
                    break

            if oxygen and atom.hydrogens >= 1:
                molecule.remove_bond(atom.id, oxygen.id)
                molecule.add_bond(atom.id, oxygen.id, 1)
                oxygen.hydrogens = 1
                atom.hydrogens += 1
                return

    def convert_ketone_to_secondary_alcohol(self, molecule):
        for atom in molecule.atoms.values():
            if atom.element != "C":
                continue

            oxygen = None

            for n, order in molecule.neighbors(atom.id):
                if molecule.atoms[n].element == "O" and order == 2:
                    oxygen = molecule.atoms[n]
                    break

            if oxygen:
                molecule.remove_bond(atom.id, oxygen.id)
                molecule.add_bond(atom.id, oxygen.id, 1)
                oxygen.hydrogens = 1
                atom.hydrogens += 1
                return

    def aqueous_substitution(self, molecule):
        product = molecule.clone()

        for atom in list(product.atoms.values()):
            if atom.element not in {"F", "Cl", "Br", "I"}:
                continue

            attached_carbon = None

            for n, order in product.neighbors(atom.id):
                if product.atoms[n].element == "C":
                    attached_carbon = n
                    break

            if attached_carbon is not None:
                product.atoms[atom.id].element = "O"
                product.atoms[atom.id].hydrogens = 1

                return ReactionResult(
                    molecule,
                    [product, "KBr / NaBr"],
                    f"{molecule.condensed()} + KOH(aq) → "
                    f"{product.condensed()}",
                    "Aqueous hydroxide substitutes the halogen.",
                )

        raise ValueError("No C–halogen bond found.")

    def elimination(self, molecule):
        product = molecule.clone()

        for atom in product.atoms.values():
            if atom.element not in {"F", "Cl", "Br", "I"}:
                continue

            carbon = None

            for n, order in product.neighbors(atom.id):
                if product.atoms[n].element == "C":
                    carbon = n
                    break

            if carbon is None:
                continue

            beta_carbon = None

            for n, order in product.neighbors(carbon):
                if n == atom.id:
                    continue

                if product.atoms[n].element == "C":
                    beta_carbon = n
                    break

            if beta_carbon is not None:
                product.remove_bond(carbon, beta_carbon)
                product.add_bond(carbon, beta_carbon, 2)
                del product.atoms[atom.id]
                product.bonds = [
                    bond for bond in product.bonds
                    if atom.id not in {bond.a, bond.b}
                ]

                return ReactionResult(
                    molecule,
                    [product, "KCl / KBr / KI"],
                    f"{molecule.condensed()} + KOH(alc), Δ → "
                    f"{product.condensed()}",
                    "Alcoholic hydroxide causes β-elimination.",
                )

        raise ValueError("No suitable haloalkane for elimination.")

    def bromination(self, molecule):
        product = molecule.clone()

        for bond in product.bonds:
            if bond.order == 2:
                a = product.atoms[bond.a]
                b = product.atoms[bond.b]

                if a.element == "C" and b.element == "C":
                    bond.order = 1

                    br1 = product.add_atom("Br")
                    br2 = product.add_atom("Br")

                    product.add_bond(a.id, br1, 1)
                    product.add_bond(b.id, br2, 1)

                    return ReactionResult(
                        molecule,
                        [product],
                        f"{molecule.condensed()} + Br₂ → "
                        f"{product.condensed()}",
                        "Bromine adds across the C=C bond.",
                    )

        raise ValueError("No C=C bond found.")

    def hydrogenation(self, molecule):
        product = molecule.clone()

        for bond in product.bonds:
            if bond.order not in {2, 3}:
                continue

            a = product.atoms[bond.a]
            b = product.atoms[bond.b]

            if a.element != "C" or b.element != "C":
                continue

            if bond.order == 2:
                # Alkene → alkane
                bond.order = 1
                a.hydrogens += 1
                b.hydrogens += 1

            elif bond.order == 3:
                # Alkyne → alkene
                bond.order = 2
                a.hydrogens += 1
                b.hydrogens += 1

            return ReactionResult(
                molecule,
                [product],
                f"{molecule.condensed()} + H₂/Ni → "
                f"{product.condensed()}",
                "Catalytic hydrogenation reduces the carbon-carbon multiple bond.",
            )

        raise ValueError(
            "No C=C or C≡C bond found. "
            "Use notation such as CH2=CH2 or CH3-C#CH."
        )

    def sodium_reaction(self, molecule):
        groups = detect_functional_groups(molecule)

        if not any(
            group in groups
            for group in {
                "primary alcohol",
                "secondary alcohol",
                "tertiary alcohol",
                "alcohol",
            }
        ):
            raise ValueError("Sodium reaction requires an alcohol.")

        return ReactionResult(
            molecule,
            [molecule.clone(), "H₂ gas", "sodium alkoxide"],
            f"{molecule.condensed()} + Na → "
            f"alkoxide + ½H₂",
            "Alcohol reacts with sodium, releasing hydrogen gas.",
        )

    def bicarbonate_reaction(self, molecule):
        groups = detect_functional_groups(molecule)

        if "carboxylic acid" not in groups:
            raise ValueError(
                "NaHCO₃ reaction requires a carboxylic acid."
            )

        return ReactionResult(
            molecule,
            [molecule.clone(), "CO₂ gas", "H₂O", "sodium carboxylate"],
            f"{molecule.condensed()} + NaHCO₃ → "
            f"sodium carboxylate + CO₂ + H₂O",
            "Carboxylic acid gives brisk effervescence with sodium bicarbonate.",
        )

    def hydrohalogenation(self, molecule, reagent):
        product = molecule.clone()

        for bond in product.bonds:
            if bond.order == 2:
                a = product.atoms[bond.a]
                b = product.atoms[bond.b]

                if a.element == "C" and b.element == "C":
                    bond.order = 1

                    halogen = {
                        "hbr": "Br",
                        "hcl": "Cl",
                        "hi": "I",
                    }[reagent]

                    halogen_id = product.add_atom(halogen)
                    product.add_bond(a.id, halogen_id, 1)

                    a.hydrogens += 1
                    b.hydrogens += 1

                    return ReactionResult(
                        molecule,
                        [product],
                        f"{molecule.condensed()} + {reagent.upper()} → "
                        f"{product.condensed()}",
                        "Hydrohalogenation of the alkene.",
                    )

        raise ValueError("No alkene found.")

    def hydration(self, molecule):
        product = molecule.clone()

        for bond in product.bonds:
            if bond.order == 2:
                a = product.atoms[bond.a]
                b = product.atoms[bond.b]

                if a.element == "C" and b.element == "C":
                    bond.order = 1

                    oxygen = product.add_atom("O", 1)
                    product.add_bond(a.id, oxygen, 1)

                    a.hydrogens += 1
                    b.hydrogens += 1

                    return ReactionResult(
                        molecule,
                        [product],
                        f"{molecule.condensed()} + H₂O → "
                        f"{product.condensed()}",
                        "Acid-catalysed hydration of the alkene.",
                    )

        raise ValueError("No alkene found.")

    def dehydration(self, molecule):
        product = molecule.clone()

        for atom in list(product.atoms.values()):
            if atom.element != "O" or atom.hydrogens != 1:
                continue

            carbon = None

            for n, order in product.neighbors(atom.id):
                if product.atoms[n].element == "C":
                    carbon = n
                    break

            if carbon is not None:
                beta_carbon = None

                for n, order in product.neighbors(carbon):
                    if n != atom.id and product.atoms[n].element == "C":
                        beta_carbon = n
                        break

                if beta_carbon is not None:
                    product.remove_bond(carbon, beta_carbon)
                    product.add_bond(carbon, beta_carbon, 2)

                    del product.atoms[atom.id]
                    product.bonds = [
                        bond for bond in product.bonds
                        if atom.id not in {bond.a, bond.b}
                    ]

                    return ReactionResult(
                        molecule,
                        [product, "H₂O"],
                        f"{molecule.condensed()} "
                        f"→ {product.condensed()} + H₂O",
                        "Alcohol dehydrated to alkene.",
                    )

        raise ValueError("No alcohol suitable for dehydration.")