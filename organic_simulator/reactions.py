from analyzer import MoleculeAnalyzer
from reagents import ReagentDatabase


class ReactionSite:
    def __init__(
        self,
        site_type,
        atom_ids,
        description,
    ):
        self.site_type = site_type
        self.atom_ids = atom_ids
        self.description = description

    def __str__(self):
        atoms = ", ".join(
            str(atom_id)
            for atom_id in self.atom_ids
        )

        return (
            f"{self.site_type} "
            f"at atoms [{atoms}]: "
            f"{self.description}"
        )


class ReactionResult:
    def __init__(self):
        self.products = []
        self.byproducts = []
        self.reacted_sites = 0
        self.unreacted_sites = 0
        self.reagent_consumed = 0.0
        self.notes = []
        self.product_molecule = None

    def display(self):
        print("\n--- Reaction Result ---")

        print(
            f"Reacted sites: {self.reacted_sites}"
        )

        print(
            f"Unreacted sites: {self.unreacted_sites}"
        )

        print(
            f"Reagent consumed: "
            f"{self.reagent_consumed:g} mol"
        )

        if self.products:
            print("\nProducts:")

            for product in self.products:
                print(f"  - {product}")

        if self.byproducts:
            print("\nBy-products:")

            for byproduct in self.byproducts:
                print(f"  - {byproduct}")

        if self.notes:
            print("\nNotes:")

            for note in self.notes:
                print(f"  - {note}")

        if self.product_molecule is not None:
            print("\nProduct molecular formula:")
            print(
                f"  {self.product_molecule.formula}"
            )

            print("\nProduct molar mass:")
            print(
                f"  {self.product_molecule.molar_mass:.3f} g/mol"
            )


class ReactionEngine:
    def __init__(self):
        self.database = ReagentDatabase()
        self.analyzer = MoleculeAnalyzer()

    def calculate_reacted_sites(
        self,
        available_sites,
        reagent_moles,
        coefficient=1,
    ):
        if reagent_moles is None:
            return available_sites

        if reagent_moles < 0:
            raise ValueError(
                "Reagent amount cannot be negative."
            )

        possible = int(
            reagent_moles // coefficient
        )

        return min(
            available_sites,
            possible,
        )

    def detect_sites(self, molecule):
        self.analyzer.analyze(molecule)

        sites = []

        for group in molecule.functional_groups:
            name = group.name.lower()

            if "primary alcohol" in name:
                sites.append(
                    ReactionSite(
                        "primary_alcohol",
                        group.atom_ids,
                        "Primary alcohol group",
                    )
                )

            elif "secondary alcohol" in name:
                sites.append(
                    ReactionSite(
                        "secondary_alcohol",
                        group.atom_ids,
                        "Secondary alcohol group",
                    )
                )

            elif "tertiary alcohol" in name:
                sites.append(
                    ReactionSite(
                        "tertiary_alcohol",
                        group.atom_ids,
                        "Tertiary alcohol group",
                    )
                )

            elif name == "alcohol":
                sites.append(
                    ReactionSite(
                        "alcohol",
                        group.atom_ids,
                        "Alcohol group",
                    )
                )

            elif "carboxylic acid" in name:
                sites.append(
                    ReactionSite(
                        "carboxylic_acid",
                        group.atom_ids,
                        "Carboxylic acid group",
                    )
                )

            elif "aldehyde" in name:
                sites.append(
                    ReactionSite(
                        "aldehyde",
                        group.atom_ids,
                        "Aldehyde group",
                    )
                )

            elif "ketone" in name:
                sites.append(
                    ReactionSite(
                        "ketone",
                        group.atom_ids,
                        "Ketone group",
                    )
                )

            elif "halo" in name:
                sites.append(
                    ReactionSite(
                        "haloalkane",
                        group.atom_ids,
                        "Carbon-halogen bond",
                    )
                )

        for bond in molecule.bonds:
            if bond.order == 2:
                atom1 = molecule.atom_by_id(
                    bond.atom1
                )

                atom2 = molecule.atom_by_id(
                    bond.atom2
                )

                if (
                    atom1.element == "C"
                    and atom2.element == "C"
                ):
                    sites.append(
                        ReactionSite(
                            "alkene",
                            [
                                bond.atom1,
                                bond.atom2,
                            ],
                            "Carbon-carbon double bond",
                        )
                    )

            elif bond.order == 3:
                atom1 = molecule.atom_by_id(
                    bond.atom1
                )

                atom2 = molecule.atom_by_id(
                    bond.atom2
                )

                if (
                    atom1.element == "C"
                    and atom2.element == "C"
                ):
                    sites.append(
                        ReactionSite(
                            "alkyne",
                            [
                                bond.atom1,
                                bond.atom2,
                            ],
                            "Carbon-carbon triple bond",
                        )
                    )

        return sites

    def choose_sites(
        self,
        sites,
        allowed_types,
        reagent_moles,
        excess,
    ):
        eligible = [
            site
            for site in sites
            if site.site_type in allowed_types
        ]

        if excess:
            reacted = len(eligible)

        else:
            reacted = self.calculate_reacted_sites(
                len(eligible),
                reagent_moles,
            )

        return eligible[:reacted], len(eligible) - reacted

    def react(
        self,
        molecule,
        reagent_text,
        reagent_moles=None,
        excess=False,
    ):
        reagent = self.database.find(
            reagent_text
        )

        if reagent is None:
            raise ValueError(
                f"Unknown reagent: {reagent_text}"
            )

        sites = self.detect_sites(molecule)

        reaction_type = reagent.reaction_type

        if reaction_type == "sodium":
            return self.react_sodium(
                molecule,
                sites,
                reagent_moles,
                excess,
            )

        if reaction_type == "bicarbonate":
            return self.react_bicarbonate(
                molecule,
                sites,
                reagent_moles,
                excess,
            )

        if reaction_type == "hydrogenation":
            return self.react_hydrogenation(
                molecule,
                sites,
                reagent_moles,
                excess,
            )

        if reaction_type == "bromination":
            return self.react_bromination(
                molecule,
                sites,
                reagent_moles,
                excess,
            )

        if reaction_type == "aqueous_koh":
            return self.react_aqueous_koh(
                molecule,
                sites,
                reagent_moles,
                excess,
            )

        if reaction_type == "alcoholic_koh":
            return self.react_alcoholic_koh(
                molecule,
                sites,
                reagent_moles,
                excess,
            )

        if reaction_type in {
            "permanganate",
            "dichromate",
            "pcc",
            "tollens",
            "fehling",
        }:
            return self.react_oxidation(
                molecule,
                sites,
                reagent,
                reagent_moles,
                excess,
            )

        if reaction_type == "esterification":
            return self.react_esterification(
                molecule,
                sites,
                reagent_moles,
                excess,
            )

        raise ValueError(
            f"Reaction not implemented: {reaction_type}"
        )

    def react_sodium(
        self,
        molecule,
        sites,
        reagent_moles,
        excess,
    ):
        result = ReactionResult()

        reacted_sites, unreacted = self.choose_sites(
            sites,
            [
                "alcohol",
                "primary_alcohol",
                "secondary_alcohol",
                "tertiary_alcohol",
                "carboxylic_acid",
            ],
            reagent_moles,
            excess,
        )

        result.reacted_sites = len(
            reacted_sites
        )

        result.unreacted_sites = unreacted
        result.reagent_consumed = len(
            reacted_sites
        )

        alcohol_count = sum(
            1
            for site in reacted_sites
            if "alcohol" in site.site_type
        )

        acid_count = sum(
            1
            for site in reacted_sites
            if site.site_type == "carboxylic_acid"
        )

        if alcohol_count:
            result.products.append(
                f"{alcohol_count} mol sodium alkoxide"
            )

        if acid_count:
            result.products.append(
                f"{acid_count} mol sodium carboxylate"
            )

        if result.reacted_sites:
            result.byproducts.append(
                f"{result.reacted_sites / 2:g} mol H2"
            )

        result.notes.append(
            "Each acidic O-H site consumes one mol Na."
        )

        return result

    def react_bicarbonate(
        self,
        molecule,
        sites,
        reagent_moles,
        excess,
    ):
        result = ReactionResult()

        reacted_sites, unreacted = self.choose_sites(
            sites,
            ["carboxylic_acid"],
            reagent_moles,
            excess,
        )

        result.reacted_sites = len(
            reacted_sites
        )

        result.unreacted_sites = unreacted
        result.reagent_consumed = len(
            reacted_sites
        )

        if reacted_sites:
            result.products.append(
                f"{len(reacted_sites)} mol sodium carboxylate"
            )

            result.byproducts.append(
                f"{len(reacted_sites)} mol CO2"
            )

            result.byproducts.append(
                f"{len(reacted_sites)} mol H2O"
            )

        result.notes.append(
            "NaHCO3 reacts with carboxylic acids, "
            "not ordinary alcohols."
        )

        return result

    def react_hydrogenation(
        self,
        molecule,
        sites,
        reagent_moles,
        excess,
    ):
        result = ReactionResult()

        reacted_sites, unreacted = self.choose_sites(
            sites,
            ["alkene", "alkyne"],
            reagent_moles,
            excess,
        )

        product = molecule.clone()

        for site in reacted_sites:
            atom1, atom2 = site.atom_ids

            bond = product.bond_between(
                atom1,
                atom2,
            )

            if bond is None:
                continue

            if site.site_type == "alkene":
                bond.order = 1

                product.add_hydrogen(atom1)
                product.add_hydrogen(atom2)

            elif site.site_type == "alkyne":
                bond.order = 2

                product.add_hydrogen(atom1)
                product.add_hydrogen(atom2)

        result.product_molecule = product
        result.reacted_sites = len(
            reacted_sites
        )

        result.unreacted_sites = unreacted
        result.reagent_consumed = len(
            reacted_sites
        )

        if reacted_sites:
            result.products.append(
                "Hydrogenated hydrocarbon"
            )

        result.notes.append(
            "C=C consumes one mol H2; "
            "C#C consumes two mol H2 for complete hydrogenation."
        )

        return result

    def react_bromination(
        self,
        molecule,
        sites,
        reagent_moles,
        excess,
    ):
        result = ReactionResult()

        reacted_sites, unreacted = self.choose_sites(
            sites,
            ["alkene"],
            reagent_moles,
            excess,
        )

        product = molecule.clone()

        for site in reacted_sites:
            atom1, atom2 = site.atom_ids

            bond = product.bond_between(
                atom1,
                atom2,
            )

            if bond is None:
                continue

            bond.order = 1

            bromine1 = product.add_atom("Br")
            bromine2 = product.add_atom("Br")

            product.add_bond(
                atom1,
                bromine1,
                1,
            )

            product.add_bond(
                atom2,
                bromine2,
                1,
            )

        result.product_molecule = product
        result.reacted_sites = len(
            reacted_sites
        )

        result.unreacted_sites = unreacted
        result.reagent_consumed = len(
            reacted_sites
        )

        if reacted_sites:
            result.products.append(
                "Vicinal dibromoalkane"
            )

        result.notes.append(
            "Each C=C bond consumes one mol Br2."
        )

        return result

    def react_aqueous_koh(
        self,
        molecule,
        sites,
        reagent_moles,
        excess,
    ):
        result = ReactionResult()

        reacted_sites, unreacted = self.choose_sites(
            sites,
            ["haloalkane"],
            reagent_moles,
            excess,
        )

        product = molecule.clone()

        for site in reacted_sites:
            carbon_id, halogen_id = site.atom_ids

            halogen = product.atom_by_id(
                halogen_id
            )

            product.replace_element(
                halogen_id,
                "O",
            )

            halogen.hydrogens = 1

        result.product_molecule = product
        result.reacted_sites = len(
            reacted_sites
        )

        result.unreacted_sites = unreacted
        result.reagent_consumed = len(
            reacted_sites
        )

        if reacted_sites:
            result.products.append(
                "Alcohol"
            )

            result.byproducts.append(
                "KX"
            )

        result.notes.append(
            "Aqueous KOH gives nucleophilic substitution."
        )

        return result

    def react_alcoholic_koh(
        self,
        molecule,
        sites,
        reagent_moles,
        excess,
    ):
        result = ReactionResult()

        reacted_sites, unreacted = self.choose_sites(
            sites,
            ["haloalkane"],
            reagent_moles,
            excess,
        )

        product = molecule.clone()

        for site in reacted_sites:
            carbon_id, halogen_id = site.atom_ids

            carbon = product.atom_by_id(
                carbon_id
            )

            halogen = product.atom_by_id(
                halogen_id
            )

            carbon_neighbours = [
                neighbour_id
                for neighbour_id in product.neighbours(
                    carbon_id
                )
                if neighbour_id != halogen_id
            ]

            beta_carbon_id = None

            for neighbour_id in carbon_neighbours:
                neighbour = product.atom_by_id(
                    neighbour_id
                )

                if neighbour.element == "C":
                    beta_carbon_id = neighbour_id
                    break

            if beta_carbon_id is None:
                continue

            product.remove_bond(
                carbon_id,
                halogen_id,
            )

            product.atoms = [
                atom
                for atom in product.atoms
                if atom.id != halogen_id
            ]

            product.set_bond_order(
                carbon_id,
                beta_carbon_id,
                2,
            )

            beta_carbon = product.atom_by_id(
                beta_carbon_id
            )

            if beta_carbon.hydrogens > 0:
                beta_carbon.hydrogens -= 1

        result.product_molecule = product
        result.reacted_sites = len(
            reacted_sites
        )

        result.unreacted_sites = unreacted
        result.reagent_consumed = len(
            reacted_sites
        )

        if reacted_sites:
            result.products.append(
                "Alkene"
            )

            result.byproducts.append(
                "KX + H2O"
            )

        result.notes.append(
            "Alcoholic KOH promotes beta-elimination."
        )

        return result

    def react_oxidation(
        self,
        molecule,
        sites,
        reagent,
        reagent_moles,
        excess,
    ):
        result = ReactionResult()

        allowed = [
            "primary_alcohol",
            "secondary_alcohol",
            "aldehyde",
        ]

        reacted_sites, unreacted = self.choose_sites(
            sites,
            allowed,
            reagent_moles,
            excess,
        )

        product = molecule.clone()

        for site in reacted_sites:
            carbon_id, oxygen_id = site.atom_ids

            carbon = product.atom_by_id(
                carbon_id
            )

            oxygen = product.atom_by_id(
                oxygen_id
            )

            if site.site_type == "primary_alcohol":
                carbon.hydrogens = max(
                    0,
                    carbon.hydrogens - 1,
                )

                product.set_bond_order(
                    carbon_id,
                    oxygen_id,
                    2,
                )

                oxygen.hydrogens = 0

            elif site.site_type == "secondary_alcohol":
                carbon.hydrogens = max(
                    0,
                    carbon.hydrogens - 1,
                )

                product.set_bond_order(
                    carbon_id,
                    oxygen_id,
                    2,
                )

                oxygen.hydrogens = 0

            elif site.site_type == "aldehyde":
                oxygen.hydrogens = 0

                carbon.hydrogens = 0

                product.set_bond_order(
                    carbon_id,
                    oxygen_id,
                    2,
                )

        result.product_molecule = product
        result.reacted_sites = len(
            reacted_sites
        )

        result.unreacted_sites = unreacted
        result.reagent_consumed = len(
            reacted_sites
        )

        if reacted_sites:
            result.products.append(
                "Oxidized product"
            )

        result.notes.append(
            f"Oxidation performed using {reagent.name}."
        )

        return result

    def react_esterification(
        self,
        molecule,
        sites,
        reagent_moles,
        excess,
    ):
        result = ReactionResult()

        alcohol_sites = [
            site
            for site in sites
            if "alcohol" in site.site_type
        ]

        acid_sites = [
            site
            for site in sites
            if site.site_type == "carboxylic_acid"
        ]

        possible = min(
            len(alcohol_sites),
            len(acid_sites),
        )

        if excess:
            reacted = possible

        else:
            reacted = self.calculate_reacted_sites(
                possible,
                reagent_moles,
            )

        result.reacted_sites = reacted
        result.unreacted_sites = possible - reacted
        result.reagent_consumed = reacted

        if reacted:
            result.products.append(
                f"{reacted} mol ester"
            )

            result.byproducts.append(
                f"{reacted} mol H2O"
            )

        result.notes.append(
            "Esterification requires a carboxylic acid, "
            "an alcohol, and acidic conditions."
        )

        return result