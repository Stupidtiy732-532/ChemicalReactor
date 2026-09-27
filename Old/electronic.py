from molecule import Molecule


class ElectronicAnalyzer:

    # ==================================================
    # VALENCE ELECTRONS
    # ==================================================

    VALENCE_ELECTRONS = {

        # Period 1
        "H": 1,
        "He": 2,

        # Period 2
        "Li": 1,
        "Be": 2,
        "B": 3,
        "C": 4,
        "N": 5,
        "O": 6,
        "F": 7,
        "Ne": 8,

        # Period 3
        "Na": 1,
        "Mg": 2,
        "Al": 3,
        "Si": 4,
        "P": 5,
        "S": 6,
        "Cl": 7,
        "Ar": 8,

        # Common heavier main-group atoms
        "Ga": 3,
        "Ge": 4,
        "As": 5,
        "Se": 6,
        "Br": 7,
        "Kr": 8,

        "In": 3,
        "Sn": 4,
        "Sb": 5,
        "Te": 6,
        "I": 7,
        "Xe": 8,

        "Tl": 3,
        "Pb": 4,
        "Bi": 5,
        "Po": 6,
        "At": 7,
        "Rn": 8
    }


    # ==================================================
    # COMMON VALENCE STATES
    #
    # These are allowed bonding capacities,
    # not strict chemical laws.
    # ==================================================

    COMMON_VALENCES = {

        "H": [1],

        "Li": [1],
        "Be": [2],
        "B": [3],

        "C": [4],

        "N": [3, 5],
        "P": [3, 5],

        "O": [2],
        "S": [2, 4, 6],

        "F": [1],
        "Cl": [1, 3, 5, 7],
        "Br": [1, 3, 5, 7],
        "I": [1, 3, 5, 7],

        "Si": [4],
        "Ge": [4],
        "Sn": [2, 4],
        "Pb": [2, 4],

        "Al": [3],
        "Ga": [3],
        "In": [1, 3],
        "Tl": [1, 3],

        "As": [3, 5],
        "Sb": [3, 5],
        "Bi": [3, 5],

        "Se": [2, 4, 6],
        "Te": [2, 4, 6],
        "Po": [2, 4, 6]
    }


    def __init__(self, molecule):

        self.molecule = molecule


    # ==================================================
    # BASIC BOND INFORMATION
    # ==================================================

    def bonds_of(self, atom_id):

        return self.molecule.bonds_of(
            atom_id
        )


    def sigma_bond_count(self, atom_id):

        # Every ordinary covalent bond has
        # one sigma component.

        return len(
            self.bonds_of(atom_id)
        )


    def pi_bond_count(self, atom_id):

        count = 0

        for bond in self.bonds_of(atom_id):

            count += max(
                0,
                bond.order - 1
            )

        return count


    def bond_order_sum(self, atom_id):

        return (
            self.molecule.bond_order_sum(
                atom_id
            )
        )


    # ==================================================
    # HYDROGEN INFORMATION
    # ==================================================

    def explicit_hydrogens(self, atom_id):

        atom = self.molecule.get_atom(
            atom_id
        )

        return atom.hydrogens


    def preferred_valence(self, atom_id):

        atom = self.molecule.get_atom(
            atom_id
        )

        element = atom.element

        possible = self.COMMON_VALENCES.get(
            element
        )

        if not possible:

            return None


        current = (
            self.bond_order_sum(atom_id)
            + atom.hydrogens
        )


        # Choose the smallest allowed valence
        # capable of accommodating the current
        # bonding state.

        for valence in possible:

            if valence >= current:

                return valence


        # If hypervalent or unusual,
        # return largest known state.

        return possible[-1]


    def implicit_hydrogens(self, atom_id):

        atom = self.molecule.get_atom(
            atom_id
        )

        # Explicitly specified hydrogen notation
        # is already included separately.

        valence = self.preferred_valence(
            atom_id
        )

        if valence is None:

            return 0


        used_valence = (
            self.bond_order_sum(atom_id)
            + atom.hydrogens
        )


        remaining = (
            valence - used_valence
        )


        # Charged atoms need special treatment.
        # For now we avoid blindly adding H.

        if atom.charge != 0:

            return 0


        return max(
            0,
            remaining
        )


    def total_hydrogens(self, atom_id):

        return (
            self.explicit_hydrogens(atom_id)
            + self.implicit_hydrogens(atom_id)
        )


    # ==================================================
    # VALENCE ELECTRONS
    # ==================================================

    def valence_electrons(self, atom_id):

        atom = self.molecule.get_atom(
            atom_id
        )

        base = self.VALENCE_ELECTRONS.get(
            atom.element
        )

        if base is None:

            return None


        # Positive charge means electrons removed.
        # Negative charge means electrons added.

        return (
            base - atom.charge
        )


    # ==================================================
    # LONE PAIRS
    # ==================================================

    def explicit_lone_pairs(self, atom_id):

        atom = self.molecule.get_atom(
            atom_id
        )

        return atom.lone_pairs


    def inferred_lone_pairs(self, atom_id):

        atom = self.molecule.get_atom(
            atom_id
        )


        # Explicit specification overrides
        # automatic inference.

        if atom.lone_pairs > 0:

            return 0


        valence_electrons = (
            self.valence_electrons(
                atom_id
            )
        )


        if valence_electrons is None:

            return 0


        # Electrons assigned to bonds.
        #
        # For formal electron bookkeeping,
        # each bond contributes one electron
        # to the atom.

        bonding_electrons = (
            self.bond_order_sum(atom_id)
            + atom.hydrogens
        )


        unpaired = (
            atom.radical_electrons
        )


        remaining_electrons = (
            valence_electrons
            - bonding_electrons
            - unpaired
        )


        if remaining_electrons < 0:

            return 0


        return (
            remaining_electrons // 2
        )


    def total_lone_pairs(self, atom_id):

        atom = self.molecule.get_atom(
            atom_id
        )


        if atom.lone_pairs > 0:

            return atom.lone_pairs


        return self.inferred_lone_pairs(
            atom_id
        )


    # ==================================================
    # RADICAL INFORMATION
    # ==================================================

    def radical_electrons(self, atom_id):

        atom = self.molecule.get_atom(
            atom_id
        )

        return atom.radical_electrons


    # ==================================================
    # FORMAL CHARGE CALCULATION
    # ==================================================

    def calculated_formal_charge(
        self,
        atom_id
    ):

        atom = self.molecule.get_atom(
            atom_id
        )

        valence = self.VALENCE_ELECTRONS.get(
            atom.element
        )

        if valence is None:

            return atom.charge


        nonbonding = (

            2 * self.total_lone_pairs(
                atom_id
            )

            + atom.radical_electrons
        )


        bonding = (
            self.bond_order_sum(atom_id)
            + atom.hydrogens
        )


        return (
            valence
            - nonbonding
            - bonding
        )


    # ==================================================
    # VALENCE VALIDATION
    # ==================================================

    def valence_usage(self, atom_id):

        atom = self.molecule.get_atom(
            atom_id
        )

        return (
            self.bond_order_sum(atom_id)
            + atom.hydrogens
        )


    def valence_status(self, atom_id):

        atom = self.molecule.get_atom(
            atom_id
        )

        possible = self.COMMON_VALENCES.get(
            atom.element
        )


        if possible is None:

            return "UNKNOWN"


        usage = self.valence_usage(
            atom_id
        )


        if usage in possible:

            return "VALID"


        if usage < min(possible):

            return "UNDERBONDED"


        if usage > max(possible):

            return "OVERBONDED"


        return "UNUSUAL"


    # ==================================================
    # HYBRIDIZATION
    # ==================================================

    def steric_number(self, atom_id):

        atom = self.molecule.get_atom(
            atom_id
        )

        sigma = (
            self.sigma_bond_count(atom_id)
            + atom.hydrogens
        )


        lone_pairs = (
            self.total_lone_pairs(
                atom_id
            )
        )


        return (
            sigma
            + lone_pairs
        )


    def hybridization(self, atom_id):

        atom = self.molecule.get_atom(
            atom_id
        )


        # Atoms directly participating in
        # multiple bonding require special
        # handling.

        pi_bonds = self.pi_bond_count(
            atom_id
        )


        steric = self.steric_number(
            atom_id
        )


        # Triple bond

        if pi_bonds >= 2:

            return "sp"


        # Double bond

        if pi_bonds == 1:

            return "sp2"


        # No multiple bond

        if steric == 2:

            return "sp"


        if steric == 3:

            return "sp2"


        if steric == 4:

            return "sp3"


        if steric == 5:

            return "sp3d"


        if steric == 6:

            return "sp3d2"


        return "unknown"


    # ==================================================
    # ORBITAL INFORMATION
    # ==================================================

    def available_p_orbital(self, atom_id):

        hybridization = (
            self.hybridization(
                atom_id
            )
        )


        return hybridization in (
            "sp",
            "sp2"
        )


    # ==================================================
    # PI PARTICIPATION
    # ==================================================

    def has_pi_participation(
        self,
        atom_id
    ):

        return (
            self.pi_bond_count(atom_id)
            > 0
        )


    def pi_electron_contribution(
        self,
        atom_id
    ):

        atom = self.molecule.get_atom(
            atom_id
        )


        # Atom in a pi bond contributes
        # one electron to that pi system.

        if self.has_pi_participation(
            atom_id
        ):

            return 1


        # Radical contribution

        if atom.radical_electrons > 0:

            return (
                atom.radical_electrons
            )


        # Lone pair contribution

        if (
            self.total_lone_pairs(
                atom_id
            ) > 0
            and self.available_p_orbital(
                atom_id
            )
        ):

            return 2


        return 0


    # ==================================================
    # COMPLETE ATOM ANALYSIS
    # ==================================================

    def analyze_atom(self, atom_id):

        atom = self.molecule.get_atom(
            atom_id
        )


        return {

            "id": atom.id,

            "element": atom.element,

            "explicit_hydrogens":
                self.explicit_hydrogens(
                    atom_id
                ),

            "implicit_hydrogens":
                self.implicit_hydrogens(
                    atom_id
                ),

            "total_hydrogens":
                self.total_hydrogens(
                    atom_id
                ),

            "bond_order_sum":
                self.bond_order_sum(
                    atom_id
                ),

            "valence_electrons":
                self.valence_electrons(
                    atom_id
                ),

            "explicit_lone_pairs":
                self.explicit_lone_pairs(
                    atom_id
                ),

            "total_lone_pairs":
                self.total_lone_pairs(
                    atom_id
                ),

            "radical_electrons":
                self.radical_electrons(
                    atom_id
                ),

            "formal_charge":
                self.calculated_formal_charge(
                    atom_id
                ),

            "stored_charge":
                atom.charge,

            "valence_status":
                self.valence_status(
                    atom_id
                ),

            "hybridization":
                self.hybridization(
                    atom_id
                ),

            "p_orbital_available":
                self.available_p_orbital(
                    atom_id
                )
        }


    # ==================================================
    # COMPLETE MOLECULAR ANALYSIS
    # ==================================================

    def analyze_molecule(self):

        result = {}

        for atom_id in self.molecule.atoms:

            result[atom_id] = (
                self.analyze_atom(
                    atom_id
                )
            )


        return result