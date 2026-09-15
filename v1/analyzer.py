from collections import defaultdict, deque

from electronic import ElectronicAnalyzer


class MolecularAnalyzer:
    """
    Graph-based molecular structure analyzer.

    Handles:
        - Molecular graph
        - Neighbors and degrees
        - Bond order information
        - Sigma/pi bond counts
        - Cycle detection
        - Ring analysis
        - Conjugation candidates
        - Pi-electron counting
        - Aromaticity estimation
        - Approximate hybridization
        - Approximate geometry
        - Valence warnings
        - Overall molecular summary

    This is a symbolic first-pass analyzer.
    It is not a quantum-mechanical calculation.
    """

    # Approximate common valences.
    COMMON_VALENCES = {
        "H": {1},
        "B": {3},
        "C": {4},
        "N": {3, 5},
        "O": {2},
        "F": {1},
        "P": {3, 5},
        "S": {2, 4, 6},
        "Cl": {1},
        "Br": {1},
        "I": {1},
    }

    def __init__(self, molecule):
        self.molecule = molecule

        # Atom ID -> list of (neighbor atom ID, Bond object)
        self.graph = defaultdict(list)

        for bond in molecule.bonds:
            self.graph[bond.atom1].append(
                (bond.atom2, bond)
            )

            self.graph[bond.atom2].append(
                (bond.atom1, bond)
            )

    # ==========================================================
    # BASIC GRAPH INFORMATION
    # ==========================================================

    def neighbors(self, atom_id):
        """
        Return neighboring atom IDs.
        """

        return [
            atom_id
            for atom_id, bond in self.graph[atom_id]
        ]

    def degree(self, atom_id):
        """
        Number of directly bonded neighboring atoms.
        """

        return len(self.graph[atom_id])

    def bonds_of(self, atom_id):
        """
        Return all Bond objects connected to an atom.
        """

        return [
            bond
            for _, bond in self.graph[atom_id]
        ]

    def bond_order_sum(self, atom_id):
        """
        Sum of all bond orders connected to an atom.

        Single bond  = 1
        Double bond  = 2
        Triple bond  = 3
        """

        return sum(
            bond.order
            for _, bond in self.graph[atom_id]
        )

    def sigma_bond_count(self, atom_id):
        """
        Every connected bond contributes one sigma bond.
        """

        return len(self.graph[atom_id])

    def pi_bond_count(self, atom_id):
        """
        Count pi bonds connected to an atom.

        Single bond:
            0 pi bonds

        Double bond:
            1 pi bond

        Triple bond:
            2 pi bonds
        """

        return sum(
            max(0, bond.order - 1)
            for _, bond in self.graph[atom_id]
        )

    # ==========================================================
    # CONNECTED COMPONENTS
    # ==========================================================

    def connected_components(self):
        """
        Find connected components of the molecular graph.
        """

        visited = set()
        components = []

        for atom_id in self.molecule.atoms:

            if atom_id in visited:
                continue

            component = []
            queue = deque([atom_id])
            visited.add(atom_id)

            while queue:

                current = queue.popleft()
                component.append(current)

                for neighbour in self.neighbors(current):

                    if neighbour not in visited:
                        visited.add(neighbour)
                        queue.append(neighbour)

            components.append(component)

        return components

    # ==========================================================
    # CYCLE DETECTION
    # ==========================================================

    def find_cycles(self):
        """
        Find simple cycles in the molecular graph.

        This uses depth-first search and cycle normalization.
        """

        cycles = []

        def dfs(
            current,
            parent,
            path,
            visited
        ):

            visited.add(current)
            path.append(current)

            for neighbour, _ in self.graph[current]:

                if neighbour == parent:
                    continue

                if neighbour in path:

                    start = path.index(neighbour)

                    cycle = path[start:]

                    if len(cycle) >= 3:

                        normalized = self.normalize_cycle(
                            cycle
                        )

                        if normalized not in cycles:
                            cycles.append(normalized)

                    continue

                if neighbour not in visited:

                    dfs(
                        neighbour,
                        current,
                        path.copy(),
                        visited.copy()
                    )

        for atom_id in self.molecule.atoms:

            dfs(
                atom_id,
                None,
                [],
                set()
            )

        return cycles

    @staticmethod
    def normalize_cycle(cycle):
        """
        Normalize a cycle so that rotations and reverse
        traversal directions produce the same result.
        """

        cycle = list(cycle)
        rotations = []

        for i in range(len(cycle)):

            rotations.append(
                tuple(
                    cycle[i:] + cycle[:i]
                )
            )

        reversed_cycle = list(
            reversed(cycle)
        )

        for i in range(len(reversed_cycle)):

            rotations.append(
                tuple(
                    reversed_cycle[i:]
                    + reversed_cycle[:i]
                )
            )

        return min(rotations)

    # ==========================================================
    # CYCLE BOND INFORMATION
    # ==========================================================

    def cycle_bonds(self, cycle):
        """
        Return bonds forming the perimeter of a cycle.
        """

        cycle = list(cycle)
        result = []

        for i in range(len(cycle)):

            atom1 = cycle[i]

            atom2 = cycle[
                (i + 1) % len(cycle)
            ]

            for neighbour, bond in self.graph[atom1]:

                if neighbour == atom2:
                    result.append(bond)
                    break

        return result

    def cycle_bond_map(self, cycle):
        """
        Return a map of atom-pair keys to Bond objects
        for bonds belonging to a cycle.
        """

        cycle = list(cycle)
        bond_map = {}

        for i in range(len(cycle)):

            atom1 = cycle[i]

            atom2 = cycle[
                (i + 1) % len(cycle)
            ]

            key = frozenset(
                (atom1, atom2)
            )

            for neighbour, bond in self.graph[atom1]:

                if neighbour == atom2:
                    bond_map[key] = bond
                    break

        return bond_map

    def double_bonds_in_cycle(self, cycle):
        """
        Count double bonds in a cycle.
        """

        return sum(
            bond.order == 2
            for bond in self.cycle_bonds(cycle)
        )

    def triple_bonds_in_cycle(self, cycle):
        """
        Count triple bonds in a cycle.
        """

        return sum(
            bond.order == 3
            for bond in self.cycle_bonds(cycle)
        )

    def cycle_pi_bond_count(self, atom_id, cycle):
        """
        Count pi bonds connected to atom_id that belong
        specifically to the supplied cycle.
        """

        cycle = list(cycle)
        count = 0

        for i in range(len(cycle)):

            atom1 = cycle[i]

            atom2 = cycle[
                (i + 1) % len(cycle)
            ]

            if (
                atom_id != atom1
                and atom_id != atom2
            ):
                continue

            for neighbour, bond in self.graph[atom1]:

                if neighbour == atom2:

                    count += max(
                        0,
                        bond.order - 1
                    )

                    break

        return count

    # ==========================================================
    # RING ANALYSIS
    # ==========================================================

    def analyze_rings(self):
        """
        Return information about every detected cycle.
        """

        cycles = self.find_cycles()
        results = []

        for cycle in cycles:

            size = len(cycle)

            double_bonds = (
                self.double_bonds_in_cycle(
                    cycle
                )
            )

            triple_bonds = (
                self.triple_bonds_in_cycle(
                    cycle
                )
            )

            results.append({
                "atoms": cycle,
                "size": size,
                "double_bonds": double_bonds,
                "triple_bonds": triple_bonds,
            })

        return results

    def cycle_rank(self):
        """
        Calculate the cyclomatic number:

            μ = E - V + C

        where:

            E = number of bonds
            V = number of atoms
            C = number of connected components

        A cycle rank greater than zero indicates rings/cycles.
        """

        E = len(self.molecule.bonds)
        V = len(self.molecule.atoms)
        C = len(self.connected_components())

        return E - V + C

    def has_rings(self):
        return self.cycle_rank() > 0

    # ==========================================================
    # ELECTRONIC INFORMATION
    # ==========================================================

    def electronic_analyzer(self):
        """
        Return an ElectronicAnalyzer for the molecule.
        """

        return ElectronicAnalyzer(
            self.molecule
        )

    def atom_has_lone_pair(self, atom_id):
        """
        Check whether an atom has at least one lone pair.
        """

        electronic = self.electronic_analyzer()

        return (
            electronic.total_lone_pairs(
                atom_id
            ) > 0
        )

    # ==========================================================
    # CONJUGATION
    # ==========================================================

    def cycle_is_conjugated(self, cycle):
        """
        Estimate whether every atom in a cycle can
        participate in a conjugated cyclic system.

        This is a first-pass symbolic rule.
        """

        electronic = self.electronic_analyzer()

        for atom_id in cycle:

            atom = self.molecule.get_atom(
                atom_id
            )

            # Atom already participates in a pi bond
            # inside this cycle.
            cycle_pi = (
                self.cycle_pi_bond_count(
                    atom_id,
                    cycle
                )
            )

            if cycle_pi > 0:
                continue

            # Radical may occupy a p orbital.
            if atom.radical_electrons > 0:
                continue

            # Positive charge may provide an empty p orbital.
            if atom.charge > 0:
                continue

            # A lone pair may participate.
            if (
                electronic.total_lone_pairs(
                    atom_id
                ) > 0
            ):
                continue

            # Otherwise this atom breaks conjugation.
            return False

        return True

    def conjugated_bonds(self):
        """
        Find simple conjugation candidates.

        A double bond is considered a candidate when another
        double bond is directly connected to either endpoint.
        """

        result = []

        for bond in self.molecule.bonds:

            if bond.order != 2:
                continue

            atom1 = self.molecule.get_atom(
                bond.atom1
            )

            atom2 = self.molecule.get_atom(
                bond.atom2
            )

            connected_bonds = (
                self.molecule.bonds_of(
                    atom1.id
                )
                + self.molecule.bonds_of(
                    atom2.id
                )
            )

            for neighbouring_bond in connected_bonds:

                if neighbouring_bond.order == 2:
                    result.append(bond)
                    break

        unique = []
        seen = set()

        for bond in result:

            key = tuple(
                sorted(
                    (bond.atom1, bond.atom2)
                )
            )

            if key not in seen:
                seen.add(key)
                unique.append(bond)

        return unique

    # ==========================================================
    # PI ELECTRON COUNTING
    # ==========================================================

    def cycle_pi_electron_contribution(
        self,
        atom_id,
        cycle
    ):
        """
        Estimate the pi-electron contribution of one atom
        to one cyclic conjugated system.
        """

        electronic = self.electronic_analyzer()

        atom = self.molecule.get_atom(
            atom_id
        )

        # Atom involved in a ring pi bond.
        cycle_pi = (
            self.cycle_pi_bond_count(
                atom_id,
                cycle
            )
        )

        if cycle_pi > 0:
            return 1

        # Radical contribution.
        if atom.radical_electrons > 0:
            return atom.radical_electrons

        # Positive charge: empty p orbital.
        if atom.charge > 0:
            return 0

        # Lone pair contribution.
        if (
            electronic.total_lone_pairs(
                atom_id
            ) > 0
        ):
            return 2

        return 0

    def pi_electron_count(self, cycle):
        """
        Count estimated pi electrons in a cycle.
        """

        electrons = 0

        for atom_id in cycle:

            electrons += (
                self.cycle_pi_electron_contribution(
                    atom_id,
                    cycle
                )
            )

        return electrons

    # ==========================================================
    # AROMATICITY
    # ==========================================================

    def is_aromatic_ring(self, cycle):
        """
        Estimate aromaticity using:

            1. Ring size >= 3
            2. Continuous conjugation
            3. Hückel rule: 4n + 2 pi electrons
        """

        if len(cycle) < 3:
            return False

        if not self.cycle_is_conjugated(cycle):
            return False

        pi_electrons = (
            self.pi_electron_count(
                cycle
            )
        )

        if pi_electrons < 2:
            return False

        # Hückel rule:
        #
        # π electrons = 4n + 2
        #
        # Therefore:
        #
        # (π - 2) mod 4 = 0

        if (
            (pi_electrons - 2) % 4 != 0
        ):
            return False

        return True

    def aromatic_rings(self):
        """
        Return only rings classified as aromatic.
        """

        result = []

        for ring in self.analyze_rings():

            if self.is_aromatic_ring(
                ring["atoms"]
            ):
                result.append(ring)

        return result

    # ==========================================================
    # HYBRIDIZATION
    # ==========================================================

    def hybridization(self, atom_id):
        """
        Approximate hybridization.

        Rules:

            Two or more pi bonds -> sp
            One pi bond         -> sp2
            Otherwise            -> sp3 for common atoms

        This is graph-based and approximate.
        """

        atom = self.molecule.get_atom(
            atom_id
        )

        pi_bonds = self.pi_bond_count(
            atom_id
        )

        sigma_bonds = self.sigma_bond_count(
            atom_id
        )

        if pi_bonds >= 2:
            return "sp"

        if pi_bonds == 1:
            return "sp2"

        if atom.element in {
            "C",
            "Si",
        }:
            if sigma_bonds <= 4:
                return "sp3"

        if atom.element in {
            "N",
            "P",
        }:
            if sigma_bonds <= 3:
                return "sp3"

        if atom.element in {
            "O",
            "S",
        }:
            if sigma_bonds <= 2:
                return "sp3"

        if atom.element == "H":
            return "s"

        return "unknown"

    # ==========================================================
    # LOCAL GEOMETRY
    # ==========================================================

    def geometry_type(self, atom_id):
        """
        Approximate electron-domain geometry.

        Note:
            sp3 does not always mean the visible molecular
            shape is tetrahedral. Lone pairs can make the
            molecular shape bent or pyramidal.
        """

        atom = self.molecule.get_atom(
            atom_id
        )

        hybridization = self.hybridization(
            atom_id
        )

        if hybridization == "sp":
            return "linear"

        if hybridization == "sp2":
            return "trigonal planar"

        if hybridization == "sp3":

            lone_pairs = (
                self.electronic_analyzer()
                .total_lone_pairs(atom_id)
            )

            if atom.element == "O" and lone_pairs >= 2:
                return "bent / tetrahedral electron-domain"

            if atom.element == "N" and lone_pairs >= 1:
                return "trigonal pyramidal / tetrahedral electron-domain"

            return "tetrahedral electron-domain"

        if hybridization == "s":
            return "spherical"

        return "unknown"

    # ==========================================================
    # VALENCE ANALYSIS
    # ==========================================================

    def valence(self, atom_id):
        """
        Calculate simple valence:

            valence = bond-order sum + explicit hydrogens
        """

        atom = self.molecule.get_atom(
            atom_id
        )

        return (
            self.bond_order_sum(atom_id)
            + atom.hydrogens
        )

    def valence_status(self, atom_id):
        """
        Compare calculated valence against common valence rules.
        """

        atom = self.molecule.get_atom(
            atom_id
        )

        actual = self.valence(
            atom_id
        )

        allowed = self.COMMON_VALENCES.get(
            atom.element
        )

        if allowed is None:

            return {
                "valid": True,
                "status": "no rule",
                "actual": actual,
                "allowed": None,
            }

        valid = actual in allowed

        return {
            "valid": valid,
            "status": (
                "valid"
                if valid
                else "warning"
            ),
            "actual": actual,
            "allowed": sorted(allowed),
        }

    def all_valence_status(self):
        """
        Return valence information for every atom.
        """

        return {
            atom_id: self.valence_status(
                atom_id
            )
            for atom_id in self.molecule.atoms
        }

    # ==========================================================
    # ATOM REPORT
    # ==========================================================

    def atom_summary(self, atom_id):
        """
        Return a detailed symbolic report for one atom.
        """

        atom = self.molecule.get_atom(
            atom_id
        )

        return {
            "id": atom.id,
            "element": atom.element,
            "hydrogens": atom.hydrogens,
            "charge": atom.charge,
            "radical_electrons": (
                atom.radical_electrons
            ),
            "lone_pairs": (
                atom.lone_pairs
            ),
            "degree": self.degree(
                atom_id
            ),
            "sigma_bonds": self.sigma_bond_count(
                atom_id
            ),
            "pi_bonds": self.pi_bond_count(
                atom_id
            ),
            "bond_order_sum": self.bond_order_sum(
                atom_id
            ),
            "valence": self.valence(
                atom_id
            ),
            "hybridization": self.hybridization(
                atom_id
            ),
            "geometry": self.geometry_type(
                atom_id
            ),
            "valence_status": self.valence_status(
                atom_id
            ),
        }

    # ==========================================================
    # COMPLETE SUMMARY
    # ==========================================================

    def summary(self):
        """
        Return a complete molecular analysis dictionary.
        """

        rings = self.analyze_rings()
        aromatic = self.aromatic_rings()

        atoms = {}

        for atom_id in self.molecule.atoms:

            atoms[atom_id] = self.atom_summary(
                atom_id
            )

        return {
            "formula": self.molecule.formula(),
            "total_charge": self.molecule.total_charge(),
            "atom_count": len(
                self.molecule.atoms
            ),
            "bond_count": len(
                self.molecule.bonds
            ),
            "connected_components": (
                self.connected_components()
            ),
            "cycle_rank": self.cycle_rank(),
            "has_rings": self.has_rings(),
            "rings": rings,
            "aromatic_rings": aromatic,
            "atoms": atoms,
        }

    # ==========================================================
    # PRINT SUMMARY
    # ==========================================================

    def print_summary(self):
        """
        Print a readable structural report.
        """

        data = self.summary()

        print("STRUCTURAL ANALYSIS")
        print("=" * 60)

        print(
            f"Formula: {data['formula']}"
        )

        print(
            f"Total charge: {data['total_charge']}"
        )

        print(
            f"Atoms: {data['atom_count']}"
        )

        print(
            f"Bonds: {data['bond_count']}"
        )

        print(
            f"Connected components: "
            f"{data['connected_components']}"
        )

        print(
            f"Cycle rank: {data['cycle_rank']}"
        )

        print(
            f"Has rings: {data['has_rings']}"
        )

        print(
            f"Rings: {len(data['rings'])}"
        )

        print(
            f"Aromatic rings: "
            f"{len(data['aromatic_rings'])}"
        )

        print()

        for atom_id, atom_data in data["atoms"].items():

            valence_data = (
                atom_data["valence_status"]
            )

            print(
                f"Atom {atom_id}: "
                f"{atom_data['element']} | "
                f"degree={atom_data['degree']} | "
                f"sigma={atom_data['sigma_bonds']} | "
                f"pi={atom_data['pi_bonds']} | "
                f"valence={atom_data['valence']} | "
                f"hybridization="
                f"{atom_data['hybridization']} | "
                f"geometry="
                f"{atom_data['geometry']} | "
                f"status="
                f"{valence_data['status']}"
            )

        if data["rings"]:

            print()
            print("RING DETAILS")
            print("-" * 60)

            for index, ring in enumerate(
                data["rings"],
                start=1
            ):

                print(
                    f"Ring {index}: "
                    f"atoms={ring['atoms']} | "
                    f"size={ring['size']} | "
                    f"double bonds="
                    f"{ring['double_bonds']} | "
                    f"triple bonds="
                    f"{ring['triple_bonds']}"
                )

        if data["aromatic_rings"]:

            print()
            print("AROMATIC RINGS")
            print("-" * 60)

            for index, ring in enumerate(
                data["aromatic_rings"],
                start=1
            ):

                pi_electrons = (
                    self.pi_electron_count(
                        ring["atoms"]
                    )
                )

                print(
                    f"Aromatic ring {index}: "
                    f"atoms={ring['atoms']} | "
                    f"pi electrons="
                    f"{pi_electrons}"
                )

    # Compatibility alias.
    def print_report(self):
        self.print_summary()