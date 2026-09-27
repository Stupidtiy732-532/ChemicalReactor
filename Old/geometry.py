import math
from molecule import Vector3D
from collections import deque


class GeometryEngine:

    # SI unit: metre
    DEFAULT_BOND_LENGTH = 1.5e-10

    TETRAHEDRAL_ANGLE = math.radians(109.5)
    TRIGONAL_ANGLE = math.radians(120.0)
    LINEAR_ANGLE = math.radians(180.0)

    def __init__(self, molecule):
        self.molecule = molecule

    # ==================================================
    # BASIC VECTOR OPERATIONS
    # ==================================================

    @staticmethod
    def add(a, b):
        return Vector3D(
            a.x + b.x,
            a.y + b.y,
            a.z + b.z
        )

    @staticmethod
    def scale(vector, factor):
        return Vector3D(
            vector.x * factor,
            vector.y * factor,
            vector.z * factor
        )

    # ==================================================
    # CONNECTED COMPONENTS
    # ==================================================

    def connected_components(self):

        remaining = set(
            self.molecule.atoms.keys()
        )

        components = []

        while remaining:

            start = next(iter(remaining))
            component = set()
            stack = [start]

            while stack:

                atom_id = stack.pop()

                if atom_id in component:
                    continue

                component.add(atom_id)
                remaining.discard(atom_id)

                for neighbor in self.molecule.neighbors(
                    atom_id
                ):
                    if neighbor not in component:
                        stack.append(neighbor)

            components.append(
                sorted(component)
            )

        return components

    # ==================================================
    # APPROXIMATE 3D LAYOUT
    # ==================================================

    def generate_layout(self):

        if not self.molecule.atoms:
            return

        components = self.connected_components()

        component_offset = 0.0

        for component in components:

            self._layout_component(
                component,
                component_offset
            )

            component_offset += (
                len(component)
                * self.DEFAULT_BOND_LENGTH
                * 2.5
            )

    # ==================================================

    def _layout_component(
        self,
        component,
        offset
    ):

        if not component:
            return

        root = component[0]

        placed = {
            root
        }

        self.molecule.get_atom(
            root
        ).position = Vector3D(
            offset,
            0.0,
            0.0
        )

        queue = deque([root])

        while queue:
            current = queue.popleft()

            current_atom = self.molecule.get_atom(
                current
            )

            neighbors = [
                neighbor
                for neighbor in self.molecule.neighbors(
                    current
                )
                if neighbor in component
            ]

            unplaced = [
                neighbor
                for neighbor in neighbors
                if neighbor not in placed
            ]

            if not unplaced:
                continue

            directions = self._directions(
                len(unplaced)
            )

            for index, neighbor in enumerate(
                unplaced
            ):

                direction = directions[
                    index % len(directions)
                ]

                neighbor_atom = (
                    self.molecule.get_atom(
                        neighbor
                    )
                )

                neighbor_atom.position = self.add(
                    current_atom.position,
                    self.scale(
                        direction,
                        self.DEFAULT_BOND_LENGTH
                    )
                )

                placed.add(neighbor)
                queue.append(neighbor)

    # ==================================================

    def _directions(self, count):
        if count == 1:
            return [Vector3D(1.0, 0.0, 0.0)]

        if count == 2:
            angle = math.radians(120.0)
            return [
                Vector3D(1.0, 0.0, 0.0),
                Vector3D(math.cos(angle), math.sin(angle), 0.0),
            ]

        if count == 3:
            y = math.sqrt(3) / 2
            return [
                Vector3D(1.0, 0.0, 0.0),
                Vector3D(-0.5, y, 0.0),
                Vector3D(-0.5, -y, 0.0),
            ]

        return [
            Vector3D(1.0, 1.0, 1.0),
            Vector3D(1.0, -1.0, -1.0),
            Vector3D(-1.0, 1.0, -1.0),
            Vector3D(-1.0, -1.0, 1.0),
        ]

    # ==================================================
    # MEASUREMENTS
    # ==================================================

    def bond_length(self, atom1, atom2):

        a = self.molecule.get_atom(atom1)
        b = self.molecule.get_atom(atom2)

        return a.position.distance_to(
            b.position
        )

    # ==================================================

    def angle(self, atom1, atom2, atom3):

        a = self.molecule.get_atom(atom1).position
        b = self.molecule.get_atom(atom2).position
        c = self.molecule.get_atom(atom3).position

        v1 = (
            a.x - b.x,
            a.y - b.y,
            a.z - b.z
        )

        v2 = (
            c.x - b.x,
            c.y - b.y,
            c.z - b.z
        )

        dot = sum(
            x * y
            for x, y in zip(v1, v2)
        )

        magnitude_1 = math.sqrt(
            sum(x * x for x in v1)
        )

        magnitude_2 = math.sqrt(
            sum(x * x for x in v2)
        )

        if magnitude_1 == 0 or magnitude_2 == 0:
            raise ValueError(
                "Cannot calculate angle "
                "using zero-length vector."
            )

        cosine = dot / (
            magnitude_1 * magnitude_2
        )

        cosine = max(
            -1.0,
            min(1.0, cosine)
        )

        return math.acos(cosine)

    # ==================================================

    def angle_degrees(
        self,
        atom1,
        atom2,
        atom3
    ):

        return math.degrees(
            self.angle(
                atom1,
                atom2,
                atom3
            )
        )

    # ==================================================
    # OUTPUT
    # ==================================================

    def coordinates(self):

        return {
            atom.id: atom.position.as_tuple()
            for atom in self.molecule.atoms.values()
        }

    # ==================================================

    def export_xyz(self):

        # XYZ uses Å by convention.
        angstrom = 1e10

        lines = [
            str(len(self.molecule.atoms)),
            self.molecule.formula()
        ]

        for atom in self.molecule.atoms.values():

            x, y, z = atom.position.as_tuple()

            lines.append(
                f"{atom.element} "
                f"{x * angstrom:.6f} "
                f"{y * angstrom:.6f} "
                f"{z * angstrom:.6f}"
            )

        return "\n".join(lines)