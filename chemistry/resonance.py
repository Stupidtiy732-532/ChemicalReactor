from dataclasses import dataclass


@dataclass(frozen=True)
class ResonanceArrangement:
    bond_orders: tuple[float, ...]
    lone_pairs: tuple[int, ...]
    charges: tuple[int, ...]
    unpaired_electrons: tuple[int, ...]


@dataclass
class ResonanceSystem:

    atom_ids: tuple[int, ...]
    bond_ids: tuple[int, ...]
    arrangements: list[ResonanceArrangement]

    def add_arrangement(
        self,
        arrangement: ResonanceArrangement,
    ):

        if len(arrangement.bond_orders) != len(self.bond_ids):
            raise ValueError(
                "Bond-order count does not match connectivity"
            )

        if len(arrangement.lone_pairs) != len(self.atom_ids):
            raise ValueError(
                "Lone-pair count does not match atom count"
            )

        if len(arrangement.charges) != len(self.atom_ids):
            raise ValueError(
                "Charge count does not match atom count"
            )

        if len(arrangement.unpaired_electrons) != len(self.atom_ids):
            raise ValueError(
                "Unpaired-electron count does not match atom count"
            )

        self.arrangements.append(arrangement)