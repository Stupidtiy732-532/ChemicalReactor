from dataclasses import dataclass


@dataclass(frozen=True)
class ElectronState:
    lone_pairs: int = 0
    unpaired_electrons: int = 0
    charge: int = 0

    @property
    def nonbonding_electrons(self) -> int:
        return (
            2 * self.lone_pairs
            + self.unpaired_electrons
        )