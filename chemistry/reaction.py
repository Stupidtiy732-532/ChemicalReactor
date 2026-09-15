from dataclasses import dataclass
from typing import Callable

from chemistry.molecule import Molecule


@dataclass(frozen=True)
class ReactionResult:
    reaction_id: str
    reactants: tuple[Molecule, ...]
    products: tuple[Molecule, ...]
    observations: tuple[str, ...] = ()


@dataclass(frozen=True)
class Reaction:
    id: str
    name: str
    reagent: str
    conditions: tuple[str, ...]
    apply: Callable[..., ReactionResult]