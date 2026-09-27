from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class SimulationState:
    time: float = 0.0
    particles: list[Any] = field(default_factory=list)
    fields: dict[str, Any] = field(default_factory=dict)
    systems: dict[str, Any] = field(default_factory=dict)

    def copy(self) -> "SimulationState":
        return SimulationState(
            time=self.time,
            particles=list(self.particles),
            fields=dict(self.fields),
            systems=dict(self.systems),
        )