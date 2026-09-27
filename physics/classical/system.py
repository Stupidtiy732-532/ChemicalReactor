from __future__ import annotations

from engine.state import SimulationState
from physics.classical.newton import update_particle


def newtonian_step(state: SimulationState, dt: float) -> None:
    for particle in state.particles:
        update_particle(particle, dt)