from __future__ import annotations

from engine.state import SimulationState
from physics.classical.newton import update_particle
from physics.electromagnetism.coulomb import force as coulomb_force


def calculate_forces(state: SimulationState) -> None:
    for particle in state.particles:
        particle.clear_force()

    particles = state.particles

    for i, particle_a in enumerate(particles):
        for particle_b in particles[i + 1:]:
            force_on_b = coulomb_force(
                source=particle_a,
                target=particle_b,
            )

            particle_b.apply_force(force_on_b)
            particle_a.apply_force(-force_on_b)


def newtonian_electromagnetic_step(
    state: SimulationState,
    dt: float,
) -> None:
    calculate_forces(state)

    for particle in state.particles:
        update_particle(particle, dt)