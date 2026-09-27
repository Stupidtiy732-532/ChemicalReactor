from __future__ import annotations

from collections.abc import Callable

from engine.state import SimulationState


ForceCalculator = Callable[[SimulationState], None]


def integrate(
    state: SimulationState,
    dt: float,
    calculate_forces: ForceCalculator,
) -> None:
    """
    Velocity Verlet integration for the complete physical system.

    The force calculator must update particle.force for every particle.

    Algorithm:

        v(t + dt/2) = v(t) + a(t) dt/2

        x(t + dt) = x(t) + v(t + dt/2) dt

        calculate F(t + dt)

        v(t + dt) = v(t + dt/2) + a(t + dt) dt/2
    """

    if dt <= 0:
        raise ValueError("Timestep must be positive.")

    # First half-step in velocity.
    for particle in state.particles:
        acceleration = particle.force / particle.mass

        particle.velocity = (
            particle.velocity
            + acceleration * (0.5 * dt)
        )

    # Full position step.
    for particle in state.particles:
        particle.position = (
            particle.position
            + particle.velocity * dt
        )

    # Forces must be recalculated for the entire system.
    calculate_forces(state)

    # Second half-step in velocity.
    for particle in state.particles:
        acceleration = particle.force / particle.mass

        particle.velocity = (
            particle.velocity
            + acceleration * (0.5 * dt)
        )