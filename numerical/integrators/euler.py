from __future__ import annotations

from physics.particles import Particle


def integrate_particle(
    particle: Particle,
    dt: float,
) -> None:
    if dt <= 0:
        raise ValueError("Timestep must be positive.")

    acceleration = particle.force / particle.mass

    particle.velocity = (
        particle.velocity
        + acceleration * dt
    )

    particle.position = (
        particle.position
        + particle.velocity * dt
    )