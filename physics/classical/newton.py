"""
Newtonian mechanics.
"""

from __future__ import annotations

from physics.particles import Particle


def acceleration(particle: Particle):
    """
    Calculate acceleration from Newton's second law.

        F = ma

        a = F / m
    """

    return particle.force / particle.mass


def update_particle(particle: Particle, dt: float) -> None:
    """
    Advance a particle by one timestep using explicit Euler integration.
    """

    if dt <= 0:
        raise ValueError("Timestep must be positive.")

    a = acceleration(particle)

    particle.velocity = particle.velocity + a * dt
    particle.position = particle.position + particle.velocity * dt