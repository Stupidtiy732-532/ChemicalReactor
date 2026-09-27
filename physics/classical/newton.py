from __future__ import annotations

from physics.particles import Particle


def acceleration(particle: Particle):
    """
    Newton's second law:

        F = ma

        a = F / m
    """

    return particle.force / particle.mass