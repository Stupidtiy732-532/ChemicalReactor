"""
Electrostatic interaction between point charges.
"""

from __future__ import annotations

from physics.constants import VACUUM_PERMITTIVITY, PI
from physics.particles import Particle
from physics.vectors import Vector3


COULOMB_CONSTANT = 1.0 / (4.0 * PI * VACUUM_PERMITTIVITY)


def force(source: Particle, target: Particle) -> Vector3:
    """
    Calculate the electrostatic force exerted by source on target.

    F = k q1 q2 / r^3 * r_vector

    where

        r_vector = r_target - r_source
    """

    displacement = target.position - source.position
    distance_squared = displacement.magnitude_squared()

    if distance_squared == 0.0:
        raise ValueError(
            "Coulomb force is undefined for coincident point charges."
        )

    distance = distance_squared ** 0.5

    magnitude = (
        COULOMB_CONSTANT
        * source.charge
        * target.charge
        / distance_squared
    )

    return displacement.normalized() * magnitude