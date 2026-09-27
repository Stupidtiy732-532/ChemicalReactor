"""
Physical particle representations.
"""

from __future__ import annotations

from physics.vectors import Vector3


class Particle:
    def __init__(
        self,
        mass: float,
        charge: float = 0.0,
        position: Vector3 | None = None,
        velocity: Vector3 | None = None,
    ) -> None:
        if mass <= 0:
            raise ValueError("Particle mass must be positive.")

        self.mass = float(mass)
        self.charge = float(charge)

        self.position = (
            position.copy()
            if position is not None
            else Vector3()
        )

        self.velocity = (
            velocity.copy()
            if velocity is not None
            else Vector3()
        )

        self.force = Vector3()

    @property
    def momentum(self) -> Vector3:
        return self.velocity * self.mass

    def clear_force(self) -> None:
        self.force = Vector3()

    def apply_force(self, force: Vector3) -> None:
        self.force = self.force + force

    def __repr__(self) -> str:
        return (
            f"Particle("
            f"mass={self.mass}, "
            f"charge={self.charge}, "
            f"position={self.position}, "
            f"velocity={self.velocity}"
            f")"
        )