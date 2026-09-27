"""
Basic three-dimensional vector mathematics.

Vectors represent physical spatial quantities such as:
position, velocity, momentum, force, electric field, and magnetic field.
"""

from __future__ import annotations

from math import sqrt


class Vector3:
    __slots__ = ("x", "y", "z")

    def __init__(
        self,
        x: float = 0.0,
        y: float = 0.0,
        z: float = 0.0,
    ) -> None:
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def __add__(self, other: Vector3) -> Vector3:
        return Vector3(
            self.x + other.x,
            self.y + other.y,
            self.z + other.z,
        )

    def __sub__(self, other: Vector3) -> Vector3:
        return Vector3(
            self.x - other.x,
            self.y - other.y,
            self.z - other.z,
        )

    def __mul__(self, scalar: float) -> Vector3:
        return Vector3(
            self.x * scalar,
            self.y * scalar,
            self.z * scalar,
        )

    def __rmul__(self, scalar: float) -> Vector3:
        return self * scalar

    def __truediv__(self, scalar: float) -> Vector3:
        if scalar == 0:
            raise ZeroDivisionError("Cannot divide a vector by zero.")

        return Vector3(
            self.x / scalar,
            self.y / scalar,
            self.z / scalar,
        )

    def __neg__(self) -> Vector3:
        return Vector3(
            -self.x,
            -self.y,
            -self.z,
        )

    def dot(self, other: Vector3) -> float:
        return (
            self.x * other.x
            + self.y * other.y
            + self.z * other.z
        )

    def cross(self, other: Vector3) -> Vector3:
        return Vector3(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )

    def magnitude_squared(self) -> float:
        return self.dot(self)

    def magnitude(self) -> float:
        return sqrt(self.magnitude_squared())

    def normalized(self) -> Vector3:
        magnitude = self.magnitude()

        if magnitude == 0:
            raise ValueError("Cannot normalize the zero vector.")

        return self / magnitude

    def copy(self) -> Vector3:
        return Vector3(self.x, self.y, self.z)

    def __repr__(self) -> str:
        return (
            f"Vector3("
            f"x={self.x}, "
            f"y={self.y}, "
            f"z={self.z}"
            f")"
        )