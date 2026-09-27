from __future__ import annotations

from collections.abc import Callable

from .clock import SimulationClock
from .state import SimulationState


PhysicsStep = Callable[[SimulationState, float], None]


class Simulation:
    def __init__(
        self,
        timestep: float = 1e-15,
        speed: float = 1.0,
    ) -> None:
        self.clock = SimulationClock(
            timestep=timestep,
            speed=speed,
        )

        self.state = SimulationState()
        self.physics: list[PhysicsStep] = []

    def add_physics(self, physics: PhysicsStep) -> None:
        self.physics.append(physics)

    def step(self) -> None:
        dt = self.clock.step()

        if dt == 0.0:
            return

        for physics in self.physics:
            physics(self.state, dt)

        self.state.time = self.clock.time

    def skip(self, dt: float) -> None:
        self.clock.skip(dt)
        self.state.time = self.clock.time

    def pause(self) -> None:
        self.clock.pause()

    def resume(self) -> None:
        self.clock.resume()

    def set_speed(self, speed: float) -> None:
        self.clock.set_speed(speed)

    def run_steps(self, steps: int) -> None:
        for _ in range(steps):
            if self.clock.paused:
                break

            self.step()