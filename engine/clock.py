from __future__ import annotations

from dataclasses import dataclass


@dataclass
class SimulationClock:
    time: float = 0.0
    timestep: float = 1e-15
    speed: float = 1.0
    paused: bool = False

    def advance(self, real_dt: float) -> float:
        if self.paused:
            return 0.0

        dt = real_dt * self.speed
        self.time += dt
        return dt

    def step(self) -> float:
        if self.paused:
            return 0.0

        self.time += self.timestep
        return self.timestep

    def skip(self, dt: float) -> None:
        self.time += dt

    def pause(self) -> None:
        self.paused = True

    def resume(self) -> None:
        self.paused = False

    def set_speed(self, speed: float) -> None:
        if speed < 0:
            raise ValueError("Simulation speed cannot be negative.")

        self.speed = speed

    def set_timestep(self, timestep: float) -> None:
        if timestep <= 0:
            raise ValueError("Timestep must be positive.")

        self.timestep = timestep