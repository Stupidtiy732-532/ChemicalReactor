from engine.simulation import Simulation
from physics.constants import ELECTRON_MASS, ELEMENTARY_CHARGE
from physics.particles import Particle
from physics.vectors import Vector3


def main() -> None:
    simulation = Simulation(
        timestep=1e-18,
        speed=1.0,
    )

    electron = Particle(
        mass=ELECTRON_MASS,
        charge=-ELEMENTARY_CHARGE,
        position=Vector3(0.0, 0.0, 0.0),
        velocity=Vector3(1.0e6, 0.0, 0.0),
    )

    simulation.state.particles.append(electron)

    print("ChemicalReactor physics core")
    print("=" * 40)

    print(f"Initial position: {electron.position}")
    print(f"Initial velocity: {electron.velocity}")
    print(f"Initial momentum: {electron.momentum}")

    print(f"Simulation time: {simulation.state.time:.3e} s")


if __name__ == "__main__":
    main()