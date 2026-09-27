from engine.simulation import Simulation
from physics.classical.system import newtonian_electromagnetic_step
from physics.constants import (
    ELECTRON_MASS,
    PROTON_MASS,
    ELEMENTARY_CHARGE,
)
from physics.particles import Particle
from physics.vectors import Vector3


def main() -> None:
    simulation = Simulation(
        timestep=1e-20,
        speed=1.0,
    )

    simulation.add_physics(
        newtonian_electromagnetic_step
    )

    electron = Particle(
        mass=ELECTRON_MASS,
        charge=-ELEMENTARY_CHARGE,
        position=Vector3(0.0, 0.0, 0.0),
        velocity=Vector3(),
    )

    proton = Particle(
        mass=PROTON_MASS,
        charge=ELEMENTARY_CHARGE,
        position=Vector3(1e-9, 0.0, 0.0),
        velocity=Vector3(),
    )

    simulation.state.particles.extend(
        [electron, proton]
    )

    print("ChemicalReactor physics core")
    print("=" * 40)

    print(f"Initial electron position: {electron.position}")
    print(f"Initial proton position:   {proton.position}")

    simulation.run_steps(1000)

    print()
    print(f"Simulation time: {simulation.state.time:.3e} s")
    print(f"Electron position: {electron.position}")
    print(f"Electron velocity: {electron.velocity}")
    print(f"Proton position:   {proton.position}")
    print(f"Proton velocity:   {proton.velocity}")


if __name__ == "__main__":
    main()