from engine.simulation import Simulation


def main() -> None:
    simulation = Simulation(
        timestep=1e-15,
        speed=1.0,
    )

    print(f"Initial time: {simulation.state.time:.3e} s")

    simulation.run_steps(10)

    print(f"After 10 steps: {simulation.state.time:.3e} s")

    simulation.skip(1e-12)

    print(f"After skip: {simulation.state.time:.3e} s")

    simulation.pause()

    simulation.step()

    print(f"After pause: {simulation.state.time:.3e} s")

    simulation.resume()
    simulation.step()

    print(f"After resume: {simulation.state.time:.3e} s")


if __name__ == "__main__":
    main()