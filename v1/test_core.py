from parser import parse_structure
from geometry import GeometryEngine
from analyzer import MolecularAnalyzer


def test():
    molecule = parse_structure("CH3-CH2-OH")

    geometry = GeometryEngine(molecule)
    geometry.generate_layout()

    print(molecule)
    print()

    print(geometry.export_xyz())
    print()

    analyzer = MolecularAnalyzer(molecule)
    analyzer.print_summary()


if __name__ == "__main__":
    test()