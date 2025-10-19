from typing import Dict, List, NamedTuple

import argparse
import math


class Isotope(NamedTuple):
    mass: float
    abundance: float


class Atom(NamedTuple):
    symbol: str
    isotopes: List[Isotope]


# Isotope data from NIST (https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl)
ATOMS: Dict[str, Atom] = {
    "H": Atom(
        symbol="H",
        isotopes=[
            Isotope(mass=1.00782503223, abundance=0.999885),
            Isotope(mass=2.01410177812, abundance=0.000115),
        ]
    ),
    
    "C": Atom(
        symbol="C",
        isotopes=[
            Isotope(mass=12.0000000000, abundance=0.9893),
            Isotope(mass=13.0033548351, abundance=0.0107),
        ]
    ),
    
    "N": Atom(
        symbol="N",
        isotopes=[
            Isotope(mass=14.0030740044, abundance=0.99636),
            Isotope(mass=15.0001088989, abundance=0.00364),
        ]
    ),
    
    "O": Atom(
        symbol="O",
        isotopes=[
            Isotope(mass=15.9949146196, abundance=0.99757),
            Isotope(mass=16.9991317565, abundance=0.00038),
            Isotope(mass=17.9991596129, abundance=0.00205),
        ]
    ),
    
    "S": Atom(
        symbol="S",
        isotopes=[
            Isotope(mass=31.9720711744, abundance=0.9499),
            Isotope(mass=32.9714589098, abundance=0.0075),
            Isotope(mass=33.9678670040, abundance=0.0425),
            # We ignore the very rare isotope S-36
        ]
    ),
    
    "P": Atom(
        symbol="P",
        isotopes=[
            Isotope(mass=30.9737619984, abundance=1.0000),
        ]
    ),
}

AVERAGINE_COMPOSITION: Dict[str, float] = {
    "C": 4.9384,
    "H": 7.7583,
    "N": 1.3577,
    "O": 1.4773,
    "S": 0.0417,
}


class BinomialCalculator:
    def __init__(self, n: int):
        self.log_factorials = []

        # log(0!) = log(1) = 0
        self.log_factorials.append(0.0)

        for i in range(1, n + 1):
            # log(n!) = log((n - 1)! * n) = log((n - 1)!) + log(n)
            self.log_factorials.append(self.log_factorials[-1] + math.log(i))

    def cnk(self, n: int, k: int) -> float:
        if n > len(self.log_factorials) - 1:
            raise ValueError("n is too large for precomputed factorials")
        
        if k < 0 or k > n:
            return 0.0
        
        log_cnk = self.log_factorials[n] - self.log_factorials[k] - self.log_factorials[n - k]
        return math.exp(log_cnk)
    
    def binomial_probability(self, n: int, k: int, p: float) -> float:
        cnk = self.cnk(n, k)
        return cnk * (p ** k) * ((1 - p) ** (n - k))
    

class Isotoper:
    def __init__(self, n: int):
        self.binomial_calculator = BinomialCalculator(n)

    def get_element_distribution(self, element: str, n: int, max_neutrons: int) -> List[float]:
        atom = ATOMS[element]

        distribution = [0.0 for _ in range(max_neutrons + 1)]

        # Here we manually unroll recursion for up to 3 isotopes
        assert len(atom.isotopes) <= 3
        for plus_1 in range(0, min(n, max_neutrons) + 1):
            p1 = atom.isotopes[1].abundance
            prob1 = self.binomial_calculator.binomial_probability(n, plus_1, p1)

            if len(atom.isotopes) == 2:
                distribution[plus_1] += prob1
                continue

            for plus_2 in range(0, min(n - plus_1, (max_neutrons - plus_1) // 2) + 1):
                p2 = atom.isotopes[2].abundance / (1 - p1)
                prob2 = self.binomial_calculator.binomial_probability(n - plus_1, plus_2, p2)

                distribution[plus_1 + 2 * plus_2] += prob1 * prob2

        return distribution
    
    def get_molecule_distribution(self, composition: Dict[str, int], max_neutrons: int) -> List[float]:
        total_distribution = [1.0] + [0.0] * max_neutrons

        for element, count in composition.items():
            element_distribution = self.get_element_distribution(element, count, max_neutrons)
            total_distribution = self._convolve(total_distribution, element_distribution)

        return total_distribution

    def _convolve(self, dist1: List[float], dist2: List[float]) -> List[float]:
        n = len(dist1)
        assert len(dist2) == n

        result = [0.0] * n

        for i in range(len(dist1)):
            for j in range(len(dist2)):
                if i + j < n:
                    result[i + j] += dist1[i] * dist2[j]

        return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate averagine-related isotopic distributions")
    parser.add_argument("--min_aas", type=int, required=True, help="Minimum number of averagines")
    parser.add_argument("--max_aas", type=int, required=True, help="Maximum number of averagines")
    parser.add_argument("--max_neutrons", type=int, required=True, help="Maximum number of additional neutrons")
    parser.add_argument("--output", type=str, required=True, help="Output file path (JSON)")

    args = parser.parse_args()

    # The maximum number of atoms is 8 hydrogens per averagine
    isotoper = Isotoper(n=args.max_aas * 8)

    result = []

    for aas in range(args.min_aas, args.max_aas + 1):
        composition = {atom: round(count * aas) for atom, count in AVERAGINE_COMPOSITION.items()}
        mass = sum(
            ATOMS[atom].isotopes[0].mass * count for atom, count in composition.items()
        )

        distribution = isotoper.get_molecule_distribution(composition, args.max_neutrons)
        result.append({"mass": mass, "distribution": distribution})

    with open(args.output, "w") as f:
        import json
        f.write('[\n')
        for idx, item in enumerate(result):
            mass_json = json.dumps(item["mass"])
            dist_json = json.dumps(item["distribution"])
            f.write('  {\n')
            f.write(f'    "mass": {mass_json},\n')
            f.write(f'    "distribution": {dist_json}\n')
            f.write('  }' + (',\n' if idx < len(result) - 1 else '\n'))
        f.write(']\n')
