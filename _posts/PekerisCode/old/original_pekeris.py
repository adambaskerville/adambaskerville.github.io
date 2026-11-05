from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy.linalg import eigvalsh


ROOT = Path(__file__).resolve().parent
RR_DIR = ROOT / "pekeris_rr"


@dataclass(frozen=True)
class HeliumParameters:
    z1: float = -1.0
    z2: float = -1.0
    z3: float = 2.0
    m1: float = 1.0
    m2: float = 1.0
    m3: float = 1.0e12
    t: float = 1.0
    k: float = 1.0
    ep: float = 1.0

    def namespace(self) -> dict[str, float]:
        return {
            "Z1": self.z1,
            "Z2": self.z2,
            "Z3": self.z3,
            "m1": self.m1,
            "m2": self.m2,
            "m3": self.m3,
            "t": self.t,
            "K": self.k,
            "ep": self.ep,
        }


@dataclass(frozen=True)
class SolveResult:
    omega: int
    basis_size: int
    epsilon: float
    energy_hartree: float
    build_seconds: float
    solve_seconds: float
    total_seconds: float
    state: str


def load_recurrence(name: str) -> dict[tuple[int, int, int], object]:
    path = RR_DIR / f"{name}.rr"
    recurrence: dict[tuple[int, int, int], object] = {}
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        dl, dm, dn, expression = line.split(maxsplit=3)
        expression = expression.replace("^", "**").replace("Pi", "math.pi")
        recurrence[(int(dl), int(dm), int(dn))] = compile(
            expression,
            str(path),
            "eval",
        )
    return recurrence


def basis_for_omega(omega: int, state: str) -> list[tuple[int, int, int]]:
    basis: list[tuple[int, int, int]] = []
    want_triplet = state == "triplet"
    for ww in range(omega + 1):
        for vv in range(ww + 1):
            for uu in range(vv + 1):
                l = uu
                m = vv - uu
                n = ww - vv
                if l + m + n > omega:
                    continue
                if want_triplet and l < m:
                    basis.append((l, m, n))
                if (not want_triplet) and l <= m:
                    basis.append((l, m, n))
    return basis


def eval_rr(code: object, l: int, m: int, n: int, params: HeliumParameters) -> float:
    namespace = {"math": math, "l": l, "m": m, "n": n, **params.namespace()}
    return float(eval(code, {"__builtins__": {}}, namespace))


def recurrence_entry(
    recurrence: dict[tuple[int, int, int], object],
    basis: list[tuple[int, int, int]],
    row: int,
    col: int,
    state: str,
    params: HeliumParameters,
) -> float:
    l, m, n = basis[row]
    lp, mp, np_ = basis[col]

    def lookup(dl: int, dm: int, dn: int, left_l: int, left_m: int) -> float:
        code = recurrence.get((dl, dm, dn))
        if code is None:
            return 0.0
        return eval_rr(code, left_l, left_m, n, params)

    direct = lookup(lp - l, mp - m, np_ - n, l, m)
    exchange_col = 0.0
    exchange_row = 0.0
    exchange_both = 0.0

    if lp != mp:
        exchange_col = lookup(mp - l, lp - m, np_ - n, l, m)
    if l != m:
        exchange_row = lookup(lp - m, mp - l, np_ - n, m, l)
    if l != m and lp != mp:
        exchange_both = lookup(mp - m, lp - l, np_ - n, m, l)

    if state == "singlet":
        return direct + exchange_col + exchange_row + exchange_both
    if state == "triplet":
        return direct - exchange_col - exchange_row + exchange_both
    raise ValueError(f"Unknown state {state!r}. Expected 'singlet' or 'triplet'.")


def build_matrix(
    recurrence: dict[tuple[int, int, int], object],
    basis: list[tuple[int, int, int]],
    state: str,
    params: HeliumParameters,
) -> np.ndarray:
    size = len(basis)
    matrix = np.zeros((size, size), dtype=np.float64)
    for row in range(size):
        for col in range(row, size):
            matrix[row, col] = recurrence_entry(recurrence, basis, row, col, state, params)
    return matrix + matrix.T - np.diag(matrix.diagonal())


def solve_pekeris(
    omega: int,
    *,
    state: str = "singlet",
    params: HeliumParameters | None = None,
) -> SolveResult:
    params = params or HeliumParameters()
    basis = basis_for_omega(omega, state)

    hh_rr = load_recurrence("HH")
    ss_rr = load_recurrence("SS")

    build_start = time.perf_counter()
    hamiltonian = build_matrix(hh_rr, basis, state, params)
    overlap = -build_matrix(ss_rr, basis, state, params)
    build_seconds = time.perf_counter() - build_start

    solve_start = time.perf_counter()
    eigenvalues = eigvalsh(hamiltonian, overlap, check_finite=False)
    epsilon = abs(float(eigenvalues[-1]))
    energy_hartree = -((math.sqrt(params.k) * epsilon) ** 2)
    solve_seconds = time.perf_counter() - solve_start

    return SolveResult(
        omega=omega,
        basis_size=len(basis),
        epsilon=epsilon,
        energy_hartree=energy_hartree,
        build_seconds=build_seconds,
        solve_seconds=solve_seconds,
        total_seconds=build_seconds + solve_seconds,
        state=state,
    )


def benchmark(
    omegas: list[int],
    *,
    state: str = "singlet",
    params: HeliumParameters | None = None,
) -> list[SolveResult]:
    return [solve_pekeris(omega, state=state, params=params) for omega in omegas]


def print_results(results: list[SolveResult]) -> None:
    print(
        " omega  size    epsilon        energy / Ha        build / s   solve / s   total / s"
    )
    for result in results:
        print(
            f"{result.omega:>6d} {result.basis_size:>5d} "
            f"{result.epsilon:>12.9f} {result.energy_hartree:>16.12f} "
            f"{result.build_seconds:>10.4f} {result.solve_seconds:>10.4f} "
            f"{result.total_seconds:>10.4f}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Minimal reproduction of the original Pekeris helium method."
    )
    parser.add_argument(
        "--omega",
        type=int,
        default=10,
        help="Pekeris matrix order. The helium singlet ground state is already good by omega=10.",
    )
    parser.add_argument(
        "--state",
        choices=["singlet", "triplet"],
        default="singlet",
        help="Exchange symmetry. The 1958 ground-state calculation is the singlet case.",
    )
    parser.add_argument(
        "--benchmark",
        type=int,
        nargs="*",
        help="Run several omega values and print a convergence/timing table.",
    )
    parser.add_argument(
        "--json",
        type=Path,
        help="Optional output path for machine-readable benchmark results.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.benchmark:
        results = benchmark(args.benchmark, state=args.state)
    else:
        results = [solve_pekeris(args.omega, state=args.state)]

    print_results(results)

    if args.json:
        payload = [asdict(result) for result in results]
        args.json.write_text(json.dumps(payload, indent=2))
        print(f"\nSaved results to {args.json}")


if __name__ == "__main__":
    main()
