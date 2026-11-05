from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import asdict, dataclass

import numpy as np
from scipy.linalg import eig


@dataclass(frozen=True)
class SolveResult:
    omega: int
    basis_size: int
    epsilon: float
    energy_hartree: float
    build_seconds: float
    solve_seconds: float
    total_seconds: float
    nuclear_charge: int


def symmetric_basis_for_omega(omega: int) -> list[tuple[int, int, int]]:
    """Ordering used in Table I of Pekeris (1958) for the symmetric ground state."""
    basis: list[tuple[int, int, int]] = []
    for ww in range(omega + 1):
        for vv in range(ww + 1):
            for uu in range(vv + 1):
                l = uu
                m = vv - uu
                n = ww - vv
                if l + m + n <= omega and l <= m:
                    basis.append((l, m, n))
    return basis


def canonical_singlet_state(l: int, m: int, n: int) -> tuple[int, int, int]:
    if l <= m:
        return (l, m, n)
    return (m, l, n)


def eq22_terms(l: int, m: int, n: int, z: int) -> list[tuple[tuple[int, int, int], float, float]]:
    """Return the 33 terms of Pekeris' Eq. (22) as (shift, a, b) triples.

    Each triple means: coefficient A(l+dl, m+dm, n+dn) is multiplied by
    (a + epsilon * b) in the recurrence equation for A(l, m, n).

    Parameters
    ----------
    l, m, n : int
        Indices of the current basis state.
    z : int
        Nuclear charge (Z=2 for helium).

    Returns
    -------
    list of ((dl, dm, dn), a, b)
        The 33 recurrence terms, grouped by their origin in the Laguerre identities.
        When a=0 the term arises purely from the xL_n coordinate-multiplication identity.
        When b=0 the term arises from the pure L_n'' derivative identity on w.
    """
    terms: list[tuple[tuple[int, int, int], float, float]] = []

    def add(shift: tuple[int, int, int], a: float, b: float) -> None:
        terms.append((shift, float(a), float(b)))

    # +2 shifts: from the xL_n''(x) identity on u, v, and their cross terms
    add((2, 0, 0), -4 * (l + 1) * (l + 2) * z, 4 * (l + 1) * (l + 2) * (1 + m + n))
    add((0, 2, 0), -4 * (m + 1) * (m + 2) * z, 4 * (m + 1) * (m + 2) * (1 + l + n))
    add((1, 1, 0), 4 * (l + 1) * (m + 1) * (1 - 2 * z), 4 * (l + 1) * (m + 1) * (2 + l + m))
    add((1, 0, 1), 2 * (l + 1) * (n + 1) * (1 - 2 * z), 2 * (l + 1) * (n + 1) * (2 + 2 * m + n))
    add((0, 1, 1), 2 * (m + 1) * (n + 1) * (1 - 2 * z), 2 * (m + 1) * (n + 1) * (2 + 2 * l + n))
    add((0, 0, 2), (n + 1) * (n + 2), 0)

    # +1 shifts: from the xL_n'(x) identity
    add(
        (1, 0, 0),
        (l + 1) * (4 * z * (4 * l + 4 * m + 2 * n + 7) - 8 * m - 4 * n - 6),
        -2 * (l + 1) * ((m + n) * (4 * m + 12 * l) + n**2 + 12 * l + 18 * m + 15 * n + 14),
    )
    add(
        (0, 1, 0),
        (m + 1) * (4 * z * (4 * l + 4 * m + 2 * n + 7) - 8 * l - 4 * n - 6),
        -2 * (m + 1) * ((l + n) * (4 * l + 12 * m) + n**2 + 12 * m + 18 * l + 15 * n + 14),
    )
    add(
        (0, 0, 1),
        4 * (n + 1) * (z * (2 * l + 2 * m + 2) - l - m - n - 2),
        4 * (n + 1) * (l**2 + m**2 - 4 * l * m - 2 * l * n - 2 * m * n - 3 * l - 3 * m - 2 * n - 2),
    )

    # mixed +2/-1 shifts: from cross-derivative terms in the perimetric PDE
    add((0, 2, -1), 0, 4 * (m + 1) * (m + 2) * n)
    add((2, 0, -1), 0, 4 * (l + 1) * (l + 2) * n)
    add((-1, 0, 2), 0, 2 * l * (n + 1) * (n + 2))
    add((0, -1, 2), 0, 2 * m * (n + 1) * (n + 2))

    # diagonal (0,0,0): collects all middle-index contributions from every xL_n identity
    add(
        (0, 0, 0),
        4 * (2 * l + 1) * (2 * m + 1)
        + 4 * (2 * n + 1) * (l + m + 1)
        + 6 * n**2
        + 6 * n
        + 2
        - 4 * z * ((l + m) * (6 * l + 6 * m + 4 * n + 12) - 4 * l * m + 4 * n + 8),
        4
        * (
            (l + m) * (10 * l * m + 10 * m * n + 10 * l * n + 10 * l + 10 * m + 18 * n + 4 * n**2 + 16)
            + l * m * (4 - 12 * n)
            + 8
            + 12 * n
            + 4 * n**2
        ),
    )

    # mixed +1/-1 shifts: off-diagonal parts of the cross-derivative PDE terms
    add((-1, 1, 0), 4 * l * (m + 1) * (1 - 2 * z), 4 * l * (m + 1) * (1 + l + m))
    add((1, -1, 0), 4 * (l + 1) * m * (1 - 2 * z), 4 * (l + 1) * m * (1 + l + m))
    add((-1, 0, 1), 2 * l * (n + 1) * (1 - 2 * z), 2 * l * (n + 1) * (2 * m - 4 * l - n))
    add((0, -1, 1), 2 * m * (n + 1) * (1 - 2 * z), 2 * m * (n + 1) * (2 * l - 4 * m - n))
    add((1, 0, -1), 2 * (l + 1) * n * (1 - 2 * z), 2 * (l + 1) * n * (2 * m - 4 * l - n - 3))
    add((0, 1, -1), 2 * (m + 1) * n * (1 - 2 * z), 2 * (m + 1) * n * (2 * l - 4 * m - n - 3))

    # -1 shifts: lowering terms from the xL_n'(x) identity
    add(
        (-1, 0, 0),
        2 * l * (-(4 * m + 2 * n + 3) + z * (8 * l + 8 * m + 4 * n + 6)),
        -2 * l * ((m + n + 1) * (12 * l + 4 * m + 2) + n + n**2),
    )
    add(
        (0, -1, 0),
        2 * m * (-(4 * l + 2 * n + 3) + z * (8 * l + 8 * m + 4 * n + 6)),
        -2 * m * ((l + n + 1) * (12 * m + 4 * l + 2) + n + n**2),
    )
    add(
        (0, 0, -1),
        4 * n * (-(l + m + n + 1) + z * (2 * l + 2 * m + 2)),
        -4 * n * ((l + m) * (1 + 2 * n - l - m) + 6 * l * m + 2 * n),
    )

    # ±2 shifts: from L_n''(x) lowering and mixed cross terms
    add((1, 0, -2), 0, 2 * n * (n - 1) * (l + 1))
    add((0, 1, -2), 0, 2 * n * (n - 1) * (m + 1))
    add((-2, 0, 1), 0, 4 * l * (l - 1) * (n + 1))
    add((0, -2, 1), 0, 4 * m * (m - 1) * (n + 1))
    add((-2, 0, 0), -4 * l * (l - 1) * z, 4 * l * (l - 1) * (1 + m + n))
    add((0, -2, 0), -4 * m * (m - 1) * z, 4 * m * (m - 1) * (1 + l + n))
    add((0, 0, -2), n * (n - 1), 0)

    # double -1 shifts: coupling to coefficients lower in two indices simultaneously
    add((-1, -1, 0), 4 * l * m * (1 - 2 * z), 4 * l * m * (l + m))
    add((-1, 0, -1), 2 * l * n * (1 - 2 * z), 2 * l * n * (2 * m + n + 1))
    add((0, -1, -1), 2 * m * n * (1 - 2 * z), 2 * m * n * (2 * l + n + 1))

    return terms


def build_paper_matrices(omega: int, z: int = 2) -> tuple[np.ndarray, np.ndarray, list[tuple[int, int, int]]]:
    basis = symmetric_basis_for_omega(omega)
    index = {state: i for i, state in enumerate(basis)}
    size = len(basis)
    a = np.zeros((size, size), dtype=np.float64)
    b = np.zeros((size, size), dtype=np.float64)

    for row, (l, m, n) in enumerate(basis):
        for (dl, dm, dn), const_part, linear_part in eq22_terms(l, m, n, z):
            lp = l + dl
            mp = m + dm
            np_ = n + dn
            if lp < 0 or mp < 0 or np_ < 0:
                continue
            if lp + mp + np_ > omega:
                continue

            state = canonical_singlet_state(lp, mp, np_)
            col = index.get(state)
            if col is None:
                continue

            a[row, col] += const_part
            b[row, col] += linear_part

    return a, b, basis


def solve_paper_determinant(omega: int, z: int = 2) -> SolveResult:
    build_start = time.perf_counter()
    a, b, basis = build_paper_matrices(omega, z=z)
    build_seconds = time.perf_counter() - build_start

    solve_start = time.perf_counter()
    eigenvalues = eig(-a, b, left=False, right=False, check_finite=False)
    solve_seconds = time.perf_counter() - solve_start

    finite = eigenvalues[np.isfinite(eigenvalues)]
    real = finite[np.abs(finite.imag) < 1.0e-9].real
    positive = real[real > 0]
    if positive.size == 0:
        raise RuntimeError("No positive real epsilon eigenvalue was found.")

    epsilon = float(np.max(positive))
    energy_hartree = -(epsilon**2)

    return SolveResult(
        omega=omega,
        basis_size=len(basis),
        epsilon=epsilon,
        energy_hartree=energy_hartree,
        build_seconds=build_seconds,
        solve_seconds=solve_seconds,
        total_seconds=build_seconds + solve_seconds,
        nuclear_charge=z,
    )


def benchmark(omegas: list[int], z: int = 2) -> list[SolveResult]:
    return [solve_paper_determinant(omega, z=z) for omega in omegas]


def print_results(results: list[SolveResult]) -> None:
    print(" omega  size    epsilon        energy / Ha        build / s   solve / s   total / s")
    for result in results:
        print(
            f"{result.omega:>6d} {result.basis_size:>5d} "
            f"{result.epsilon:>12.9f} {result.energy_hartree:>16.12f} "
            f"{result.build_seconds:>10.4f} {result.solve_seconds:>10.4f} "
            f"{result.total_seconds:>10.4f}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Direct implementation of Pekeris' 1958 Eq. (22) recurrence for the helium ground state."
    )
    parser.add_argument("--omega", type=int, default=10, help="Polynomial order omega from Pekeris Table III.")
    parser.add_argument("--Z", type=int, default=2, help="Nuclear charge. Helium is Z=2.")
    parser.add_argument("--benchmark", type=int, nargs="*", help="Run several omega values and print a convergence table.")
    parser.add_argument("--json", type=str, help="Optional output path for machine-readable results.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.benchmark:
        results = benchmark(args.benchmark, z=args.Z)
    else:
        results = [solve_paper_determinant(args.omega, z=args.Z)]

    print_results(results)

    if args.json:
        with open(args.json, "w", encoding="utf-8") as handle:
            json.dump([asdict(result) for result in results], handle, indent=2)
        print(f"\nSaved results to {args.json}")


if __name__ == "__main__":
    main()
