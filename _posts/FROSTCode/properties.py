"""Orbital properties of the trial orbitals that certify.py bounds (certify.trial), in exact rational arithmetic.

    <r^k>_i = int f_i^2 r^(k+2) dr / int f_i^2 r^2 dr     k >= -2, per electron in orbital i
    rho(0)  = sum_i 2 f_i(0)^2 / (4 pi N_i)                 total electron density at the nucleus
    cusp_i  = f_i'(0) / f_i(0)                              -Z exactly: the Frobenius recurrence builds it in

Each orbital is normalised on its own. For n > 1 only density properties (rho(0), sum_i 2 <r^k>_i) are accurate
to the wall error: two independent Be runs agree on them to ~200 s.f. at 200 digits, but on the individual
<r^k>_i only to ~|S_12| (~1e-130), a rotation within the occupied space that leaves the density unchanged
to first order.

    python properties.py he500a.json [--k -2 -1 1 2 3]
"""
import argparse
import json
from math import comb

import gmpy2
from gmpy2 import mpfr, mpq, mpz

from certify import pint, pmul, trial


def properties(data, ks=(-2, -1, 1, 2)):
    """{'cusp': [mpq per orbital], 'rho0': mpfr, 'T': [mpq per orbital], 'r^k': {k: [mpq per orbital]}};
    T_i = 1/2 int f_i'^2 r^2 / N_i, the kinetic energy of one electron in orbital i."""
    els, F, A = trial(data)
    n, ks = len(F[0]), sorted(set(ks) | {0})
    M = [{k: mpq(0) for k in ks} for _ in range(n)]  # int f_i^2 r^(k+2) dr, all scaled by 2^(2B)
    T = [mpq(0)] * n
    for (a, h, _), Fe in zip(els, F):
        al, et = mpz(a * 2 ** A), mpz(h * 2 ** A)  # r = (al + et s) / 2^A
        for i, f in enumerate(Fe):
            if f:
                df = [(k + 1) * c for k, c in enumerate(f[1:])]
                T[i] += pint(pmul(pmul(df, df), [al * al, 2 * al * et, et * et])) / (2 * et * 2 ** A)  # as in certify
                rho = pmul(f, f)
                mu = [pint([0] * j + rho) for j in range(ks[-1] + 3)]  # int_0^1 rho s^j ds
                for k in ks:
                    m = k + 2
                    M[i][k] += mpq(et, 2 ** (A * (m + 1))) * sum(
                        comb(m, j) * al ** (m - j) * et ** j * mu[j] for j in range(m + 1))
    gmpy2.get_context().precision = int(data["precision_bits"])
    return {
        "cusp": [mpq(Fe0[1]) / (els[0][1] * Fe0[0]) for Fe0 in F[0]],
        "rho0": sum(2 * mpfr(F[0][i][0] ** 2 / M[i][0]) for i in range(n)) / (4 * gmpy2.const_pi()),
        "T": [T[i] / M[i][0] for i in range(n)],
        "r^k": {k: [M[i][k] / M[i][0] for i in range(n)] for k in ks if k},
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("orbital")
    ap.add_argument("--k", type=int, nargs="+", default=[-2, -1, 1, 2])
    ap.add_argument("--print-digits", type=int)
    a = ap.parse_args()
    data = json.load(open(a.orbital))
    d = a.print_digits or int(data["digits"]) + 10
    p = properties(data, a.k)
    print("cusp f'(0)/f(0) =", ", ".join(str(c) for c in p["cusp"]))
    print(f"rho(0) = {p['rho0']:.{d}e}")
    T = 2 * sum(p["T"])
    print(f"T = {mpfr(T):.{d}e}   virial ratio V/T + 2 = 1 + E/T = {mpfr(1 + mpq(data['E']) / T):.3e}")
    for k, v in p["r^k"].items():
        print(f"<r^{k}> = " + ", ".join(f"{mpfr(x):.{d}e}" for x in v) + f"   sum_i 2 <r^{k}>_i = {mpfr(2 * sum(v)):.{d}e}")
