"""Hartree-Fock limit of a two-electron atom (closed-shell 1s^2) by Taylor-series shooting.

RHF for 1s^2 is a 1-D nonlinear eigenproblem in the radial function f (phi = f / sqrt(4 pi),
int f^2 r^2 dr = 1) and Y = r V_H(r), V_H the Hartree potential of one electron:

    r f'' + 2 f' = 2 (Y - Z - eps r) f,        Y'' = -r f^2
    f(0) = a,  Y(0) = 0,  Y'(0) = b;            f(R) = 0,  Y'(R) = 0,  Y(R) = 1

Frobenius series at r = 0 (the cusp f'(0) = -Z f(0) comes out of the recurrence), then Taylor
steps sized from each series' own radius of convergence. For He the orbital's nearest singularity
sits near r = -1.1, so the radius at r0 is ~ r0 + 1.1 and the steps grow geometrically.
Since Y = b r + Y~ exactly, f depends on b and eps only through mu = eps - b: nested secant
solves find mu (f(R) = 0) and a (normalisation), and b = -Y~'(R) follows. Continuation in R;
working precision and order grow with R.

With the Dirichlet wall this is the exact HF minimiser on [0, R], so E_R >= E_HF and
E_R - E_HF ~ exp(-2 kappa R), kappa = sqrt(-2 eps). R is chosen from --digits.

Writes the orbital as exact dyadic Taylor coefficients per step, for certify.py.

    python oracle.py --digits 100 --out he_100.json
"""
import argparse
import json
import math
import time

import gmpy2
from gmpy2 import mpfr

LN10 = math.log(10)


def dot(a, b):
    """sum_i a[i] b[n-1-i]: one Cauchy-product coefficient."""
    return sum(map(gmpy2.mul, a, reversed(b)), mpfr(0))


def series(r0, s, eps, Z, p, energy=False):
    """Order-p Taylor coefficients about r0 of f and Y, plus J = int r f^2 Y and
    T = 1/2 int r^2 f'^2 when energy=True. s = (f, f', Y, Y', J, T) at r0;
    at r0 = 0 the Frobenius recurrence is used and s[1] is ignored."""
    f = [s[0]] if r0 == 0 else [s[0], s[1]]
    Y, J, T = [s[2], s[3]], [s[4]], [s[5]]
    ff, ffY, df, dd = [], [], [], []
    for k in range(p):
        ff.append(dot(f[: k + 1], f[: k + 1]))
        Y.append(-(r0 * ff[k] + (ff[k - 1] if k else 0)) / ((k + 2) * (k + 1)))
        rhs = 2 * (dot(Y[: k + 1], f[: k + 1]) - Z * f[k] - eps * (r0 * f[k] + (f[k - 1] if k else 0)))
        if r0 == 0:
            f.append(rhs / ((k + 1) * (k + 2)))
        else:
            f.append((rhs - (k + 1) * (k + 2) * f[k + 1]) / (r0 * (k + 2) * (k + 1)))
        if energy:
            ffY.append(dot(ff, Y[: k + 1]))
            J.append((r0 * ffY[k] + (ffY[k - 1] if k else 0)) / (k + 1))
            df.append((k + 1) * f[k + 1])
            dd.append(dot(df, df))
            T.append((r0 * r0 * dd[k] + 2 * r0 * (dd[k - 1] if k else 0) + (dd[k - 2] if k > 1 else 0)) / (2 * (k + 1)))
    return f, Y, J, T


def value(c, h):
    v = mpfr(0)
    for x in reversed(c):
        v = v * h + x
    return v


def slope(c, h):
    v = mpfr(0)
    for i in range(len(c) - 1, 0, -1):
        v = v * h + i * c[i]
    return v


def step_size(f, c, r0, R):
    """c times the radius of convergence estimated from the last coefficients (root test),
    rounded down to a 2^-16 grid so knots stay exact dyadics."""
    k0 = len(f) - 9
    rho = min((abs(f[0]) / abs(f[k])) ** (mpfr(1) / k) for k in range(k0, len(f)) if f[k])
    h = mpfr(math.floor(float(c * rho) * 65536)) / 65536
    return min(max(h, mpfr(1) / 65536), R - r0)


def shoot(x, R, Z, p, c, mesh=None, energy=False):
    """Integrate 0 -> R with Y'(0) = 0. Y = b r + Y~ exactly (Y'' does not involve b), and b r
    only shifts Y/r by a constant, so f depends on b and eps only through mu = eps - b.
    x = (a, mu). Returns (state at R, mesh, [(r0, h, f-coefficients)])."""
    a, mu = x
    s = (a, mpfr(0), mpfr(0), mpfr(0), mpfr(0), mpfr(0))
    r0, hs, elems = mpfr(0), [], []
    while r0 < R:
        f, Y, J, T = series(r0, s, mu, Z, p, energy)
        h = step_size(f, c, r0, R) if mesh is None else mesh[len(hs)]
        hs.append(h)
        if energy:
            elems.append((r0, h, f))
        s = (value(f, h), slope(f, h), value(Y, h), slope(Y, h), value(J, h), value(T, h))
        r0 += h
    return s, hs, elems


def secant(g, x0, x1, tol):
    """Root of the scalar function g by the secant method from x0, x1."""
    g0, g1 = g(x0), g(x1)
    for _ in range(60):
        if g1 == g0:  # at the noise floor
            return x1
        x0, x1, g0 = x1, x1 - g1 * (x1 - x0) / (g1 - g0), g1
        if abs(x1 - x0) <= tol * abs(x1):
            return x1
        g1 = g(x1)
    raise RuntimeError("secant did not converge")


def stage(x, R, Z, p, c, tol, rel):
    """Solve f(R) = 0 and N = Y~(R) - R Y~'(R) = int_0^R f^2 r^2 = 1 for x = (a, mu) as two nested
    scalar problems. f(R) is linear in the growing-mode amplitude, so mu(a) takes ~2 secant steps
    and leaves no growing mode behind; N(a) is then smooth. (A joint 2x2 Newton fails at large R:
    N carries ~exp(2 kappa R) c^2, which makes its Jacobian row parallel to f(R)'s.)
    The mesh is frozen per stage so both residuals are smooth in x. Returns (x, mesh, shots)."""
    _, mesh, _ = shoot(x, R, Z, p, c)
    shots, mu = [0], x[1]

    def end(a, m):
        shots[0] += 1
        return shoot((a, m), R, Z, p, c, mesh)[0]

    def norm(a):
        nonlocal mu
        mu = secant(lambda m: end(a, m)[0], mu, mu * (1 + rel), tol)
        s = end(a, mu)
        return s[2] - R * s[3] - 1

    a = secant(norm, x[0], x[0] * (1 + rel), tol)
    norm(a)  # mu for the final a
    return [a, mu], mesh, shots[0]


def solve(digits, Z=2, R=None, c=0.15, guard=40, log=print):
    """Continuation in R up to the wall giving ~`digits` significant figures.
    Returns a dict ready for json (orbital as exact dyadic [mantissa, exponent] pairs)."""
    t0 = time.time()
    zeta = Z - 5 / 16  # screened hydrogenic start
    eps = zeta * zeta / 2 - Z * zeta + 5 * zeta / 8
    up = lambda v: math.ceil(v * 64) / 64  # walls on a 1/64 grid: exact dyadics, and x1.2 stages at any Z
    x, r = [mpfr(2 * zeta ** 1.5), mpfr(eps - zeta)], up(7 / zeta)
    while True:
        kappa = math.sqrt(-2 * float(eps))
        Rf = R or up((digits + 8) * LN10 / (2 * kappa))
        Rk = min(r, Rf)
        d = min(digits, 2 * kappa * Rk / LN10) + 10  # what this wall can deliver, plus margin
        gmpy2.get_context().precision = int((d + guard) / math.log10(2))
        p = math.ceil((d + 5) / math.log10(1 / c))
        x = [mpfr(v) for v in x]
        x, mesh, its = stage(x, mpfr(Rk), Z, p, mpfr(c), mpfr(10) ** -d, mpfr(10) ** (-(d + guard) / 2))
        eps = x[1] - shoot(x, mpfr(Rk), Z, p, mpfr(c), mesh)[0][3]  # eps = mu + b, b = -Y~'(R)
        log(f"R={Rk:8.3f} digits={d:5.1f} order={p:4d} steps={len(mesh):3d} shots={its:3d} eps={float(eps):.15f}  {time.time() - t0:7.1f}s")
        if Rk == Rf:
            break
        r = up(Rk * 1.2)
    s, mesh, elems = shoot(x, mpfr(Rk), Z, p, mpfr(c), mesh, energy=True)
    b, N = -s[3], s[2] - Rk * s[3]
    E = 2 * (x[1] + b) - (s[4] + b * N)  # 2 eps - J, with J = J~ + b N
    nd = int(d) + 5
    exact = lambda v: [str(m) for m in mpfr(v).as_mantissa_exp()]
    return {
        "Z": Z, "R": Rk, "digits": digits, "order": p, "step": c, "precision_bits": gmpy2.get_context().precision,
        "seconds": round(time.time() - t0, 1),
        "E": format(E, f".{nd}f"), "E_virial": format(-2 * s[5], f".{nd}f"), "eps": format(x[1] + b, f".{nd}f"),
        "f0": format(x[0], f".{nd}f"), "VH0": format(b, f".{nd}f"),
        "elements": [[exact(r0), exact(h), [exact(v) for v in f]] for r0, h, f in elems],
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--digits", type=int, required=True, help="target significant figures (sets the wall R)")
    ap.add_argument("--Z", type=int, default=2)
    ap.add_argument("--R", type=float, help="override the wall radius")
    ap.add_argument("--step", type=float, default=0.15, help="step / radius of convergence (sets the order)")
    ap.add_argument("--guard", type=int, default=40, help="extra working digits")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = solve(a.digits, a.Z, a.R, a.step, a.guard)
    with open(a.out, "w") as fh:
        json.dump(out, fh)
    print(f"E      = {out['E']}\nE_vir  = {out['E_virial']}\n({out['seconds']}s, wrote {a.out})")
