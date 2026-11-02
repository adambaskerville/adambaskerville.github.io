"""Hartree-Fock limit of closed-shell atoms whose occupied orbitals are all s, by Taylor-series shooting.
n = 1 is the He-like 1s^2 problem of oracle.py, n = 2 the Be-like 1s^2 2s^2.

Radial functions f_i (phi_i = f_i / sqrt(4 pi)), pair potentials Y_ij = r V_ij with
V_ij(r) = int f_i f_j s^2 / max(r, s) ds. The canonical Fock equations are

    r f_i'' + 2 f_i' = 2 [(sum_j 2 Y_jj - Y_ii - Z - eps_i r) f_i - sum_{j != i} Y_ij f_j]
    Y_ij'' = -r f_i f_j,      Y_ij(0) = 0,  Y_ij'(R) = 0

Y_ij = b_ij r + Y~_ij exactly, so the f_i depend on eps and b only through
mu_i = eps_i - sum_j 2 b_jj + b_ii and nu_ij = b_ij (i != j).

Each orbital has its own wall R_i, at the same exp(-2 kappa_i R_i): integrating 1s out to the 2s wall would
let its growing mode exp(kappa_1 r) demand kappa_1 / kappa_2 (~4 for Be) times the working precision.
Orbital i is a variational function on [0, R_i], so the energy is still an upper bound. With unequal walls
exact orthogonality would need a multiplier lambda_ij, and enforcing it is badly conditioned (S_12 responds
to lambda at ~1e-4, and the step that asks for forces a 1s tail far above its natural exp(-kappa_1 R_1),
which the exchange term feeds into the 2s growing mode). We keep the canonical equations instead
(lambda = 0, nu_ij = b_ij); the orbitals are then orthogonal only to S_ij ~ exp(-kappa_1 R_1), which costs
O(S_ij^2) ~ the wall error in energy, and the energy is evaluated with the non-orthogonal determinant
formula (P = S^-1), exact for the span of the f_i, here and in certify.py.

Shooting: nested scalar solves for mu (f_i(R_i) = 0), outer Broyden for (a, nu) (S_ii = 1, nu = b). The
start is a homotopy in the electron-electron coupling t from hydrogenic orbitals, exact at t = 0.

    python oracle_ns.py --Z 4 --digits 100 --out be100.json
"""
import argparse
import json
import math
import time

import gmpy2
from gmpy2 import mpfr

from oracle import LN10, dot, slope, value


def pairs(n):
    return [(i, j) for i in range(n) for j in range(i, n)]


def series(r0, st, mu, nu, Z, p, active, t, energy=False):
    """Order-p Taylor coefficients about r0 of the f_i and Y~_ij, and if energy of the pair integrals
    T_ij = 1/2 int r^2 f_i' f_j', V_ij = -Z int r f_i f_j and J_(ij,kl) = int r f_i f_j Y~_kl.
    st = state at r0; inactive orbitals (beyond their wall) are identically zero."""
    n, P = len(mu), pairs(len(mu))
    f = [None if not active[i] else [st["f"][i][0]] if r0 == 0 else list(st["f"][i]) for i in range(n)]
    Y = {q: list(st["Y"][q]) for q in P}
    E = {key: {q: [v] for q, v in st[key].items()} for key in ("T", "V", "J")}
    ff, C, df = {q: [] for q in P}, [[] for _ in range(n)], [[] for _ in range(n)]
    dd, fY = {q: [] for q in P}, {qq: [] for qq in st["J"]}
    tail = lambda c, k: r0 * c[k] + (c[k - 1] if k else 0)  # noqa: E731  (r c)_k about r0
    for k in range(p):
        for q in P:
            fi, fj = f[q[0]], f[q[1]]
            ff[q].append(dot(fi[: k + 1], fj[: k + 1]) if fi and fj else mpfr(0))
            Y[q].append(-tail(ff[q], k) / ((k + 1) * (k + 2)))
        new = []
        for i in range(n):
            C[i].append(t * (2 * sum(Y[(j, j)][k] for j in range(n)) - Y[(i, i)][k]))
            if not f[i]:
                new.append(None)
                continue
            g = dot(C[i], f[i][: k + 1]) - Z * f[i][k] - mu[i] * tail(f[i], k)
            for j in range(n):
                if j != i and f[j]:
                    q = (min(i, j), max(i, j))
                    g -= t * dot(Y[q][: k + 1], f[j][: k + 1]) + nu[q] * tail(f[j], k)
            g, kk = 2 * g, (k + 1) * (k + 2)
            new.append(g / kk if r0 == 0 else (g - kk * f[i][k + 1]) / (r0 * kk))
        for i in range(n):
            if f[i]:
                f[i].append(new[i])
                if energy:
                    df[i].append((k + 1) * f[i][k + 1])
        if energy:
            for q in P:
                both = f[q[0]] and f[q[1]]
                dd[q].append(dot(df[q[0]], df[q[1]]) if both else mpfr(0))
                E["T"][q].append((r0 * r0 * dd[q][k] + 2 * r0 * (dd[q][k - 1] if k else 0) + (dd[q][k - 2] if k > 1 else 0)) / (2 * (k + 1)))
                E["V"][q].append(-Z * tail(ff[q], k) / (k + 1))
            for qq in fY:
                fY[qq].append(dot(ff[qq[0]], Y[qq[1]][: k + 1]))
                E["J"][qq].append(tail(fY[qq], k) / (k + 1))
    return f, Y, E


def step_size(fs, c, r0, R):
    """c times the root-test radius of convergence (min over orbitals), on a 2^-16 grid. The scale is
    max(|f_0|, |f_1|) so that an expansion point at the 2s node does not fake a tiny radius."""
    rho = min(
        (max(abs(f[0]), abs(f[1])) / abs(f[k])) ** (mpfr(1) / k)
        for f in fs if f for k in range(len(f) - 9, len(f)) if f[k]
    )
    h = mpfr(math.floor(float(c * rho) * 65536)) / 65536
    return min(max(h, mpfr(1) / 65536), R - r0)


def shoot(a, mu, nu, walls, Z, p, c, t=1, mesh=None, energy=False, stop=None):
    """Integrate 0 -> stop (default max(walls)). Returns (state at the end, [f_i(R_i)], mesh, elements)."""
    n = len(a)
    P = pairs(n)
    st = {"f": [(v, mpfr(0)) for v in a], "Y": {q: (mpfr(0), mpfr(0)) for q in P},
          "T": {q: mpfr(0) for q in P}, "V": {q: mpfr(0) for q in P}, "J": {(q, r): mpfr(0) for q in P for r in P}}
    r0, hs, elems, fR = mpfr(0), [], [], [None] * n
    while r0 < (stop or walls[-1]):
        active = [r0 < w for w in walls]
        f, Y, E = series(r0, st, mu, nu, Z, p, active, t, energy)
        h = step_size(f, c, r0, min(w for w in walls if w > r0)) if mesh is None else mesh[len(hs)]
        hs.append(h)
        if energy:
            elems.append((r0, h, f))
        st = {"f": [(value(fi, h), slope(fi, h)) if fi else (mpfr(0), mpfr(0)) for fi in f],
              "Y": {q: (value(y, h), slope(y, h)) for q, y in Y.items()},
              **{key: {q: value(x, h) for q, x in E[key].items()} for key in E}}
        r0 += h
        for i in range(n):
            if active[i] and r0 == walls[i]:
                fR[i], st["f"][i] = st["f"][i][0], (mpfr(0), mpfr(0))  # f_i = 0 beyond its wall
    return st, fR, hs, elems


def solve_lin(A, b):
    """Gaussian elimination with row equilibration and partial pivoting (n <= 3)."""
    n = len(b)
    M = []
    for row, v in zip(A, b):
        s = max(abs(x) for x in row)
        M.append([x / s for x in row] + [v / s])
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(M[r][col]))
        M[col], M[piv] = M[piv], M[col]
        for r in range(col + 1, n):
            m = M[r][col] / M[col][col]
            M[r] = [x - m * y for x, y in zip(M[r], M[col])]
    x = [mpfr(0)] * n
    for r in reversed(range(n)):
        x[r] = (M[r][n] - sum(M[r][k] * x[k] for k in range(r + 1, n))) / M[r][r]
    return x


def newton(F, x, tol, rel, J=None):
    """Broyden's method on F(x) = 0 (the secant method in 1-D), from J or a forward-difference Jacobian.
    Stops when the step is below tol, or when max|F| has not improved for three iterations: the residual
    is at its noise floor (see Stage.inner). A stale J that fails to improve is rebuilt once.
    Returns (x, J) so J can seed the next call."""
    fx = F(x)
    best, stall, fresh = (max(map(abs, fx)), x, fx), 0, J is None
    for _ in range(60):
        if J is None:
            cols = []
            for k in range(len(x)):
                h = rel * max(abs(x[k]), 1)
                xk = list(x)
                xk[k] += h
                cols.append([(u - v) / h for u, v in zip(F(xk), fx)])
            J, fresh = [list(r) for r in zip(*cols)], True
        if any(not any(row) for row in J) or (len(fx) == 1 and fx[0] == 0):
            return best[1], J  # residual exactly unchanged or zero: noise floor
        dx = solve_lin(J, [-v for v in fx])
        xn = [u + d for u, d in zip(x, dx)]
        if max(abs(d) / max(abs(u), 1) for d, u in zip(dx, xn)) <= tol:
            return xn, J
        fn = F(xn)
        if max(map(abs, fn)) < best[0]:
            best, stall = (max(map(abs, fn)), xn, fn), 0
        else:
            stall += 1
            if stall == 3:
                return best[1], J
            if not fresh:
                J, x, fx = None, best[1], best[2]
                continue
        r = [u - v - sum(Jr[k] * dx[k] for k in range(len(dx))) for u, v, Jr in zip(fn, fx, J)]
        dd = sum(d * d for d in dx)
        J = [[Jr[k] + ri * dx[k] / dd for k in range(len(dx))] for Jr, ri in zip(J, r)]
        x, fx = xn, fn
    raise RuntimeError("newton did not converge")


class Stage:
    """Nested solve at fixed walls, mesh, precision and coupling t."""

    def __init__(self, Z, walls, p, c, t, tol, tol_in, rel):
        self.Z, self.walls, self.p, self.c, self.t, self.tol, self.tol_in, self.rel = Z, walls, p, c, t, tol, tol_in, rel
        self.mesh, self.shots = None, 0

    def shoot(self, a, mu, nu, energy=False, stop=None):
        self.shots += 1
        return shoot(a, mu, nu, self.walls, self.Z, self.p, self.c, self.t, self.mesh, energy, stop)

    def inner(self, a, nu, mu):
        """mu with f_i(R_i) = 0 by nested scalar solves, innermost orbital first: every trial mu_i re-solves
        mu_0..mu_{i-1}, so no inner orbital carries a growing mode into the coupling terms (a joint solve
        lands on spurious roots that this pollution creates). The re-solved inner orbitals still differ by
        ~exp(kappa_i R_i) x rounding near R_i, which the exchange coupling passes to f_j(R_j) as noise: mu_j
        is then fixed only to ~10^(-0.6 d), which the energy (quadratic in orbital errors) does not feel,
        hence newton's noise-floor stop. Returns (mu, state of the final shot)."""
        mu, last, top = list(mu), {}, len(mu) - 1

        def res(i, m):
            mu[i] = m
            if i:
                (mu[i - 1],), self.J[i - 1] = newton(lambda v: [res(i - 1, v[0])], [mu[i - 1]], self.tol_in, self.rel, self.J[i - 1])
                if i > 1:
                    res(i - 1, mu[i - 1])
            last["shot"] = self.shoot(a, mu, nu, stop=None if i == top else self.walls[i])  # f_i(R_i) needs only [0, R_i]
            return last["shot"][1][i]

        (mu[top],), self.J[top] = newton(lambda v: [res(top, v[0])], [mu[top]], self.tol_in, self.rel, self.J[top])
        res(top, mu[top])
        return mu, last["shot"][0]

    def run(self, a, mu, nu):
        """The outer Jacobian is rebuilt every stage: the solution moves by ~10^(-previous digits), so a fresh
        one converges in one or two steps, while a stale one only linearly, ~0.2 d orders short per stage."""
        n, off = len(a), [q for q in pairs(len(a)) if q[0] != q[1]]
        self.mesh = shoot(a, mu, nu, self.walls, self.Z, self.p, self.c, self.t)[2]
        self.J = [None] * n  # per-orbital d f_i(R_i) / d mu_i, kept across inner solves of this stage
        state = {"mu": mu}

        def outer(x):
            aa, nn = x[:n], dict(zip(off, x[n:]))
            state["mu"], st = self.inner(aa, nn, state["mu"])
            R = self.walls[-1]
            S = {q: st["Y"][q][0] - R * st["Y"][q][1] for q in pairs(n)}
            return [S[(i, i)] - 1 for i in range(n)] + [nn[q] + self.t * st["Y"][q][1] for q in off]

        x = newton(outer, list(a) + [nu[q] for q in off], self.tol, self.rel)[0]
        a, nu = x[:n], dict(zip(off, x[n:]))
        return a, self.inner(a, nu, state["mu"])[0], nu


def solve(digits, Z=4, n=2, c=0.15, guard=40, log=print):
    t0 = time.time()
    up = lambda v: math.ceil(v * 64) / 64  # noqa: E731
    ctx = gmpy2.get_context()

    def setup(d, g=guard):
        ctx.precision = int((d + g) / math.log10(2))
        # inner tolerance below the outer one: mu is solved to near the noise floor
        return (math.ceil((d + 5) / math.log10(1 / c)), mpfr(10) ** -d, mpfr(10) ** (-d - g / 2),
                mpfr(10) ** (-(d + g) / 2))

    # homotopy t: 0 -> 1 from hydrogenic orbitals (exact at t = 0) inside walls at 12 / kappa_i
    a = [mpfr(2) * (mpfr(Z) / (i + 1)) ** 1.5 for i in range(n)]
    mu = [-mpfr(Z) ** 2 / (2 * (i + 1) ** 2) for i in range(n)]
    nu = {q: mpfr(0) for q in pairs(n) if q[0] != q[1]}
    walls = [up(12 * (i + 1) / Z) for i in range(n)]
    d = 10  # the homotopy only has to deliver a start
    p, tol, tol_in, rel = setup(d, 16)
    K, prev = 40, None  # f(R_i) is ~exponential in mu: small t-steps with a linear predictor
    for k in range(1, K + 1):
        cur = (a, mu, nu)
        if prev:
            a, mu = [2 * u - v for u, v in zip(a, prev[0])], [2 * u - v for u, v in zip(mu, prev[1])]
            nu = {q: 2 * nu[q] - prev[2][q] for q in nu}
        stg = Stage(Z, [mpfr(w) for w in walls], p, mpfr(c), mpfr(k) / K, tol, tol_in, rel)
        (a, mu, nu), prev = stg.run(a, mu, nu), cur
    log(f"homotopy done: walls={walls} mu={[float(m) for m in mu]}  {time.time() - t0:.1f}s")

    def energies(st):
        b = {q: -st["Y"][q][1] for q in pairs(n)}
        eps = [mu[i] + 2 * sum(b[(j, j)] for j in range(n)) - b[(i, i)] for i in range(n)]
        return b, eps

    while True:
        b, eps = energies(stg.shoot(a, mu, nu)[0])
        kap = [math.sqrt(-2 * float(e)) for e in eps]
        Rf = up((digits + 8) * LN10 / (2 * kap[-1]))
        Rn = min(up(walls[-1] * 1.2), Rf)
        walls = [max(walls[i], up(Rn * kap[-1] / kap[i])) for i in range(n - 1)] + [Rn]
        d = min(digits, 2 * kap[-1] * Rn / LN10) + 10
        p, tol, tol_in, rel = setup(d)
        a, mu, nu = [mpfr(v) for v in a], [mpfr(v) for v in mu], {q: mpfr(v) for q, v in nu.items()}
        stg = Stage(Z, [mpfr(w) for w in walls], p, mpfr(c), mpfr(1), tol, tol_in, rel)
        a, mu, nu = stg.run(a, mu, nu)
        log(f"R={walls} digits={d:5.1f} order={p:4d} steps={len(stg.mesh):3d} shots={stg.shots:3d} "
            f"eps={[round(float(e), 12) for e in eps]}  {time.time() - t0:7.1f}s")
        if Rn == Rf:
            break
    st, fR, mesh, elems = stg.shoot(a, mu, nu, energy=True)
    b, eps = energies(st)
    R, P = mpfr(walls[-1]), pairs(n)
    key = lambda i, j: (min(i, j), max(i, j))  # noqa: E731
    S = {q: st["Y"][q][0] - R * st["Y"][q][1] for q in P}
    G = {(q, r): st["J"][(q, r)] + b[r] * S[q] for q in P for r in P}  # (q|r), b r part of Y_r added back
    Pm = [solve_lin([[S[key(i, j)] for j in range(n)] for i in range(n)], [mpfr(int(i == k)) for i in range(n)]) for k in range(n)]
    ij = [(i, j) for i in range(n) for j in range(n)]
    E = sum(2 * Pm[i][j] * (st["T"][key(i, j)] + st["V"][key(i, j)]) for i, j in ij) + sum(
        Pm[i][j] * Pm[k][l] * (2 * G[(key(i, j), key(k, l))] - G[(key(i, k), key(j, l))]) for i, j in ij for k, l in ij)
    Evir = -2 * sum(Pm[i][j] * st["T"][key(i, j)] for i, j in ij)
    nd = int(d) + 5
    exact = lambda v: [str(m) for m in mpfr(v).as_mantissa_exp()]  # noqa: E731
    fmt = lambda v: format(v, f".{nd}e")  # noqa: E731
    return {
        "Z": Z, "n": n, "walls": walls, "digits": digits, "order": p, "step": c,
        "precision_bits": ctx.precision, "seconds": round(time.time() - t0, 1),
        "E": format(E, f".{nd}f"), "E_virial": format(Evir, f".{nd}f"),
        "eps": [fmt(e) for e in eps], "f0": [fmt(v) for v in a], "S_offdiag": [fmt(S[q]) for q in P if q[0] != q[1]],
        "elements": [[exact(r0), exact(h), [[exact(v) for v in fi] if fi else [] for fi in f]] for r0, h, f in elems],
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--digits", type=int, required=True)
    ap.add_argument("--Z", type=int, default=4)
    ap.add_argument("--n", type=int, default=2, help="number of doubly occupied s orbitals")
    ap.add_argument("--step", type=float, default=0.15)
    ap.add_argument("--guard", type=int, default=40)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = solve(a.digits, a.Z, a.n, a.step, a.guard)
    with open(a.out, "w") as fh:
        json.dump(out, fh)
    print(f"E      = {out['E']}\nE_vir  = {out['E_virial']}\neps    = {out['eps']}\nS_ij   = {out['S_offdiag']}\n"
          f"({out['seconds']}s, wrote {a.out})")
