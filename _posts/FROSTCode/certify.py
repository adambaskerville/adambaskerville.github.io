"""Rigorous upper bound on the closed-shell RHF energy of an atom with s orbitals, in exact rational arithmetic.

Reads the piecewise Taylor polynomials written by oracle.py (one orbital) or oracle_ns.py (several, each
with its own wall) and turns them into admissible trial orbitals *exactly*: dyadic coefficients,
continuity forced at every knot, f_i(R_i) = 0 and f_i = 0 beyond R_i. For one orbital the energy is

    E[f] = 2 (T + V) / N + J / N^2
    N = int f^2 r^2,  T = 1/2 int f'^2 r^2,  V = -Z int f^2 r,
    J = 2 int f^2 r Q(r),  Q(r) = int_0^r f^2 s^2 ds        (= (phi phi | phi phi))

and in general the determinant energy with P = S^-1 (see certify()), using integers and fractions
only: every integrand is a polynomial on a dyadic interval.
By the variational principle E_HF <= E[f], with no rounding anywhere. The decimal printed
is rounded towards +infinity, so it is itself an upper bound. Shares no code with oracle.py.

    python certify.py he_100.json [--order q] [--print-digits n]
"""
import argparse
import json

from gmpy2 import lcm, log10, mpfr, mpq, mpz

_L = {}


def lcm_to(n):
    """lcm(1..n), cached."""
    if n not in _L:
        v = mpz(1)
        for k in range(2, n + 1):
            v = lcm(v, k)
        _L[n] = v
    return _L[n]


def pmul_school(a, b):
    out = [mpz(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def pmul(a, b):
    """Exact product of integer polynomials by Kronecker substitution: each sign part is packed into one
    integer (k bytes per coefficient, enough that no slot of the product overflows), multiplied once by
    GMP, and unpacked; a = a+ - a-."""
    if min(len(a), len(b)) < 16:
        return pmul_school(a, b)
    n = len(a) + len(b) - 1
    k = (max(map(abs, a)).bit_length() + max(map(abs, b)).bit_length() + min(len(a), len(b)).bit_length() + 8) // 8

    def pack(c):
        return mpz(int.from_bytes(b"".join(int(v).to_bytes(k, "little") for v in c), "little"))

    def parts(c):
        return pack([max(v, 0) for v in c]), pack([max(-v, 0) for v in c])

    def unpack(x):
        d = int(x).to_bytes(n * k, "little")
        return [int.from_bytes(d[i * k:(i + 1) * k], "little") for i in range(n)]

    (ap, an), (bp, bn) = parts(a), parts(b)
    pos, neg = unpack(ap * bp + an * bn), unpack(ap * bn + an * bp)
    return [mpz(x - y) for x, y in zip(pos, neg)]


def pint(a):
    """int_0^1 a(s) ds for an integer polynomial, exactly."""
    L = lcm_to(len(a))
    return mpq(sum(c * (L // (k + 1)) for k, c in enumerate(a)), L)


def dyadic(pair):
    m, e = int(pair[0]), int(pair[1])
    return mpq(m * 2 ** e) if e >= 0 else mpq(m, 2 ** -e)


def inverse(S):
    """Exact inverse of a small positive-definite rational matrix (Gauss-Jordan; pivots are positive)."""
    n = len(S)
    M = [list(row) + [mpq(int(i == j)) for j in range(n)] for i, row in enumerate(S)]
    for c in range(n):
        M[c] = [v / M[c][c] for v in M[c]]
        for r in range(n):
            if r != c:
                M[r] = [x - M[r][c] * y for x, y in zip(M[r], M[c])]
    return [row[n:] for row in M]


def trial(data, order=None):
    """The trial orbitals of `data` exactly: oracle.py writes one per element, oracle_ns.py ("walls" present)
    a list, each zero beyond its own wall. Returns (els, F, A): knots (r0, h) as mpq, and per element the
    integer coefficients F[e][i] in s = (r - r0) / h (all scaled by one 2^B), with continuity and
    f_i(R_i) = 0 imposed exactly; 2^A makes every knot an integer."""
    els = []
    for r0, h, fs in data["elements"]:
        r0, h = dyadic(r0), dyadic(h)
        fs = fs if "walls" in data else [fs]
        fs = [[dyadic(c) for c in f][: None if order is None else order + 1] for f in fs]
        els.append((r0, h, [[c * h ** k for k, c in enumerate(f)] for f in fs]))  # coefficients in s = (r - r0) / h
    if not els or els[0][0] != 0:
        raise ValueError("first element must start at r = 0")
    if any(h <= 0 for _, h, _ in els) or any(a0 + h0 != a1 for (a0, h0, _), (a1, _, _) in zip(els, els[1:])):
        raise ValueError("elements must tile [0, R] exactly")  # a reversed interval would fake a lower 'bound'
    n = len(els[0][2])
    last = [max(e for e, (_, _, fs) in enumerate(els) if fs[i]) for i in range(n)]
    if any(not fs[i] for i in range(n) for _, _, fs in els[: last[i] + 1]):
        raise ValueError("each orbital must be supported on [0, R_i]")
    # one common power of two makes every coefficient and knot an integer (E is scale invariant)
    B = max(c.denominator.bit_length() - 1 for _, _, fs in els for f in fs for c in f)
    A = max(v.denominator.bit_length() - 1 for a, h, _ in els for v in (a, h))
    F = [[[mpz(c * 2 ** B) for c in f] for f in fs] for _, _, fs in els]
    for i in range(n):
        for e in range(1, last[i] + 1):  # exact continuity: next value at s=0 := previous value at s=1
            F[e][i][0] = sum(F[e - 1][i])
        F[last[i]][i][1] -= sum(F[last[i]][i])  # exact f_i(R_i) = 0
    return els, F, A


def certify(data, order=None):
    """Exact E[f] (an mpq) for the trial orbitals of `data`. E is the closed-shell determinant's,
        E = 2 sum_ij P_ij h_ij + sum_ijkl P_ij P_kl [2 (ij|kl) - (ik|jl)],   P = S^-1,
    exact for the span of the orbitals whether or not they are orthonormal; one orbital gives
    2 (T + V) / N + J / N^2."""
    Z = int(data["Z"])
    els, F, A = trial(data, order)
    n = len(F[0])
    pairs = [(i, j) for i in range(n) for j in range(i, n)]
    S, T, V, Q = ({q: mpq(0) for q in pairs} for _ in range(4))
    C = {(q, p): mpq(0) for q in pairs for p in pairs}  # int rho_q r Q_p, Q_p(r) = int_0^r rho_p s^2
    for (a, h, _), Fe in zip(els, F):
        al, et = mpz(a * 2 ** A), mpz(h * 2 ** A)  # r = (al + et s) / 2^A
        u1, u2 = [al, et], [al * al, 2 * al * et, et * et]
        live = [q for q in pairs if Fe[q[0]] and Fe[q[1]]]
        rho = {q: pmul(Fe[q[0]], Fe[q[1]]) for q in live}
        w, g = {q: pmul(rho[q], u2) for q in live}, {q: pmul(rho[q], u1) for q in live}
        dP = [[(k + 1) * P[k + 1] for k in range(len(P) - 1)] for P in Fe]
        Lw = {q: lcm_to(len(w[q])) for q in live}
        W = {q: [mpz(0)] + [c * (Lw[q] // (k + 1)) for k, c in enumerate(w[q])] for q in live}  # Lw * int_0^s w
        gi = {q: pint(g[q]) for q in live}
        for q in live:
            for p in pairs:
                C[(q, p)] += mpq(et, 2 ** (2 * A)) * Q[p] * gi[q]
                if p in live:
                    C[(q, p)] += mpq(et * et, 2 ** (5 * A) * Lw[p]) * pint(pmul(g[q], W[p]))
            T[q] += pint(pmul(pmul(dP[q[0]], dP[q[1]]), u2)) / (2 * et * 2 ** A)
            V[q] -= Z * mpq(et, 2 ** (2 * A)) * gi[q]
        for q in live:
            dn = mpq(et, 2 ** (3 * A)) * pint(w[q])
            S[q] += dn
            Q[q] += dn
    key = lambda i, j: (min(i, j), max(i, j))  # noqa: E731
    G = lambda q, p: C[(q, p)] + C[(p, q)]  # noqa: E731  (q|p)
    P = inverse([[S[key(i, j)] for j in range(n)] for i in range(n)])
    ij = [(i, j) for i in range(n) for j in range(n)]
    return sum(2 * P[i][j] * (T[key(i, j)] + V[key(i, j)]) for i, j in ij) + sum(
        P[i][j] * P[k][l] * (2 * G(key(i, j), key(k, l)) - G(key(i, k), key(j, l))) for i, j in ij for k, l in ij)


def decimal_up(q, n):
    """q rounded towards +infinity to n decimals, as a string."""
    s = 10 ** n
    v = -((-q.numerator * s) // q.denominator)
    sign, v = ("-", -v) if v < 0 else ("", v)
    return f"{sign}{v // s}.{str(v % s).zfill(n)}"


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("orbital")
    ap.add_argument("--order", type=int, help="truncate each step's polynomial to this order")
    ap.add_argument("--print-digits", type=int)
    a = ap.parse_args()
    data = json.load(open(a.orbital))
    E = certify(data, a.order)
    n = a.print_digits or int(data["digits"]) + 10
    gap = mpfr(E - mpq(data["E"]), 64)
    print(f"E_cert  <= {decimal_up(E, n)}   (exact rational upper bound, rounded up)")
    print(f"E_oracle = {data['E'][: n + 3]}")
    sf = -log10(abs(gap) / abs(mpfr(E, 64))) if gap else float("inf")
    print(f"E_cert - E_oracle = {gap:.3e}  ->  agree to {float(sf):.1f} significant figures")
