"""Runnable checks for perimetric4.tex (four-body perimetric coordinates).

    python3 checks.py        # every claim used in the note is an assert

Conventions: particles 1..4, pair order PAIRS, s = r**2 componentwise,
Gram matrix with particle 4 at the origin: P_ij = (s_i4 + s_j4 - s_ij)/2.
"""
import itertools
import numpy as np
import sympy as sp

rng = np.random.default_rng(1)
PAIRS = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
IDX = {p: k for k, p in enumerate(PAIRS)}


def sq(s, i, j):
    return s[IDX[(min(i, j), max(i, j))]]


def gram(s):
    return [[(sq(s, i, 4) + sq(s, j, 4) - (sq(s, i, j) if i != j else 0)) / 2
             for j in (1, 2, 3)] for i in (1, 2, 3)]


def dists(x):  # x: 4 points (rows), any dimension / complex ok
    return np.array([np.sum((x[i - 1] - x[j - 1]) ** 2) for i, j in PAIRS])


# ---- 1. Schoenberg: realizable <=> Gram PSD, Gram linear in s -------------
for _ in range(100):
    x = rng.normal(size=(4, 3))
    X = x[:3] - x[3]
    assert np.allclose(np.array(gram(dists(x))), X @ X.T)

# ---- 2. Rung L1: Heron factors linearly, Cayley-Menger does not -----------
R = sp.symbols('r12 r13 r14 r23 r24 r34', positive=True)
Pr = sp.Matrix(gram([r ** 2 for r in R]))
face = sp.factor_list(Pr[:2, :2].det())[1]
assert len(face) == 4 and all(sp.Poly(g, *R).total_degree() == 1 for g, _ in face)
cm = sp.factor_list(sp.expand(Pr.det()))[1]
assert len(cm) == 1 and sp.Poly(cm[0][0], *R).total_degree() == 6

# Kreuter: rectangle and its relabelling are realizable, their mean is not.
# order (r12,r13,r14,r23,r24,r34); nucleus = particle 4
a, b = 1.0, 2.0
c = np.hypot(a, b)
RA = np.array([c, b, a, a, b, c])
RB = np.array([c, a, b, b, a, c])
mineig = lambda r: np.linalg.eigvalsh(np.array(gram(r ** 2))).min()
assert mineig(RA) > -1e-12 and mineig(RB) > -1e-12 and mineig((RA + RB) / 2) < -0.1

# ---- 3. Rung L2: tangent cone at a collinear configuration is round -------
# Boundary points near collinear r0 have unit normals sweeping a continuum.
x0 = np.array([[0., 0, 0], [1, 0, 0], [2.5, 0, 0], [4, 0, 0]])
detP = sp.lambdify(R, Pr.det())
gradP = sp.lambdify(R, [sp.diff(Pr.det(), r) for r in R])
normals = []
for th in np.linspace(0, np.pi, 40, endpoint=False):  # th, th+pi are mirror images
    x = x0.copy()  # planar (xy) offsets of particles 2,3 in relative direction th
    x[1, 1], x[2, 1] = 1e-4 * np.cos(th), 1e-4 * np.sin(th)
    r = np.sqrt(dists(x))
    assert abs(detP(*r)) < 1e-9  # still on the boundary (planar)
    g = np.array(gradP(*r))
    normals.append(g / np.linalg.norm(g))
normals = np.array(normals)
gap = min(min(np.linalg.norm(normals[i] - normals[j]), np.linalg.norm(normals[i] + normals[j]))
          for i, j in itertools.combinations(range(40), 2))
assert gap > 1e-3  # 40 distinct boundary normals within 1e-4 of r0: not a polyhedral corner

# ---- 4. Rung L3+: polynomial squared distances on the closed orthant -------
# Cholesky P = L L^T on [0,inf)^3 x R^3 (degree 2); each (off-diag, diag) pair of L
# is (Re, Im) of z^2, z = q_a + i q_b in the closed quadrant (degree 4 overall).
IU = np.tril_indices(3)


def L_of_q(q):
    z = [complex(q[2 * k], q[2 * k + 1]) ** 2 for k in range(3)]
    L = np.zeros((3, 3))
    L[0, 0], L[1, 0] = z[0].imag, z[0].real
    L[1, 1], L[2, 0] = z[1].imag, z[1].real
    L[2, 2], L[2, 1] = z[2].imag, z[2].real
    return L


def q_of_P(P):  # inverse on the interior (and a preimage on the boundary)
    try:
        L = np.linalg.cholesky(P)
    except np.linalg.LinAlgError:  # singular (boundary) P: any limit factor will do
        L = np.linalg.cholesky(P + 1e-13 * np.trace(P) * np.eye(3))
    w = [np.sqrt(complex(o, max(d, 0.0))) for d, o in
         ((L[0, 0], L[1, 0]), (L[1, 1], L[2, 0]), (L[2, 2], L[2, 1]))]
    return np.array([c for x in w for c in (x.real, x.imag)])


def jac(f, x, h=1e-6):
    return np.array([(f(x + h * e) - f(x - h * e)) / (2 * h) for e in np.eye(len(x))]).T


for _ in range(300):
    q = rng.exponential(size=6)
    L = L_of_q(q)
    P = L @ L.T
    assert np.linalg.eigvalsh(P).min() > 0                     # interior -> interior
    assert np.allclose(q_of_P(P), q, rtol=1e-6)                # bijective on the interior
    qq = q.copy(); qq[rng.integers(6)] = 0                     # facet -> boundary
    ev = np.linalg.eigvalsh(L_of_q(qq) @ L_of_q(qq).T)
    assert abs(ev[0]) < 1e-12 * ev[-1]
    X = rng.normal(size=(3, rng.integers(1, 4)))               # every PSD (rank 1..3) is hit
    qX = q_of_P(X @ X.T)
    assert qX.min() >= -1e-12 and np.allclose(L_of_q(qX) @ L_of_q(qX).T, X @ X.T, atol=1e-5)
for _ in range(20):  # Jacobians entering the polynomial measure dtau ~ 8 l11^2 l22 dL
    q = rng.exponential(size=6)
    L = L_of_q(q)
    assert np.isclose(abs(np.linalg.det(jac(lambda x: L_of_q(x)[IU], q))),
                      np.prod([4 * (q[2 * k] ** 2 + q[2 * k + 1] ** 2) for k in range(3)]), rtol=1e-6)
    def P_of_l(v):
        M = np.zeros((3, 3)); M[IU] = v
        return (M @ M.T)[IU]
    assert np.isclose(abs(np.linalg.det(jac(P_of_l, L[IU]))), 8 * L[0, 0] ** 3 * L[1, 1] ** 2 * L[2, 2], rtol=1e-6)

# ---- 5. Rung L4: boundary variety  Zbar ~ (X x X)/G,  X = K3 -------------
# X subset P^5: alpha_ij^2 = z_j - z_i  <=>  three diagonal quadrics.
Q = np.array([  # coefficient rows over alpha^2 in PAIRS order
    [1, -1, 0, 1, 0, 0],   # a12^2 + a23^2 - a13^2
    [1, 0, -1, 0, 1, 0],   # a12^2 + a24^2 - a14^2
    [0, 1, -1, 0, 0, 1],   # a13^2 + a34^2 - a14^2
])


def kummer_pt(z):  # point of the affine cone over X from z in C^4
    return np.sqrt(np.array([z[j - 1] - z[i - 1] for i, j in PAIRS], dtype=complex))


for _ in range(100):
    al = kummer_pt(rng.normal(size=4) + 1j * rng.normal(size=4))
    be = kummer_pt(rng.normal(size=4) + 1j * rng.normal(size=4))
    assert np.allclose(Q @ al ** 2, 0) and np.allclose(Q @ be ** 2, 0)
    r = al * be                                     # r_ij = alpha_ij beta_ij
    P = np.array(gram(r ** 2))
    assert abs(np.linalg.det(P)) < 1e-10 * np.abs(r).max() ** 6   # lies on CM = 0
    # dominance: d(s)/d(z2..z4, w2..w4) has rank 5 = dim of the cone Z
    J = np.zeros((6, 6), complex)
    zeta, omega = al ** 2, be ** 2
    for k, (i, j) in enumerate(PAIRS):
        for m, coef in ((j, 1), (i, -1)):
            if m != 1:
                J[k, m - 2] += coef * omega[k]
                J[k, 3 + m - 2] += coef * zeta[k]
    sv = np.linalg.svd(J, compute_uv=False)
    assert sv[4] > 1e-8 * sv[0] and sv[5] < 1e-10 * sv[0]

# Singular locus of X: Jacobian rank of the quadrics on the support of alpha.
# Supports come from z in P^2 = C^4/(C + C*): intersections of the 6 lines.
E = sp.Matrix([[(m == j) - (m == i) for m in (1, 2, 3, 4)] for i, j in PAIRS])  # zeta = E z
ones = sp.Matrix([1, 1, 1, 1])
supports = {tuple(range(6))} | {tuple(k for k in range(6) if k != m) for m in range(6)}
for f, g in itertools.combinations(range(6), 2):  # intersection point of two lines in P^2
    z = [v for v in sp.Matrix.vstack(E[f, :], E[g, :]).nullspace()
         if sp.Matrix.hstack(v, ones).rank() == 2][0]
    vals = E * z
    supports.add(tuple(k for k in range(6) if vals[k] != 0))
Qs = sp.Matrix(Q)
sing = [S for S in supports if Qs[:, list(S)].rank() < 3]
assert sorted(len(S) for S in sing) == [3, 3, 3, 3]          # the 4 triple points only
assert sum(2 ** len(S) // 2 for S in sing) == 16             # 16 nodes (Kummer)
for S in sing:  # each node is A1: a rank-3 quadric in the 3 vanishing coords
    zero = [k for k in range(6) if k not in S]
    M = Qs[:, list(S)].T     # combos y with y^T Q vanishing on S
    y = M.nullspace()
    assert len(y) == 1
    quad = (y[0].T * Qs)
    assert all(quad[k] == 0 for k in S) and all(quad[k] != 0 for k in zero)

# Independent check that X is not rational(ly connected): Esnault's theorem
# says a smooth proj. rationally connected variety over F_p has #X = 1 mod p.
# Resolving the 16 nodes adds a multiple of p, so #X(F_p) mod p decides.
def count_cover(p, n=4):  # #points of {alpha_ij^2 = z_j - z_i, 1<=i<j<=n} in P
    chi = np.zeros(p, int)
    chi[(np.arange(1, p) ** 2) % p] = 1
    chi[1:][chi[1:] == 0] = -1
    pairs = list(itertools.combinations(range(n), 2))
    N = 0
    for zt in itertools.product(range(p), repeat=n - 1):
        z = (0,) + zt
        prod = 1
        for i, j in pairs:
            prod *= 1 + chi[(z[j] - z[i]) % p]
            if prod == 0:
                break
        N += prod
    return (N - 1) // (p - 1)


assert all(count_cover(p, 3) == p + 1 for p in (13, 17, 29))  # control: 3-body cover is a conic
# Fingerprint of a singular K3 with CM by Z[i]: #X = (p+1)^2 for p = 3 mod 4
# (supersingular), and #X = 1 + 2(a^2 - b^2) mod p for p = a^2 + b^2, a odd,
# i.e. the trace of the weight-3 CM newform eta(4 tau)^6 (h^{2,0} = 1).
residues = {}
for p in (7, 11, 13, 17, 19, 29, 37, 41):
    n = count_cover(p)
    residues[p] = n % p
    if p % 4 == 3:
        assert n == (p + 1) ** 2
    else:
        a, b = next((a, b) for a in range(1, p, 2) for b in range(0, p, 2) if a * a + b * b == p)
        assert residues[p] == (1 + 2 * (a * a - b * b)) % p != 1

# ---- 6. Minimal degree: no closed rank-1 triangle set carries an RP^2 class ----
# Over the 20 triangles of the 5-simplex: T is "closed" if every triangle whose
# three edges lie in E(T) belongs to T (forced by positivity for quadratic F).
# We check that no closed T admits s in Z^1(T;Z2), Z in Z_2(T;Z2) with <s u s, Z> = 1.
def _nullspace2(M):
    M = M.copy() % 2; rows, cols = M.shape; piv = []; r = 0
    for c in range(cols):
        p = next((i for i in range(r, rows) if M[i, c]), None)
        if p is None:
            continue
        M[[r, p]] = M[[p, r]]
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        piv.append(c); r += 1
    out = []
    for f in (c for c in range(cols) if c not in piv):
        x = np.zeros(cols, np.int64); x[f] = 1
        for i, pc in enumerate(piv):
            x[pc] = M[i, f]
        out.append(x % 2)
    return out


_E = list(itertools.combinations(range(6), 2)); _ei = {e: i for i, e in enumerate(_E)}
_T = list(itertools.combinations(range(6), 3))
_te = [[_ei[(a, b)], _ei[(b, c)], _ei[(a, c)]] for a, b, c in _T]
_tm = [sum(1 << i for i in te) for te in _te]


def rp2_class(tris):
    edges = sorted({e for t in tris for e in _te[t]}); ec = {e: k for k, e in enumerate(edges)}
    D = np.zeros((len(tris), len(edges)), np.int64)
    for r, t in enumerate(tris):
        for e in _te[t]:
            D[r, ec[e]] = 1
    for s in _nullspace2(D):
        sq = np.array([s[ec[_ei[_T[t][:2]]]] * s[ec[_ei[_T[t][1:]]]] for t in tris]) % 2
        if any((sq @ z) % 2 for z in _nullspace2(D.T.copy())):
            return True
    return False


hemi = [(0, 1, 2), (0, 2, 3), (0, 3, 4), (0, 4, 5), (0, 1, 5), (1, 2, 4), (2, 3, 5), (1, 3, 4), (1, 3, 5), (2, 4, 5)]
assert rp2_class([_T.index(f) for f in hemi])                  # positive control: RP^2_6
n_closed = 0
for mask in range(1, 1 << 20):
    em = 0
    for t in range(20):
        if mask >> t & 1:
            em |= _tm[t]
    if all((_tm[t] & ~em) or (mask >> t & 1) for t in range(20)):
        n_closed += 1
        assert not rp2_class([t for t in range(20) if mask >> t & 1])
assert n_closed == 6114

# ---- 7. Degree 2, global step: every RP^2-carrying cycle is a hemi-icosahedron, and the ----
# forced coplanarity of its 10 complementary triangles collapses all directions to a plane.
_B = np.zeros((15, 20), np.int64)
for _t in range(20):
    for _e in _te[_t]:
        _B[_e, _t] = 1
_cb = _nullspace2(_B)
assert len(_cb) == 10                                               # Z_2 of the 2-skeleton
pairs = []
for bits in itertools.product([0, 1], repeat=10):
    if not any(bits):
        continue
    z = np.zeros(20, np.int64)
    for b, v in zip(bits, _cb):
        if b:
            z ^= v
    tris = [t for t in range(20) if z[t]]
    edges = sorted({e for t in tris for e in _te[t]}); ec = {e: k for k, e in enumerate(edges)}
    D = np.zeros((len(tris), len(edges)), np.int64)
    for r, t in enumerate(tris):
        for e in _te[t]:
            D[r, ec[e]] = 1
    Z1 = _nullspace2(D)
    for sb in itertools.product([0, 1], repeat=len(Z1)):
        s = np.zeros(len(edges), np.int64)
        for b, v in zip(sb, Z1):
            if b:
                s ^= v
        if sum(s[ec[_ei[_T[t][:2]]]] * s[ec[_ei[_T[t][1:]]]] for t in tris) % 2 != 1:
            continue
        forced = [_T[t] for t in range(20) if not z[t] and all(e in ec for e in _te[t])
                  and sum(s[ec[e]] for e in _te[t]) % 2]
        pairs.append((tuple(_T[t] for t in tris), tuple(forced)))
assert len(pairs) == 384 and len({p[0] for p in pairs}) == 12
assert all(len(z) == 10 and set(f) == set(_T) - set(z) for z, f in pairs)


def _partitions(s):
    if not s:
        yield []
        return
    for p in _partitions(s[1:]):
        for i in range(len(p)):
            yield p[:i] + [[s[0]] + p[i]] + p[i + 1:]
        yield [[s[0]]] + p


for faces in {p[0] for p in pairs}:
    comp = set(_T) - set(faces)
    for P in _partitions(list(range(6))):
        cls = {v: i for i, b in enumerate(P) for v in b}
        if len(P) < 3:
            continue
        lines = [set(cls[v] for v in t) for t in comp if len({cls[v] for v in t}) == 3]
        merged = True
        while merged:
            merged = False
            for i, j in itertools.combinations(range(len(lines)), 2):
                if len(lines[i] & lines[j]) >= 2:
                    lines[i] |= lines.pop(j); merged = True
                    break
        assert any(len(L) == len(P) for L in lines)               # all directions in one plane

print("all checks passed;  #X(F_p) mod p =", residues)
