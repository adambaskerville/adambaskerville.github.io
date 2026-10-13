# %% [markdown]
# # Do four-body perimetric coordinates exist?
# Companion code for the T>T blog post. Each cell matches a section of the post.
# Requirements: numpy, scipy, sympy (all pip-installable, nothing exotic).

# %%
import itertools
import numpy as np
import sympy as sp

rng = np.random.default_rng(1)

# Six interparticle distances of four particles, always in this order
PAIRS = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
IDX = {p: k for k, p in enumerate(PAIRS)}


def sq(s, i, j):
    """Squared distance between particles i and j from the vector s = r**2."""
    return s[IDX[(min(i, j), max(i, j))]]


def gram(s):
    """Gram matrix of the three vectors X_i = x_i - x_4 (particle 4 at the origin).
    It is LINEAR in the squared distances s."""
    return np.array([[(sq(s, i, 4) + sq(s, j, 4) - (sq(s, i, j) if i != j else 0)) / 2
                      for j in (1, 2, 3)] for i in (1, 2, 3)])


def dists2(x):
    """Squared distances of four points x (4 x dim array)."""
    return np.array([np.sum((x[i - 1] - x[j - 1]) ** 2) for i, j in PAIRS])


def realizable(r, tol=1e-12):
    """Six distances describe four points in 3D  <=>  the Gram matrix is PSD."""
    return np.linalg.eigvalsh(gram(np.asarray(r) ** 2)).min() > -tol


# Schoenberg's theorem, checked on random tetrahedra
for _ in range(1000):
    x = rng.normal(size=(4, 3))
    X = x[:3] - x[3]
    assert np.allclose(gram(dists2(x)), X @ X.T)
print("Gram matrix = X X^T for 1000 random tetrahedra: OK")

# %% [markdown]
# ## Rung 1: no linear coordinates (the domain is not convex)

# %%
a, b = 1.0, 2.0
c = np.hypot(a, b)                       # diagonal of an a x b rectangle
RA = np.array([c, b, a, a, b, c])        # rectangle, particle 4 at a corner
RB = np.array([c, a, b, b, a, c])        # the same rectangle, relabelled
mid = (RA + RB) / 2

print("rectangle A realizable:", realizable(RA))
print("rectangle B realizable:", realizable(RB))
print("their average realizable:", realizable(mid))
print("smallest Gram eigenvalue of the average:", np.linalg.eigvalsh(gram(mid ** 2)).min())

# %% [markdown]
# ## Why three bodies are easy and four are not: Heron vs Cayley-Menger

# %%
r12, r13, r14, r23, r24, r34 = R = sp.symbols('r12 r13 r14 r23 r24 r34', positive=True)
P = sp.Matrix([[r14**2, (r14**2 + r24**2 - r12**2) / 2, (r14**2 + r34**2 - r13**2) / 2],
               [(r14**2 + r24**2 - r12**2) / 2, r24**2, (r24**2 + r34**2 - r23**2) / 2],
               [(r14**2 + r34**2 - r13**2) / 2, (r24**2 + r34**2 - r23**2) / 2, r34**2]])

triangle = sp.factor(P[:2, :2].det())    # one face of the tetrahedron (a triangle)
print("triangle (Heron):", triangle)

tet = sp.factor_list(sp.expand(P.det()))  # the whole tetrahedron (Cayley-Menger)
print("tetrahedron: number of factors =", len(tet[1]),
      "| degree =", sp.Poly(tet[1][0][0], *R).total_degree())

# %% [markdown]
# ## Rung 2: no coordinates that are smooth up to the boundary

# %%
detP = sp.lambdify(R, P.det())
grad = sp.lambdify(R, [sp.diff(P.det(), v) for v in R])

x0 = np.array([[0., 0, 0], [1, 0, 0], [2.5, 0, 0], [4, 0, 0]])   # four collinear particles
normals = []
for th in np.linspace(0, np.pi, 40, endpoint=False):
    x = x0.copy()
    x[1, 1], x[2, 1] = 1e-4 * np.cos(th), 1e-4 * np.sin(th)        # nudge into a plane
    r = np.sqrt(dists2(x))
    g = np.array(grad(*r))
    normals.append(g / np.linalg.norm(g))
normals = np.array(normals)
gap = min(min(np.linalg.norm(u - v), np.linalg.norm(u + v))
          for u, v in itertools.combinations(normals, 2))
print("40 boundary normals within 1e-4 of a collinear configuration")
print("smallest angle-gap between any two of them:", round(gap, 4), "(all distinct)")

# %% [markdown]
# ## Rung 2, with data: near four collinear particles the boundary is a round cone

# %%
def transverse_block(G, vhat):
    """2x2 Schur complement of the Gram matrix across the collinear direction vhat.
    PSD exactly when G is, singular exactly when the configuration is flat."""
    Q = np.linalg.qr(np.column_stack([vhat, np.eye(3)]))[0]
    Gq = Q.T @ G @ Q
    return Gq[1:, 1:] - np.outer(Gq[1:, 0], Gq[0, 1:]) / Gq[0, 0]


def cone_slice(x, vhat):
    """PSD 2x2 matrices [[a, b], [b, c]] are the round cone z >= sqrt(x^2 + y^2) with
    (x, y, z) = ((a - c)/2, b, (a + c)/2). Return the point on the slice z = 1."""
    C = transverse_block(gram(dists2(x)), vhat)
    a, b, c = C[0, 0], C[0, 1], C[1, 1]
    return np.array([(a - c) / 2, b]) / ((a + c) / 2)


eps = 0.02
vhat = (x0[:3] - x0[3])[:, 0] / np.linalg.norm((x0[:3] - x0[3])[:, 0])
genuine = np.array([cone_slice(x0 + rng.normal(size=(4, 3)) * eps, vhat) for _ in range(2000)])
flat = []
for _ in range(2000):
    phi = rng.uniform(0, np.pi)                       # a random plane through the line
    d = np.outer(rng.normal(size=4) * eps, [0, np.cos(phi), np.sin(phi)])
    d[:, 0] += rng.normal(size=4) * eps
    flat.append(cone_slice(x0 + d, vhat))
flat = np.array(flat)
radius_flat = np.linalg.norm(flat, axis=1)
angles_hit = np.unique((np.degrees(np.arctan2(flat[:, 1], flat[:, 0])) // 10).astype(int))
print("genuine tetrahedra strictly inside the circle:", bool(np.linalg.norm(genuine, axis=1).max() < 1))
print("flat tetrahedra: largest |radius - 1| =", f"{np.abs(radius_flat - 1).max():.1e}")
print("10-degree sectors of the circle hit by flat tetrahedra:", len(angles_hit), "of 36")
assert np.linalg.norm(genuine, axis=1).max() < 1 and np.abs(radius_flat - 1).max() < 1e-6

# %% [markdown]
# ## Rung 3: independent coordinates DO exist (squared distances polynomial)

# %%
def L_of_q(q):
    """Lower-triangular Cholesky factor built from six coordinates q >= 0.
    Each (off-diagonal, diagonal) pair is (Re, Im) of z^2 with z = q_a + i q_b."""
    z = [complex(q[2 * k], q[2 * k + 1]) ** 2 for k in range(3)]
    L = np.zeros((3, 3))
    L[0, 0], L[1, 0] = z[0].imag, z[0].real
    L[1, 1], L[2, 0] = z[1].imag, z[1].real
    L[2, 2], L[2, 1] = z[2].imag, z[2].real
    return L


def s_of_q(q):
    """Squared distances (degree-4 polynomials in q) from the Gram matrix L L^T."""
    G = L_of_q(q) @ L_of_q(q).T
    s = np.empty(6)
    s[IDX[(1, 4)]], s[IDX[(2, 4)]], s[IDX[(3, 4)]] = np.diag(G)
    for i, j in ((1, 2), (1, 3), (2, 3)):
        s[IDX[(i, j)]] = G[i - 1, i - 1] + G[j - 1, j - 1] - 2 * G[i - 1, j - 1]
    return s


def q_of_s(s):
    """Inverse map: Cholesky factor, then the principal square root of each pair."""
    L = np.linalg.cholesky(gram(s))
    w = [np.sqrt(complex(o, d)) for d, o in ((L[0, 0], L[1, 0]), (L[1, 1], L[2, 0]), (L[2, 2], L[2, 1]))]
    return np.array([c for x in w for c in (x.real, x.imag)])


worst = 0
for _ in range(1000):
    q = rng.uniform(0.2, 3.0, size=6)
    s = s_of_q(q)
    assert realizable(np.sqrt(s))
    worst = max(worst, np.max(np.abs(q_of_s(s) / q - 1)))
print("1000 random q in the open orthant -> valid tetrahedra -> back to q")
print("worst relative round-trip error:", f"{worst:.1e}")

# %% [markdown]
# ## Rung 4: the hidden K3 surface

# %%
# X: three diagonal quadrics in P^5  (alpha_ij^2 = z_j - z_i, written without z)
QX = np.array([[1, -1, 0, 1, 0, 0],    # a12^2 + a23^2 - a13^2 = 0
               [1, 0, -1, 0, 1, 0],    # a12^2 + a24^2 - a14^2 = 0
               [0, 1, -1, 0, 0, 1]])   # a13^2 + a34^2 - a14^2 = 0


def point_on_X(z):
    """A point of X from four complex numbers z_1..z_4 (null coordinates)."""
    return np.sqrt(np.array([z[j - 1] - z[i - 1] for i, j in PAIRS], dtype=complex))


for _ in range(1000):
    al = point_on_X(rng.normal(size=4) + 1j * rng.normal(size=4))
    be = point_on_X(rng.normal(size=4) + 1j * rng.normal(size=4))
    assert np.allclose(QX @ al ** 2, 0) and np.allclose(QX @ be ** 2, 0)
    r = al * be                                        # r_ij = alpha_ij * beta_ij
    assert abs(np.linalg.det(gram(r ** 2))) < 1e-10 * np.abs(r).max() ** 6
print("1000 random pairs (alpha, beta) on X x X give r with Cayley-Menger = 0: OK")

# %%
# Singular points of X: rank of the quadric coefficients on the support of alpha
E = sp.Matrix([[(m == j) - (m == i) for m in (1, 2, 3, 4)] for i, j in PAIRS])
ones = sp.Matrix([1, 1, 1, 1])
supports = {tuple(range(6))} | {tuple(k for k in range(6) if k != m) for m in range(6)}
for f, g in itertools.combinations(range(6), 2):
    z = [v for v in sp.Matrix.vstack(E[f, :], E[g, :]).nullspace()
         if sp.Matrix.hstack(v, ones).rank() == 2][0]
    vals = E * z
    supports.add(tuple(k for k in range(6) if vals[k] != 0))
singular = [S for S in supports if sp.Matrix(QX)[:, list(S)].rank() < 3]
print("singular supports:", singular)
print("number of singular points:", sum(2 ** len(S) // 2 for S in singular))

# %%
# Arithmetic fingerprint: count points of X over the finite field F_p
def legendre_symbol(p):
    chi = -np.ones(p, int)
    chi[0] = 0
    chi[(np.arange(1, p) ** 2) % p] = 1
    return chi


def count_points(p, n=4):
    """Projective points of X over F_p: fix z_1 = 0, run over the other z's, and weight
    each choice by the number of square roots alpha_ij = sqrt(z_j - z_i) it allows."""
    chi = legendre_symbol(p)
    z = [np.zeros((p,) * (n - 1), int)] + list(np.indices((p,) * (n - 1)))
    prod = np.ones((p,) * (n - 1), int)
    for i, j in itertools.combinations(range(n), 2):
        prod *= 1 + chi[(z[j] - z[i]) % p]          # number of square roots
    return (int(prod.sum()) - 1) // (p - 1)         # projective points


print(" p   #X(F_p)   mod p   1 + 2(a^2 - b^2) mod p")
for p in (7, 11, 13, 17, 19, 29, 37, 41):
    n = count_points(p)
    if p % 4 == 3:
        print(f"{p:3d} {n:8d} {n % p:7d}   (p = 3 mod 4: #X = (p+1)^2 = {(p + 1) ** 2})")
    else:
        a_, b_ = next((a_, b_) for a_ in range(1, p, 2) for b_ in range(0, p, 2) if a_ * a_ + b_ * b_ == p)
        print(f"{p:3d} {n:8d} {n % p:7d}   {(1 + 2 * (a_ * a_ - b_ * b_)) % p:5d}")
print("control (three bodies, a conic):", [bool(count_points(p, 3) == p + 1) for p in (13, 17, 29)])

# %%
# The exact formula: |X(F_p)| = (p+1)^2 + a_p^2, with a_p read off the doughnut y^2 = x^3 - x
def a_p(p):
    """#E(F_p) = p + 1 - a_p for the curve E: y^2 = x^3 - x."""
    x = np.arange(p)
    return -int(legendre_symbol(p)[(x ** 3 - x) % p].sum())


primes = list(sp.primerange(5, 200))
assert all(count_points(p) == (p + 1) ** 2 + a_p(p) ** 2 for p in primes)
assert all(count_points(p, 3) == p + 1 for p in primes)
print(f"|X(F_p)| = (p+1)^2 + a_p^2 for all {len(primes)} primes 5 <= p < 200")
print("first few a_p:", {p: a_p(p) for p in primes[:8]})

# %% [markdown]
# ## The degree-2 hunt: every RP^2-carrying cycle is a hemi-icosahedron

# %%
def nullspace_mod2(M):
    M = M.copy() % 2
    rows, cols = M.shape
    piv, r = [], 0
    for c in range(cols):
        p = next((i for i in range(r, rows) if M[i, c]), None)
        if p is None:
            continue
        M[[r, p]] = M[[p, r]]
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        piv.append(c)
        r += 1
    out = []
    for f in (c for c in range(cols) if c not in piv):
        x = np.zeros(cols, np.int64)
        x[f] = 1
        for i, pc in enumerate(piv):
            x[pc] = M[i, f]
        out.append(x % 2)
    return out


EDGES = list(itertools.combinations(range(6), 2))
EI = {e: i for i, e in enumerate(EDGES)}
TRIS = list(itertools.combinations(range(6), 3))
TE = [[EI[(a, b)], EI[(b, c)], EI[(a, c)]] for a, b, c in TRIS]

boundary = np.zeros((15, 20), np.int64)
for t in range(20):
    for e in TE[t]:
        boundary[e, t] = 1
cycle_basis = nullspace_mod2(boundary)
print("independent mod-2 2-cycles of the 5-simplex:", len(cycle_basis))

hits = []
for bits in itertools.product([0, 1], repeat=10):
    if not any(bits):
        continue
    z = np.zeros(20, np.int64)
    for bit, v in zip(bits, cycle_basis):
        if bit:
            z ^= v
    tris = [t for t in range(20) if z[t]]
    edges = sorted({e for t in tris for e in TE[t]})
    col = {e: k for k, e in enumerate(edges)}
    D = np.zeros((len(tris), len(edges)), np.int64)
    for row, t in enumerate(tris):
        for e in TE[t]:
            D[row, col[e]] = 1
    for sbits in itertools.product([0, 1], repeat=len(nullspace_mod2(D))):
        sigma = np.zeros(len(edges), np.int64)
        for bit, v in zip(sbits, nullspace_mod2(D)):
            if bit:
                sigma ^= v
        cup = sum(sigma[col[EI[TRIS[t][:2]]]] * sigma[col[EI[TRIS[t][1:]]]] for t in tris) % 2
        if cup == 1:
            hits.append(tuple(TRIS[t] for t in tris))
print("(cycle, sign) pairs carrying the RP^2 class:", len(hits))
print("distinct supporting cycles:", len(set(hits)), "| triangles in each:", {len(h) for h in hits})

# %% [markdown]
# ## The Galois obstruction: a sparse helium scheme that cannot see the sign of 1/r12
# (uses he_graded.py from the same folder)

# %%
import scipy.linalg as sla
from he_graded import operators

for lam in (+1.0, -1.0):
    A, B, _ = operators(6, 3, Z=2.0, kappa=1.0, lam_ee=lam)
    w = sla.eigvals(A.toarray(), B.toarray()) * 4
    w = np.sort(w[np.isfinite(w) & (np.abs(w.imag) < 1e-9)].real)
    print(f"e-e coupling {lam:+.0f}: lowest real eigenvalues", np.round(w[:4], 6))
