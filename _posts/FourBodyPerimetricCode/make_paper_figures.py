"""Figures for the paper: the round tangent cone (Sec. 4) and the exact point count (Remark rem:arith).

Writes figures/fig_cone.pdf and figures/fig_counts.pdf. Greyscale-safe; scatter layers rasterised.
"""
import os
import numpy as np
import sympy as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

plt.rcParams.update({"font.size": 9, "mathtext.fontset": "cm", "font.family": "serif",
                     "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6})
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(3)
GREY = "0.6"


def gram(x):
    """Gram matrix of x_i - x_last."""
    d = x[:-1] - x[-1]
    return d @ d.T


def transverse_block(P, vhat):
    """Schur complement of P onto the orthogonal complement of the collinear direction vhat,
    i.e. the block U^T H U of Lemma lem:tcone after eliminating the range direction."""
    Q = np.linalg.qr(np.column_stack([vhat, np.eye(len(vhat))]))[0]
    Pq = Q.T @ P @ Q
    return Pq[1:, 1:] - np.outer(Pq[1:, 0], Pq[0, 1:]) / Pq[0, 0]


def unit(x0):
    p = (x0[:-1] - x0[-1])[:, 0]
    return p / np.linalg.norm(p)


eps = 0.02

# (a) three bodies near a collinear triangle: the transverse block is a scalar >= 0 (a flat wall)
x3 = np.array([[0.0, 0, 0], [1.0, 0, 0], [2.6, 0, 0]])
v3 = unit(x3)
def tri_point(d):
    P = gram(x3 + d)
    return [(P - np.outer(v3, v3) * (v3 @ P @ v3))[0, 1] / eps, transverse_block(P, v3)[0, 0] / eps ** 2]
tri_in = np.array([tri_point(np.c_[rng.normal(size=(3, 2)) * eps, np.zeros(3)]) for _ in range(1500)])
tri_bd = np.array([tri_point(np.c_[rng.normal(size=3) * eps, np.zeros((3, 2))]) for _ in range(300)])
assert np.abs(tri_bd[:, 1]).max() < 1e-6 and tri_in[:, 1].min() > 0

# (b) four bodies near a collinear configuration: the 2x2 block [[a,b],[b,c]] lies in the round cone
#     z >= sqrt(x^2 + y^2), (x, y, z) = ((a-c)/2, b, (a+c)/2); plot the slice z = 1
x4 = np.array([[0.0, 0, 0], [1.0, 0, 0], [2.5, 0, 0], [4.0, 0, 0]])
v4 = unit(x4)
def slice_point(x):
    C = transverse_block(gram(x), v4)
    a, b, c = C[0, 0], C[0, 1], C[1, 1]
    return np.array([(a - c) / 2, b]) / ((a + c) / 2)
cone_in = np.array([slice_point(x4 + rng.normal(size=(4, 3)) * eps) for _ in range(2500)])
cone_bd = []
for _ in range(1500):                                   # coplanar: all four in a random plane through the line
    phi = rng.uniform(0, np.pi)
    d = np.outer(rng.normal(size=4) * eps, [0, np.cos(phi), np.sin(phi)])
    d[:, 0] += rng.normal(size=4) * eps
    cone_bd.append(slice_point(x4 + d))
cone_bd = np.array(cone_bd)
rad_bd = np.linalg.norm(cone_bd, axis=1)
sectors = np.unique((np.degrees(np.arctan2(cone_bd[:, 1], cone_bd[:, 0])) // 10).astype(int))
assert np.linalg.norm(cone_in, axis=1).max() < 1 and np.abs(rad_bd - 1).max() < 1e-6 and len(sectors) == 36
print(f"cone: max |radius-1| (coplanar) = {np.abs(rad_bd - 1).max():.1e}, sectors hit = {len(sectors)}/36")

fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.75), gridspec_kw={"width_ratios": [1.35, 1]})
ax = axs[0]
ax.scatter(*tri_in.T, s=2, color=GREY, lw=0, rasterized=True, label="non-degenerate")
ax.scatter(*tri_bd.T, s=7, color="k", marker="s", lw=0, rasterized=True, label="collinear")
ax.axhline(0, color="k", lw=0.6)
ax.set_xlim(*np.percentile(tri_in[:, 0], [1, 99])); ax.set_ylim(-0.6, np.percentile(tri_in[:, 1], 97))
ax.set_xlabel("along-line component"); ax.set_ylabel("transverse block")
ax.legend(frameon=False, loc="upper left", fontsize=8, markerscale=2)
ax.text(-0.02, 1.02, "(a)", transform=ax.transAxes, ha="right", va="bottom")
ax = axs[1]; ax.set_aspect("equal")
ax.scatter(*cone_in.T, s=2, color=GREY, lw=0, rasterized=True, label="non-degenerate")
ax.scatter(*cone_bd.T, s=4, color="k", lw=0, rasterized=True, label="coplanar")
ax.add_patch(Circle((0, 0), 1, fill=False, ec="k", lw=0.5, ls=(0, (3, 2))))
ax.set_xlim(-1.15, 1.15); ax.set_ylim(-1.15, 1.15)
ax.set_xlabel(r"$(a-c)/(a+c)$"); ax.set_ylabel(r"$2b/(a+c)$")
ax.text(-0.02, 1.02, "(b)", transform=ax.transAxes, ha="right", va="bottom")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_cone.pdf"), dpi=400)
plt.close(fig)


# ---- exact point count  #X(F_p) = (p+1)^2 + a_p^2 ---------------------------------------
def legendre_symbol(p):
    chi = -np.ones(p, int); chi[0] = 0
    chi[(np.arange(1, p) ** 2) % p] = 1
    return chi


def count_points(p, n=4):
    """Projective points of {alpha_ij^2 = z_j - z_i} over F_p (n = 4: X; n = 3: the three-body conic)."""
    chi = legendre_symbol(p)
    z = [np.zeros((p,) * (n - 1), int)] + list(np.indices((p,) * (n - 1)))
    prod = np.ones((p,) * (n - 1), int)
    for i in range(n):
        for j in range(i + 1, n):
            prod *= 1 + chi[(z[j] - z[i]) % p]
    return (int(prod.sum()) - 1) // (p - 1)


def a_p(p):
    """#E(F_p) = p + 1 - a_p for E: y^2 = x^3 - x."""
    x = np.arange(p)
    return -int(legendre_symbol(p)[(x ** 3 - x) % p].sum())


primes = [int(q) for q in sp.primerange(5, 200)]
dev_X = np.array([count_points(p) - (p + 1) ** 2 for p in primes])
dev_C = np.array([count_points(p, 3) - (p + 1) for p in primes])
ap2 = np.array([a_p(p) ** 2 for p in primes])
assert np.array_equal(dev_X, ap2) and not dev_C.any()
print(f"counts: #X = (p+1)^2 + a_p^2 for all {len(primes)} primes 5 <= p < 200; conic deviation all zero")

fig, axs = plt.subplots(1, 2, figsize=(6.3, 2.5), gridspec_kw={"width_ratios": [1.5, 1]})
ax = axs[0]
ax.plot(primes, dev_X, "o", ms=6.5, mfc="none", mec="k", mew=0.7, label=r"$|X(\mathbb{F}_p)|-(p+1)^2$ (count on $X$)")
ax.plot(primes, ap2, "+", ms=5, color="0.35", mew=0.9, label=r"$a_p^2$ (count on $E_i$)")
ax.plot(primes, dev_C, "x", ms=3.5, color=GREY, mew=0.8, label=r"$|C(\mathbb{F}_p)|-(p+1)$, three bodies")
ax.set_xlabel(r"$p$"); ax.legend(frameon=False, loc="upper left", fontsize=7.5)
ax.text(-0.02, 1.02, "(a)", transform=ax.transAxes, ha="right", va="bottom")
ax = axs[1]
m = ap2.max() * 1.05
ax.plot([0, m], [0, m], color=GREY, lw=0.6, ls=(0, (4, 3)))
ax.plot(ap2, dev_X, "o", ms=5, mfc="none", mec="k", mew=0.7)
ax.set_xlabel(r"$a_p^2$"); ax.set_ylabel(r"$|X(\mathbb{F}_p)|-(p+1)^2$")
ax.set_xlim(-25, m); ax.set_ylim(-25, m)
ax.text(-0.02, 1.02, "(b)", transform=ax.transAxes, ha="right", va="bottom")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_counts.pdf"))
plt.close(fig)
print("figures written to", OUT)
