"""Coordinate diagrams and data visualisations for the four-body perimetric post."""
import numpy as np
import sympy as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon, FancyArrowPatch, Arc

SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE = "#2a78d6", "#eb6834"
UC = ["#2a78d6", "#1baf7a", "#4a3aa7"]                    # u1, u2, u3 (validated set)
OUT = "../../assets/img/four_body_perimetric/"
plt.rcParams.update({"font.size": 12, "text.color": INK, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2})
rng = np.random.default_rng(3)


def save(fig, name):
    fig.savefig(OUT + name, dpi=160, facecolor=SURF, bbox_inches="tight")
    plt.close(fig)


def incircle(A, B, C):
    a, b, c = np.linalg.norm(B - C), np.linalg.norm(C - A), np.linalg.norm(A - B)
    I = (a * A + b * B + c * C) / (a + b + c)
    s = (a + b + c) / 2
    return I, np.sqrt((s - a) * (s - b) * (s - c) / s)


def tangent_point(P, Q, I):
    d = (Q - P) / np.linalg.norm(Q - P)
    return P + d * np.dot(I - P, d)


def perimetric_triangle(ax, A, B, C, labels, lw=6, vlabels=None, fs=15):
    I, r = incircle(A, B, C)
    ax.add_patch(Circle(I, r, fill=False, ec=INK2, lw=1.4, ls=(0, (3, 2))))
    V = [A, B, C]
    for k, P in enumerate(V):
        for Q in (V[(k + 1) % 3], V[(k + 2) % 3]):
            T = tangent_point(P, Q, I)
            ax.plot(*np.array([P, T]).T, color=UC[k], lw=lw, solid_capstyle="butt")
            ax.plot(*T, "o", ms=7, color=SURF, mec=INK2, mew=1.4, zorder=5)
    for k, P in enumerate(V):
        ax.plot(*P, "o", ms=10, color=INK, mec=SURF, mew=2, zorder=6)
        off = (P - I) / np.linalg.norm(P - I)
        if vlabels:
            ax.text(*(P + off * 0.62), vlabels[k], ha="center", va="center", fontsize=fs, color=INK)
    return I


# ===== A1: three-body perimetric coordinates =====================================
fig, ax = plt.subplots(figsize=(7.5, 5.2), facecolor=SURF)
ax.set_facecolor(SURF); ax.axis("off"); ax.set_aspect("equal")
A, B, C = np.array([0.0, 0.0]), np.array([4.4, 0.0]), np.array([1.5, 3.1])
I = perimetric_triangle(ax, A, B, C, None, vlabels=["particle 1", "particle 2", "particle 3"])
for P, Q, lab in [(A, B, "$r_{12}=u_1+u_2$"), (B, C, "$r_{23}=u_2+u_3$"), (C, A, "$r_{13}=u_1+u_3$")]:
    mid = (P + Q) / 2; nrm = mid - I; nrm /= np.linalg.norm(nrm)
    ha = "center" if abs(nrm[0]) < 0.3 else ("left" if nrm[0] > 0 else "right")
    ax.text(*(mid + nrm * 0.3), lab, ha=ha, va="center", fontsize=13, color=INK)
for k, P in enumerate([A, B, C]):                           # u_i labels along a tangent
    T = tangent_point(P, [B, C, A][k], I); m = (P + T) / 2
    nrm = np.array([-(T - P)[1], (T - P)[0]]); nrm /= np.linalg.norm(nrm)
    if np.dot(nrm, m - I) > 0: nrm = -nrm
    ax.text(*(m + nrm * 0.24), f"$u_{k + 1}$", ha="center", va="center", fontsize=14, color=INK)
ax.set_xlim(-1.6, 6.0); ax.set_ylim(-0.95, 3.75)
save(fig, "coords_three_body.png")

# ===== A2: four-body interparticle distances and the Gram vectors ===================
fig = plt.figure(figsize=(12, 5.2), facecolor=SURF)
V = np.array([[1.4, -0.2, 0.0], [0.3, 1.2, 0.0], [0.9, 0.5, 1.15], [0.0, 0.0, 0.0]])   # particles 1,2,3,4
cen = V.mean(0)
def frame(ax):
    ax.view_init(20, -70); ax.axis("off"); ax.set_box_aspect((1, 1, 0.85))
    ax.set_xlim(-0.1, 1.5); ax.set_ylim(-0.3, 1.3); ax.set_zlim(-0.05, 1.25)
ax = fig.add_subplot(1, 2, 1, projection="3d"); ax.set_facecolor(SURF)
for (i, j) in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]:
    ax.plot(*V[[i, j]].T, color=ORANGE, lw=2.5)
    m = V[[i, j]].mean(0); m = m + (m - cen) * 0.25
    ax.text(*m, f"$r_{{{min(i, j) + 1}{max(i, j) + 1}}}$", fontsize=13, color=INK, ha="center", va="center")
ax.scatter(*V.T, s=110, color=INK, depthshade=False)
for k, v in enumerate(V):
    off = (v - cen) * 0.28 if k != 1 else np.array([-0.05, 0.15, 0.12])
    ax.text(*(v + off), str(k + 1), fontsize=14, color=INK, ha="center", va="center")
frame(ax); ax.set_title("six distances: the edges of a tetrahedron", fontsize=13, color=INK)
ax = fig.add_subplot(1, 2, 2, projection="3d"); ax.set_facecolor(SURF)
for (i, j) in [(0, 1), (0, 2), (1, 2)]:
    ax.plot(*V[[i, j]].T, color=GRID, lw=1.5)
for k in range(3):
    ax.quiver(0, 0, 0, *V[k], color=UC[k], lw=3, arrow_length_ratio=0.1)
    off = (V[k] - cen) * 0.28 if k != 1 else np.array([-0.05, 0.15, 0.12])
    ax.text(*(V[k] + off), f"$X_{k + 1}$", fontsize=15, color=INK, ha="center", va="center")
ax.scatter([0], [0], [0], s=110, color=INK, depthshade=False)
ax.text(*(V[3] + (V[3] - cen) * 0.28), "4", fontsize=14, color=INK, ha="center", va="center")
frame(ax); ax.set_title("Gram matrix $P_{ij}=X_i\\cdot X_j$, with particle 4 at the origin", fontsize=13, color=INK)
save(fig, "coords_four_body.png")

# ===== A3: Frolov's coordinates: perimetric coordinates on every face ================
def apex(Pi, Pj, Pk, di, dj):
    """Unfold: the point at distances di, dj from Pi, Pj, on the far side of line PiPj from Pk."""
    L = np.linalg.norm(Pj - Pi); u = (Pj - Pi) / L; n = np.array([-u[1], u[0]])
    x = (di ** 2 - dj ** 2 + L ** 2) / (2 * L); y = np.sqrt(max(di ** 2 - x ** 2, 0))
    cand = Pi + x * u + y * n
    return cand if np.dot(cand - Pi, n) * np.dot(Pk - Pi, n) < 0 else Pi + x * u - y * n
T3 = np.array([[0.0, 0, 0], [3.0, 0.0, 0], [1.4, 2.4, 0], [1.5, 0.8, 1.6]])     # a real tetrahedron 1,2,3,4
d = lambda i, j: np.linalg.norm(T3[i] - T3[j])
P1, P2, P3 = T3[0, :2], T3[1, :2], T3[2, :2]
net = [(P1, P2, P3)]
for (i, j, k) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
    Pi, Pj, Pk = T3[i, :2], T3[j, :2], T3[k, :2]
    net.append((Pi, Pj, apex(Pi, Pj, Pk, d(i, 3), d(j, 3))))
fig, ax = plt.subplots(figsize=(8.5, 7.2), facecolor=SURF)
ax.set_facecolor(SURF); ax.axis("off"); ax.set_aspect("equal")
for k, (P, Q, Rr) in enumerate(net):
    ax.add_patch(Polygon([P, Q, Rr], closed=True, fill=True, fc=SURF if k == 0 else "#f3f2ee", ec=INK2, lw=1.5))
    I, r = incircle(P, Q, Rr)
    ax.add_patch(Circle(I, r, fill=False, ec=INK2, lw=1.1, ls=(0, (3, 2))))
    for X_, Y_ in [(P, Q), (Q, Rr), (Rr, P)]:
        ax.plot(*tangent_point(X_, Y_, I), "o", ms=6, color=SURF, mec=ORANGE, mew=1.6, zorder=4)
c0 = (P1 + P2 + P3) / 3
for P, lab in [(P1, "1"), (P2, "2"), (P3, "3")]:
    ax.plot(*P, "o", ms=10, color=INK, mec=SURF, mew=2, zorder=6)
    ax.text(*(P + (P - c0) / np.linalg.norm(P - c0) * 0.3), lab, fontsize=14, ha="center", va="center")
for (_, _, Rr) in net[1:]:
    ax.plot(*Rr, "o", ms=10, color=INK, mec=SURF, mew=2, zorder=6)
    ax.text(*(Rr + (Rr - c0) / np.linalg.norm(Rr - c0) * 0.3), "4", fontsize=14, ha="center", va="center")
allpts = np.vstack([np.array(f) for f in net])
ax.text(allpts[:, 0].mean(), allpts[:, 1].min() - 0.5,
        "the tetrahedron unfolded flat: each face gets its own incircle and three tangent lengths, 12 numbers.\n"
        "Every edge is shared by two faces, giving 6 linear constraints. Keep 6 of the 12: Frolov's linear coordinates.",
        ha="center", va="top", fontsize=11.5, color=INK)
ax.set_xlim(allpts[:, 0].min() - 0.6, allpts[:, 0].max() + 0.6)
ax.set_ylim(allpts[:, 1].min() - 1.6, allpts[:, 1].max() + 0.6)
save(fig, "coords_frolov.png")

# ===== A4: Cholesky / normal coordinates, and the z^2 trick ===========================
fig = plt.figure(figsize=(15, 5), facecolor=SURF)
ax = fig.add_subplot(1, 3, 1, projection="3d"); ax.set_facecolor(SURF)
L = np.array([[1.3, 0, 0], [0.45, 1.0, 0], [0.7, 0.35, 0.9]])
for k in range(3):
    ax.quiver(0, 0, 0, *L[k], color=UC[k], lw=3, arrow_length_ratio=0.08)
ax.plot([L[1, 0], L[1, 0]], [0, L[1, 1]], [0, 0], color=UC[1], lw=1, ls=(0, (3, 2)))
ax.plot([L[2, 0], L[2, 0]], [L[2, 1], L[2, 1]], [0, L[2, 2]], color=UC[2], lw=1, ls=(0, (3, 2)))
ax.plot([0, L[2, 0]], [L[2, 1], L[2, 1]], [0, 0], color=UC[2], lw=1, ls=(0, (3, 2)))
for v, lab in [(L[0] * 1.08, "$X_1=(l_{11},0,0)$"), (L[1] * 1.12, "$X_2=(l_{21},l_{22},0)$"),
               (L[2] * 1.1, "$X_3=(l_{31},l_{32},l_{33})$")]:
    ax.text(*v, lab, fontsize=12, color=INK)
for v in np.eye(3) * 1.5:
    ax.plot(*np.array([[0, 0, 0], v]).T, color=GRID, lw=1)
ax.view_init(20, -60); ax.axis("off"); ax.set_box_aspect((1, 1, 0.8))
ax.set_title("pin $X_1$ to an axis and $X_2$ to a plane: six numbers $l_{ij}$\n"
             "$l_{11},l_{22},l_{33}\\geq 0$, the rest any real number", fontsize=12, color=INK)
t = np.linspace(0, 1.6, 9); a_ = np.linspace(0, 1.6, 100)
for col, f, title in [(2, lambda z: z, "before: the quarter-plane $q_a, q_b \\geq 0$"),
                      (3, lambda z: z ** 2, "after $z \\mapsto z^2$: the whole upper half-plane\n"
                                             "(real, imag) $= (l_{21}, l_{11})$, so $l_{21}$ can be negative")]:
    ax = fig.add_subplot(1, 3, col); ax.set_facecolor(SURF); ax.set_aspect("equal")
    for c0 in t:                                             # grid of the quarter-plane and its image
        for zz in (a_ + 1j * c0, c0 + 1j * a_):
            w = f(zz)
            ax.plot(w.real, w.imag, color=BLUE, lw=0.9, alpha=0.55)
    ax.axhline(0, color=INK2, lw=1); ax.axvline(0, color=INK2, lw=1)
    ax.set_xlim(-2.7, 2.7); ax.set_ylim(-0.3, 2.9); ax.axis("off")
    ax.set_title(title, fontsize=12, color=INK)
save(fig, "coords_cholesky.png")

# ===== A5: null coordinates of a flat configuration ==================================
fig, ax = plt.subplots(figsize=(8.5, 5.4), facecolor=SURF); ax.set_facecolor(SURF)
Z = np.array([0.4 + 1.2j, 2.3 + 1.9j, 3.1 + 0.6j, 1.4 + 0.35j])
ax.axhline(0, color=INK2, lw=1); ax.axvline(0, color=INK2, lw=1)
for k, zz in enumerate(Z):
    ax.plot(zz.real, zz.imag, "o", ms=11, color=BLUE, mec=SURF, mew=2, zorder=5)
    ax.text(zz.real, zz.imag + 0.22, f"$z_{k + 1}$", ha="center", fontsize=14)
    ax.plot(zz.real, -zz.imag, "o", ms=9, color=SURF, mec=ORANGE, mew=2, zorder=5)
    ax.text(zz.real, -zz.imag - 0.32, f"$w_{k + 1}$", ha="center", fontsize=14)
    ax.plot([zz.real, zz.real], [zz.imag, -zz.imag], color=GRID, lw=1, ls=(0, (3, 3)))
ax.plot([Z[0].real, Z[1].real], [Z[0].imag, Z[1].imag], color=BLUE, lw=2)
ax.plot([Z[0].real, Z[1].real], [-Z[0].imag, -Z[1].imag], color=ORANGE, lw=2)
ax.text(4.0, 1.4, "$z=x+iy$  (the plane)\n$w=x-iy$  (its mirror image)\n\n"
        "$r_{12}^2=(z_1-z_2)(w_1-w_2)$\n$r_{12}=\\sqrt{z_1-z_2}\\cdot\\sqrt{w_1-w_2}$",
        va="center", fontsize=13, color=INK)
ax.set_xlim(-0.4, 7.0); ax.set_ylim(-2.6, 2.6); ax.set_aspect("equal"); ax.axis("off")
save(fig, "coords_null.png")

# ===== A6: helium in (r1, r2, cos theta12) ==========================================
fig, ax = plt.subplots(figsize=(6.5, 4.6), facecolor=SURF); ax.set_facecolor(SURF)
N, E1, E2 = np.array([0, 0]), np.array([3.0, 0.6]), np.array([1.0, 2.6])
for P, Q, col, lab, off in [(N, E1, UC[0], "$r_1$", (0.05, -0.3)), (N, E2, UC[1], "$r_2$", (-0.35, 0.05)),
                            (E1, E2, ORANGE, "$r_{12}=\\sqrt{r_1^2+r_2^2-2r_1r_2\\cos\\theta_{12}}$", (1.55, 0.25))]:
    ax.plot(*np.array([P, Q]).T, color=col, lw=3)
    ax.text(*((P + Q) / 2 + np.array(off)), lab, fontsize=13, ha="center")
ang1, ang2 = np.degrees(np.arctan2(E1[1], E1[0])), np.degrees(np.arctan2(E2[1], E2[0]))
ax.add_patch(Arc(N, 1.4, 1.4, theta1=ang1, theta2=ang2, color=INK2, lw=1.5))
ax.text(0.85, 0.75, "$\\theta_{12}$", fontsize=14)
for P, lab in [(N, "nucleus"), (E1, "electron 1"), (E2, "electron 2")]:
    ax.plot(*P, "o", ms=13 if lab == "nucleus" else 10, color=INK, mec=SURF, mew=2, zorder=5)
    ax.text(P[0], P[1] - 0.38 if lab != "electron 2" else P[1] + 0.22, lab, ha="center", fontsize=12, color=INK2)
ax.set_xlim(-0.8, 4.6); ax.set_ylim(-0.7, 3.2); ax.set_aspect("equal"); ax.axis("off")
save(fig, "coords_helium.png")


# ===== B1: the round cone, from real configurations ===================================
def gram(X):        # X: positions (n, 3); Gram of x_i - x_last
    D = X[:-1] - X[-1]
    return D @ D.T


def transverse_block(P, vhat):
    """Schur complement of P onto the directions orthogonal to the collinear direction vhat."""
    Q = np.linalg.qr(np.column_stack([vhat, np.eye(len(vhat))]))[0][:, :len(vhat)]
    if Q[:, 0] @ vhat < 0: Q[:, 0] *= -1
    Pp = Q.T @ P @ Q
    return Pp[1:, 1:] - np.outer(Pp[1:, 0], Pp[0, 1:]) / Pp[0, 0]


eps = 0.02
fig = plt.figure(figsize=(16, 5.2), facecolor=SURF)
# --- (a) 3 bodies: near a collinear triangle the boundary is a flat floor
x0 = np.array([[0.0, 0, 0], [1.0, 0, 0], [2.6, 0, 0]])
p = (x0[:-1] - x0[-1])[:, 0]; vhat = p / np.linalg.norm(p)
def tri_point(d):
    P = gram(x0 + d); s = transverse_block(P, vhat)[0, 0]
    return [(P - np.outer(p, p))[0, 1] / eps, s / eps ** 2]
pts_in = np.array([tri_point(np.c_[rng.normal(size=(3, 2)) * eps, np.zeros(3)]) for _ in range(1500)])
pts_bd = np.array([tri_point(np.c_[rng.normal(size=3) * eps, np.zeros((3, 2))]) for _ in range(300)])
ax = fig.add_subplot(1, 3, 1); ax.set_facecolor(SURF)
ax.scatter(*pts_in.T, s=6, color=BLUE, alpha=0.4, lw=0, label="genuine triangles")
ax.scatter(*pts_bd.T, s=12, color=ORANGE, lw=0, label="flat (collinear) triangles")
ax.axhline(0, color=INK2, lw=1)
ax.set_xlabel("change along the line"); ax.set_ylabel("transverse direction")
ax.set_ylim(-1, np.percentile(pts_in[:, 1], 98)); ax.set_xlim(*np.percentile(pts_in[:, 0], [1, 99]))
ax.legend(frameon=False, loc="upper left", fontsize=10.5)
ax.set_title("3 bodies near a line: a flat floor", fontsize=13, color=INK)
for s_ in ax.spines.values(): s_.set_visible(False)
ax.set_xticks([]); ax.set_yticks([])
# --- (b) 4 bodies: near four collinear particles the boundary is a round cone
x0 = np.array([[0.0, 0, 0], [1.0, 0, 0], [2.5, 0, 0], [4.0, 0, 0]])
p = (x0[:-1] - x0[-1])[:, 0]; vhat = p / np.linalg.norm(p)
def cone_xyz(Cb):
    a, b, c = Cb[0, 0], Cb[0, 1], Cb[1, 1]
    return [(a - c) / 2, b, (a + c) / 2]
pts_in = np.array([cone_xyz(transverse_block(gram(x0 + rng.normal(size=(4, 3)) * eps), vhat))
                   for _ in range(2500)]) / eps ** 2
pts_bd = []
for _ in range(1500):                                  # flat tetrahedra: all four in one random plane
    phi = rng.uniform(0, np.pi); n = np.array([0, np.cos(phi), np.sin(phi)])
    d = np.outer(rng.normal(size=4) * eps, n); d[:, 0] += rng.normal(size=4) * eps
    pts_bd.append(cone_xyz(transverse_block(gram(x0 + d), vhat)))
pts_bd = np.array(pts_bd) / eps ** 2
zmax = np.percentile(pts_bd[:, 2], 90)
keep_in, keep_bd = pts_in[:, 2] < zmax, pts_bd[:, 2] < zmax
ax = fig.add_subplot(1, 3, 2, projection="3d"); ax.set_facecolor(SURF)
th, hh = np.meshgrid(np.linspace(0, 2 * np.pi, 60), np.linspace(0, zmax, 12))
ax.plot_surface(hh * np.cos(th), hh * np.sin(th), hh, color=ORANGE, alpha=0.08, lw=0, shade=False)
ax.scatter(*pts_in[keep_in].T, s=3, color=BLUE, alpha=0.3, depthshade=False)
ax.scatter(*pts_bd[keep_bd].T, s=6, color=ORANGE, alpha=0.95, depthshade=False)
ax.view_init(18, -55); ax.axis("off"); ax.set_box_aspect((1, 1, 1.1), zoom=1.25)
ax.set_title("4 bodies near a line: a round cone", fontsize=13, color=INK, y=1.0)
# --- (c) slice the cone at height 1: a disk, with flat tetrahedra exactly on the rim
ax = fig.add_subplot(1, 3, 3); ax.set_facecolor(SURF); ax.set_aspect("equal")
sl_in = pts_in[:, :2] / pts_in[:, [2]]; sl_bd = pts_bd[:, :2] / pts_bd[:, [2]]
ax.add_patch(Circle((0, 0), 1, fill=False, ec=INK2, lw=1, ls=(0, (4, 3))))
ax.scatter(*sl_in.T, s=6, color=BLUE, alpha=0.4, lw=0, label="genuine 3D tetrahedra")
ax.scatter(*sl_bd.T, s=12, color=ORANGE, lw=0, label="flat (coplanar) tetrahedra")
ax.set_xlim(-1.25, 1.25); ax.set_ylim(-1.25, 1.45); ax.axis("off")
ax.legend(frameon=False, loc="upper center", fontsize=10.5, ncol=2, bbox_to_anchor=(0.5, 1.02))
ax.set_title("a horizontal slice of the cone: a perfect circle", fontsize=13, color=INK)
print("max |radius - 1| of flat tetrahedra on the slice:", np.abs(np.linalg.norm(sl_bd, axis=1) - 1).max())
print("max radius of genuine tetrahedra:", np.linalg.norm(sl_in, axis=1).max())
save(fig, "data_cone.png")


# ===== B2: the doughnut's fingerprint in the four-body surface =========================
def legendre_table(p):
    chi = -np.ones(p, int); chi[0] = 0; chi[(np.arange(1, p) ** 2) % p] = 1
    return chi


def count_X(p):
    chi = legendre_table(p); r = np.arange(p)
    z2, z3, z4 = np.meshgrid(r, r, r, indexing="ij")
    z = [np.zeros_like(z2), z2, z3, z4]; prod = np.ones_like(z2)
    for i in range(4):
        for j in range(i + 1, 4):
            prod *= 1 + chi[(z[j] - z[i]) % p]
    return (int(prod.sum()) - 1) // (p - 1)


def count_conic(p):                                   # three-body analogue
    chi = legendre_table(p); r = np.arange(p)
    z2, z3 = np.meshgrid(r, r, indexing="ij")
    prod = (1 + chi[z2 % p]) * (1 + chi[z3 % p]) * (1 + chi[(z3 - z2) % p])
    return (int(prod.sum()) - 1) // (p - 1)


def a_E(p):
    chi = legendre_table(p); x = np.arange(p)
    return -int(chi[(x ** 3 - x) % p].sum())


primes = [int(q) for q in sp.primerange(5, 200)]
dough = np.array([a_E(q) ** 2 for q in primes])
four = np.array([count_X(q) - (q + 1) ** 2 for q in primes])
three = np.array([count_conic(q) - (q + 1) for q in primes])
assert np.array_equal(dough, four) and not three.any()

fig, axs = plt.subplots(1, 2, figsize=(13, 5.2), facecolor=SURF, gridspec_kw={"width_ratios": [1.25, 1]})
ax = axs[0]; ax.set_facecolor(SURF)
ax.plot(primes, four, "o", ms=8, color=ORANGE, mec=SURF, mew=1.5, label="four-body surface $X$:  $|X(\\mathbb{F}_p)|-(p+1)^2$")
ax.plot(primes, dough, "o", ms=3.5, color=INK, label="doughnut $y^2=x^3-x$:  $a_p^2$")
ax.plot(primes, three, "s", ms=5, color=BLUE, mec=SURF, mew=1, label="three-body conic:  $|C(\\mathbb{F}_p)|-(p+1)$")
ax.set_xlabel("prime $p$"); ax.set_ylabel("deviation from a rational shape")
ax.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
for s_ in ax.spines.values(): s_.set_visible(False)
ax.legend(frameon=False, loc="upper left", fontsize=10.5)
ax.set_title("counting points: the four-body surface carries the doughnut's signal", fontsize=12.5, color=INK)
ax = axs[1]; ax.set_facecolor(SURF)
ax.plot([0, dough.max()], [0, dough.max()], color=INK2, lw=1, ls=(0, (5, 4)))
ax.plot(dough, four, "o", ms=8, color=ORANGE, mec=SURF, mew=1.5)
ax.set_xlabel("from the doughnut:  $a_p^2$")
ax.set_ylabel("from the four-body surface:  $|X(\\mathbb{F}_p)|-(p+1)^2$")
ax.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
for s_ in ax.spines.values(): s_.set_visible(False)
ax.set_title(f"exact agreement for all {len(primes)} primes $5\\leq p<200$", fontsize=12.5, color=INK)
save(fig, "data_doughnut.png")
print("all diagrams saved")
