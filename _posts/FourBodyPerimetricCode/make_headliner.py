"""Headliner image: why perimetric coordinates stop at three bodies."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle, Polygon

SURF, INK, INK2 = "#fcfcfb", "#0b0b0b", "#52514e"
BLUE, ORANGE = "#2a78d6", "#eb6834"
UC = ["#2a78d6", "#1baf7a", "#4a3aa7"]          # u1, u2, u3
OUT = "../../assets/img/four_body_perimetric/headliner.png"

fig = plt.figure(figsize=(16, 8.4), facecolor=SURF)
fig.text(0.5, 0.955, "Why perimetric coordinates stop at three bodies", ha="center", va="center",
         fontsize=24, color=INK, weight="bold")

# --- layout boxes [left, bottom, width, height]
rowy = {"top": 0.50, "bot": 0.05}
H = 0.38

def arrow(y, ok, text):
    ax = fig.add_axes([0.355, y, 0.17, H]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    col = BLUE if ok else ORANGE
    ax.add_patch(FancyArrowPatch((0.05, 0.5), (0.95, 0.5), arrowstyle="-|>", mutation_scale=28,
                                 lw=3, color=col, linestyle="-" if ok else (0, (6, 4))))
    ax.text(0.5, 0.64, text, ha="center", va="bottom", fontsize=13, color=INK)
    ax.text(0.5, 0.36, "✓" if ok else "✗", ha="center", va="top", fontsize=34, color=col, weight="bold")

def caption(x, y, title, sub):
    fig.text(x, y, title, fontsize=15, color=INK, weight="bold", va="top")
    fig.text(x, y - 0.035, sub, fontsize=12, color=INK2, va="top")

# ===== TOP ROW: triangle -> flat-walled corner =====================================
ax = fig.add_axes([0.04, rowy["top"], 0.30, H]); ax.axis("off"); ax.set_aspect("equal")
A, B, C = np.array([0.0, 0.0]), np.array([4.0, 0.0]), np.array([1.3, 3.0])
a, b, c = np.linalg.norm(B - C), np.linalg.norm(C - A), np.linalg.norm(A - B)
I = (a * A + b * B + c * C) / (a + b + c)
s = (a + b + c) / 2; rin = np.sqrt((s - a) * (s - b) * (s - c) / s)
ax.add_patch(Polygon([A, B, C], closed=True, fill=False, ec=INK2, lw=2))
ax.add_patch(Circle(I, rin, fill=False, ec=BLUE, lw=2))
verts = [A, B, C]
for k, P in enumerate(verts):                           # two equal tangents from each vertex
    for Q in (verts[(k + 1) % 3], verts[(k + 2) % 3]):
        d = (Q - P) / np.linalg.norm(Q - P); T = P + d * np.dot(I - P, d)
        ax.plot(*np.array([P, T]).T, color=UC[k], lw=6, solid_capstyle="butt")
        ax.plot(*T, "o", ms=8, color=SURF, mec=INK2, mew=1.5, zorder=5)
for k, (P, lab) in enumerate([(A, "$u_1$"), (B, "$u_2$"), (C, "$u_3$")]):
    ax.plot(*P, "o", ms=11, color=INK, mec=SURF, mew=2, zorder=6)
    off = (P - I) / np.linalg.norm(P - I) * 0.45
    ax.text(*(P + off), lab, ha="center", va="center", fontsize=18, color=INK)
ax.set_xlim(-0.8, 4.8); ax.set_ylim(-0.7, 3.6)
caption(0.045, rowy["top"] + 0.025, "3 bodies: a triangle",
        "the two tangents from each corner are equal: that length is $u_i \\geq 0$")

arrow(rowy["top"], True, "flat walls straighten\ninto a box")

ax3 = fig.add_axes([0.53, rowy["top"] - 0.03, 0.30, H], projection="3d"); ax3.set_facecolor(SURF)
g = np.linspace(0, 1, 2); G1, G2 = np.meshgrid(g, g); Z0 = np.zeros_like(G1)
for X, Y, Z in [(G1, G2, Z0), (G1, Z0, G2), (Z0, G1, G2)]:
    ax3.plot_surface(X, Y, Z, color=BLUE, alpha=0.18, edgecolor=BLUE, lw=0.6, shade=False)
for k, v in enumerate(np.eye(3)):
    ax3.quiver(0, 0, 0, *(1.2 * v), color=UC[k], lw=3, arrow_length_ratio=0.1)
for v, lab in zip(np.eye(3) * 1.36, ["$u_1$", "$u_2$", "$u_3$"]):
    ax3.text(*v, lab, fontsize=16, color=INK, ha="center")
ax3.view_init(22, 35); ax3.set_box_aspect((1, 1, 0.9)); ax3.axis("off"); ax3.set_zlim(0, 1.45)
caption(0.84, rowy["top"] + 0.30, "a box corner", "every point is a triangle;\ndistances are linear:\n$r_{12}=u_1+u_2$")

# ===== BOTTOM ROW: tetrahedron -> curved wall with a hidden doughnut ==============
ax4 = fig.add_axes([0.05, rowy["bot"] - 0.01, 0.28, H + 0.04], projection="3d"); ax4.set_facecolor(SURF)
V = np.array([[0, 0, 0], [1.2, 0.1, 0], [0.35, 1.1, 0], [0.55, 0.42, 1.0]])
for i, j in [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 3)]:
    ax4.plot(*V[[i, j]].T, color=INK2, lw=2)
ax4.scatter(*V.T, s=90, color=INK, depthshade=False)
cen = V.mean(axis=0)
for k, v in enumerate(V):
    off = (v - cen) / np.linalg.norm(v - cen) * 0.16
    if k == 2:                                   # vertex 3 projects near the centre: nudge it clear
        off = np.array([-0.2, 0.0, 0.1])
    ax4.text(*(v + off), str(k + 1), fontsize=14, color=INK, ha="center", va="center")
ax4.view_init(20, -75); ax4.set_box_aspect((1, 1, 0.85)); ax4.axis("off")
caption(0.045, rowy["bot"] + 0.025, "4 bodies: a tetrahedron", "six edges, tied together by a curved volume condition")

arrow(rowy["bot"], False, "no box can unroll\na curved wall")

ax5 = fig.add_axes([0.53, rowy["bot"] - 0.03, 0.30, H + 0.06], projection="3d"); ax5.set_facecolor(SURF)
t = np.linspace(0, 2 * np.pi, 60); h = np.linspace(0, 1.3, 20); TT, HH = np.meshgrid(t, h)
ax5.plot_surface(HH * np.cos(TT), HH * np.sin(TT), 1.15 * HH, color=ORANGE, alpha=0.16,
                 edgecolor=ORANGE, lw=0.25, shade=False)          # round "ice-cream cone" wall
u, v = np.meshgrid(np.linspace(0, 2 * np.pi, 50), np.linspace(0, 2 * np.pi, 30))
Rt, rt, zc = 0.42, 0.15, 0.95
ax5.plot_surface((Rt + rt * np.cos(v)) * np.cos(u), (Rt + rt * np.cos(v)) * np.sin(u),
                 zc + rt * np.sin(v), color=ORANGE, alpha=0.75, lw=0, shade=True)   # hidden doughnut
ax5.view_init(18, 30); ax5.set_box_aspect((1, 1, 1)); ax5.axis("off")
ax5.set_xlim(-1.3, 1.3); ax5.set_ylim(-1.3, 1.3); ax5.set_zlim(0, 1.5)
caption(0.84, rowy["bot"] + 0.30, "a round wall", "with a doughnut inside it\n(a K3 surface):\ndistances can't be rational")

fig.savefig(OUT, dpi=150, facecolor=SURF)
print("saved", OUT)
