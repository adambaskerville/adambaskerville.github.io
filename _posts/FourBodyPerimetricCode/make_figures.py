"""Figures for the four-body perimetric post (static PNGs)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE = "#2a78d6", "#eb6834"
OUT = "../../assets/img/four_body_perimetric/"
plt.rcParams.update({"font.size": 11, "text.color": INK, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.edgecolor": GRID})

# ---- Figure 1: the rectangle that breaks linear coordinates -------------------------
a, b = 1.0, 2.0
c = np.hypot(a, b); s = (a + b) / 2
fig, axs = plt.subplots(1, 3, figsize=(11, 3.6), facecolor=SURF)
def rect(ax, pts, labels, title):
    ax.set_facecolor(SURF); ax.set_aspect("equal"); ax.axis("off")
    P = np.array(pts)
    for i, j in [(0, 1), (1, 2), (2, 3), (3, 0)]:
        ax.plot(*P[[i, j]].T, color=BLUE, lw=2)
    for i, j in [(0, 2), (1, 3)]:
        ax.plot(*P[[i, j]].T, color=BLUE, lw=1.2, ls=(0, (4, 3)))
    for (x, y), lab in zip(P, labels):
        ax.plot(x, y, "o", ms=9, color=BLUE, mec=SURF, mew=2)
        ax.text(x, y + (0.15 if y > 0.5 else -0.15), lab, ha="center",
                va="bottom" if y > 0.5 else "top", fontsize=12)
    ax.set_title(title, fontsize=11, color=INK)
    ax.set_xlim(-0.5, 2.6); ax.set_ylim(-0.4, 1.6)
rect(axs[0], [(0, 0), (b, 0), (b, a), (0, a)], ["4", "1", "3", "2"],
     "Configuration A\n$r_{14}=b,\\ r_{24}=a,\\ r_{34}=c$")
rect(axs[1], [(0, 0), (b, 0), (b, a), (0, a)], ["4", "2", "3", "1"],
     "Configuration B (same rectangle, relabelled)\n$r_{14}=a,\\ r_{24}=b,\\ r_{34}=c$")
ax = axs[2]; ax.set_facecolor(SURF); ax.axis("off"); ax.set_aspect("equal")
ax.set_xlim(-0.5, 2.6); ax.set_ylim(-0.4, 1.6)
ax.text(1.05, 0.95, "Average of A and B:\nall four sides $=\\frac{a+b}{2}=1.5$\nboth diagonals $=c\\approx2.236$",
        ha="center", va="center", fontsize=11)
ax.text(1.05, 0.05, "but any four points with four equal sides $\\ell$\nhave $d_1^2+d_2^2\\leq 4\\ell^2=9$, while $2c^2=10$",
        ha="center", va="center", fontsize=10, color=ORANGE)
ax.set_title("The average: not a tetrahedron\n(no arrangement has these distances)", fontsize=11)
fig.tight_layout(); fig.savefig(OUT + "rectangle_counterexample.png", dpi=170, facecolor=SURF)

# ---- Figure 2: Pekeris converges, the graded scheme drifts -------------------------
pek_n = [3, 7, 13, 22, 34, 50, 70, 95, 125, 161, 203, 252]
pek_E = [-2.890543433666, -2.903269662751, -2.903596041997, -2.903688986121, -2.903712109027,
         -2.903719562451, -2.903722283080, -2.903723389072, -2.903723878621, -2.903724111149,
         -2.903724228322, -2.903724290411]
gr_n = [216, 512, 1000, 1728, 2744, 4096]
gr_E = [-2.896608160, -2.907032893, -2.911892758, -2.914526149, -2.916153033, -2.917246621]
E0 = -2.903724377034

fig, ax = plt.subplots(figsize=(8, 4.6), facecolor=SURF); ax.set_facecolor(SURF)
ax.axhline(E0, color=INK2, lw=1, ls=(0, (5, 4)))
ax.text(3.2, E0 + 0.0006, "exact helium energy  $-2.9037244$", color=INK2, fontsize=10, va="bottom")
ax.plot(pek_n, pek_E, "-o", color=BLUE, lw=2, ms=6, mec=SURF, mew=1.5, label="Pekeris (perimetric coordinates)")
ax.plot(gr_n, gr_E, "-o", color=ORANGE, lw=2, ms=6, mec=SURF, mew=1.5, label=r"graded scheme in $(r_1,r_2,\cos\theta_{12})$")
ax.text(30, -2.8985, "Pekeris: converges from above", color=INK, fontsize=10, ha="center", va="bottom")
ax.text(4096 * 0.95, gr_E[-1] - 0.0003, "sparse but wrong:\ndrifts straight past", color=INK, fontsize=10, ha="right", va="top")
ax.set_xscale("log"); ax.set_xlabel("number of basis functions / matrix size")
ax.set_ylabel("ground-state energy (hartree)")
ax.set_ylim(-2.920, -2.8895)
ax.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
for sp_ in ax.spines.values(): sp_.set_visible(False)
ax.legend(frameon=False, loc="upper right", labelcolor=INK)
fig.tight_layout(); fig.savefig(OUT + "helium_graded_drift.png", dpi=170, facecolor=SURF)
print("saved")
