"""Figures for the FROST III post (SVG). Run from this folder: python make_figures.py   (~2 s)

data/amplitudes.json is max over xi of |F_al(xi)|, the size of each Legendre channel of each orbital, for runs B of the
paper (from the frost repository, paper3/data/runs.json). The significant figures are those of Table I of the paper.
"""
import cmath
import json
import math

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402

BG, AX, GRID, INK, INK2 = "#f8fafc", "#ffffff", "#e2e8f0", "#0f172a", "#475569"
BLUE, RED, GREY, GREEN, ORANGE = "#2563eb", "#e11d48", "#94a3b8", "#16a34a", "#eb6834"
OUT = "../../assets/img/frost3_hf/"
plt.rcParams.update({"font.family": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans", "sans-serif"],
                     "font.size": 12, "svg.fonttype": "none", "figure.facecolor": BG, "axes.facecolor": AX,
                     "axes.edgecolor": GRID, "axes.grid": True, "grid.color": GRID, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "legend.frameon": False,
                     "savefig.bbox": "tight"})


def save(fig, name):
    fig.savefig(OUT + name)
    plt.close(fig)


# ------------------------------------------------------------------ records: certified figures against the previous best
# molecule, certified significant figures, previous best (Jensen 1999; Tasi & Csaszar 2007; Beck et al. 2003 and x2dhf)
REC = [("H$_2$", 103, 11), ("HeH$^+$", 101, 9), ("LiH", 26, 12)]
x = range(len(REC))
fig, ax = plt.subplots(figsize=(7, 4.6))
for i, (_, new, old) in zip(x, REC):
    ax.plot([i, i], [old, new], color=GRID, lw=2, zorder=1)
ax.plot(x, [r[2] for r in REC], "o", color=GREY, ms=9, label="previous best", zorder=2)
ax.plot(x, [r[1] for r in REC], "*", color=RED, ms=18, label="FROST III, certified upper bound", zorder=3)
for i, (_, new, old) in zip(x, REC):
    ax.annotate(f"{new}", (i, new), xytext=(0, 11), textcoords="offset points", ha="center", fontsize=10, color=INK2)
    ax.annotate(f"{old}", (i, old), xytext=(0, -20), textcoords="offset points", ha="center", fontsize=10, color=INK2)
ax.set_xticks(list(x), [r[0] for r in REC])
ax.set_xlim(-0.5, len(REC) - 0.5)
ax.set_ylabel("correct significant figures")
ax.set_ylim(0, 150)
ax.legend(loc="upper right")
save(fig, "frost3_records.svg")

# ------------------------------------------------------------------ prolate spheroidal coordinates
# z = (R/2) xi eta along the axis, x = (R/2) sqrt((xi^2 - 1)(1 - eta^2)) across it; nuclei A (eta = -1), B (eta = 1).
R, XW = 2.0, 2.2
fig, ax = plt.subplots(figsize=(8, 4.6))
ts = [i / 400 for i in range(401)]
for xi in (1.05, 1.2, 1.45, 1.8):
    zs = [R / 2 * xi * math.cos(math.pi * t) for t in ts]
    ys = [R / 2 * math.sqrt(xi * xi - 1) * math.sin(math.pi * t) for t in ts]
    ax.plot(zs, ys, color=BLUE, lw=1.2)
    ax.plot(zs, [-y for y in ys], color=BLUE, lw=1.2)
for eta in (-0.8, -0.5, -0.2, 0.2, 0.5, 0.8):
    xis = [1 + (XW - 1) * t for t in ts]
    zs = [R / 2 * xi * eta for xi in xis]
    ys = [R / 2 * math.sqrt((xi * xi - 1) * (1 - eta * eta)) for xi in xis]
    ax.plot(zs, ys, color=GREY, lw=1, ls=(0, (3, 2)))
    ax.plot(zs, [-y for y in ys], color=GREY, lw=1, ls=(0, (3, 2)))
zs = [R / 2 * XW * math.cos(2 * math.pi * t) for t in ts]
ys = [R / 2 * math.sqrt(XW * XW - 1) * math.sin(2 * math.pi * t) for t in ts]
ax.plot(zs, ys, color=RED, lw=2.5)
ax.plot([-R / 2, R / 2], [0, 0], color=INK, lw=3, solid_capstyle="butt")
ax.plot([-R / 2, R / 2], [0, 0], "o", color=INK, ms=9)
ax.annotate("nucleus $A$\n$\\eta=-1$", (-R / 2, 0), xytext=(-3.6, 0), ha="center", va="center", fontsize=11,
            arrowprops=dict(arrowstyle="-", color=GREY))
ax.annotate("nucleus $B$\n$\\eta=1$", (R / 2, 0), xytext=(3.6, 0), ha="center", va="center", fontsize=11,
            arrowprops=dict(arrowstyle="-", color=GREY))
ax.annotate(r"$\xi=1$: the line between the nuclei", (0, 0), xytext=(0, -2.3), ha="center", va="top", fontsize=11,
            arrowprops=dict(arrowstyle="-", color=GREY))
ax.annotate(r"the wall, $\xi=\xi_w$", (0, R / 2 * math.sqrt(XW * XW - 1)), xytext=(0.5, 2.25), fontsize=11,
            color=RED, arrowprops=dict(arrowstyle="-", color=RED))
ax.annotate(r"constant $\xi$: ellipses", (R / 2 * 1.45 * math.cos(1.0), R / 2 * math.sqrt(1.45 ** 2 - 1) * math.sin(1.0)),
            xytext=(1.75, 1.45), fontsize=11, color=BLUE, arrowprops=dict(arrowstyle="-", color=BLUE))
ax.annotate(r"constant $\eta$: hyperbolas", (-0.5 * 1.6, math.sqrt((1.6 ** 2 - 1) * 0.75)), xytext=(-3.1, 1.6),
            fontsize=11, color=INK2, arrowprops=dict(arrowstyle="-", color=GREY))
ax.set_aspect("equal")
ax.set_xlim(-4.2, 4.2)
ax.set_ylim(-2.7, 2.6)
ax.axis("off")
save(fig, "frost3_coords.svg")

# ------------------------------------------------------------------ Legendre channels: measured against predicted decay
# Paper I's helium singularity r* scaled to charge Z, mapped to eta* = -1 + 2 r*/R on the axis; the Bernstein ellipse
# through eta* has rho = |eta* +- sqrt(eta*^2 - 1)| > 1 and the channels fall as rho^-l.
RSTAR_HE = complex(-1.099, 0.577)


def log10_rho(Z, R):
    eta = -1 + 2 * RSTAR_HE * 2 / Z / R
    w = eta + cmath.sqrt(eta * eta - 1)
    return math.log10(max(abs(w), 1 / abs(w)))


amp = json.load(open("data/amplitudes.json"))
fig, ax = plt.subplots(figsize=(10, 5.2))
for mol, label, Z, R, c in (("h2", r"H$_2$", 1, 1.4, BLUE), ("heh", r"HeH$^+$", 2, 1.455, ORANGE),
                            ("lih", "LiH", 3, 3.015, GREEN)):
    a = amp[mol]
    for k, mx in enumerate(a["max"]):
        pts = [(l, math.log10(v)) for l, v in zip(a["ls"], mx) if v]
        name = label + ("" if len(a["max"]) == 1 else [r" $1\sigma$", r" $2\sigma$"][k])
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-" if k == 0 else "s--", color=c, ms=4, lw=1.2,
                mfc=c if k == 0 else AX, label=name)
    l0, v0 = pts[len(pts) // 2]
    s = log10_rho(Z, R)
    ax.plot([l0, pts[-1][0] + 4], [v0 + 1, v0 + 1 - s * (pts[-1][0] + 4 - l0)], color=INK, lw=1.3, ls=(0, (1, 1.5)))
ax.plot([], [], color=INK, lw=1.3, ls=(0, (1, 1.5)), label="predicted")
ax.set_xlabel(r"Legendre channel $l$")
ax.set_ylabel(r"$\log_{10}\max_\xi|F_l(\xi)|$")
ax.legend(loc="lower left")
save(fig, "frost3_channels.svg")
print("wrote", OUT)
