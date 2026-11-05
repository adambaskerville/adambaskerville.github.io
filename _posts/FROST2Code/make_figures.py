"""Figures for the FROST II post (SVG). Run from this folder: python make_figures.py   (~2 s)

data/tails.json is log10|P_a(r)| of the krypton orbitals of the certified 100-digit run (from the frost repository,
paper2/data/). The significant figures below are those printed in Appendix B of the paper.
"""
import json
import math

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402

BG, AX, GRID, INK, INK2 = "#f8fafc", "#ffffff", "#e2e8f0", "#0f172a", "#475569"
BLUE, RED, GREY = "#2563eb", "#e11d48", "#94a3b8"
OUT = "../../assets/img/frost2_hf/"
plt.rcParams.update({"font.family": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans", "sans-serif"],
                     "font.size": 12, "svg.fonttype": "none", "figure.facecolor": BG, "axes.facecolor": AX,
                     "axes.edgecolor": GRID, "axes.grid": True, "grid.color": GRID, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "legend.frameon": False,
                     "savefig.bbox": "tight"})


def save(fig, name):
    fig.savefig(OUT + name)
    plt.close(fig)


# ------------------------------------------------------------------ records: certified figures against the previous best
# atom, Z, certified significant figures, previous best (Yamamoto & Hatano 2026 up to Ar, Cinal 2020 beyond)
REC = [("Ne", 10, 504, 30), ("Mg", 12, 103, 30), ("Ar", 18, 105, 30), ("Ca", 20, 102, 25), ("Zn", 30, 104, 25),
       ("Kr", 36, 105, 25), ("Sr", 38, 103, 25), ("Pd", 46, 105, 25), ("Cd", 48, 104, 25), ("Xe", 54, 104, 25)]
x = range(len(REC))
fig, ax = plt.subplots(figsize=(10, 5.2))
for i, (_, _, new, old) in zip(x, REC):
    ax.plot([i, i], [old, new], color=GRID, lw=2, zorder=1)
ax.plot(x, [r[3] for r in REC], "o", color=GREY, ms=9, label="previous best, quadruple precision", zorder=2)
ax.plot(x, [r[2] for r in REC], "*", color=RED, ms=18, label="FROST II, certified upper bound", zorder=3)
for i, (_, _, new, _) in zip(x, REC):
    ax.annotate(f"{new}", (i, new), xytext=(0, 11), textcoords="offset points", ha="center", fontsize=10,
                color=INK2)
ax.set_xticks(list(x), [f"{a}\n$Z$ = {z}" for a, z, _, _ in REC])
ax.set_ylabel("correct significant figures")
ax.set_ylim(0, 560)
ax.legend(loc="upper right")
save(fig, "frost2_records.svg")

# ------------------------------------------------------------------ krypton: inner orbitals ride on the 4p tail
# Solid: computed |P_a|; dotted: the orbital's own decay exp(-kappa_a r), kappa_a = sqrt(-2 eps_a), from its maximum.
t = json.load(open("data/tails.json"))
rows, X = t["rows"], 40
fig, ax = plt.subplots(figsize=(10, 5.2))
for name, c in (("1s", "#2a78d6"), ("3d", "#eb6834"), ("4p", "#1baf7a")):
    i = t["orbitals"].index(name)
    pts = [(r[0], r[i + 1]) for r in rows if r[i + 1] is not None and r[0] <= X]
    ax.plot([p[0] for p in pts], [p[1] for p in pts], color=c, lw=2, solid_capstyle="round")
    ax.annotate(rf"${name}$", pts[-1], xytext=(5, 0), textcoords="offset points", va="center", fontsize=12)
    if name != "4p":
        k = math.sqrt(-2 * float(t["eps"][i]))
        rp, yp = max(pts, key=lambda p: p[1])
        r_end = min(X, rp + (yp + 45) * math.log(10) / k)
        ax.plot([rp, r_end], [yp, yp - k * (r_end - rp) / math.log(10)], color=c, lw=1.5, ls=(0, (1, 1.5)))
        ax.annotate(rf"own decay $e^{{-\kappa_{{{name}}}r}}$", (r_end, -45), xytext=(4, 8), textcoords="offset points",
                    fontsize=10, color=INK2)
i1 = t["orbitals"].index("1s")
node = min((r for r in rows if 1.0 < r[0] < 1.6 and r[i1 + 1] is not None), key=lambda r: r[i1 + 1])
ax.annotate("exchange-induced\nnode of the 1s", (node[0], node[i1 + 1]), xytext=(6, -30), fontsize=10, color=INK2,
            arrowprops=dict(arrowstyle="-", color=GREY))
ax.set_xlabel(r"$r$ (bohr)")
ax.set_ylabel(r"$\log_{10}|P_a(r)|$")
ax.set_xlim(0, X + 3)
ax.set_ylim(-45, 2)
save(fig, "frost2_tails.svg")
print("wrote", OUT)
