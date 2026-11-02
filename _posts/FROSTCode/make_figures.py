"""Figures for the FROST post (SVG). Run from this folder: python make_figures.py   (~20 s)

data/ holds summaries of the paper's runs: wall.json (He energies for walls R = 4..100 at 200 digits),
he120.json (the 120-digit run, with the root-test radius of every step) and he500.json (the two 500-digit
runs B and C, the certificate of B, and the pointwise comparison of their orbitals). The shooting figure is
computed here with oracle.py.
"""
import json
import math
from decimal import Decimal as D, getcontext

import gmpy2
import matplotlib
from gmpy2 import mpfr

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

from oracle import shoot, solve, value  # noqa: E402

BG, AX, GRID, INK, INK2 = "#f8fafc", "#ffffff", "#e2e8f0", "#0f172a", "#475569"
BLUE, RED, GREY, GREEN, AMBER = "#2563eb", "#e11d48", "#94a3b8", "#16a34a", "#d97706"
OUT = "../../assets/img/frost_hf/"
plt.rcParams.update({"font.family": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans", "sans-serif"],
                     "font.size": 12, "svg.fonttype": "none", "figure.facecolor": BG, "axes.facecolor": AX,
                     "axes.edgecolor": GRID, "axes.grid": True, "grid.color": GRID, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "legend.frameon": False,
                     "savefig.bbox": "tight"})
getcontext().prec = 700
LN10 = math.log(10)
HE = json.load(open("data/he500.json"))
REF = D(HE["b"]["E"])  # run C (R = 445): wall error ~1e-519, the reference
kappa = math.sqrt(-2 * float(HE["a"]["eps"]))


def sf(e, ref=REF):
    """Correct significant figures of e against ref."""
    return float(-((D(e) - ref).copy_abs() / ref.copy_abs()).log10())


def save(fig, name):
    fig.savefig(OUT + name)
    plt.close(fig)


# ------------------------------------------------------------------ history of the He and Be HF energies
# Values before 2026 as compiled by Hatano & Yamamoto (2025), Table 4 (significant figures)
HIST = {"He": [(1977, 14), (1981, 13), (1995, 10), (2001, 10), (2009, 12), (2018, 13), (2019, 10), (2020, 25),
               (2020, 30), (2025, 66), (2026, 503)],
        "Be": [(1995, 10), (2001, 10), (2009, 12), (2019, 11), (2020, 25), (2022, 30), (2025, 43), (2026, 501)]}
fig, ax = plt.subplots(figsize=(10, 5.2))
for (atom, pts), c, dx in zip(HIST.items(), (BLUE, AMBER), (-0.4, 0.4)):
    ax.plot([y + dx for y, _ in pts[:-1]], [s for _, s in pts[:-1]], "o", color=c, ms=8, label=f"{atom}, basis sets")
    ax.plot([pts[-1][0] + dx], [pts[-1][1]], "*", color=RED if atom == "He" else "#9f1239", ms=18,
            label=f"{atom}, FROST")
ax.annotate("2018: our Laguerre basis\n(King, Baskerville & Cox), 13", (2018, 13), (2003, 40), color=INK2, fontsize=10,
            arrowprops={"arrowstyle": "-", "color": GREY})
ax.annotate("2025: Hatano & Yamamoto, 66", (2025, 66), (2005, 110), color=INK2, fontsize=10,
            arrowprops={"arrowstyle": "-", "color": GREY})
ax.annotate("2026: FROST, He 503 / Be 501", (2026, 503), (2003, 380), color=RED, fontsize=11, fontweight="bold",
            arrowprops={"arrowstyle": "-", "color": RED})
ax.set_yscale("log")
ax.set_ylim(7, 900)
ax.set_xlim(1975, 2029)
ax.set_xlabel("year")
ax.set_ylabel("correct significant figures")
ax.legend(loc="center left")
save(fig, "frost_history.svg")

# ------------------------------------------------------------------ the pipeline
fig, ax = plt.subplots(figsize=(10, 3.4))
ax.set_axis_off()
ax.set_xlim(0, 10)
ax.set_ylim(0, 3.4)
boxes = [("Hartree–Fock\nequations", "2 coupled radial ODEs\nper orbital + pair", GREY),
         ("Taylor shooting", "MPFR floats, any\nheuristic allowed", BLUE),
         ("Piecewise\npolynomial", "dyadic rational\ncoefficients (exact)", AMBER),
         ("Exact energy", "integers only\n(GMP), no rounding", AMBER),
         ("Upper bound", "variational\nprinciple", GREEN)]
for k, (title, sub, c) in enumerate(boxes):
    x = 0.1 + 2 * k
    ax.add_patch(FancyBboxPatch((x, 1.15), 1.6, 1.5, boxstyle="round,pad=0.05,rounding_size=0.12", fc=c + "22", ec=c, lw=2))
    ax.text(x + 0.8, 2.25, title, ha="center", va="center", fontsize=12, fontweight="bold", color=INK)
    ax.text(x + 0.8, 1.55, sub, ha="center", va="center", fontsize=9.5, color=INK2)
    if k:
        ax.annotate("", (x - 0.05, 1.9), (x - 0.35, 1.9), arrowprops={"arrowstyle": "-|>", "color": INK2, "lw": 1.5})
ax.text(3.0, 0.65, "solver (fast, untrusted)", ha="center", color=BLUE, fontsize=11)
ax.text(7.0, 0.65, "certificate (separate program, exact)", ha="center", color=AMBER, fontsize=11)
ax.plot([0.1, 5.7], [0.95, 0.95], color=BLUE, lw=2)
ax.plot([4.1, 9.7], [0.85, 0.85], color=AMBER, lw=2)
save(fig, "frost_pipeline.svg")

# ------------------------------------------------------------------ shooting: a wrong mu lets the growing mode in
d = solve(20, log=lambda *_: None)
gmpy2.get_context().precision = 200
a, mu = mpfr(d["f0"]), mpfr(d["eps"]) - mpfr(d["VH0"])


def curve(m, R=12):
    _, _, el = shoot((a, m), mpfr(R), 2, 60, mpfr("0.15"), energy=True)
    rs, ys = [], []
    for r0, h, f in el:
        for j in range(20):
            x = h * j / 20
            rs.append(float(r0 + x))
            ys.append(abs(float((r0 + x) * value(f, x))))
    return rs, ys


fig, ax = plt.subplots(figsize=(10, 5.2))
for dm, c, ls in ((1e-3, AMBER, "--"), (-1e-3, AMBER, ":"), (1e-8, BLUE, "--"), (-1e-8, BLUE, ":")):
    ax.plot(*curve(mu * (1 + dm)), ls, color=c, lw=1.8, label=rf"$\mu$ off by {dm:+.0e}")
ax.plot(*curve(mu), "-", color=RED, lw=2.5, label=r"$\mu$ correct to 20 digits")
ax.set_yscale("log")
ax.set_ylim(1e-9, 100)
ax.set_xlim(0, 12)
ax.set_xlabel(r"$r$ (bohr)")
ax.set_ylabel(r"$|r\,f(r)|$")
ax.legend(loc="lower left")
save(fig, "frost_shooting.svg")

# ------------------------------------------------------------------ radius of convergence along the orbital
he120 = json.load(open("data/he120.json"))
rstar = complex(-1.099, -0.577)  # nearest singularities r* = -1.099 +- 0.577i (fitted in the paper's make_paper.py)
rad = he120["radius"]
fig, ax = plt.subplots(figsize=(10, 5.2))
ax.plot([r[0] for r in rad], [r[2] for r in rad], "o", color=BLUE, ms=6, label="root-test estimate (sets each step)")
xr = [k / 10 for k in range(0, 10 * int(he120["R"]) + 1)]
ax.plot(xr, [abs(x - rstar) for x in xr], "-", color=RED, lw=2, label=r"distance to the singularities, $|r_0-r_*|$")
ax.axhline(he120["order"] / (math.e * kappa), color=GREY, ls="--", lw=1.5,
           label=r"order-$p$ ceiling of the $e^{-\kappa r}$ tail, $p/(e\kappa)$")
ax.set_xscale("symlog", linthresh=1)
ax.set_yscale("log")
ax.set_xlim(0, 1.1 * he120["R"])
ax.set_ylim(0.8, 800)
ax.set_xlabel(r"expansion point $r_0$ (bohr)")
ax.set_ylabel("radius of convergence (bohr)")
ax.legend(loc="upper left")
save(fig, "frost_radius.svg")

# ------------------------------------------------------------------ the wall law  E_R - E_HF ~ C R^(2 sigma) exp(-2 kappa R)
wall = json.load(open("data/wall.json"))
WALLS = sorted(int(R) for R in wall)
sigma = 1 / kappa  # (Z - 1) / kappa for He
Rs = [R for R in WALLS if R >= 20]
lnC = sum(float((D(wall[str(R)]) - REF).ln()) - 2 * sigma * math.log(R) + 2 * kappa * R for R in Rs) / len(Rs)
model = lambda R: -((lnC + 2 * sigma * math.log(R) - 2 * kappa * R) / LN10 - math.log10(float(-REF)))  # noqa: E731
Ra = HE["a"]["R"]
fig, ax = plt.subplots(figsize=(10, 5.2))
ax.plot(range(4, 450), [model(R) for R in range(4, 450)], "-", color=GREY, lw=2,
        label=r"$C\,R^{2\sigma}e^{-2\kappa R}$, fitted to $20\leq R\leq100$ only")
ax.plot(WALLS, [sf(wall[str(R)]) for R in WALLS], "o", color=BLUE, ms=7, label="solver energy, wall at $R$")
ax.plot([Ra], [sf(HE["cert"]["E_upper"])], "*", color=RED, ms=18, label=f"certified bound, R = {Ra:.0f}")
ax.set_xlabel(r"wall radius $R$ (bohr)")
ax.set_ylabel("correct significant figures")
ax.legend(loc="upper left")
save(fig, "frost_wall.svg")

# ------------------------------------------------------------------ pointwise accuracy of the 500-digit orbital
LOC = HE["local"]
rows = [r for r in LOC["rows"] if r[1] is not None]
fig, ax = plt.subplots(figsize=(10, 5.2))
ax.plot([r[0] for r in rows], [-r[1] for r in rows], "o", color=BLUE, ms=6, label="run B against run C (measured)")
ax.plot([r[0] for r in LOC["rows"] if r[0] >= 10], [-r[3] for r in LOC["rows"] if r[0] >= 10], "-", color=GREY, lw=2,
        label="run B's wall error, predicted by Whittaker functions")
ax.plot([r[0] for r in rows], [-r[2] for r in rows], "s", color=RED, ms=5, label=r"run C against the Coulomb tail $A\,W(2\kappa r)/r$")
ax.set_xlabel(r"$r$ (bohr)")
ax.set_ylabel("correct significant figures of $f(r)$")
ax.set_xlim(0, 1.03 * LOC["RB"])
ax.set_ylim(0, 540)
ax.legend(loc="upper right")
save(fig, "frost_local.svg")
print("wrote", OUT)
