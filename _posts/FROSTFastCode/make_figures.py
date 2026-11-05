"""Figures for the FROST Fast post (SVG). Run from this folder: python make_figures.py   (~2 s)

data/ holds the benchmark output (frost-data/fast): bench_atoms.json (HelFEM at its fastest setting for each target
against hfatom, medians of repeats), helfem_scan.jsonl (every HelFEM run of the settings scan), other_pyscf.json
(Gaussian basis sets, PySCF RHF) and other_x2dhf.json (x2dhf on its own atomic test inputs). The FROST MPFR times
are the 12-digit runs of the certified code on the same machine (FROST papers I-II code, final version).
"""
import json
import math

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402

BG, AX, GRID, INK, INK2 = "#f8fafc", "#ffffff", "#e2e8f0", "#0f172a", "#475569"
# one hue per method, fixed order, validated for colour-blind separation; each also has its own marker
C = {"hfatom": "#e11d48", "HelFEM": "#2563eb", "Gaussian (PySCF)": "#d97706", "x2dhf": "#7c3aed",
     "FROST MPFR": "#0d9488"}
MK = {"hfatom": "*", "HelFEM": "o", "Gaussian (PySCF)": "^", "x2dhf": "s", "FROST MPFR": "D"}
OUT = "../../assets/img/frost_fast/"
plt.rcParams.update({"font.family": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans", "sans-serif"],
                     "font.size": 12, "svg.fonttype": "none", "figure.facecolor": BG, "axes.facecolor": AX,
                     "axes.edgecolor": GRID, "axes.grid": True, "grid.color": GRID, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "legend.frameon": False,
                     "savefig.bbox": "tight"})
ATOMS = ["He", "Be", "Ne", "Mg", "Ar", "Ca", "Zn", "Kr", "Sr", "Pd", "Cd", "Xe"]
MPFR = {"Ne": (1.3, 2.5e-17), "Ar": (11.0, 9.8e-17), "Kr": (93.0, 1.5e-16)}   # seconds, error (12-digit target)

bench = json.load(open("data/bench_atoms.json"))
scan = [json.loads(line) for line in open("data/helfem_scan.jsonl")]
# contracted ANO-RCC is built for relativistic Hamiltonians (0.5 Eh off for Ar in non-relativistic HF): left out
gauss = [r for r in json.load(open("data/other_pyscf.json")) if r["basis"] != "ano-rcc"]
# x2dhf: runs that finished on the ground state (two test sets are other states of Kr and Zn)
x2 = [r for r in json.load(open("data/other_x2dhf.json")) if r["E"] and abs(r["err"]) < 1e-3]


def save(fig, name):
    fig.savefig(OUT + name)
    plt.close(fig)


def row(atom, tau):
    return next(r for r in bench if r["atom"] == atom and r["tau"] == tau)


# ------------------------------------------------------------------ hero: wall time for each atom, each method at its best
fig, ax = plt.subplots(figsize=(10, 5.4))
x = range(len(ATOMS))
ours = [row(a, 1e-12)["ours_wall"] for a in ATOMS]
# HelFEM's fastest setting below 1e-9; where no setting got there (Sr), its 1e-6 time, drawn hollow
hel = [row(a, 1e-9)["helfem_wall"] or row(a, 1e-6)["helfem_wall"] for a in ATOMS]
miss = [i for i, a in enumerate(ATOMS) if row(a, 1e-9)["helfem_wall"] is None]
for i, (o, h) in enumerate(zip(ours, hel)):
    ax.plot([i, i], [o, h], color=GRID, lw=2, zorder=1)
ax.plot([i for i in x if i not in miss], [h for i, h in enumerate(hel) if i not in miss], MK["HelFEM"],
        color=C["HelFEM"], ms=9, label="HelFEM, fastest setting below 1e-9 Eh", zorder=3)
for i in miss:
    ax.plot([i], [hel[i]], MK["HelFEM"], mfc=AX, mec=C["HelFEM"], mew=2, ms=9, zorder=3)
    ax.annotate("1e-6 only", (i, hel[i]), xytext=(0, 9), textcoords="offset points", ha="center", fontsize=8,
                color=INK2)
best = {}
for r in x2:   # x2dhf's fastest test input for each atom
    if r["atom"] in ATOMS and (r["atom"] not in best or r["seconds"] < best[r["atom"]]["seconds"]):
        best[r["atom"]] = r
ax.plot([ATOMS.index(a) for a in best], [r["seconds"] for r in best.values()], MK["x2dhf"], color=C["x2dhf"], ms=9,
        label="x2dhf, fastest of its own test inputs", zorder=3)
ax.plot([ATOMS.index(a) for a in MPFR], [v[0] for v in MPFR.values()], MK["FROST MPFR"], color=C["FROST MPFR"], ms=8,
        label="FROST MPFR, 12-digit target", zorder=3)
ax.plot(x, ours, MK["hfatom"], color=C["hfatom"], ms=16, label="hfatom, 1e-12 target: correctly rounded", zorder=4)
for i, (o, h) in enumerate(zip(ours, hel)):
    ax.annotate(f"{h / o:.0f}×", (i, math.sqrt(o * h)), xytext=(7, 0), textcoords="offset points", fontsize=9,
                color=INK2, va="center")
ax.set_yscale("log")
ax.set_xticks(list(x), ATOMS)
ax.set_ylabel("wall time of the whole calculation (s)")
ax.set_title("Hartree–Fock limits of closed-shell atoms, Apple M4 Pro", loc="left", fontsize=13, color=INK)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=2, fontsize=10)
save(fig, "frost_fast_times.svg")

# ------------------------------------------------------------------ error against time: every method, Ne and Kr
fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), sharey=True)
for ax, atom in zip(axes, ("Ne", "Kr")):
    h = [r for r in scan if r["atom"] == atom and r["E"] and abs(r["err"]) > 0]
    ax.plot([r["wall"] for r in h], [abs(r["err"]) for r in h], MK["HelFEM"], color=C["HelFEM"], ms=4, alpha=.35,
            mec="none", label="HelFEM, every setting scanned")
    g = [r for r in gauss if r["atom"] == atom]
    ax.plot([r["seconds"] for r in g], [abs(r["err"]) for r in g], MK["Gaussian (PySCF)"], color=C["Gaussian (PySCF)"],
            ms=9, label="Gaussian basis (PySCF RHF)")
    for r in g:   # label offsets chosen so the neighbouring 5Z and unc-ANO-RCC labels do not collide
        off = {"cc-pvtz": ((-7, -4), "right"), "cc-pv5z": ((6, -11), "left"), "unc-ano-rcc": ((6, 6), "left")}
        if r["basis"] in off:
            ax.annotate(r["basis"], (r["seconds"], abs(r["err"])), xytext=off[r["basis"]][0],
                        textcoords="offset points", ha=off[r["basis"]][1], fontsize=8, color=INK2)
    xr = [r for r in x2 if r["atom"] == atom]
    ax.plot([r["seconds"] for r in xr], [abs(r["err"]) for r in xr], MK["x2dhf"], color=C["x2dhf"], ms=9, label="x2dhf")
    if atom in MPFR:
        ax.plot([MPFR[atom][0]], [MPFR[atom][1]], MK["FROST MPFR"], color=C["FROST MPFR"], ms=8, label="FROST MPFR")
    o = [row(atom, t) for t in (1e-6, 1e-9, 1e-12)]
    ax.plot([r["ours_wall"] for r in o], [max(abs(r["ours_err"]), 1e-16) for r in o], MK["hfatom"], color=C["hfatom"],
            ms=15, label="hfatom, --tol 1e-6, 1e-9, 1e-12")
    half = math.ulp(abs(float(json.load(open("data/refs.json"))[atom]))) / 2
    ax.axhline(half, color=INK2, lw=1, ls=(0, (4, 3)))
    ax.annotate("half an ulp of E: correctly rounded below", (0.98, half), xycoords=("axes fraction", "data"),
                xytext=(0, -11), textcoords="offset points", ha="right", fontsize=8, color=INK2)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_title(atom, loc="left", fontsize=13, color=INK)
    ax.set_xlabel("wall time (s)")
axes[0].set_ylabel("|E − E_HF| (hartree)")
fig.legend(*axes[1].get_legend_handles_labels(), loc="upper center", bbox_to_anchor=(0.5, -0.01), ncol=3,
           fontsize=9)
save(fig, "frost_fast_pareto.svg")
print("wrote", OUT)
