"""Time-to-error benchmark, closed-shell atoms: hfatom against HelFEM's `atomic` (HF), same machine, all cores.

scan  : HelFEM over primitive basis, nodes, elements and Rmax; every run (energy, wall time) appended to
        helfem_scan.jsonl. Its best setting for each threshold is the fastest run whose true error (against the
        certified limit) is below it -- an oracle choice, generous to HelFEM.
time  : the chosen HelFEM settings and hfatom --tol <threshold> (no tuning) timed again, median of repeats,
        nothing else running. Wall time of the whole process for both.

    python3 bench_atoms.py scan He Ne Ar ...
    python3 bench_atoms.py time He Ne Ar ...      -> bench_atoms.md, bench_atoms.json
"""
import json
import os
import re
import statistics
import subprocess
import sys
import tempfile
import time
from decimal import Decimal

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from refs import REF  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
HELFEM = os.path.expanduser("~/Documents/frost-data/competitors/HelFEM/build4/src/atomic")
HFATOM = os.path.expanduser("~/Documents/GitHub/frost/fast/hfatom")
SCAN = os.path.join(HERE, "helfem_scan.jsonl")
TAUS = (1e-6, 1e-9, 1e-12)
LMAX = {o: max("spdf".index(x[-1]) for x in REF[o][1].split()) for o in REF}


def helfem(name, s, timeout=900):
    """One HelFEM run: (energy or None, wall seconds). Run in a scratch dir (it writes checkpoints)."""
    Z = REF[name][0]
    L = LMAX[name]
    cmd = [HELFEM, f"--Z={Z}", f"--lmax={L}", f"--mmax={L}", f"--nelem={s['nelem']}", f"--nnodes={s['nnodes']}",
           f"--Rmax={s['Rmax']}", f"--primbas={s['primbas']}", "--method=HF", "--verbosity=1"]
    with tempfile.TemporaryDirectory() as d:
        t = time.perf_counter()
        try:
            out = subprocess.run(cmd, capture_output=True, text=True, cwd=d, timeout=timeout).stdout
        except subprocess.TimeoutExpired:
            return None, timeout
        dt = time.perf_counter() - t
    m = re.search(r"Converged to energy\s+(\S+)!", out)
    return (m.group(1) if m else None), dt


def ours(name, tau):
    t = time.perf_counter()
    out = subprocess.run([HFATOM, name, "--tol", repr(tau)], capture_output=True, text=True, check=True).stdout
    dt = time.perf_counter() - t
    return float(re.search(r"E\s+= (\S+)", out).group(1)), dt     # %.17g round-trips: Decimal(float) is exact


def err(name, E):
    return float(Decimal(E) - Decimal(REF[name][2])) if E else float("nan")


def scan(names):
    done = set()
    if os.path.exists(SCAN):
        for line in open(SCAN):
            r = json.loads(line)
            done.add((r["atom"], r["primbas"], r["nnodes"], r["Rmax"], r["nelem"]))
    with open(SCAN, "a") as fh:
        for name in names:
            light = REF[name][0] <= 12   # the full grid for He-Mg; heavier atoms only where the light ones' best runs were
            for primbas in ((4, 3, 5) if light else (4, 3)):
                for nnodes in ((15,) if primbas == 5 else (8, 10, 12, 15, 20) if light else (10, 15, 20)):
                    for Rmax in ((20, 30, 40, 50) if light else (30, 40, 50)):
                        worse = 0
                        best = float("inf")
                        for nelem in range(2, 16):
                            s = dict(primbas=primbas, nnodes=nnodes, Rmax=Rmax, nelem=nelem)
                            if (name, primbas, nnodes, Rmax, nelem) in done:
                                continue
                            E, dt = helfem(name, s)
                            e = err(name, E)
                            fh.write(json.dumps(dict(atom=name, **s, E=E, err=e, wall=dt)) + "\n")
                            fh.flush()
                            print(name, s, E, f"{e:+.1e}", f"{dt:.2f}s", flush=True)
                            a = abs(e) if E else float("inf")
                            # stop once at the floor, or after two refinements without improvement
                            if a < 1e-13:
                                break
                            worse = worse + 1 if a >= best else 0
                            best = min(best, a)
                            if worse >= 2 or dt > 300:
                                break


def best_settings(name, tau):
    runs = [json.loads(line) for line in open(SCAN)]
    ok = [r for r in runs if r["atom"] == name and r["E"] and abs(r["err"]) < tau]
    return min(ok, key=lambda r: r["wall"]) if ok else None


def timing(names, reps=5):
    rows = []
    for name in names:
        for tau in TAUS:
            b = best_settings(name, tau)
            hw, he = None, None
            if b:
                s = {k: b[k] for k in ("primbas", "nnodes", "Rmax", "nelem")}
                ts = []
                for _ in range(max(2, reps // 2)):
                    E, dt = helfem(name, s)
                    ts.append(dt)
                hw, he = statistics.median(ts), err(name, E)
            ts = []
            for _ in range(2 * reps + 1):     # milliseconds each, so more repeats against background noise
                E, dt = ours(name, tau)
                ts.append(dt)
            row = dict(atom=name, tau=tau, helfem=b and {k: b[k] for k in ("primbas", "nnodes", "Rmax", "nelem")},
                       helfem_wall=hw, helfem_err=he, ours_wall=statistics.median(ts), ours_err=err(name, E))
            rows.append(row)
            print(json.dumps(row), flush=True)
    json.dump(rows, open(os.path.join(HERE, "bench_atoms.json"), "w"), indent=1)
    with open(os.path.join(HERE, "bench_atoms.md"), "w") as fh:
        fh.write("| atom | target | HelFEM best (primbas, nnodes, Rmax, nelem) | HelFEM wall | HelFEM error | "
                 "hfatom wall | hfatom error | speedup |\n|---|---|---|---|---|---|---|---|\n")
        for r in rows:
            h = r["helfem"]
            hs = f"{h['primbas']}, {h['nnodes']}, {h['Rmax']}, {h['nelem']}" if h else "not reached"
            hw = f"{r['helfem_wall']:.3f} s" if h else "-"
            he = f"{r['helfem_err']:+.1e}" if h else "-"
            sp = f"{r['helfem_wall'] / r['ours_wall']:.0f}x" if h else "-"
            fh.write(f"| {r['atom']} | {r['tau']:.0e} | {hs} | {hw} | {he} | {1e3 * r['ours_wall']:.1f} ms | "
                     f"{r['ours_err']:+.1e} | {sp} |\n")


if __name__ == "__main__":
    {"scan": scan, "time": timing}[sys.argv[1]](sys.argv[2:])
