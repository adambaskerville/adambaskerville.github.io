"""The other methods on the same atoms: Gaussian basis sets (PySCF RHF) and x2dhf (2D finite differences, its own
atomic test-set inputs). True errors against the certified limits. Run with nothing else on the machine:

    python3 bench_other.py pyscf Ne Ar Kr Xe      -> other_pyscf.json
    python3 bench_other.py x2dhf he ne kr/set-01  -> other_x2dhf.json (an atom's sets, or one set)
x2dhf is the OpenMP build (x2dhfctl -b -o -B gfortran), rebuilt with -std=legacy: gfortran 15 rejects its legacy
FORMAT strings at run time otherwise. No source change.
PySCF times are of mf.kernel() alone (integrals and SCF), without the Python import: generous to PySCF.
"""
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from decimal import Decimal

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from refs import REF  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
X2 = os.path.expanduser("~/Documents/frost-data/competitors/x2dhf")
BASES = ["cc-pvdz", "cc-pvtz", "cc-pvqz", "cc-pv5z", "cc-pcv5z", "ano-rcc", "x2c-qzvpall", "unc-ano-rcc"]


def pyscf(names):
    from pyscf import gto, scf
    rows = []
    for name in names:
        for basis in BASES:
            try:
                mol = gto.M(atom=f"{name} 0 0 0", basis=basis, verbose=0, cart=False)
            except Exception:   # basis not defined for this element
                continue
            mf = scf.RHF(mol)
            mf.conv_tol = 1e-12
            t = time.perf_counter()
            E = float(mf.kernel())
            dt = time.perf_counter() - t
            row = dict(atom=name, basis=basis, nbf=mol.nao, E=repr(E), err=float(Decimal(repr(E)) - Decimal(REF[name][2])),
                       seconds=dt, converged=bool(mf.converged))
            rows.append(row)
            print(json.dumps(row), flush=True)
    json.dump(rows, open(os.path.join(HERE, "other_pyscf.json"), "w"), indent=1)


def x2dhf(names):
    rows = []
    sym = {r[0]: k for k, r in REF.items()}
    env = dict(os.environ, PATH=f"/opt/homebrew/opt/grep/libexec/gnubin:{X2}/bin:" + os.environ["PATH"],
               X2DHF_DIRECTORY=X2, TERM="dumb")
    for name in names:
        for d in sorted(glob.glob(f"{X2}/test-sets/{name}/run" if "/" in name else f"{X2}/test-sets/{name}/set-*/run")):
            d = os.path.dirname(d)
            inps = [open(f).read() for f in sorted(glob.glob(f"{d}/input*.data"))]
            # sets whose final stage is HF (earlier stages may be LDA starts, which count towards the cost)
            if not inps or not re.search(r"^method hf", inps[-1], re.M) or "nuclei" not in inps[0]:
                continue
            atom = sym.get(int(float(inps[0].split("nuclei")[1].split()[0])))
            if atom is None:
                continue
            # the set's own run script (all its stages, as the authors run it), in a scratch copy
            with tempfile.TemporaryDirectory() as w:
                for f in glob.glob(f"{d}/input*.data") + [f"{d}/run"]:
                    shutil.copy(f, w)
                t = time.perf_counter()
                subprocess.run(["bash", "run"], cwd=w, env=env, capture_output=True, timeout=7200)
                wall = time.perf_counter() - t
                out = open(f"{w}/current.lst").read() if os.path.exists(f"{w}/current.lst") else ""
            Es = re.findall(r"total energy:\s{6,}(\S+)", out)
            E = Es[-1].replace("D", "E") if Es else None
            grid = re.search(r"^grid\s+(.*)$", inps[-1], re.M)
            # x2dhf's own wall clock summed over the stages: leaves out the driver script, generous to x2dhf
            own = sum(float(x) for x in re.findall(r"start-stop:\s+(\S+)", out))
            row = dict(atom=atom, set=os.path.basename(d), stages=len(inps), grid=grid.group(1) if grid else "", E=E,
                       err=float(Decimal(E) - Decimal(REF[atom][2])) if E else None, seconds=own, wall=wall)
            rows.append(row)
            print(json.dumps(row), flush=True)
    json.dump(rows, open(os.path.join(HERE, "other_x2dhf.json"), "w"), indent=1)


if __name__ == "__main__":
    {"pyscf": pyscf, "x2dhf": x2dhf}[sys.argv[1]](sys.argv[2:])
