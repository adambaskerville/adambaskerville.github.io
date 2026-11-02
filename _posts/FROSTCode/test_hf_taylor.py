"""End-to-end check (~50 s): oracle runs reproduce reference HF energies and orbital properties, and their exact-rational
certificates are valid (above the reference) and tight upper bounds."""
from gmpy2 import mpq

import oracle_ns
from certify import certify
from oracle import solve
from properties import properties

# From he_like_series.txt: two independent oracle runs per Z agree to >= 103 s.f.
# He also matches all 66 digits of Hatano & Yamamoto, J. Comput. Chem. Jpn. Int. Ed. 11 (2025).
E_REF = {
    1: mpq("-0.48792973437085174869283837275680759069853910785709626903828113406656"),
    2: mpq("-2.861679995612238878775543740022563429412879776752447053906864697434188362627955"),
    10: mpq("-93.861113519231338154978163881035306429111285436383821944783379879614"),
}


def test_isoelectronic_25_digits():
    """H- (largest wall, no Coulomb tail), He, Ne8+ (smallest wall: continuation must scale with Z)."""
    for Z, ref in E_REF.items():
        data = solve(25, Z=Z, log=lambda *_: None)
        tol = abs(ref) / 10**24
        assert abs(mpq(data["E"]) - ref) < tol, (Z, data["E"])
        E = certify(data)
        assert ref < E < ref + tol, (Z, float(E - ref))

# Hatano & Yamamoto, J. Comput. Chem. Jpn. Int. Ed. 11 (2025): 43 s.f., truncated
E_BE = mpq("-14.5730231683163995824114473598700144538777")


def test_beryllium_20_digits():
    """Two orbitals with separate walls, canonical (lambda = 0) shooting, non-orthogonal certificate."""
    data = oracle_ns.solve(20, Z=4, n=2, log=lambda *_: None)
    tol = abs(E_BE) / 10**19
    assert abs(mpq(data["E"]) - E_BE) < tol, data["E"]
    E = certify(data)
    assert E_BE < E < E_BE + tol, float(E - E_BE)


def test_he_properties_30_digits():
    """Exact cusp; rho(0), <r^-2>, <r> against Hatano & Yamamoto (2025) (32, 35, 63 s.f.); <1/r> = V_H(0)."""
    data = solve(30, log=lambda *_: None)
    p = properties(data, (-2, -1, 1))
    assert p["cusp"] == [-2]
    for got, ref in ((p["rho0"], "3.5959182637622489207069875630017"), (p["r^k"][-2][0], "5.9955036673125245972735587422203366"),
                     (p["r^k"][1][0], "0.927273404731493197418927707229715980644270061322365802985452499"),
                     (p["r^k"][-1][0], data["VH0"])):
        assert abs(mpq(got) / mpq(ref) - 1) < mpq(1, 10**29), (got, ref)


if __name__ == "__main__":
    test_isoelectronic_25_digits()
    test_beryllium_20_digits()
    test_he_properties_30_digits()
    print("ok")
