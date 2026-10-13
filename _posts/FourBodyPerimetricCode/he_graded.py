"""Pekeris-type sparse recursion for helium in (r1, r2, cos th12) with a Z/2-graded
expansion   psi = exp(-(x1+x2)/2) [ f0(x1,x2,c) + rho * f1(x1,x2,c) ],
x_i = 2 kappa r_i,  rho = 2 kappa r12 = sqrt(s),  s = x1^2 + x2^2 - 2 x1 x2 c.
f0, f1 are expanded in L_i(x1) L_j(x2) P_l(c).  The even equation is multiplied
by x1^2 x2^2 and the odd one by x1^2 x2^2 s; then every term is a polynomial times
x d/dx, x^2 d^2/dx^2, Legendre's operator or (1-c^2) d/dc, all of which act on
the basis with a finite stencil -> sparse  A c = eps B c,  eps = E / (4 kappa^2).

    python3 he_graded.py            # self-check + convergence table
"""
import sys
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla

PAD = 8


def lag(Np):  # columns = image of L_n under x, x d/dx, x^2 d^2/dx^2
    n = np.arange(Np, dtype=float)
    X = sps.diags([-n[1:], 2 * n + 1, -n[1:]], [1, 0, -1], format='csr')
    XD = sps.diags([-n[1:], n], [1, 0], format='csr')
    X2D2 = X @ XD - XD - X @ sps.diags(n)
    return X, XD, X2D2


def leg(Lp):  # columns = image of P_l under c, Legendre op, (1-c^2) d/dc
    l = np.arange(Lp, dtype=float)
    C = sps.diags([l[1:] / (2 * l[1:] + 1), (l[:-1] + 1) / (2 * l[:-1] + 1)], [1, -1], format='csr')
    Lam = sps.diags(-l * (l + 1))
    w = l * (l + 1) / (2 * l + 1)
    W = sps.diags([w[1:], -w[:-1]], [1, -1], format='csr')
    return C, Lam, W


def operators(N, L, Z=2.0, kappa=1.0, lam_ee=1.0):
    """Sparse (A, B) on the box i,j < N, l < L for both components."""
    Np, Lp = N + PAD, L + PAD
    X, XD, X2D2 = lag(Np)
    C, Lam, W = leg(Lp)
    In, Il = sps.identity(Np), sps.identity(Lp)
    K = lambda a, b, c: sps.kron(sps.kron(a, b), c, format='csr')
    x1, x2, c = K(X, In, Il), K(In, X, Il), K(In, In, C)
    d1, d2 = K(XD, In, Il), K(In, XD, Il)          # x_i d/dx_i
    dd1, dd2 = K(X2D2, In, Il), K(In, X2D2, Il)    # x_i^2 d^2/dx_i^2
    lam3, w3 = K(In, In, Lam), K(In, In, W)
    I3 = K(In, In, Il)
    z, lam = Z / (2 * kappa), lam_ee / (2 * kappa)

    x1s, x2s = x1 @ x1, x2 @ x2
    Mx = x1s @ x2s
    s = x1s + x2s - 2 * x1 @ x2 @ c
    kin = -0.5 * (x2s @ (dd1 - x1 @ d1 + 0.25 * x1s + 2 * d1 - x1)
                  + x1s @ (dd2 - x2 @ d2 + 0.25 * x2s + 2 * d2 - x2)
                  + (x1s + x2s) @ lam3)
    nuc = -z * (x1 @ x2s + x1s @ x2)
    extra = -0.5 * (2 * Mx @ d1 - 2 * x1 @ x2s @ x2 @ c @ d1
                    + 2 * Mx @ d2 - 2 * x1s @ x1 @ x2 @ c @ d2
                    - 2 * x1 @ x2 @ (x1s + x2s) @ w3
                    - Mx @ (x1 + x2) @ (I3 - c) + 4 * Mx)
    A00, A11 = kin + nuc, s @ (kin + nuc) + extra
    B00, B11 = Mx, s @ Mx

    ii, jj, ll = np.meshgrid(np.arange(N), np.arange(N), np.arange(L), indexing='ij')
    box = ((ii * Np + jj) * Lp + ll).ravel()
    r = lambda M: M[box][:, box]
    A = sps.bmat([[r(A00), lam * r(Mx)], [lam * r(Mx), r(A11)]], format='csr')
    B = sps.bmat([[r(B00), None], [None, r(B11)]], format='csr')
    full = dict(A00=A00, A01=lam * Mx, A11=A11, B00=B00, B11=B11, Np=Np, Lp=Lp)
    return A, B, full


def self_check():
    """Apply the true Hamiltonian symbolically to random f0, f1 and compare."""
    import sympy as sp
    from scipy.special import eval_genlaguerre, eval_legendre
    rng = np.random.default_rng(0)
    N, L, Z, kappa = 4, 3, 2.0, 1.3
    _, _, F = operators(N, L, Z, kappa)
    Np, Lp = F['Np'], F['Lp']
    c0, c1 = np.zeros(Np * Np * Lp), np.zeros(Np * Np * Lp)
    for vec in (c0, c1):
        for i in range(N):
            for j in range(N):
                for l in range(L):
                    vec[(i * Np + j) * Lp + l] = rng.normal()
    eps = -0.7  # arbitrary: checks (H - E) psi, including the E-terms

    def series(vec, X1, X2, Cc):
        tot = 0.0
        for k in np.nonzero(vec)[0]:
            i, rem = divmod(k, Np * Lp)
            j, l = divmod(rem, Lp)
            tot += vec[k] * eval_genlaguerre(i, 0, X1) * eval_genlaguerre(j, 0, X2) * eval_legendre(l, Cc)
        return tot

    g0 = (F['A00'] - eps * F['B00']) @ c0 + F['A01'] @ c1      # = x1^2 x2^2 * even part
    g1 = F['A01'] @ c0 + (F['A11'] - eps * F['B11']) @ c1      # = x1^2 x2^2 s * odd part

    a1, a2, cc = sp.symbols('x1 x2 c')
    f = []
    for vec in (c0, c1):
        e = 0
        for k in np.nonzero(vec)[0]:
            i, rem = divmod(k, Np * Lp)
            j, l = divmod(rem, Lp)
            e += float(vec[k]) * sp.assoc_laguerre(i, 0, a1) * sp.assoc_laguerre(j, 0, a2) * sp.legendre(l, cc)
        f.append(e)
    rho = sp.sqrt(a1 ** 2 + a2 ** 2 - 2 * a1 * a2 * cc)
    psi = sp.exp(-(a1 + a2) / 2) * (f[0] + rho * f[1])
    Lc = lambda g: (1 - cc ** 2) * sp.diff(g, cc, 2) - 2 * cc * sp.diff(g, cc)
    H = (-sp.Rational(1, 2) * sum(sp.diff(psi, v, 2) + 2 / v * sp.diff(psi, v) for v in (a1, a2))
         - sp.Rational(1, 2) * (1 / a1 ** 2 + 1 / a2 ** 2) * Lc(psi)
         - Z / (2 * kappa) * (1 / a1 + 1 / a2) * psi + 1 / (2 * kappa) / rho * psi - eps * psi)
    Hn = sp.lambdify((a1, a2, cc), H * sp.exp((a1 + a2) / 2), 'numpy')
    for _ in range(20):
        X1, X2, Cc = rng.uniform(0.2, 6), rng.uniform(0.2, 6), rng.uniform(-0.95, 0.95)
        S = X1 ** 2 + X2 ** 2 - 2 * X1 * X2 * Cc
        pred = series(g0, X1, X2, Cc) / (X1 * X2) ** 2 + np.sqrt(S) * series(g1, X1, X2, Cc) / ((X1 * X2) ** 2 * S)
        exact = Hn(X1, X2, Cc)
        assert abs(pred - exact) < 1e-8 * max(1, abs(exact)), (pred, exact)
    # hydrogenic limit: no e-e term, kappa = Z, f0 = 1 is exact with eps = -1/4
    A, B, _ = operators(3, 2, Z=2.0, kappa=2.0, lam_ee=0.0)
    e000 = np.zeros(A.shape[0]); e000[0] = 1
    assert np.allclose(A @ e000, -0.25 * (B @ e000))
    # Galois obstruction: the sparse system cannot tell e-e repulsion from attraction.
    # diag(I,-I) A(+lam) diag(I,-I) = A(-lam), B unchanged -> identical spectra,
    # although +1/r12 (He) and -1/r12 have different spectra.
    import scipy.linalg as sla
    Ap, Bp, _ = operators(6, 3, lam_ee=+1.0)
    Am, Bm, _ = operators(6, 3, lam_ee=-1.0)
    wp = np.sort_complex(sla.eigvals(Ap.toarray(), Bp.toarray()))
    wm = np.sort_complex(sla.eigvals(Am.toarray(), Bm.toarray()))
    ok = np.isfinite(wp) & np.isfinite(wm)
    assert np.allclose(wp[ok], wm[ok], rtol=1e-8, atol=1e-8)
    print("self-check passed (incl. repulsion/attraction spectral identity)")


def ground(N, L, Z=2.0, kappa=1.0, guess=-2.9037):
    A, B, _ = operators(N, L, Z, kappa)
    sig = guess / (4 * kappa ** 2)
    vals = spla.eigs(A.tocsc(), k=4, M=B.tocsc(), sigma=sig, which='LM', return_eigenvectors=False)
    vals = vals[np.abs(vals.imag) < 1e-8].real
    return 4 * kappa ** 2 * vals[np.argmin(np.abs(vals - sig))], A


if __name__ == '__main__':
    self_check()
    if len(sys.argv) > 1 and sys.argv[1] == 'check':
        sys.exit()
    E_ref = -2.903724377034119598
    kappa = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    for N, L in [(6, 3), (8, 4), (10, 5), (12, 6), (14, 7), (16, 8)]:
        E, A = ground(N, L, kappa=kappa)
        nnz_row = A.getnnz(axis=1)
        print(f"N={N:2d} L={L} dim={A.shape[0]:6d}  E={E:.12f}  err={E - E_ref:+.2e}  "
              f"nnz/row: mean {nnz_row.mean():.0f} max {nnz_row.max()}")
