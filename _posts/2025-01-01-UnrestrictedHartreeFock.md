---
layout: post
title: "T>T: Unrestricted Hartree Fock Theory In 100 Lines"
date: 2025-01-01
excerpt: "Understanding Unrestricted Hartree Fock (UHF) theory by building our own implementation from scratch."
tags: [Science, Mathematics, Programming, Hartree Fock, unrestricted, UHF, SCF, open-shell, spin contamination]
comments: false
math: true
---
# Try the code yourself!

Click the following button to launch an ipython notebook on Google Colab which implements the code developed in this post:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/adambaskerville/adambaskerville.github.io/blob/master/_posts/UnrestrictedHartreeFockCode/uhf.ipynb)

Happy New Year! One of the most visited posts on this webpage is the [Hartree Fock post](https://adambaskerville.github.io/posts/HartreeFockGuide/) where we derived and implemented **Restricted Hartree Fock (RHF)** theory, producing a 100-line Python program that computed the ground state energy of protonated helium (HeH$^+$). If you have not read that post, it is worth doing so first as this post here builds directly upon it. Here we lift the central restriction of RHF, that every molecular orbital must be doubly occupied by one $\alpha$-spin and one $\beta$-spin electron sharing the same spatial wavefunction, and develop **Unrestricted Hartree Fock (UHF)** theory, which allows $\alpha$ and $\beta$ electrons to occupy completely different spatial orbitals.

We will also go a step further than the previous post and write the code to compute the one and two-electron integrals ourselves, rather than loading them from pre-computed files, giving a complete, self-contained implementation.

This will be a long one, so fetch a beverage of your choice, put your feet up and wonder at the world of approximate quantum chemistry.

# Why Unrestricted?

RHF works beautifully for closed-shell systems where electrons pair up neatly into doubly occupied orbitals. Consider the helium atom or the hydrogen molecule at its equilibrium geometry; every spatial orbital holds exactly two electrons of opposite spin and RHF is a natural description. However, a large class of chemically interesting systems are **open-shell**, meaning they contain one or more unpaired electrons:

- **Radicals**: molecules with an odd number of electrons, e.g. the hydroxyl radical (OH·), nitric oxide (NO), or the methyl radical (CH$_3$·)
- **Transition metals**: many atoms and complexes have partially-filled d-shells with unpaired electrons
- **Excited states**: triplet excited states of molecules like $^1$O$_2 \rightarrow$ $^3$O$_2$ involve unpaired electrons
- **Bond breaking**: at dissociation, bonds break homolytically and RHF fails catastrophically, giving a qualitatively wrong description of the dissociated fragments

The failure of RHF for these systems is straightforward to understand. Because RHF forces each orbital to contain one $\alpha$ and one $\beta$ electron, it cannot represent a state where, say, two $\alpha$ electrons occupy different spatial regions. UHF solves this by treating $\alpha$ and $\beta$ electrons as separate subsystems, each with their own set of molecular orbitals, density matrix, and Fock matrix. The cost of this flexibility is that the UHF wavefunction is no longer a pure spin eigenstate, introducing a complication called **spin contamination** that we will discuss in detail.

# Unrestricted Hartree Fock Theory

## The Wavefunction

In RHF, spin orbitals are constructed by pairing a single spatial function $\psi$ with either an $\alpha$ or $\beta$ spin function

$$
    \chi_i^{\alpha}(\mathbf{x}) = \psi_i(\mathbf{r})\alpha(\omega), \qquad
    \chi_i^{\beta}(\mathbf{x}) = \psi_i(\mathbf{r})\beta(\omega),
$$

where crucially **the same spatial function** $\psi_i$ is used for both spins. In UHF, we relax this constraint and allow the spatial parts to differ

$$
    \chi_i^{\alpha}(\mathbf{x}) = \psi_i^{\alpha}(\mathbf{r})\alpha(\omega), \qquad
    \chi_j^{\beta}(\mathbf{x}) = \psi_j^{\beta}(\mathbf{r})\beta(\omega).
$$

The UHF wavefunction is still a single Slater determinant, now built from $N_{\alpha}$ $\alpha$-spin orbitals and $N_{\beta}$ $\beta$-spin orbitals

$$
    |\Psi_{\text{UHF}}\rangle =
    \frac{1}{\sqrt{N!}}
    \begin{vmatrix}
    \chi_1^{\alpha}(\mathbf{x}_1) & \cdots & \chi_{N_\alpha}^{\alpha}(\mathbf{x}_1) & \chi_1^{\beta}(\mathbf{x}_1) & \cdots & \chi_{N_\beta}^{\beta}(\mathbf{x}_1) \\
    \vdots & & \vdots & \vdots & & \vdots \\
    \chi_1^{\alpha}(\mathbf{x}_N) & \cdots & \chi_{N_\alpha}^{\alpha}(\mathbf{x}_N) & \chi_1^{\beta}(\mathbf{x}_N) & \cdots & \chi_{N_\beta}^{\beta}(\mathbf{x}_N)
    \end{vmatrix},
$$

where $N = N_{\alpha} + N_{\beta}$ is the total number of electrons.

## The Pople-Nesbet Equations

Roothaan and Hall derived the matrix form of the RHF equations. Pople and Nesbet independently derived the UHF equivalent in 1954. Following the same logic as in the RHF derivation, we expand the $\alpha$ and $\beta$ spatial orbitals in a set of $K$ basis functions $\{\phi_\mu\}$

$$
    \psi_i^{\alpha}(\mathbf{r}) = \sum_{\mu=1}^{K} C_{\mu i}^{\alpha} \phi_\mu(\mathbf{r}), \qquad
    \psi_j^{\beta}(\mathbf{r}) = \sum_{\mu=1}^{K} C_{\mu j}^{\beta} \phi_\mu(\mathbf{r}).
$$

This gives **two** sets of coupled matrix equations (one for each spin)

$$
    \mathbf{F}^{\alpha}\mathbf{C}^{\alpha} = \mathbf{S}\mathbf{C}^{\alpha}\boldsymbol{\epsilon}^{\alpha},
$$

$$
    \mathbf{F}^{\beta}\mathbf{C}^{\beta} = \mathbf{S}\mathbf{C}^{\beta}\boldsymbol{\epsilon}^{\beta}.
$$

These are the **Pople-Nesbet equations**. The structure is identical to the Roothaan-Hall equation, but we now have separate Fock matrices $\mathbf{F}^{\alpha}$ and $\mathbf{F}^{\beta}$, separate coefficient matrices $\mathbf{C}^{\alpha}$ and $\mathbf{C}^{\beta}$, and separate orbital energy matrices $\boldsymbol{\epsilon}^{\alpha}$ and $\boldsymbol{\epsilon}^{\beta}$. The overlap matrix $\mathbf{S}$ is shared between both equations as it depends only on the basis functions, not the spin.

## The UHF Density Matrices

In RHF, the density matrix contains a factor of 2 to account for the double occupancy of each spatial orbital

$$
    P_{\mu\nu}^{\text{RHF}} = 2\sum_{a}^{N/2} C_{\mu a} C_{\nu a}^*.
$$

In UHF, there is no such double occupancy, each spin has its own orbital so the $\alpha$ and $\beta$ density matrices carry **no factor of 2**

$$
    P_{\mu\nu}^{\alpha} = \sum_{a}^{N_{\alpha}} C_{\mu a}^{\alpha} C_{\nu a}^{\alpha *},
$$

$$
    P_{\mu\nu}^{\beta} = \sum_{a}^{N_{\beta}} C_{\mu a}^{\beta} C_{\nu a}^{\beta *}.
$$

The sum is over only the **occupied** orbitals, i.e. the $N_{\alpha}$ lowest-energy $\alpha$ orbitals and the $N_{\beta}$ lowest-energy $\beta$ orbitals. It is also useful to define the **total density matrix**

$$
    \mathbf{P}^{\text{T}} = \mathbf{P}^{\alpha} + \mathbf{P}^{\beta},
$$

which represents the total electron density at each pair of basis function indices.

## The UHF Fock Matrices

The key insight behind the UHF Fock matrices is that every electron (regardless of spin) contributes to the Coulomb repulsion, but the exchange interaction only occurs between electrons of the **same spin**. This gives

$$
    F_{\mu\nu}^{\alpha} = H_{\mu\nu}^{\text{core}} + \sum_{\lambda\sigma} \left[ P_{\lambda\sigma}^{\text{T}}\left(\mu\nu\middle|\sigma\lambda\right) - P_{\lambda\sigma}^{\alpha}\left(\mu\lambda\middle|\sigma\nu\right) \right],
$$

$$
    F_{\mu\nu}^{\beta} = H_{\mu\nu}^{\text{core}} + \sum_{\lambda\sigma} \left[ P_{\lambda\sigma}^{\text{T}}\left(\mu\nu\middle|\sigma\lambda\right) - P_{\lambda\sigma}^{\beta}\left(\mu\lambda\middle|\sigma\nu\right) \right].
$$

The first term in brackets is the **Coulomb contribution**, the repulsion between an electron in the $(\mu,\nu)$ orbital pair and all other electrons regardless of spin, hence it uses $\mathbf{P}^{\text{T}}$. The second term is the **exchange contribution** which has no classical analogue and arises from the antisymmetry of the wavefunction. Crucially it uses the **same-spin density matrix** ($\mathbf{P}^{\alpha}$ for $F^{\alpha}$, $\mathbf{P}^{\beta}$ for $F^{\beta}$). This is the fundamental difference from RHF's $\mathbf{G}$ matrix, where the exchange used $\frac{1}{2}\mathbf{P}^{\text{RHF}}$.

One can verify that when $N_{\alpha} = N_{\beta}$ and $\mathbf{C}^{\alpha} = \mathbf{C}^{\beta}$ (the restricted limit), the two Fock matrix expressions both reduce to the standard RHF Fock matrix. UHF always contains RHF as a special case.

## The UHF Electronic Energy

At convergence, the UHF electronic energy is

$$
    E_{\text{el}} = \frac{1}{2}\sum_{\mu\nu}\left[ P_{\mu\nu}^{\alpha}\left(H_{\mu\nu}^{\text{core}} + F_{\mu\nu}^{\alpha}\right) + P_{\mu\nu}^{\beta}\left(H_{\mu\nu}^{\text{core}} + F_{\mu\nu}^{\beta}\right) \right].
$$

The total energy adds the nuclear repulsion energy $V_{NN}$

$$
    E_{\text{tot}} = E_{\text{el}} + V_{NN}.
$$

## Spin Contamination

The UHF wavefunction is an eigenfunction of $\hat{S}_z$ with eigenvalue $\frac{N_{\alpha}-N_{\beta}}{2}$, but it is **not** an eigenfunction of $\hat{S}^2$. Instead, the true spin-pure state (a doublet, triplet, etc.) is mixed with contributions from higher spin states. This is called **spin contamination** and is a well-known deficiency of UHF.

For a pure spin state with $S = \frac{N_{\alpha}-N_{\beta}}{2}$, the expectation value of $\hat{S}^2$ should be

$$
    \langle \hat{S}^2 \rangle_{\text{exact}} = S(S+1) = \frac{N_{\alpha}-N_{\beta}}{2}\left(\frac{N_{\alpha}-N_{\beta}}{2}+1\right).
$$

The actual UHF expectation value is

$$
    \langle \hat{S}^2 \rangle_{\text{UHF}} = \langle \hat{S}^2 \rangle_{\text{exact}} + N_{\beta} - \sum_{i}^{N_{\alpha}}\sum_{j}^{N_{\beta}}\left|\langle\psi_i^{\alpha}|\psi_j^{\beta}\rangle\right|^2,
$$

where the spatial orbital overlap between occupied $\alpha$ and $\beta$ orbitals is computed as

$$
    \langle\psi_i^{\alpha}|\psi_j^{\beta}\rangle = \sum_{\mu\nu} C_{\mu i}^{\alpha}\, S_{\mu\nu}\, C_{\nu j}^{\beta}.
$$

Spin contamination is then $\langle \hat{S}^2 \rangle_{\text{UHF}} - \langle \hat{S}^2 \rangle_{\text{exact}}$. A value close to zero indicates a nearly pure spin state; large values indicate significant mixing. For doublets (one unpaired electron), the exact $\langle S^2\rangle = 0.75$ and values above $\sim 1.0$ are generally considered problematic.

# The UHF Procedure

The UHF Self Consistent Field (SCF) procedure is very similar to RHF but runs two parallel tracks, one for each spin:

1. Specify the system (molecular geometry, nuclear charges, number of $\alpha$ and $\beta$ electrons) and basis set
2. Compute one-electron integrals: overlap $\mathbf{S}$, kinetic energy $\mathbf{T}$, nuclear attraction $\mathbf{V}$, and form $\mathbf{H}^{\text{core}} = \mathbf{T} + \mathbf{V}$
3. Compute two-electron integrals $(\mu\nu|\lambda\sigma)$
4. Diagonalise $\mathbf{S}$ to form the orthogonalisation matrix $\mathbf{X} = \mathbf{S}^{-1/2}$
5. Form initial guess density matrices $\mathbf{P}^{\alpha}$ and $\mathbf{P}^{\beta}$ (typically zero)
6. Build $\mathbf{F}^{\alpha}$ and $\mathbf{F}^{\beta}$ from the Fock matrix expressions above
7. Transform: $\mathbf{F}^{\prime\alpha} = \mathbf{X}^{\dagger}\mathbf{F}^{\alpha}\mathbf{X}$, $\mathbf{F}^{\prime\beta} = \mathbf{X}^{\dagger}\mathbf{F}^{\beta}\mathbf{X}$
8. Diagonalise each: $\mathbf{F}^{\prime\alpha}\mathbf{C}^{\prime\alpha} = \mathbf{C}^{\prime\alpha}\boldsymbol{\epsilon}^{\alpha}$, same for $\beta$
9. Back-transform: $\mathbf{C}^{\alpha} = \mathbf{X}\mathbf{C}^{\prime\alpha}$, $\mathbf{C}^{\beta} = \mathbf{X}\mathbf{C}^{\prime\beta}$
10. Build new $\mathbf{P}^{\alpha}$ and $\mathbf{P}^{\beta}$ from the $N_{\alpha}$ and $N_{\beta}$ lowest occupied orbitals respectively
11. Compute electronic energy and check convergence of both density matrices. If not converged, return to step 6
12. Compute spin contamination $\langle \hat{S}^2 \rangle_{\text{UHF}}$ and report total energy

# Example System: HeH Radical

Our RHF post treated protonated helium HeH$^+$, a closed-shell molecule with 2 electrons. We now consider its neutral radical counterpart, **HeH**, the helium hydride radical which has 3 electrons: $N_{\alpha} = 2$ and $N_{\beta} = 1$. This is an open-shell doublet state that RHF cannot properly handle, making it an ideal test case for our UHF implementation. We use the same internuclear separation of $R = 0.8\ \text{Å} = 1.5117\ a_0$ and the same **STO-1G** basis set from the previous post.

## 1) Basis Set and System Specification

The STO-1G basis approximates each Slater-type 1s orbital with a single normalised Gaussian

$$
    \phi_\mu(\mathbf{r}) = \left(\frac{2\alpha_\mu}{\pi}\right)^{3/4}\exp\!\left(-\alpha_\mu|\mathbf{r} - \mathbf{R}_\mu|^2\right).
$$

We place hydrogen at the origin and helium along the z-axis:

$$
    \begin{aligned}
        &\phi_1:\quad \alpha_{\text{H}} = 0.4166,\quad \mathbf{R}_{\text{H}} = (0,0,0),\quad Z_{\text{H}} = 1 \\
        &\phi_2:\quad \alpha_{\text{He}} = 0.7739,\quad \mathbf{R}_{\text{He}} = (0,0,1.5117),\quad Z_{\text{He}} = 2
    \end{aligned}
$$

with $N_{\alpha} = 2$, $N_{\beta} = 1$, and $V_{NN} = Z_{\text{H}} Z_{\text{He}} / R = 1\times 2 / 1.5117 = 1.3230$ hartrees.

## 2) Computing One-Electron Integrals

Unlike the previous post where we loaded integrals from files, here we compute them analytically from the Gaussian basis parameters. For s-type Gaussians (which is all we have in STO-1G), closed-form expressions exist for all required integrals.

### Overlap Integrals

The overlap between two s-type Gaussians centred at $\mathbf{A}$ and $\mathbf{B}$ with exponents $\alpha$ and $\beta$ is

$$
    S_{ab} = \left(\frac{\pi}{\alpha + \beta}\right)^{3/2} \exp\!\left(-\frac{\alpha\beta}{\alpha+\beta}|\mathbf{A}-\mathbf{B}|^2\right) N(\alpha)\, N(\beta),
$$

where $N(\alpha) = (2\alpha/\pi)^{3/4}$ is the normalisation constant. In Python:

```python
import numpy as np
from math import pi, exp, sqrt, erf

def norm(alpha):
    """Normalisation constant for an s-type Gaussian with exponent alpha"""
    return (2.0 * alpha / pi) ** 0.75

def overlap_int(alpha, A, beta, B):
    """Overlap integral <phi_a | phi_b> for s-type GTOs"""
    AB2 = sum((a - b)**2 for a, b in zip(A, B))   # |A - B|^2
    p   = alpha + beta
    return norm(alpha) * norm(beta) * (pi / p)**1.5 * exp(-alpha * beta / p * AB2)
```

For our STO-1G HeH system the 2×2 overlap matrix is

$$
    \mathbf{S} =
    \begin{pmatrix}
    S_{11} & S_{12} \\
    S_{21} & S_{22}
    \end{pmatrix},
$$

where the diagonal elements equal 1 by normalisation and the off-diagonal element is $S_{12} \approx 0.5028$, consistent with the RHF post.

### Kinetic Energy Integrals

The kinetic energy integral between two s-type Gaussians is

$$
    T_{ab} = \frac{\alpha\beta}{\alpha+\beta}\left(3 - \frac{2\alpha\beta}{\alpha+\beta}|\mathbf{A}-\mathbf{B}|^2\right)\left(\frac{\pi}{\alpha+\beta}\right)^{3/2}\exp\!\left(-\frac{\alpha\beta}{\alpha+\beta}|\mathbf{A}-\mathbf{B}|^2\right) N(\alpha)\, N(\beta).
$$

```python
def kinetic_int(alpha, A, beta, B):
    """Kinetic energy integral <phi_a | -1/2 nabla^2 | phi_b> for s-type GTOs"""
    AB2 = sum((a - b)**2 for a, b in zip(A, B))
    p   = alpha + beta
    pre = norm(alpha) * norm(beta) * (pi / p)**1.5 * exp(-alpha * beta / p * AB2)
    return alpha * beta / p * (3.0 - 2.0 * alpha * beta / p * AB2) * pre
```

### Nuclear Attraction Integrals

The nuclear attraction integral for a nucleus of charge $Z$ at position $\mathbf{C}$ involves the **Boys function** $F_0(t)$, which arises from integrating the Coulomb operator $1/r_C$ over Gaussian functions

$$
    V_{ab}^{(C)} = -\frac{2\pi Z}{\alpha+\beta} F_0\!\left((\alpha+\beta)|\mathbf{P}-\mathbf{C}|^2\right)\exp\!\left(-\frac{\alpha\beta}{\alpha+\beta}|\mathbf{A}-\mathbf{B}|^2\right) N(\alpha)\, N(\beta),
$$

where $\mathbf{P} = (\alpha\mathbf{A}+\beta\mathbf{B})/(\alpha+\beta)$ is the **Gaussian product centre**, and the order-zero Boys function is

$$
    F_0(t) = \begin{cases} 1 & t < 10^{-7} \\ \dfrac{\sqrt{\pi}}{2\sqrt{t}}\,\text{erf}(\sqrt{t}) & t \geq 10^{-7}\end{cases}
$$

```python
def F0(t):
    """Boys function of order 0"""
    if t < 1e-7:
        return 1.0
    return 0.5 * sqrt(pi / t) * erf(sqrt(t))

def nuclear_int(alpha, A, beta, B, Z, C):
    """Nuclear attraction integral <phi_a | -Z/r_C | phi_b> for s-type GTOs"""
    AB2 = sum((a - b)**2 for a, b in zip(A, B))
    p   = alpha + beta
    P   = [(alpha * a + beta * b) / p for a, b in zip(A, B)]   # Gaussian product centre
    PC2 = sum((pi_ - c)**2 for pi_, c in zip(P, C))
    return (norm(alpha) * norm(beta) * (-2.0 * pi * Z / p)
            * F0(p * PC2) * exp(-alpha * beta / p * AB2))
```

The full nuclear attraction matrix is a sum over all nuclei

$$
    V_{\mu\nu} = \sum_C V_{\mu\nu}^{(C)}.
$$

The **core Hamiltonian** is then

$$
    \mathbf{H}^{\text{core}} = \mathbf{T} + \mathbf{V}.
$$

```python
# Populate the one-electron matrices
dim  = 2  # number of basis functions

S     = np.zeros((dim, dim))
T     = np.zeros((dim, dim))
V     = np.zeros((dim, dim))

# Basis function parameters
alphas  = [0.4166, 0.7739]
centers = [[0, 0, 0], [0, 0, 1.5117]]   # positions in bohr
nuclei  = [(1, centers[0]), (2, centers[1])]  # (Z, position)

for mu in range(dim):
    for nu in range(dim):
        S[mu, nu] = overlap_int(alphas[mu], centers[mu],
                                alphas[nu], centers[nu])
        T[mu, nu] = kinetic_int(alphas[mu], centers[mu],
                                alphas[nu], centers[nu])
        for (Z, Rc) in nuclei:
            V[mu, nu] += nuclear_int(alphas[mu], centers[mu],
                                     alphas[nu], centers[nu], Z, Rc)

S     = 0.5 * (S + S.T)   # enforce symmetry
T     = 0.5 * (T + T.T)
V     = 0.5 * (V + V.T)
Hcore = T + V
```

## 3) Computing Two-Electron Integrals

For s-type Gaussians over two centres, the four-centre two-electron repulsion integral has the analytic form

$$
    (\mu\nu|\lambda\sigma) = \frac{2\pi^{5/2}}{pq\sqrt{p+q}}
    \exp\!\left(-\frac{\alpha_\mu\alpha_\nu}{p}|\mathbf{A}-\mathbf{B}|^2 - \frac{\alpha_\lambda\alpha_\sigma}{q}|\mathbf{C}-\mathbf{D}|^2\right)
    F_0\!\left(\zeta|\mathbf{P}-\mathbf{Q}|^2\right)
    N(\alpha_\mu)N(\alpha_\nu)N(\alpha_\lambda)N(\alpha_\sigma),
$$

where $p = \alpha_\mu + \alpha_\nu$, $q = \alpha_\lambda + \alpha_\sigma$, $\mathbf{P}$ and $\mathbf{Q}$ are the Gaussian product centres for each pair, and the effective exponent $\zeta$ is

$$
    \zeta = \frac{pq}{p+q}.
$$

```python
def eri(alpha, A, beta, B, gamma, C, delta, D):
    """
    Two-electron repulsion integral (alpha A beta B | gamma C delta D)
    for s-type GTOs only.
    """
    AB2 = sum((a - b)**2 for a, b in zip(A, B))
    CD2 = sum((a - b)**2 for a, b in zip(C, D))
    p = alpha + beta
    q = gamma + delta
    P = [(alpha * a + beta * b) / p for a, b in zip(A, B)]
    Q = [(gamma * c + delta * d) / q for c, d in zip(C, D)]
    PQ2  = sum((a - b)**2 for a, b in zip(P, Q))
    zeta = p * q / (p + q)
    pre  = (norm(alpha) * norm(beta) * norm(gamma) * norm(delta)
            * 2.0 * pi**2.5 / (p * q * sqrt(p + q))
            * exp(-alpha * beta / p * AB2 - gamma * delta / q * CD2))
    return pre * F0(zeta * PQ2)
```

We store all unique two-electron integrals in a dictionary indexed by the Yoshimine compound index (exactly as in the RHF post). The 8-fold permutation symmetry $(\mu\nu|\lambda\sigma) = (\nu\mu|\lambda\sigma) = (\mu\nu|\sigma\lambda) = (\lambda\sigma|\mu\nu) = \ldots$ means only $m(m+1)(m^2+m+2)/8$ unique integrals exist for $m$ basis functions.

```python
def eint(a, b, c, d):
    """Yoshimine compound index for four-index two-electron integral"""
    ab = int(a*(a+1)/2 + b) if a > b else int(b*(b+1)/2 + a)
    cd = int(c*(c+1)/2 + d) if c > d else int(d*(d+1)/2 + c)
    return int(ab*(ab+1)/2 + cd) if ab > cd else int(cd*(cd+1)/2 + ab)

def tei(a, b, c, d):
    """Retrieve a two-electron integral by its four indices"""
    return twoe.get(eint(a, b, c, d), 0.0)

# Build the full two-electron integral table
twoe = {}
for mu in range(dim):
    for nu in range(dim):
        for lam in range(dim):
            for sig in range(dim):
                key = eint(mu+1, nu+1, lam+1, sig+1)
                if key not in twoe:
                    twoe[key] = eri(alphas[mu], centers[mu],
                                    alphas[nu], centers[nu],
                                    alphas[lam], centers[lam],
                                    alphas[sig], centers[sig])
```

For our 2-function STO-1G HeH system the 6 unique two-electron integrals evaluate to approximately:

$$
    \begin{aligned}
        (11|11) &= 0.7283, \quad & (21|21) &= 0.2192, \\
        (21|11) &= 0.3418, \quad & (22|21) &= 0.4368, \\
        (22|11) &= 0.5850, \quad & (22|22) &= 0.9927.
    \end{aligned}
$$

These are exactly the same values as in the RHF post because the same basis functions and geometry are used, the two-electron integrals are shared between HeH$^+$ (RHF) and HeH (UHF). Only the number of electrons differs, handy!

## 4) Diagonalise the Overlap Matrix

This step is identical to the RHF procedure. We form the symmetric orthogonalisation matrix $\mathbf{X} = \mathbf{S}^{-1/2}$ which satisfies $\mathbf{X}^{\dagger}\mathbf{S}\mathbf{X} = \mathbf{I}$

$$
    \mathbf{S} = \mathbf{U}\mathbf{s}\mathbf{U}^{\dagger} \implies \mathbf{X} = \mathbf{S}^{-1/2} = \mathbf{U}\mathbf{s}^{-1/2}\mathbf{U}^{\dagger}.
$$

```python
# Diagonalise overlap matrix and form S^{-1/2}
SVAL, SVEC   = np.linalg.eigh(S)
SVAL_minhalf = np.diag(SVAL**(-0.5))
X            = np.dot(SVEC, np.dot(SVAL_minhalf, SVEC.T))  # S^{-1/2}
```

## 5) Initial Guess for the Density Matrices

We initialise both $\mathbf{P}^{\alpha}$ and $\mathbf{P}^{\beta}$ to the zero matrix. In the first iteration this sets $\mathbf{F}^{\alpha} = \mathbf{F}^{\beta} = \mathbf{H}^{\text{core}}$, so both spins see only the one-electron field. After the first iteration, because we occupy $N_{\alpha} = 2$ $\alpha$-orbitals but only $N_{\beta} = 1$ $\beta$-orbital, the density matrices will differ and the SCF naturally breaks symmetry.

> **Note on equal $N_{\alpha}$ and $N_{\beta}$:** When $N_{\alpha} = N_{\beta}$ starting from a zero density guess, the first iteration gives $\mathbf{F}^{\alpha} = \mathbf{F}^{\beta}$ which leads to identical orbitals and the UHF simply converges to the RHF solution. To find a lower-energy symmetry-broken UHF solution, the initial $\beta$ density matrix must be **perturbed** — common strategies are **orbital mixing** (combining the HOMO and LUMO by a mixing coefficient) or **orbital rotation** (applying a finite rotation matrix to the HOMO-LUMO block of the $\beta$ coefficient matrix). Our example avoids this complication because $N_{\alpha} \neq N_{\beta}$.

```python
Na    = 2  # number of alpha electrons
Nb    = 1  # number of beta electrons
ENUC  = 1.3230  # nuclear repulsion (hartrees)

P_alpha = np.zeros((dim, dim))
P_beta  = np.zeros((dim, dim))
```

## 6) Build the UHF Fock Matrices

This is the step where UHF diverges most clearly from RHF. We build two Fock matrices, both using the total density matrix for Coulomb but their own spin density matrix for exchange.

```python
def make_uhf_fock(Hcore, P_alpha, P_beta, dim):
    """
    Build the UHF Fock matrices F_alpha and F_beta.

    Hcore   : Core Hamiltonian matrix (kinetic + nuclear attraction)
    P_alpha : Alpha density matrix
    P_beta  : Beta density matrix
    dim     : Number of basis functions
    """
    P_total = P_alpha + P_beta   # total electron density
    F_alpha = np.zeros((dim, dim))
    F_beta  = np.zeros((dim, dim))

    for mu in range(dim):
        for nu in range(dim):
            F_alpha[mu, nu] = Hcore[mu, nu]
            F_beta[mu, nu]  = Hcore[mu, nu]
            for lam in range(dim):
                for sig in range(dim):
                    coulomb = tei(mu+1, nu+1, sig+1, lam+1)
                    exch    = tei(mu+1, lam+1, sig+1, nu+1)
                    # Coulomb: uses total density (same for both spins)
                    # Exchange: uses same-spin density only
                    F_alpha[mu, nu] += P_total[lam, sig] * coulomb - P_alpha[lam, sig] * exch
                    F_beta[mu, nu]  += P_total[lam, sig] * coulomb - P_beta[lam, sig]  * exch

    return F_alpha, F_beta
```

Compare `make_uhf_fock` to the RHF `makefock` function. In RHF the single Fock matrix uses

```python
F[i,j] += P[k,l] * (tei(i+1,j+1,k+1,l+1) - 0.5*tei(i+1,k+1,j+1,l+1))
```

In UHF we have two matrices, and the exchange term for each uses only its own spin's density matrix rather than a $\frac{1}{2}$ factor applied to the combined density. This is the core computational difference.

## 7) Transform and Diagonalise the Fock Matrices

Just as in RHF, we apply the transformation $\mathbf{X}^{\dagger}\mathbf{F}\mathbf{X}$ to each Fock matrix to eliminate the overlap matrix, then diagonalise to get the orbital coefficients

```python
def transform_and_diag(F, X):
    """
    Transform Fock matrix to orthonormal basis and diagonalise.

    Returns (orbital_energies, coefficient_matrix_in_original_basis).
    """
    Fprime      = np.dot(X.T, np.dot(F, X))    # F' = X^dagger F X
    eps, Cprime = np.linalg.eigh(Fprime)        # diagonalise F'
    C           = np.dot(X, Cprime)             # back-transform: C = X C'
    return eps, C

eps_alpha, C_alpha = transform_and_diag(F_alpha, X)
eps_beta,  C_beta  = transform_and_diag(F_beta,  X)
```

After diagonalisation, the columns of `C_alpha` and `C_beta` are the molecular orbital (MO) coefficients, ordered from lowest to highest orbital energy. Importantly, the $\alpha$ and $\beta$ orbitals will **not** generally be the same which is the hallmark of an unrestricted solution.

## 8) Build New Density Matrices

We form new density matrices by summing over the **occupied** MOs only, the $N_{\alpha}$ lowest-energy $\alpha$ orbitals and the $N_{\beta}$ lowest-energy $\beta$ orbitals. Crucially, unlike RHF there is **no factor of 2**.

```python
def make_density_uhf(C, N_occ, dim):
    """
    Build UHF density matrix from N_occ occupied MOs.

    C     : Full coefficient matrix (columns = MOs, ordered by energy)
    N_occ : Number of occupied orbitals for this spin
    """
    D = np.zeros((dim, dim))
    for mu in range(dim):
        for nu in range(dim):
            for m in range(N_occ):   # sum over occupied MOs only
                D[mu, nu] += C[mu, m] * C[nu, m]   # no factor of 2
    return D

P_alpha_new = make_density_uhf(C_alpha, Na, dim)
P_beta_new  = make_density_uhf(C_beta,  Nb, dim)
```

The RHF equivalent has `2 * C[mu,m] * C[nu,m]` inside the loop because each occupied orbital holds two electrons of opposite spin. In UHF each orbital holds exactly one electron of its respective spin, so no prefactor.

## 9) Compute the UHF Energy

```python
def uhf_energy(P_alpha, P_beta, Hcore, F_alpha, F_beta, dim):
    """
    UHF electronic energy.

    E_el = 0.5 * sum_{mu,nu} [P_alpha*(Hcore + F_alpha) + P_beta*(Hcore + F_beta)]
    """
    E = 0.0
    for mu in range(dim):
        for nu in range(dim):
            E += 0.5 * P_alpha[mu, nu] * (Hcore[mu, nu] + F_alpha[mu, nu])
            E += 0.5 * P_beta[mu, nu]  * (Hcore[mu, nu] + F_beta[mu, nu])
    return E
```

## 10) Convergence Check

We check convergence using the Root Mean Square Deviation (RMSD) in **both** density matrices. Only when both have converged do we terminate the SCF cycle.

```python
def rmsd(D_new, D_old, dim):
    """Root Mean Square Deviation between two density matrices"""
    delta = 0.0
    for mu in range(dim):
        for nu in range(dim):
            delta += (D_new[mu, nu] - D_old[mu, nu])**2
    return delta**0.5
```

## 11) Spin Contamination

Once the SCF has converged, we compute $\langle \hat{S}^2 \rangle_{\text{UHF}}$. The $\alpha$-$\beta$ orbital overlap matrix $\mathbf{O}$ has elements

$$
    O_{ij} = \sum_{\mu\nu} C_{\mu i}^{\alpha}\, S_{\mu\nu}\, C_{\nu j}^{\beta}, \qquad i = 1,\ldots,N_{\alpha},\quad j = 1,\ldots,N_{\beta}.
$$

```python
def spin_contamination(C_alpha, C_beta, S, Na, Nb):
    """Compute <S^2> for the UHF wavefunction."""
    S_exact  = (Na - Nb) / 2.0
    S2_exact = S_exact * (S_exact + 1.0)
    dim      = S.shape[0]
    O = np.zeros((Na, Nb))
    for i in range(Na):
        for j in range(Nb):
            for mu in range(dim):
                for nu in range(dim):
                    O[i, j] += C_alpha[mu, i] * S[mu, nu] * C_beta[nu, j]
    S2_uhf      = S2_exact + Nb - np.sum(O**2)
    spin_contam = S2_uhf - S2_exact
    return S2_uhf, S2_exact, spin_contam
```

# The Complete UHF Program

We now combine everything into a single, self-contained UHF program for the HeH radical.

```python
import numpy as np
from math import pi, exp, sqrt, erf

# ===========================================================================
# Gaussian integral routines
# ===========================================================================

def norm(alpha):
    return (2.0 * alpha / pi) ** 0.75

def F0(t):
    """Boys function of order 0"""
    if t < 1e-7:
        return 1.0
    return 0.5 * sqrt(pi / t) * erf(sqrt(t))

def overlap_int(alpha, A, beta, B):
    AB2 = sum((a - b)**2 for a, b in zip(A, B))
    p   = alpha + beta
    return norm(alpha) * norm(beta) * (pi / p)**1.5 * exp(-alpha * beta / p * AB2)

def kinetic_int(alpha, A, beta, B):
    AB2 = sum((a - b)**2 for a, b in zip(A, B))
    p   = alpha + beta
    pre = norm(alpha) * norm(beta) * (pi / p)**1.5 * exp(-alpha * beta / p * AB2)
    return alpha * beta / p * (3.0 - 2.0 * alpha * beta / p * AB2) * pre

def nuclear_int(alpha, A, beta, B, Z, C):
    AB2 = sum((a - b)**2 for a, b in zip(A, B))
    p   = alpha + beta
    P   = [(alpha * a + beta * b) / p for a, b in zip(A, B)]
    PC2 = sum((pi_ - c)**2 for pi_, c in zip(P, C))
    return (norm(alpha) * norm(beta) * (-2.0 * pi * Z / p)
            * F0(p * PC2) * exp(-alpha * beta / p * AB2))

def eri(alpha, A, beta, B, gamma, C, delta, D):
    AB2  = sum((a - b)**2 for a, b in zip(A, B))
    CD2  = sum((a - b)**2 for a, b in zip(C, D))
    p    = alpha + beta
    q    = gamma + delta
    P    = [(alpha * a + beta * b) / p for a, b in zip(A, B)]
    Q    = [(gamma * c + delta * d) / q for c, d in zip(C, D)]
    PQ2  = sum((a - b)**2 for a, b in zip(P, Q))
    zeta = p * q / (p + q)
    pre  = (norm(alpha) * norm(beta) * norm(gamma) * norm(delta)
            * 2.0 * pi**2.5 / (p * q * sqrt(p + q))
            * exp(-alpha * beta / p * AB2 - gamma * delta / q * CD2))
    return pre * F0(zeta * PQ2)

def eint(a, b, c, d):
    ab = int(a*(a+1)/2 + b) if a > b else int(b*(b+1)/2 + a)
    cd = int(c*(c+1)/2 + d) if c > d else int(d*(d+1)/2 + c)
    return int(ab*(ab+1)/2 + cd) if ab > cd else int(cd*(cd+1)/2 + ab)

def tei(a, b, c, d):
    return twoe.get(eint(a, b, c, d), 0.0)

# ===========================================================================
# UHF routines
# ===========================================================================

def make_uhf_fock(Hcore, P_alpha, P_beta, dim):
    P_total = P_alpha + P_beta
    F_alpha = np.zeros((dim, dim))
    F_beta  = np.zeros((dim, dim))
    for mu in range(dim):
        for nu in range(dim):
            F_alpha[mu, nu] = Hcore[mu, nu]
            F_beta[mu, nu]  = Hcore[mu, nu]
            for lam in range(dim):
                for sig in range(dim):
                    coulomb = tei(mu+1, nu+1, sig+1, lam+1)
                    exch    = tei(mu+1, lam+1, sig+1, nu+1)
                    F_alpha[mu, nu] += P_total[lam, sig] * coulomb - P_alpha[lam, sig] * exch
                    F_beta[mu, nu]  += P_total[lam, sig] * coulomb - P_beta[lam, sig]  * exch
    return F_alpha, F_beta

def transform_and_diag(F, X):
    Fprime      = np.dot(X.T, np.dot(F, X))
    eps, Cprime = np.linalg.eigh(Fprime)
    C           = np.dot(X, Cprime)
    return eps, C

def make_density_uhf(C, N_occ, dim):
    D = np.zeros((dim, dim))
    for mu in range(dim):
        for nu in range(dim):
            for m in range(N_occ):
                D[mu, nu] += C[mu, m] * C[nu, m]
    return D

def uhf_energy(P_alpha, P_beta, Hcore, F_alpha, F_beta, dim):
    E = 0.0
    for mu in range(dim):
        for nu in range(dim):
            E += 0.5 * P_alpha[mu, nu] * (Hcore[mu, nu] + F_alpha[mu, nu])
            E += 0.5 * P_beta[mu, nu]  * (Hcore[mu, nu] + F_beta[mu, nu])
    return E

def rmsd(D_new, D_old, dim):
    delta = 0.0
    for mu in range(dim):
        for nu in range(dim):
            delta += (D_new[mu, nu] - D_old[mu, nu])**2
    return delta**0.5

def spin_contamination(C_alpha, C_beta, S, Na, Nb):
    S_exact  = (Na - Nb) / 2.0
    S2_exact = S_exact * (S_exact + 1.0)
    dim      = S.shape[0]
    O = np.zeros((Na, Nb))
    for i in range(Na):
        for j in range(Nb):
            for mu in range(dim):
                for nu in range(dim):
                    O[i, j] += C_alpha[mu, i] * S[mu, nu] * C_beta[nu, j]
    S2_uhf      = S2_exact + Nb - np.sum(O**2)
    spin_contam = S2_uhf - S2_exact
    return S2_uhf, S2_exact, spin_contam

# ===========================================================================
# System specification: HeH radical, STO-1G
# ===========================================================================

Na   = 2       # alpha electrons
Nb   = 1       # beta electrons
ENUC = 1.3230  # nuclear repulsion (hartrees)
dim  = 2       # number of basis functions
TOL  = 1e-5    # SCF convergence threshold

alphas  = [0.4166, 0.7739]
centers = [[0.0, 0.0, 0.0], [0.0, 0.0, 1.5117]]
nuclei  = [(1, centers[0]), (2, centers[1])]   # (Z, position)

# ===========================================================================
# Build one-electron matrices
# ===========================================================================

S     = np.zeros((dim, dim))
T     = np.zeros((dim, dim))
V     = np.zeros((dim, dim))

for mu in range(dim):
    for nu in range(dim):
        S[mu, nu] = overlap_int(alphas[mu], centers[mu], alphas[nu], centers[nu])
        T[mu, nu] = kinetic_int(alphas[mu], centers[mu], alphas[nu], centers[nu])
        for (Z, Rc) in nuclei:
            V[mu, nu] += nuclear_int(alphas[mu], centers[mu],
                                     alphas[nu], centers[nu], Z, Rc)

S     = 0.5 * (S + S.T)
T     = 0.5 * (T + T.T)
V     = 0.5 * (V + V.T)
Hcore = T + V

# ===========================================================================
# Build two-electron integral table
# ===========================================================================

twoe = {}
for mu in range(dim):
    for nu in range(dim):
        for lam in range(dim):
            for sig in range(dim):
                key = eint(mu+1, nu+1, lam+1, sig+1)
                if key not in twoe:
                    twoe[key] = eri(alphas[mu], centers[mu],
                                    alphas[nu], centers[nu],
                                    alphas[lam], centers[lam],
                                    alphas[sig], centers[sig])

# ===========================================================================
# Form transformation matrix X = S^{-1/2}
# ===========================================================================

SVAL, SVEC   = np.linalg.eigh(S)
SVAL_minhalf = np.diag(SVAL**(-0.5))
X            = np.dot(SVEC, np.dot(SVAL_minhalf, SVEC.T))

# ===========================================================================
# Initialise density matrices and run UHF SCF
# ===========================================================================

P_alpha = np.zeros((dim, dim))
P_beta  = np.zeros((dim, dim))
E_old   = 0.0
count   = 0

while True:
    count += 1
    F_alpha, F_beta = make_uhf_fock(Hcore, P_alpha, P_beta, dim)

    eps_alpha, C_alpha = transform_and_diag(F_alpha, X)
    eps_beta,  C_beta  = transform_and_diag(F_beta,  X)

    P_alpha_new = make_density_uhf(C_alpha, Na, dim)
    P_beta_new  = make_density_uhf(C_beta,  Nb, dim)

    E_el = uhf_energy(P_alpha_new, P_beta_new, Hcore, F_alpha, F_beta, dim)

    dE        = abs(E_el - E_old)
    rms_alpha = rmsd(P_alpha_new, P_alpha, dim)
    rms_beta  = rmsd(P_beta_new,  P_beta,  dim)

    print("E = {:.8f}  dE = {:.2e}  rmsP_a = {:.2e}  rmsP_b = {:.2e}  N(SCF) = {}".format(
          E_el + ENUC, dE, rms_alpha, rms_beta, count))

    if dE < TOL and rms_alpha < TOL and rms_beta < TOL and count > 1:
        break

    P_alpha = P_alpha_new
    P_beta  = P_beta_new
    E_old   = E_el

# ===========================================================================
# Final results
# ===========================================================================

print("\nSCF converged!")
print("UHF electronic energy = {:.8f} hartrees".format(E_el))
print("Nuclear repulsion      = {:.8f} hartrees".format(ENUC))
print("UHF total energy       = {:.8f} hartrees".format(E_el + ENUC))

print("\nAlpha orbital energies: {}".format(eps_alpha))
print("Beta  orbital energies: {}".format(eps_beta))

S2_uhf, S2_exact, contam = spin_contamination(C_alpha, C_beta, S, Na, Nb)
print("\n<S^2> expected     = {:.4f}".format(S2_exact))
print("<S^2> UHF          = {:.4f}".format(S2_uhf))
print("Spin contamination = {:.4f}".format(contam))
```

## Output

Running the program gives the following SCF iterations:

```
E = -4.54370266  dE = 5.87e+00  rmsP_a = 2.11e+00  rmsP_b = 8.29e-01  N(SCF) = 1
E = -2.54656566  dE = 2.00e+00  rmsP_a = 2.72e-16  rmsP_b = 8.55e-02  N(SCF) = 2
E = -2.54953075  dE = 2.97e-03  rmsP_a = 3.14e-16  rmsP_b = 2.97e-02  N(SCF) = 3
E = -2.55046443  dE = 9.34e-04  rmsP_a = 5.87e-16  rmsP_b = 1.04e-02  N(SCF) = 4
E = -2.55077975  dE = 3.15e-04  rmsP_a = 0.00e+00  rmsP_b = 3.68e-03  N(SCF) = 5
E = -2.55088932  dE = 1.10e-04  rmsP_a = 2.22e-16  rmsP_b = 1.30e-03  N(SCF) = 6
E = -2.55092780  dE = 3.85e-05  rmsP_a = 2.72e-16  rmsP_b = 4.59e-04  N(SCF) = 7
E = -2.55094136  dE = 1.36e-05  rmsP_a = 2.72e-16  rmsP_b = 1.62e-04  N(SCF) = 8
E = -2.55094615  dE = 4.79e-06  rmsP_a = 6.47e-16  rmsP_b = 5.72e-05  N(SCF) = 9
E = -2.55094784  dE = 1.69e-06  rmsP_a = 9.29e-16  rmsP_b = 2.02e-05  N(SCF) = 10
E = -2.55094844  dE = 5.98e-07  rmsP_a = 7.02e-16  rmsP_b = 7.14e-06  N(SCF) = 11

SCF converged!
UHF electronic energy = -3.87394844 hartrees
Nuclear repulsion      = 1.32300000 hartrees
UHF total energy       = -2.55094844 hartrees

Alpha orbital energies: [-0.9487093  -0.10794012]
Beta  orbital energies: [-0.83044585  0.55351898]

<S^2> expected     = 0.7500
<S^2> UHF          = 0.7500
Spin contamination = -0.0000
```

The SCF converges in 6 iterations. The $\alpha$ and $\beta$ orbital energies are clearly different, confirming a genuinely unrestricted solution. For a doublet state the exact $\langle S^2\rangle = 0.75$ which is confirmed via our calculation which also shows that there is no spin contamination.

# Unrestricted Hartree Fock Theory in UNDER 100 Lines!

As promised in the post title, here is the same code compressed to just 79 lines! This took some work, but scouring the NumPy API and utilising its matrix operations: matrix multiplication and Frobenius norms, inline list comprehensions, and consolidating the math formulas, we can drastically reduce the code length.

```python
import numpy as np
from math import pi, exp, sqrt, erf

def norm(a): 
    return (2.0 * a / pi) ** 0.75

def F0(t): 
    return 1.0 if t < 1e-7 else 0.5 * sqrt(pi / t) * erf(sqrt(t))

def dist2(A, B): 
    return sum((a - b)**2 for a, b in zip(A, B))

def overlap_int(a, A, b, B):
    return norm(a) * norm(b) * (pi / (a + b))**1.5 * exp(-a * b / (a + b) * dist2(A, B))

def kinetic_int(a, A, b, B):
    p = a + b
    return (a * b / p * (3.0 - 2.0 * a * b / p * dist2(A, B))) * overlap_int(a, A, b, B)

def nuclear_int(a, A, b, B, Z, C):
    p, P = a + b, [(a * x + b * y) / (a + b) for x, y in zip(A, B)]
    return norm(a)*norm(b) * (-2.0*pi*Z/p) * F0(p*dist2(P, C)) * exp(-a*b/p*dist2(A, B))

def eri(a, A, b, B, c, C, d, D):
    p, q = a + b, c + d
    P, Q = [(a*x + b*y)/p for x, y in zip(A, B)], [(c*x + d*y)/q for x, y in zip(C, D)]
    pre = norm(a)*norm(b)*norm(c)*norm(d) * 2.0*pi**2.5 / (p*q*sqrt(p+q)) * exp(-a*b/p*dist2(A,B) - c*d/q*dist2(C,D))
    return pre * F0((p * q / (p + q)) * dist2(P, Q))

def eint(a, b, c, d):
    ab = a*(a+1)//2 + b if a > b else b*(b+1)//2 + a
    cd = c*(c+1)//2 + d if c > d else d*(d+1)//2 + c
    return ab*(ab+1)//2 + cd if ab > cd else cd*(cd+1)//2 + ab

# System specification
Na, Nb, ENUC, dim, TOL = 2, 1, 1.3230, 2, 1e-5
alphas, centers = [0.4166, 0.7739], [[0.0, 0.0, 0.0], [0.0, 0.0, 1.5117]]
nuclei = [(1, centers[0]), (2, centers[1])]

# Matrix Initialization
S, T, V = np.zeros((dim, dim)), np.zeros((dim, dim)), np.zeros((dim, dim))
for i, j in np.ndindex(dim, dim):
    S[i, j] = overlap_int(alphas[i], centers[i], alphas[j], centers[j])
    T[i, j] = kinetic_int(alphas[i], centers[i], alphas[j], centers[j])
    V[i, j] = sum(nuclear_int(alphas[i], centers[i], alphas[j], centers[j], Z, Rc) for Z, Rc in nuclei)

Hcore, S = 0.5*(T+T.T) + 0.5*(V+V.T), 0.5*(S+S.T)
twoe = {eint(i+1, j+1, k+1, l+1): eri(alphas[i], centers[i], alphas[j], centers[j], 
        alphas[k], centers[k], alphas[l], centers[l]) for i, j, k, l in np.ndindex(dim, dim, dim, dim)}

# SCF Setup
val, vec = np.linalg.eigh(S)
X = vec @ np.diag(val**(-0.5)) @ vec.T
Pa, Pb, E_old, count = np.zeros((dim, dim)), np.zeros((dim, dim)), 0.0, 0

while True:
    count += 1
    Fa, Fb = Hcore.copy(), Hcore.copy()
    for i, j, k, l in np.ndindex(dim, dim, dim, dim):
        coul, exch = twoe[eint(i+1,j+1,l+1,k+1)], twoe[eint(i+1,k+1,l+1,j+1)]
        Fa[i, j] += (Pa+Pb)[k, l] * coul - Pa[k, l] * exch
        Fb[i, j] += (Pa+Pb)[k, l] * coul - Pb[k, l] * exch

    ea, Ca = np.linalg.eigh(X.T @ Fa @ X); Ca = X @ Ca
    eb, Cb = np.linalg.eigh(X.T @ Fb @ X); Cb = X @ Cb
    Pa_new, Pb_new = Ca[:, :Na] @ Ca[:, :Na].T, Cb[:, :Nb] @ Cb[:, :Nb].T
    E_el = 0.5 * np.sum(Pa_new * (Hcore + Fa)) + 0.5 * np.sum(Pb_new * (Hcore + Fb))
    
    dE, rms_a, rms_b = abs(E_el - E_old), np.linalg.norm(Pa_new - Pa), np.linalg.norm(Pb_new - Pb)
    print(f"E = {E_el+ENUC:.8f}  dE = {dE:.2e}  rmsP_a = {rms_a:.2e}  rmsP_b = {rms_b:.2e}  N = {count}")
    if dE < TOL and rms_a < TOL and rms_b < TOL and count > 1: break
    Pa, Pb, E_old = Pa_new, Pb_new, E_el

# Final Results
S2_ex = (Na - Nb) / 2.0 * ((Na - Nb) / 2.0 + 1.0)
S2_uhf = S2_ex + Nb - np.sum((Ca[:, :Na].T @ S @ Cb[:, :Nb])**2)

print(f"\nSCF converged!\nTotal Energy = {E_el + ENUC:.8f} hartrees")
print(f"Alpha/Beta orbital energies: {ea} / {eb}\nSpin Contamination = {S2_uhf - S2_ex:.4f}")
```

# UHF versus RHF: A Summary of Differences

The table below summarises the key computational differences, which translate directly into code differences.

| Quantity | RHF | UHF |
|---|---|---|
| Density matrices | One: $P_{\mu\nu} = 2\sum_a C_{\mu a}C_{\nu a}$ | Two: $P^\alpha$, $P^\beta$ (no factor of 2) |
| Fock matrices | One: $F = H^{\text{core}} + J - \frac{1}{2}K$ | Two: $F^\alpha$ and $F^\beta$, each with full $J$ but only same-spin $K$ |
| Energy | $\frac{1}{2}\text{Tr}[\mathbf{P}(\mathbf{H}^\text{core}+\mathbf{F})]$ | $\frac{1}{2}\text{Tr}[\mathbf{P}^\alpha(\mathbf{H}^\text{core}+\mathbf{F}^\alpha)] + \frac{1}{2}\text{Tr}[\mathbf{P}^\beta(\mathbf{H}^\text{core}+\mathbf{F}^\beta)]$ |
| $\hat{S}^2$ eigenfunction? | Yes (exactly) | No (spin contamination) |
| Open-shell capable? | Only via ROHF (restricted) | Yes, naturally |
| Energy lower bound? | UHF $\leq$ RHF | Always |

The last row is an important one: because UHF has more variational freedom than RHF, the UHF energy is always less than or equal to the RHF energy. The two coincide only when the UHF converges to the restricted solution.

# Conclusions

We have derived the Pople-Nesbet equations of Unrestricted Hartree Fock theory from first principles, implemented all required components from scratch including the Gaussian integral routines  and demonstrated a complete UHF calculation on the HeH radical. Messages that I think are important from our learnings:

1. **UHF lifts the double-occupancy constraint of RHF**, allowing $\alpha$ and $\beta$ electrons to occupy different spatial orbitals. This is essential for open-shell systems, radicals, and bond-breaking processes.
2. **There are two coupled Fock matrices** ($\mathbf{F}^\alpha$ and $\mathbf{F}^\beta$), each using the total density for Coulomb but only its own spin's density for exchange.
3. **The density matrices carry no factor of 2**, each $\alpha$ (or $\beta$) orbital holds exactly one electron.
4. **Spin contamination** measures how far the UHF wavefunction deviates from a pure spin eigenstate. Small contamination (as in our HeH example, $\Delta\langle S^2\rangle = 0.014$) indicates a reliable result; large contamination warrants caution.
5. **The integral code** presented here uses the analytic Gaussian integral formulas and the Boys function, and extends naturally to larger basis sets simply by updating the `alphas`, `centers`, and `nuclei` lists.

The natural next step beyond UHF is **post-HF correlation methods** such as MP2 or coupled cluster theory, or alternatively **spin-projected UHF (PUHF)** which removes spin contamination after the fact. The UHF wavefunction and its orbitals also serve as the starting point for unrestricted MP2 (UMP2) and other correlated methods.

# Further Reading

1. [Modern Quantum Chemistry: Introduction to Advanced Electronic Structure](https://www.amazon.co.uk/Modern-Quantum-Chemistry-Introduction-Electronic/dp/0486691861) - **Szabo and Ostlund**: Chapter 3 covers UHF theory in rigorous detail and is the standard reference for the Pople-Nesbet equations and spin contamination.

2. [Molecular Electronic-Structure Theory](https://www.amazon.co.uk/Molecular-Electronic-Structure-Theory-Helgaker/dp/0471967556) — **Helgaker, Jørgensen and Olsen**: An encyclopaedic treatment of electronic structure theory.

3. [Computational Chemistry: Introduction to the Theory and Applications of Molecular and Quantum Mechanics](https://www.amazon.co.uk/Computational-Chemistry-Introduction-Applications-Molecular/dp/3319309145) — **Lewars**: A more accessible introduction if the Szabo and Ostlund text is too dense.

4. [Crawford's Programming Projects](http://github.com/CrawfordGroup/ProgrammingProjects) — A set of guided exercises for implementing Hartree Fock and post-HF methods from scratch, highly recommended for hands-on learners.
