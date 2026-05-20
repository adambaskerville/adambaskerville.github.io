---
layout: post
title: "T>T: The Integrals at the Heart of Hartree-Fock Theory"
date: 2025-04-18
excerpt: "Every Hartree-Fock calculation rests on a small set of definite integrals. This post derives each one from first principles, explains the mathematics that makes Gaussians uniquely tractable, and also provides some interactive tools for building intuition about why these integrals look the way they do."
tags: [science, mathematics, programming, python, quantum chemistry, hartree-fock, integrals, gaussian, boys function, two-electron integrals, basis set, visualisation]
comments: false
math: true
---

**This post is part of a series on Hartree-Fock theory:**
- [Restricted Hartree-Fock in 100 Lines]({% post_url 2021-04-12-HartreeFockGuide %}): The original post; recommended reading before this one
- [Unrestricted Hartree-Fock in 100 Lines]({% post_url 2025-01-01-UnrestrictedHartreeFock %}): Extends RHF to open-shell systems
- **The Integrals at the Heart of Hartree-Fock Theory**: You are here

---

If you have read the [Hartree-Fock post]({% post_url 2021-04-12-HartreeFockGuide %}) or the [UHF post]({% post_url 2025-01-01-UnrestrictedHartreeFock %}), you will have seen the Roothaan-Hall equations

$$
\mathbf{F}\mathbf{C} = \mathbf{S}\mathbf{C}\boldsymbol{\epsilon}
$$

and the Fock matrix elements

$$
F_{\mu\nu} = H_{\mu\nu}^{\text{core}} + \sum_{\lambda\sigma} P_{\lambda\sigma}\left[(\mu\nu|\sigma\lambda) - \tfrac{1}{2}(\mu\lambda|\sigma\nu)\right].
$$

Both the core Hamiltonian $H_{\mu\nu}^{\text{core}}$ and the two-electron repulsion integrals $(\mu\nu\mid\sigma\lambda)$ are definite integrals over the chosen basis functions. Without them you cannot form the Fock matrix; without the Fock matrix you cannot run the SCF cycle; without the SCF cycle you have no energy. They are not optional, they are essentially the core of the calculation.

The previous posts either loaded these integrals from pre-computed files (RHF) or stated the analytic formulae and moved on (UHF). I wanted to write a post which slows down a bit and derives every integral from scratch, explains the mathematical machinery that makes them tractable, and provides interactive tools so you can build a physical intuition for what these numbers actually represent. I found this the hardest part when learning Hartree Fock theory so hopefully this posta nd these visualisations can bring some light (and life) to the mathematical foundations of Hartree Fock integrals.

**This post covers:**
- Why Gaussian-type orbitals (GTOs) are used instead of the physically correct Slater-type orbitals
- The **Gaussian Product Theorem**: what makes HF integrals analytically tractable
- The three **one-electron integrals**: overlap $S_{\mu\nu}$, kinetic $T_{\mu\nu}$, nuclear attraction $V_{\mu\nu}$
- The **Boys function** $F_n(t)$: the special function that handles the Coulomb operator
- The **two-electron repulsion integral** (ERI) and how it inherits the same machinery
- **8-fold permutation symmetry**, the Yoshimine sort, and why $O(K^4)$ scaling dominates the cost of HF

---

## At a Glance: The Four Integral Types

All four integrals are derived from the same GTO basis functions $\phi_\mu(\mathbf{r}) = N_\mu\exp(-\alpha_\mu\vert\mathbf{r}-\mathbf{A}_\mu\vert^2)$.

| Integral | Symbol | What it encodes | Appears in |
|---|---|---|---|
| **Overlap** | $S_{\mu\nu}$ | Non-orthogonality of basis functions | Right-hand side of Roothaan-Hall ($\mathbf{S}$) |
| **Kinetic energy** | $T_{\mu\nu}$ | Mean kinetic energy of an electron | Core Hamiltonian $\mathbf{H}^{\text{core}}$ |
| **Nuclear attraction** | $V_{\mu\nu}$ | Attraction of an electron to all nuclei | Core Hamiltonian $\mathbf{H}^{\text{core}}$ |
| **Two-electron repulsion** | $(\mu\nu\mid\lambda\sigma)$ | Coulomb repulsion between two charge distributions | Fock matrix $\mathbf{F}$ (Coulomb and exchange) |

The first three are **one-electron integrals**, they depend on where a single electron is. The last is a **two-electron integral** which depends on the simultaneous positions of two electrons and is the most expensive object in a Hartree Fock calculation.

---

## Why Not Use the Physically Correct Orbitals?

The natural choice of basis function for atoms is the **Slater-type orbital** (STO): $\phi(r) \propto e^{-\zeta r}$. STOs are motivated by the exact hydrogen atom solution and decay with the right radial behaviour. The problem is that multi-centre integrals involving STOs, particularly the two-electron repulsion integrals over four different atomic centres, have no closed analytic form and must be evaluated numerically, which is extremely slow.

The practical solution, introduced by Boys in 1950, is to replace STOs with **Gaussian-type orbitals** (GTOs):

$$
\phi_\mu(\mathbf{r}) = N_\mu\, \exp\!\left(-\alpha_\mu\,|\mathbf{r} - \mathbf{A}_\mu|^2\right),
$$

where $\alpha_\mu$ is the **exponent** controlling how diffuse or compact the function is, $\mathbf{A}_\mu$ is the **centre** (typically a nucleus), and $N_\mu = (2\alpha_\mu/\pi)^{3/4}$ is the **normalisation constant** chosen so that $\int \vert\phi_\mu\vert^2\,\mathrm{d}\mathbf{r} = 1$.

Gaussian functions decay as $e^{-\alpha r^2}$ rather than $e^{-\zeta r}$, which is physically wrong, real electron densities have a cusp at the nucleus and decay exponentially, not as a Gaussian, but this physical inaccuracy is the price of a mathematical advantage: **the product of two Gaussians is another Gaussian**. This single property makes all required integrals analytically tractable in closed form,no numerics to slow us down.

In practice, a sum of several GTOs (an STO-$n$G contraction) approximates an STO well enough that the physical error is manageable, whilst the computational gain is enormous.

The interactive visualiser below lets you compare the two directly. The left panel shows the full radial profiles; the right panel zooms into the cusp region at $r=0$. Adjust $\alpha$ (GTO) and $\zeta$ (STO) to match their widths and observe how the GTO smoothly rounds off where the STO has a kink and how the GTO tail collapses too quickly at large $r$.

<iframe
  src="{{ '/assets/interactive/hf_gto_shape.html' | relative_url }}"
  title="GTO vs STO Shape Explorer"
  loading="lazy"
  style="width: 100%; height: 620px; border: 1px solid #2a2a2a; border-radius: 10px; background: #0f0f0f;">
</iframe>

---

## The Gaussian Product Theorem

In my opinion this is the fundamental result underlying every HF integral. It states that the product of two Gaussian functions centred at different points $\mathbf{A}$ and $\mathbf{B}$ is a single Gaussian centred at a weighted average point $\mathbf{P}$:

$$
\phi_\mu(\mathbf{r})\,\phi_\nu(\mathbf{r}) = K_{AB}\,\exp\!\left(-p\,|\mathbf{r} - \mathbf{P}|^2\right),
$$

where:

$$
p = \alpha_\mu + \alpha_\nu, \qquad
\mathbf{P} = \frac{\alpha_\mu\,\mathbf{A} + \alpha_\nu\,\mathbf{B}}{\alpha_\mu + \alpha_\nu}, \qquad
K_{AB} = \exp\!\left(-\frac{\alpha_\mu\,\alpha_\nu}{\alpha_\mu + \alpha_\nu}|\mathbf{A}-\mathbf{B}|^2\right).
$$

The point $\mathbf{P}$ is the **Gaussian product centre**, a weighted average of the two original centres, pulled more strongly toward whichever Gaussian has the larger exponent. The prefactor $K_{AB}$ decreases exponentially with the squared separation $\vert\mathbf{A}-\mathbf{B}\vert^2$: when the two Gaussians are far apart, their product is negligibly small everywhere, and the overlap integral is essentially zero. This is the mathematical basis for **integral screening** in large-scale HF calculations.

The theorem converts every two-centre integral into a one-centre integral. Instead of integrating a product of two Gaussians at $\mathbf{A}$ and $\mathbf{B}$, we integrate a single Gaussian at $\mathbf{P}$ which has an analytic solution.

The interactive tool below attempts to clarify this. Drag the sliders to change the exponents $\alpha_1$, $\alpha_2$ and the separation $d$. The blue and red curves are $\phi_{1}$ and $\phi_{2}$; the purple filled region is their product (the integrand for the overlap integral). The gold dashed line marks the product centre $\mathbf{P}$, and the infobox shows $p$, $\mathbf{P}$, $K_{AB}$, and the resulting overlap integral $S_{\mu\nu}$ in three dimensions.

<iframe
  src="{{ '/assets/interactive/hf_gpt.html' | relative_url }}"
  title="Gaussian Product Theorem Explorer"
  loading="lazy"
  style="width: 100%; height: 620px; border: 1px solid #2a2a2a; border-radius: 10px; background: #0f0f0f;">
</iframe>

Two things we should notice. First: as $d$ increases, $K_{AB}$ drops toward zero and the purple area shrinks. Second: for asymmetric exponents ($\alpha_{1} \gg \alpha_{2}$), $\mathbf{P}$ shifts toward the more compact Gaussian; it is not simply the midpoint unless $\alpha_{1} = \alpha_{2}$.

In Python:

```python
def norm(alpha):
    return (2.0 * alpha / pi) ** 0.75

def gaussian_product(alpha_A, A, alpha_B, B):
    """Returns (p, P, K_AB) for the Gaussian product theorem."""
    p   = alpha_A + alpha_B
    P   = [(alpha_A * a + alpha_B * b) / p for a, b in zip(A, B)]
    AB2 = sum((a - b)**2 for a, b in zip(A, B))
    K   = exp(-alpha_A * alpha_B / p * AB2)
    return p, P, K
```

---

## One-Electron Integrals

### Overlap Integrals

The overlap matrix $\mathbf{S}$ encodes the non-orthogonality of the basis functions. Its elements are:

$$
S_{\mu\nu} = \int \phi_\mu(\mathbf{r})\,\phi_\nu(\mathbf{r})\,\mathrm{d}\mathbf{r}.
$$

Applying the Gaussian product theorem, the integrand becomes a single Gaussian at $\mathbf{P}$, and the remaining integral is the standard three-dimensional Gaussian integral $\int e^{-p\vert\mathbf{r}-\mathbf{P}\vert^2}\mathrm{d}\mathbf{r} = (\pi/p)^{3/2}$:

$$
S_{\mu\nu} = N_\mu\,N_\nu\,K_{AB}\,\left(\frac{\pi}{p}\right)^{3/2}.
$$

The diagonal elements satisfy $S_{\mu\mu} = 1$ by construction (the normalisation condition). The off-diagonal elements are positive, less than 1, and decrease exponentially with separation, which is why basis functions on distant atoms do not need to be treated as overlapping. In code we an implement this as

```python
def overlap_int(alpha, A, beta, B):
    AB2 = sum((a - b)**2 for a, b in zip(A, B))
    p   = alpha + beta
    return norm(alpha) * norm(beta) * (pi / p)**1.5 * exp(-alpha * beta / p * AB2)
```

For equal exponents the formula simplifies to $S_{\mu\nu} = e^{-\alpha R^2/2}$ (a pure Gaussian decay in the separation $R$). The visualiser below shows this directly. Drag the $R$ slider to watch the two Gaussians separate: the shaded product (left panel) shrinks, and the dot traces along the $S_{\mu\nu}$ vs $R$ decay curve (right panel).

<iframe
  src="{{ '/assets/interactive/hf_overlap.html' | relative_url }}"
  title="Overlap Integral Explorer"
  loading="lazy"
  style="width: 100%; height: 620px; border: 1px solid #2a2a2a; border-radius: 10px; background: #0f0f0f;">
</iframe>

### Kinetic Energy Integrals

The kinetic energy matrix element is:

$$
T_{\mu\nu} = \int \phi_\mu(\mathbf{r})\left(-\tfrac{1}{2}\nabla^2\right)\phi_\nu(\mathbf{r})\,\mathrm{d}\mathbf{r}.
$$

Applying $\nabla^2$ to a Gaussian produces a combination of Gaussians of the same and higher order. For s-type GTOs, the result simplifies to:

$$
T_{\mu\nu} = \frac{\alpha\beta}{p}\left(3 - \frac{2\alpha\beta}{p}|\mathbf{A}-\mathbf{B}|^2\right)\left(\frac{\pi}{p}\right)^{3/2} K_{AB}\,N_\mu\,N_\nu.
$$

This has the structure of the overlap integral multiplied by a factor involving the exponents and separation. When the two Gaussians are centred on the same atom ($\vert\mathbf{A}-\mathbf{B}\vert^2 = 0$), the factor reduces to $3\alpha\beta/p$, which is a positive kinetic energy contribution. As the separation increases, the sign of the bracket can change, reflecting that spatially separated Gaussians have a more complicated kinetic contribution from the oscillatory nature of any orbital built from them

```python
def kinetic_int(alpha, A, beta, B):
    AB2 = sum((a - b)**2 for a, b in zip(A, B))
    p   = alpha + beta
    pre = norm(alpha) * norm(beta) * (pi / p)**1.5 * exp(-alpha * beta / p * AB2)
    return alpha * beta / p * (3.0 - 2.0 * alpha * beta / p * AB2) * pre
```

### Nuclear Attraction Integrals and the Boys Function

The nuclear attraction integral for a nucleus of charge $Z$ at position $\mathbf{C}$ is:

$$
V_{\mu\nu}^{(C)} = \int \phi_\mu(\mathbf{r})\left(-\frac{Z}{|\mathbf{r}-\mathbf{C}|}\right)\phi_\nu(\mathbf{r})\,\mathrm{d}\mathbf{r}.
$$

The Gaussian product theorem reduces the two-Gaussian integrand to a single Gaussian at $\mathbf{P}$. The problem is now to integrate a Gaussian against $1/\vert\mathbf{r}-\mathbf{C}\vert$. This is where we run out of Gaussian tricks: the Coulomb operator is not a Gaussian and the integral has no elementary antiderivative.

The strategy is to **convert $1/r$ into a Gaussian** by using its integral representation:

$$
\frac{1}{|\mathbf{r}-\mathbf{C}|} = \frac{2}{\sqrt{\pi}}\int_0^\infty \exp\!\left(-|\mathbf{r}-\mathbf{C}|^2\,t^2\right)\mathrm{d}t.
$$

This identity (verifiable by substituting $s = \vert\mathbf{r}-\mathbf{C}\vert t$ and using $\int_0^\infty e^{-s^2}\mathrm{d}s = \sqrt{\pi}/2$) writes the tricky $1/r$ as a Gaussian in $t$. Now the integrand in $\mathbf{r}$ is a product of two Gaussians, one from the GPT, one from the identity and the $\mathbf{r}$ integral is now tractable. The remaining integral over $t$ cannot be simplified further; it defines the **Boys function**. Substituting and evaluating the $\mathbf{r}$ integral gives:

$$
V_{\mu\nu}^{(C)} = -\frac{2\pi Z}{p}\,K_{AB}\,F_0\!\left(p|\mathbf{P}-\mathbf{C}|^2\right)\,N_\mu\,N_\nu,
$$

where $F_0$ is the **order-zero Boys function**:

$$
F_0(t) = \int_0^1 e^{-tu^2}\,\mathrm{d}u = \begin{cases} 1 & t \to 0 \\ \dfrac{\sqrt{\pi}}{2\sqrt{t}}\,\operatorname{erf}(\sqrt{t}) & t > 0. \end{cases}
$$

The Boys function appears because the Coulomb operator prevents a fully analytic resolution. It leaves behind a one-dimensional integral that happens to be tabulated. The argument $t = p\vert\mathbf{P}-\mathbf{C}\vert^2$ measures how far the Gaussian product centre is from the nucleus, scaled by the combined exponent. When the product centre sits exactly on the nucleus ($t = 0$), $F_0(0) = 1$ and the integral takes its maximum value. As the centre moves away ($t$ increases), $F_0(t)$ decays smoothly to zero.

For higher angular momentum basis functions (p, d, f shells), higher-order Boys functions $F_n(t)$ arise. They satisfy the recurrence:

$$
F_{n+1}(t) = \frac{(2n+1)\,F_n(t) - e^{-t}}{2t}, \qquad F_n(0) = \frac{1}{2n+1}.
$$

The tool below visualises the Boys function. Drag the nucleus slider to move the proton along the axis. The left panel shows the Gaussian charge density (blue), the integrand density/distance (green, clipped), and the nucleus position (gold). The right panel plots $V_{\mu\mu}$ vs nucleus position $C$ alongside the classical $-1/C$ Coulomb reference and the Boys function $F_0$ — which is exactly the shape of $V_{\mu\mu}(C)$ up to a prefactor.

<iframe
  src="{{ '/assets/interactive/hf_nuclear.html' | relative_url }}"
  title="Nuclear Attraction and Boys Function Explorer"
  loading="lazy"
  style="width: 100%; height: 620px; border: 1px solid #2a2a2a; border-radius: 10px; background: #0f0f0f;">
</iframe>

For higher angular momentum basis functions (p, d, f shells which we will not think about here), higher-order Boys functions $F_n(t)$ arise which satisfy the recurrence

$$
F_{n+1}(t) = \frac{(2n+1)\,F_n(t) - e^{-t}}{2t}, \qquad F_n(0) = \frac{1}{2n+1}.
$$

The visualiser below shows $F_n(t)$ for orders $n = 0$ to $4$. The left panel displays the integrand $f(u) = u^{2n}e^{-tu^2}$ on $[0,1]$, the shaded area is $F_n(t)$ by definition. Use the order buttons to switch $n$ and the slider to change $t$; watch the area shrink as $t$ increases.

<iframe
  src="{{ '/assets/interactive/hf_boys.html' | relative_url }}"
  title="Boys Function Explorer"
  loading="lazy"
  style="width: 100%; height: 620px; border: 1px solid #2a2a2a; border-radius: 10px; background: #0f0f0f;">
</iframe>

In code, $F_0$ is evaluated using the standard error function (available in Python's default`math` module), and higher orders follow from the recurrence:

```python
from math import erf, sqrt, pi

def F0(t):
    """Boys function of order 0."""
    if t < 1e-7:
        return 1.0
    return 0.5 * sqrt(pi / t) * erf(sqrt(t))

def nuclear_int(alpha, A, beta, B, Z, C):
    """Nuclear attraction integral <phi_a | -Z/|r-C| | phi_b>."""
    AB2 = sum((a - b)**2 for a, b in zip(A, B))
    p   = alpha + beta
    P   = [(alpha * a + beta * b) / p for a, b in zip(A, B)]
    PC2 = sum((pi_ - c)**2 for pi_, c in zip(P, C))
    return (norm(alpha) * norm(beta) * (-2.0 * pi * Z / p)
            * F0(p * PC2) * exp(-alpha * beta / p * AB2))
```

The full nuclear attraction matrix is the sum over all nuclei: $V_{\mu\nu} = \sum_C V_{\mu\nu}^{(C)}$. The core Hamiltonian is then $\mathbf{H}^{\text{core}} = \mathbf{T} + \mathbf{V}$.

---

## Two-Electron Repulsion Integrals

We now arrive at the two-electron repulsion integrals (ERI) which are the most expensive object in a Hartree-Fock calculation

$$
(\mu\nu|\lambda\sigma) = \int\!\int \phi_\mu(\mathbf{r}_1)\,\phi_\nu(\mathbf{r}_1)\,\frac{1}{|\mathbf{r}_1 - \mathbf{r}_2|}\,\phi_\lambda(\mathbf{r}_2)\,\phi_\sigma(\mathbf{r}_2)\,\mathrm{d}\mathbf{r}_1\,\mathrm{d}\mathbf{r}_2.
$$

It represents the electrostatic repulsion between the **charge distribution** $\phi_{\mu}(\mathbf{r}_{1})\phi_{\nu}(\mathbf{r}_{1})$ and the charge distribution $\phi_{\lambda}(\mathbf{r}_{2})\phi_{\sigma}(\mathbf{r}_{2})$. In the Fock matrix, Coulomb integrals of the form $(\mu\nu\mid\sigma\lambda)$ represent the average repulsion felt by an electron in the $(\mu,\nu)$ distribution from all other electrons, and exchange integrals $(\mu\lambda\mid\sigma\nu)$ arise from the antisymmetry of the wavefunction.

### Analytic Evaluation

Apply the Gaussian product theorem separately to each pair:

- Pair $(\mu,\nu)$: product Gaussian at $\mathbf{P}$ with exponent $p = \alpha_\mu + \alpha_\nu$
- Pair $(\lambda,\sigma)$: product Gaussian at $\mathbf{Q}$ with exponent $q = \alpha_\lambda + \alpha_\sigma$

The integral now involves two one-centre Gaussians interacting via $1/\vert\mathbf{r}_{1} - \mathbf{r}_{2}\vert$. Apply the same integral representation of $1/r_{12}$ used for the nuclear attraction case. Both electron integrations become Gaussian integrals and can be evaluated in closed form. The remaining $t$-integral again yields a Boys function:

$$
(\mu\nu|\lambda\sigma) = \frac{2\pi^{5/2}}{pq\sqrt{p+q}}\,K_{\mu\nu}\,K_{\lambda\sigma}\,F_0\!\left(\zeta|\mathbf{P}-\mathbf{Q}|^2\right)\,N_\mu\,N_\nu\,N_\lambda\,N_\sigma,
$$

where the **reduced exponent** is:

$$
\zeta = \frac{pq}{p+q}.
$$

The argument $\zeta\vert\mathbf{P}-\mathbf{Q}\vert^2$ measures the distance between the two product centres, scaled by the reduced exponent. When both product centres coincide ($\vert\mathbf{P}-\mathbf{Q}\vert = 0$), $F_0(0) = 1$ and the integral takes the value of the self-repulsion. As the centres separate, the ERI decays, approaching the classical Coulomb law $1/R$ at large separation (which is intuitive for me at least (I think...)).

```python
from math import exp, pi, sqrt, erf

def dist2(A, B):
    return sum((a - b)**2 for a, b in zip(A, B))

def eri(alpha, A, beta, B, gamma, C, delta, D):
    """Two-electron repulsion integral (alpha A beta B | gamma C delta D), s-type GTOs."""
    p, P, K_AB = gaussian_product(alpha, A, beta,  B)
    q, Q, K_CD = gaussian_product(gamma, C, delta, D)
    PQ2  = dist2(P, Q)
    zeta = p * q / (p + q)
    pre  = (norm(alpha) * norm(beta) * norm(gamma) * norm(delta)
            * 2.0 * pi**2.5 / (p * q * sqrt(p + q))
            * exp(-alpha*beta/p*dist2(A,B) - gamma*delta/q*dist2(C,D)))
    return pre * F0(zeta * PQ2)
```

### The Coulomb Limit

A particularly illuminating special case is the integral $(AA\mid BB)$: both electrons in spherical Gaussian charge clouds, one centred at $\mathbf{A}$ and one at $\mathbf{B}$, all four basis functions sharing the same exponent $\alpha$. In this case the analytic formula reduces to:

$$
(AA|BB) = \frac{2}{\sqrt{\pi}}\sqrt{\alpha}\,F_0\!\left(\alpha\,R^2\right),
$$

where $R = \vert\mathbf{A}-\mathbf{B}\vert$. At $R = 0$ the integral equals $2\sqrt{\alpha}/\sqrt{\pi}$: a finite, well-defined self-repulsion with no Coulomb singularity (the Gaussian smearing regularises the $1/r$ divergence). At large $R$, using the asymptotic $F_0(t) \to \sqrt{\pi}/(2\sqrt{t})$:

$$
(AA|BB) \xrightarrow{R\to\infty} \frac{1}{R}.
$$

If we view the two Gaussian charge clouds from far away they look exactly like point charges, so the ERI recovers the classical Coulomb law. The crossover between the quantum-regularised short-range behaviour and the classical long-range Coulomb behaviour is controlled entirely by $\alpha$ — larger $\alpha$ (more compact Gaussians) means the crossover happens at shorter distances.

The interactive tool below shows $(AA\mid BB)$ as a function of separation $R$ for adjustable exponent $\alpha$. The dashed red line is the classical $1/R$ Coulomb reference, and the shaded region between the two curves shows the quantum correction due to the finite extent of the charge distributions.

<iframe
  src="{{ '/assets/interactive/hf_eri.html' | relative_url }}"
  title="Two-Electron Integral Explorer"
  loading="lazy"
  style="width: 100%; height: 620px; border: 1px solid #2a2a2a; border-radius: 10px; background: #0f0f0f;">
</iframe>

Drag the $\alpha$ slider and watch how the integral at $R = 0$ changes. For small $\alpha$ (diffuse Gaussians), the charge clouds are spread out, the self-repulsion is small, and the curve approaches $1/R$ quickly. For large $\alpha$ (compact Gaussians), the self-repulsion is higher and the quantum correction persists to larger separations.

### 8-Fold Permutation Symmetry

The ERI satisfies a set of eight permutation symmetries arising from the separability of the integrand:

$$
(\mu\nu|\lambda\sigma) = (\nu\mu|\lambda\sigma) = (\mu\nu|\sigma\lambda) = (\nu\mu|\sigma\lambda) = (\lambda\sigma|\mu\nu) = (\sigma\lambda|\mu\nu) = (\lambda\sigma|\nu\mu) = (\sigma\lambda|\nu\mu).
$$

The first four follow because the integrand is symmetric under swapping $\mu \leftrightarrow \nu$ or $\lambda \leftrightarrow \sigma$. The second four follow because the Coulomb operator is symmetric under swapping the two electrons $\mathbf{r}_{1} \leftrightarrow \mathbf{r}_{2}$. For a basis of $K$ functions, the number of unique ERIs is:

$$
N_{\text{unique}} = \frac{K(K+1)}{2} \cdot \frac{\frac{K(K+1)}{2}+1}{2} = \frac{K(K+1)(K^2+K+2)}{8}.
$$

For a minimal basis of $K = 2$ (our routine HeH example), this gives 6 unique integrals rather than $2^4 = 16$. For a calculation on a medium-sized molecule with $K = 200$, the reduction is from 1.6 billion to 200 million, arguably significant, but the quartic scaling $O(K^4)$ remains the fundamental cost of Hartree-Fock.

### Storing and Retrieving ERIs: the Yoshimine Sort

With many unique ERIs, the question of how to store and retrieve them efficiently becomes important. A flat array indexed by four basis function indices $(i,j,k,l)$ would be wasteful and require $O(K^4)$ memory. The standard solution (standard to me at least), introduced by Yoshimine, is to map the four indices to a single **compound index** that automatically respects all permutation symmetries:

```python
def eint(a, b, c, d):
    """Yoshimine compound index for four-index ERI, 1-based."""
    ab = a*(a+1)//2 + b if a > b else b*(b+1)//2 + a
    cd = c*(c+1)//2 + d if c > d else d*(d+1)//2 + c
    return ab*(ab+1)//2 + cd if ab > cd else cd*(cd+1)//2 + ab
```

This compresses the four-index storage into a Python dictionary keyed by a single integer, and any permutation of the four indices maps to the same key which means lookup is automatically symmetric.

---

## Putting Everything Together

Here is the complete, self-contained integral evaluation for the two-function HeH STO-1G system, combining all components derived above:

```python
import numpy as np
from math import pi, exp, sqrt, erf

def norm(a):
    return (2.0 * a / pi) ** 0.75

def dist2(A, B):
    return sum((a - b)**2 for a, b in zip(A, B))

def F0(t):
    if t < 1e-7: return 1.0
    return 0.5 * sqrt(pi / t) * erf(sqrt(t))

def overlap_int(a, A, b, B):
    return norm(a) * norm(b) * (pi / (a+b))**1.5 * exp(-a*b/(a+b) * dist2(A,B))

def kinetic_int(a, A, b, B):
    p = a + b
    return (a*b/p * (3.0 - 2.0*a*b/p * dist2(A,B))) * overlap_int(a, A, b, B)

def nuclear_int(a, A, b, B, Z, C):
    p = a + b
    P = [(a*x + b*y) / p for x, y in zip(A, B)]
    return (norm(a)*norm(b) * (-2.0*pi*Z/p)
            * F0(p * dist2(P, C)) * exp(-a*b/p * dist2(A, B)))

def eri(a, A, b, B, c, C, d, D):
    p = a + b;  P = [(a*x+b*y)/p for x,y in zip(A,B)]
    q = c + d;  Q = [(c*x+d*y)/q for x,y in zip(C,D)]
    zeta = p*q/(p+q)
    pre = (norm(a)*norm(b)*norm(c)*norm(d)
           * 2.0*pi**2.5/(p*q*sqrt(p+q))
           * exp(-a*b/p*dist2(A,B) - c*d/q*dist2(C,D)))
    return pre * F0(zeta * dist2(P, Q))

def eint(a, b, c, d):
    ab = a*(a+1)//2+b if a>b else b*(b+1)//2+a
    cd = c*(c+1)//2+d if c>d else d*(d+1)//2+c
    return ab*(ab+1)//2+cd if ab>cd else cd*(cd+1)//2+ab

# STO-1G HeH: hydrogen at origin, helium along z-axis
alphas  = [0.4166, 0.7739]
centers = [[0.0, 0.0, 0.0], [0.0, 0.0, 1.5117]]
nuclei  = [(1, centers[0]), (2, centers[1])]
dim     = 2

# One-electron matrices
S = np.zeros((dim, dim))
T = np.zeros((dim, dim))
V = np.zeros((dim, dim))

for i in range(dim):
    for j in range(dim):
        S[i,j] = overlap_int(alphas[i], centers[i], alphas[j], centers[j])
        T[i,j] = kinetic_int(alphas[i], centers[i], alphas[j], centers[j])
        for (Z, C) in nuclei:
            V[i,j] += nuclear_int(alphas[i], centers[i], alphas[j], centers[j], Z, C)

S = 0.5*(S + S.T);  T = 0.5*(T + T.T);  V = 0.5*(V + V.T)
Hcore = T + V

# Two-electron integrals (unique only, via Yoshimine indexing)
twoe = {}
for i in range(dim):
    for j in range(dim):
        for k in range(dim):
            for l in range(dim):
                key = eint(i+1, j+1, k+1, l+1)
                if key not in twoe:
                    twoe[key] = eri(alphas[i], centers[i], alphas[j], centers[j],
                                    alphas[k], centers[k], alphas[l], centers[l])

print("Overlap matrix S:")
print(np.round(S, 4))

print("\nCore Hamiltonian Hcore = T + V:")
print(np.round(Hcore, 4))

print("\nUnique two-electron integrals:")
for key, val in sorted(twoe.items()):
    print(f"  {key}: {val:.4f}")
```

The tool below computes and displays all four matrices for the HeH basis as you drag the bond length slider. The left side shows $\mathbf{S}$, $\mathbf{T}$, $\mathbf{V}$, and $\mathbf{H}^{\text{core}}$ as colour-coded grids (blue = positive, red = negative). The right panel plots the off-diagonal cross terms and one diagonal element vs $R$. Note how the coupling elements decay to zero as the bond stretches.

<iframe
  src="{{ '/assets/interactive/hf_fock.html' | relative_url }}"
  title="Fock Matrix Builder"
  loading="lazy"
  style="width: 100%; height: 650px; border: 1px solid #2a2a2a; border-radius: 10px; background: #0f0f0f;">
</iframe>

Running this reproduces the integral values quoted in both the [RHF post]({% post_url 2021-04-12-HartreeFockGuide %}) and the [UHF post]({% post_url 2025-01-01-UnrestrictedHartreeFock %}), confirming that the analytic formulae and the pre-computed data files agree:

```
Overlap matrix S:
[[1.     0.5028]
 [0.5028 1.    ]]

Core Hamiltonian Hcore = T + V:
[[-1.6606 -1.3159]
 [-1.3159 -2.7031]]

Unique two-electron integrals:
  5: 0.7283
  12: 0.3418
  14: 0.2192
  17: 0.5850
  19: 0.4368
  20: 0.9927
```

---

## Why the Scaling Matters

Every element $F_{\mu\nu}$ of the Fock matrix requires summing over all $\lambda$ and $\sigma$, each contributing an ERI. For $K$ basis functions, the number of unique ERIs scales as $K^4/8$, and forming the Fock matrix at each SCF iteration therefore scales as $O(K^4)$.

| System | Basis functions $K$ | Unique ERIs ($K^4/8$) | Feasible? |
|---|---|---|---|
| H₂ (STO-3G) | 2 | 6 | Yes |
| Water (6-31G\*) | 19 | ~3 300 | Yes |
| Benzene (cc-pVDZ) | 114 | ~21 million | Yes (minutes) |
| Protein (STO-3G) | ~5 000 | ~39 billion | Borderline |
| Protein (cc-pVDZ) | ~50 000 | ~39 trillion | No |

The jump from benzene to a protein is not a matter of buying better hardware, it is a change in the feasibility of the calculation. This is why the entire field of linear-scaling DFT and integral-direct methods exists (for a future post!).

There are strategies which break the quartic wall.

**Integral screening**: if $\vert K_{\mu\nu}\vert \cdot \vert K_{\lambda\sigma}\vert < \epsilon$, then $\vert(\mu\nu\mid\lambda\sigma)\vert \lesssim \epsilon \cdot (\mu\nu\mid\mu\nu)^{1/2} \cdot (\lambda\sigma\mid\lambda\sigma)^{1/2}$ (Cauchy-Schwarz inequality applied to ERIs). The exponential decay of $K_{AB}$ with separation means that for large molecules, most integral pairs involve at least one nearly-zero $K$ factor and can be discarded. This brings the effective scaling closer to $O(K^2)$ for spatially extended systems.

The tool below shows this decay. The left panel plots $K_{AB}(R) = e^{-\alpha R^2/2}$ for several exponents; the dashed line is the screening threshold $\epsilon$ and the vertical markers show where each curve crosses it ($R_{\text{screen}}$). The right panel shows $R_{\text{screen}}$ vs $\alpha$, compact GTOs (large $\alpha$) screen at shorter range, diffuse GTOs need larger separations. Drag the threshold slider to watch the entire landscape shift.

<iframe
  src="{{ '/assets/interactive/hf_screening.html' | relative_url }}"
  title="Integral Screening Explorer"
  loading="lazy"
  style="width: 100%; height: 620px; border: 1px solid #2a2a2a; border-radius: 10px; background: #0f0f0f;">
</iframe>

**Resolution of identity (density fitting)**: the four-centre two-electron integrals are approximated by inserting an auxiliary basis, converting the quartic problem into a cubic one. This is the basis of methods like RI-MP2 and DF-HF.

---

## Conclusions

This ended up being more complicated than I had intended but I hope it has given a more intuitive look into the integrals which underpin HF theory. 

Every term in the Hartree-Fock Fock matrix reduces to one of four integral types: overlap $\mathbf{S}$, kinetic $\mathbf{T}$, nuclear attraction $\mathbf{V}$ (all one-electron), and the two-electron repulsion integral. All four are made analytically tractable by the same two ideas: the Gaussian product theorem, which collapses multi-centre integrals into single-centre Gaussian integrals, and the Boys function, which handles the Coulomb operator $1/r$ that would otherwise prevent a closed-form answer.

The ERI is the most important and most expensive component. Its quartic count is why Hartree-Fock is not used for large systems, why integral screening is essential even for medium-sized molecules, and why the entire enterprise of post-HF methods is built on the assumption that these integrals have already been computed and can be stored.

For a new person to the field I personally think understanding where these numbers come from not just how to call a library function that computes them gives a much clearer picture of why Hartree-Fock has the structure it does, and why every method built on top of it inherits both its strengths and its limitations.
---

## Further Reading

1. **Szabo and Ostlund**, *Modern Quantum Chemistry: Introduction to Advanced Electronic Structure* — Appendix A gives the complete derivation of GTO integrals. The standard reference.

2. **Helgaker, Jørgensen and Olsen**, *Molecular Electronic-Structure Theory* — Chapter 9 covers integral evaluation in full generality including higher angular momentum, the McMurchie-Davidson and Obara-Saika recurrence schemes.

3. **Boys, S. F.** (1950), "Electronic wave functions I. A general method of calculation for the stationary states of any molecular system," *Proc. R. Soc. A* **200**, 542. The original paper introducing GTOs and the Boys function.

4. **Yoshimine, M.** (1973), "The use of direct access devices in handling large lists of two-electron integrals and CI Hamiltonian," IBM Research Report RJ-555. The source of the Yoshimine sort.

5. [Crawford's Programming Projects](https://github.com/CrawfordGroup/ProgrammingProjects) — a guided series of HF programming exercises that covers integral evaluation in detail, highly recommended for working through the ideas by hand.
