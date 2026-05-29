---
layout: post
title: "The First Solution is Rarely the Best: A Lesson in Numerical Integration"
date: 2025-07-13
categories: [Computational Chemistry, Python, Mathematics]
tags: [Hartree-Fock, Numerical Integration, Arbitrary Precision, mpmath, Python]
math: true
---

If you spend enough time developing scientific software, you eventually learn that the first solution to a problem is rarely the best one. Science, and the code we write to execute it, is a continual learning process which is one of the things I love most about the field. I have been repeatedly guilty of brute-forcing a mathematical derivation to prove it *can* be done, only to realise later that I was doing it the hard way.

Today, I want to pull back the curtain on a specific problem I tackled a few years ago involving the exact calculation of Hartree-Fock (HF) two-electron integrals for two-electron atoms, and how I only realised today that there is a much simpler trick which I was not aware of at the time.

### The Original Approach: The 7-Nested Summation

Back in 2018, I co-authored a paper titled *Hartree-Fock implementation using a Laguerre-based wavefunction for the ground-state and correlation energies of two-electron atoms* which can be found [here](https://sussex.figshare.com/articles/journal_contribution/Hartree-Fock_implementation_using_a_Laguerre-based_wavefunction_for_the_ground-state_and_correlation_energies_of_two-electron_atoms/23451620?file=41160692). Our goal was to implement the HF method using a Laguerre-based wavefunction in perimetric coordinates. 

Because the Laguerre functions are orthogonal on the interval $[0, \infty)$, the one-electron integrals were elegantly solved using series solutions. However, the two-electron Coulomb and exchange integrals were more difficult. The integration generated terms that did not satisfy the Laguerre orthogonality condition. 

To solve them analytically, we converted the integrals to perimetric coordinates, expressed the Laguerre polynomials as power series via the binomial expansion, and solved them term-by-term. Mathematically, it was rigorous. Computationally, it was a behemoth resulting in the following 7-nested summation


$$
\begin{aligned}
&\sum_{p_i,q_i,u_i,v_i=0}^{p,q,u,v}
 \sum_{a_i=0}^{p_i+q_i}
 \sum_{b_i=0}^{p_i}
 \sum_{c_i=0}^{a_i+b_i-\phi+1}
 (-1)^{\phi+1} \pi^2 \\
&\quad \times
\Big(
a_i^2 - 2a_i b_i + b_i^2
- p^2
- 2p_i q_i
- 2p_i u_i
- 2p_i v_i
- q_i^2 \\
&\qquad
- 2q_i u_i
- 2q_i v_i
- u^2
- 2u_i v_i
- v^2
+ a_i + b_i
- 7\phi - 10
\Big) \\
&\quad \times
(\phi - a_i - b_i)!
(p_i + q_i)!
(u_i + v_i)! \\
&\quad \times
\frac{
p!q!u!v!
}{
p_i!^2 q_i!^2 u_i!^2 v_i!^2
(p_i+q_i-a_i)!
(q_i+u_i-b_i)!
(p-p_i)!
(q-q_i)!
(u-u_i)!
(v-v_i)!
}.
\end{aligned}
$$

What is the problem? Just look at it... also catastrophic cancellation. When dealing with alternating signs $(-1)^{\phi+1}$ and large factorials, standard 64-bit floating-point numbers fail. At higher basis set sizes, the polynomial terms reach magnitudes of $10^{25}$, while the final integrated volume is around $10^2$ so we lose all significant digits. 

To bypass this, I originally wrote a dedicated C++ program utilising the Arb ball arithmetic library to handle the massive factorials and enforce rigourous error bounds. It worked, but it was a heavy, memory-hungry hammer (although I did learn how to use the Arb library which I think is fantastic).

### The Epiphany: Gauss-Laguerre Quadrature

Recently, while solving an unrelated problem, I realised I had missed a more elegant mathematical mapping. Instead of expanding the polynomials algebraically and fighting the resulting factorials, we can evaluate the integral exactly using **Gauss-Laguerre Quadrature**.

Because the $1/r_{12}$ Coulomb operator perfectly cancels with the $r_{12}$ term in the perimetric Jacobian, the spatial integration becomes a pure polynomial (up to degree $4m$) multiplied by a separable exponential decay $e^{-A(r_1+r_2)}$. 

This takes the exact form required for Gauss-Laguerre quadrature

$$
\iiint_0^\infty P(z_1, z_2, z_3) e^{-\frac{A}{2}z_1} e^{-\frac{A}{2}z_2} e^{-Az_3} dz_1 dz_2 dz_3
$$

I have used Gaussian quadrature numerous times before and automatically think of it being an approximate, numerical technique but it is not so here; an $N$-point rule (where $N = 2m + 1$) integrates a polynomial of degree $\le 4m$ *exactly*. 

By mapping the integral to the quadrature nodes, the 3D integral algebraically separates into independent 1D projections. The $O(m^7)$ nightmare collapses into an $O(m^4)$ dynamic evaluation. 

### The Arbitrary Precision Solution with Python

Even with quadrature, the massive polynomial values at $m=20$ still cause precision issues in standard `NumPy` float arrays. But instead of relying on heavy C++ binaries, we can solve this cleanly in pure Python using `mpmath` for arbitrary precision.

Here is the core logic of the new implementation. It drops the C++ overhead entirely and allows us to generate the full $(pq\|uv)$ tensor with a canonical Yoshimine sort (see my original [Hartree Fock]({% post_url 2021-04-12-HartreeFockGuide %}) post to learn more about this)  in seconds, setting the environment to 50 decimal places to absorb any residual cancellation.


```python
import mpmath
import time

# Set precision to 50 decimal places to absorb cancellation
mpmath.mp.dps = 50

def get_gauss_laguerre_mp(N):
    """Generates exact Gauss-Laguerre roots and weights using Golub-Welsch."""
    T = mpmath.matrix(N, N)
    for i in range(N):
        T[i, i] = mpmath.mpf(2 * i + 1)
        if i < N - 1:
            val = mpmath.mpf(i + 1)
            T[i, i + 1] = val
            T[i + 1, i] = val
            
    E, ER = mpmath.eigsy(T)
    roots = [E[i] for i in range(N)]
    weights = [ER[0, i]**2 for i in range(N)]
    return roots, weights

def generate_canonical_integrals(m, A):
    N = 2 * m + 1
    roots, weights = get_gauss_laguerre_mp(N)
    
    # 2D Grid: X_jk = x_j + 0.5 * x_k
    X_jk = [[roots[j] + mpmath.mpf('0.5') * roots[k] for k in range(N)] for j in range(N)]
    
    # Precompute Laguerre Polynomials
    L_grid = [[[mpmath.laguerre(p, 0, X_jk[j][k]) for k in range(N)] for j in range(N)] for p in range(m)]
        
    # Evaluate inner sums (Separation of Variables)
    S = {}
    for p in range(m):
        for q in range(p, m):
            S_pq = []
            for k in range(N):
                val = mpmath.mpf(0)
                for j in range(N):
                    val += weights[j] * X_jk[j][k] * L_grid[p][j][k] * L_grid[q][j][k]
                S_pq.append(val)
            S[(p, q)] = S_pq
            S[(q, p)] = S_pq

    # Assemble unique integrals via canonical sorting
    scalar = (mpmath.mpf(8) * mpmath.pi**2) / (mpmath.mpf(A)**5)
    unique_integrals = []
    
    # Yoshimine Sort: p >= q, u >= v, pq >= uv
    for p in range(m):
        for q in range(p + 1):
            idx_pq = p * (p + 1) // 2 + q
            for u in range(m):
                for v in range(u + 1):
                    idx_uv = u * (u + 1) // 2 + v
                    
                    if idx_pq >= idx_uv:
                        val = mpmath.mpf(0)
                        for k in range(N):
                            val += weights[k] * S[(p, q)][k] * S[(u, v)][k]
                        
                        unique_integrals.append((p+1, q+1, u+1, v+1, scalar * val))
                        
    return unique_integrals

if __name__ == "__main__":
    m_basis = 20
    A_val = 1.0 
    
    print(f"Calculating and sorting unique integrals for m = {m_basis}...")
    canonical_integrals = generate_canonical_integrals(m_basis, A_val)
    print(f"Extracted {len(canonical_integrals)} unique integrals.")
```
This code is **much** simpler than the original and significantly faster whilst offering the required accuracy for each integral. Even better on my MacBook I can generate all 22155 integrals for a 20-term Hartree Fock wavefunction in just 3.7 seconds! The more you learn.