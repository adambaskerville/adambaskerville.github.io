---
layout: post
title: "T>T: The First Solution is Rarely the Best: A Lesson in Numerical Integration"
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

### Why Gaussian Quadrature is Exact in This Instance (Not an Approximation)

When I hear "numerical integration," I think of Riemann sums or the Trapezoidal rule that chop an area into tiny pieces to get a "good enough" approximation. I was mistakenly placing Gaussian quadrature in this category, but it is different. 

The fundamental theorem of Gaussian quadrature states that an $N$-point rule will integrate a polynomial of degree up to $2N - 1$ with **mathematical exactness**. There is no truncation error and no missing area; the numerical sum is analytically equivalent to the true integral.

To make this clearer for our two-electron integrals, we just need to count the polynomial degrees in our perimetric integrand. 

Once we separate the exponential decay $e^{-A(r_1+r_2)}$ (which the Gauss-Laguerre method natively absorbs into its weight functions), we are left analysing a pure polynomial. Let's look at the maximum possible degree in any single coordinate (e.g., $z_1$) for a basis set of size $m$:

1. **The Basis Functions:** We are multiplying four Laguerre polynomials together to form the $(pq\|uv)$ integral: $L_p \times L_q \times L_u \times L_v$. Since the maximum index of our basis functions is $m-1$, multiplying four of them together yields a maximum polynomial degree of $4(m-1) = 4m - 4$.
2. **The Jacobian Volume Element:** Converting the Cartesian volume to perimetric coordinates introduces the polynomial $(z_2+z_3)(z_3+z_1)(z_1+z_2)$. If you expand this out, the highest power any single coordinate will reach (e.g., $z_1^2$) is exactly 2.

Adding these together, the absolute highest polynomial degree we will ever encounter along a single integration axis is

$$(4m - 4) + 2 = \mathbf{4m - 2}$$

Now, we just apply the Gaussian quadrature rule:

$$
\begin{aligned}
    2N - 1 &\ge 4m - 2 \\
    2N &\ge 4m - 1 \\
    N &\ge 2m - 0.5
\end{aligned}
$$


This tells us that an $N = 2m$ point grid is mathematically sufficient. By setting our code to $N = 2m + 1$, we guarantee that the highest degree terms of our largest basis functions are captured exactly. 

We didn't swap an exact analytical method for an approximate numerical one. We swapped a computationally inefficient exact method (the binomial expansion) for a highly stable, computationally elegant exact method.

### The Arbitrary Precision Solution with Python

Even with quadrature, the massive polynomial values at $m=20$ still cause precision issues in standard `NumPy` float arrays. But instead of relying on heavy C++ binaries, we can solve this cleanly in pure Python using `mpmath` for arbitrary precision.

Here is the core logic of the new implementation. It drops the C++ overhead entirely and allows us to generate the full $(pq\|uv)$ tensor with a canonical Yoshimine sort (see my original [Hartree Fock]({% post_url 2021-04-12-HartreeFockGuide %}) post to learn more about this)  in seconds, setting the environment to 50 decimal places to absorb any residual cancellation.


```python
import mpmath
import time

# Set arbitrary precision to 50 decimal places
mpmath.mp.dps = 50


def get_gauss_laguerre_mp(N):
    """
    Generates Gauss-Laguerre roots and weights using the Golub-Welsch
    algorithm in arbitrary precision.
    """
    T = mpmath.matrix(N, N)
    for i in range(N):
        # Diagonal elements: alpha_i = 2i + 1
        T[i, i] = mpmath.mpf(2 * i + 1)
        if i < N - 1:
            # Off-diagonal elements: sqrt(beta_i) = i + 1
            val = mpmath.mpf(i + 1)
            T[i, i + 1] = val
            T[i + 1, i] = val

    # Diagonalize the symmetric Jacobi matrix
    E, ER = mpmath.eigsy(T)

    roots = [E[i] for i in range(N)]
    # Weights are the square of the first component of each normalised eigenvector
    weights = [ER[0, i] ** 2 for i in range(N)]

    return roots, weights


def generate_canonical_integrals_mp(m, A):
    """
    Computes and canonically sorts the unique (pq|uv) two-electron repulsion
    integrals using arbitrary precision Gauss-Laguerre Quadrature.
    """
    N = 2 * m + 1

    print(
        f"-> Generating Gauss-Laguerre grid (N={N}) at {mpmath.mp.dps} decimal places..."
    )
    roots, weights = get_gauss_laguerre_mp(N)

    print("-> Pre-computing coordinate grid and Laguerre polynomials...")
    # 2D Grid: X_jk = x_j + 0.5 * x_k
    X_jk = [
        [roots[j] + mpmath.mpf("0.5") * roots[k] for k in range(N)] for j in range(N)
    ]

    L_grid = []
    for p in range(m):
        # mpmath.laguerre(n, alpha, x)
        grid_p = [
            [mpmath.laguerre(p, 0, X_jk[j][k]) for k in range(N)] for j in range(N)
        ]
        L_grid.append(grid_p)

    print("-> Separating 3D integral and evaluating inner sums S_pq...")
    S = {}
    for p in range(m):
        for q in range(p, m):
            S_pq = []
            for k in range(N):
                val = mpmath.mpf(0)
                for j in range(N):
                    term = weights[j] * X_jk[j][k] * L_grid[p][j][k] * L_grid[q][j][k]
                    val += term
                S_pq.append(val)
            # Store symmetrically
            S[(p, q)] = S_pq
            S[(q, p)] = S_pq

    num_pairs = m * (m + 1) // 2
    total_unique = num_pairs * (num_pairs + 1) // 2
    print(f"-> Assembling {total_unique} canonically sorted unique integrals...")

    scalar = (mpmath.mpf(8) * mpmath.pi**2) / (mpmath.mpf(A) ** 5)

    unique_integrals = []
    index_1d = 0

    # Yoshimine Sort: p >= q, u >= v, pq >= uv
    for p in range(m):
        for q in range(p + 1):
            idx_pq = p * (p + 1) // 2 + q
            for u in range(m):
                for v in range(u + 1):
                    idx_uv = u * (u + 1) // 2 + v

                    if idx_pq >= idx_uv:
                        # Final contraction over k
                        val = mpmath.mpf(0)
                        for k in range(N):
                            term = weights[k] * S[(p, q)][k] * S[(u, v)][k]
                            val += term

                        final_val = scalar * val
                        unique_integrals.append((index_1d, p, q, u, v, final_val))
                        index_1d += 1

    return unique_integrals


# ==========================================
# Execution Example
# ==========================================
if __name__ == "__main__":
    m_basis = 20
    A_val = 1.0  # Set to 1.0 for unscaled reference validation

    start_time = time.time()
    canonical_integrals = generate_canonical_integrals_mp(m_basis, A_val)
    elapsed = time.time() - start_time

    total_unique = len(canonical_integrals)
    print(f"\nCompleted {total_unique} integrals in {elapsed:.2f} seconds.")

    output_filename = f"canonical_integrals_mpmath_m{m_basis}.txt"
    print(f"Writing arbitrary-precision output to {output_filename}...")

    with open(output_filename, "w") as f:
        f.write(f"# Canonical Two-Electron Integrals (m={m_basis}, A={A_val})\n")
        f.write(f"# Computed with mpmath at {mpmath.mp.dps} decimal places\n")
        f.write("# Format: 1D_Index | p  q  u  v | Integral_Value\n")
        f.write("-" * 65 + "\n")

        for row in canonical_integrals:
            idx = int(row[0])
            # +1 shift for 1-based orbital index notation
            p = int(row[1]) + 1
            q = int(row[2]) + 1
            u = int(row[3]) + 1
            v = int(row[4]) + 1

            # Print to 30 decimal places
            val_str = mpmath.nstr(row[5], 30, min_fixed=0)

            # Formatting to handle negative signs cleanly
            f.write(f"{idx:6d} | {p:2d} {q:2d} {u:2d} {v:2d} |  {val_str}\n")

    print("Done.")
```
This code is **much** simpler than the original and significantly faster whilst offering the required accuracy for each integral. Even better on my MacBook I can generate all 22155 integrals for a 20-term Hartree Fock wavefunction in just 3.7 seconds! The more you learn.