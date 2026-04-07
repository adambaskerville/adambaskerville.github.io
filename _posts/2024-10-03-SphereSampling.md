---
layout: post
title: "T>T: On the Use of Fibonacci Lattices for Spherical Point Sets"
date: 2024-10-03
excerpt: "Exploring why random point placement on a sphere creates voids even with the correct distribution, and how using a sunflower lattice uses the irrationality of the golden ratio to produce near-optimal coverage."
tags: [science, mathematics, programming, python, numpy, sphere, sampling, molecular dynamics, fibonacci, golden ratio, computational chemistry, quadrature]
comments: false
math: true
---

## The Problem

I came across an issue on a project I was working on a few years ago which, for whatever reason, popped into my head yesterday whilst building a shed. I was working on a method for reconstructing ESP surfaces using spherical harmonic coefficients, where the surface needs to be sampled at $N$ points distributed around the molecule. For small angular momentum quantum numbers $l$ the result is largely insensitive to the exact point placement, any reasonable distribution works. But as $l$ grows, the reconstruction quality degrades sharply if the sampling is uneven, and the degradation is not subtle. The RMSD between the true and reconstructed ESP surface grows by an order of magnitude past $l \approx 20$ for naive random sampling.

## The Naive Approach: Getting the Distribution Wrong

The first mistake I made is to pick $\theta$ and $\phi$ uniformly at random from their natural ranges:

$$
    \theta \sim \mathcal{U}(0, \pi), \quad \phi \sim \mathcal{U}(0, 2\pi).
$$

Both coordinates are sampled uniformly. What could be wrong?

The spherical area element is

$$
    dA = \sin\theta \, d\theta \, d\phi.
$$

Near the equator, where $\theta \approx \pi/2$, the factor $\sin\theta \approx 1$ so each increment of angle corresponds to roughly one increment of surface area. Near the poles, where $\theta \approx 0$ or $\theta \approx \pi$, $\sin\theta \approx 0$: a range of angles that subtends almost no surface area at all. Sampling $\theta$ uniformly in $[0, \pi]$ allocates the same number of points to polar regions that have almost no surface area as it allocates to equatorial regions that have a great deal. The result is severe over-crowding near the poles and under-coverage at mid-latitudes.

```python
import numpy as np

# Incorrect: uniform in angle clusters near the poles
def sample_wrong(n: int) -> tuple[np.ndarray, np.ndarray]:
    theta = np.random.uniform(0, np.pi, n)
    phi   = np.random.uniform(0, 2 * np.pi, n)
    return theta, phi
```

The over-concentration at the poles is not subtle. The fraction of total sphere area contained within a polar cap of half-angle $\theta_0$ is

$$
    f(\theta_0) = \tfrac{1}{2}(1 - \cos\theta_0).
$$

For $\theta_0 = 10°$ this is about 1.5% of the sphere, yet uniform-$\theta$ sampling allocates $10/180 \approx 5.6\%$ of points to that region, nearly four times too many near each pole.

## The Correct Random Approach: Still Has Problems

The standard fix for pole clustering is well-known (once I looked into it further). Rather than sampling $\theta$ uniformly, sample $u = \cos\theta$ uniformly on $[-1, 1]$:

$$
    u \sim \mathcal{U}(-1, 1), \quad \theta = \arccos(u).
$$

This accounts for the $\sin\theta$ factor in the area element and produces a correctly uniform distribution over the sphere surface. It appears throughout computational physics textbooks including *Numerical Recipes in C* (Press et al., 1992) as the standard method for sphere sampling, and is what the code I was working with used:

```python
import numpy as np


def generate_random_points_on_sphere(n_samples: int) -> tuple[np.ndarray, np.ndarray]:
    """Uniform random points on the sphere using the arccos transformation."""
    u      = np.random.uniform(size=n_samples)
    thetas = np.arccos(2 * u - 1)
    phis   = 2 * np.pi * np.random.uniform(size=n_samples)
    return thetas, phis
```

The transformation `arccos(2u - 1)` with $u \in [0, 1]$ is equivalent to sampling $\cos\theta$ uniformly on $[-1, 1]$. The marginal distribution is now correct. So why does the resulting surface show holes?

The answer is that the distribution is correct but the sample is still **random**, and random processes have inherent clustering and voids.

## Why Randomness Itself Is the Problem

Consider placing $N$ points on a sphere according to a uniform Poisson process (which is what independent random sampling produces). What is the probability that some particular region of area $A$ is completely empty?

If the total sphere surface area is $4\pi$, each point independently falls inside the region with probability $p = A / 4\pi$. The probability that all $N$ points miss the region is

$$
    P(\text{empty}) = (1 - p)^N \approx e^{-Np} = e^{-NA/4\pi}.
$$

Now choose the region to be a "typical cell", the equal-area patch we might hope each of the $N$ points is responsible for, with area $A = 4\pi/N$. Then

$$
    P(\text{empty cell}) \approx e^{-1} \approx 0.368.
$$

For any fixed equal-area partition of the sphere into $N$ cells, the **expected empty-cell fraction tends to about 37%** under random Poisson placement. This is not an edge case like I thought it might be or a consequence of using the wrong distribution; it is a direct consequence of Poisson statistics and does not disappear as $N$ grows.

The practical effect in my case was direct. The spherical harmonic reconstruction at order $l$ uses basis functions that oscillate at angular frequencies proportional to $l$. With $N$ random points, the effective Nyquist limit is not clean: some regions are oversampled, some are empty, and the reconstruction misinterprets the imbalance as high-frequency oscillation in the potential surface. This was the root cause of the RMSD instability past $l \approx 20$.

## The Fibonacci Lattice

The solution for me was to abandon randomness entirely and use a **deterministic** construction designed to cover the sphere as uniformly as possible. I thought the most elegant construction to be the **Fibonacci lattice**, also called the **sunflower spiral**.

The $i$-th point, for $i = 0, 1, \ldots, N-1$, is placed at

$$
\cos\theta_i = 1 - \frac{2(i + \tfrac{1}{2})}{N}, \qquad
\phi_i = \frac{2\pi \, i}{\varphi},
$$

where $\varphi = (1 + \sqrt{5})/2 \approx 1.618\ldots$ is the **golden ratio**.

The two ingredients are independent and each solves a distinct part of the coverage problem:

**The $\theta_i$ spacing.** The values $\cos\theta_i = 1 - (2i+1)/N$ are the midpoints of $N$ equal-width bins partitioning $[-1, 1]$. Since

$$
dA = \sin\theta \, d\theta \, d\phi = d(\cos\theta)\, d\phi,
$$

spacing $\cos\theta$ equally ensures that each latitude band carries exactly the same surface area, $4\pi/N$ per point. No latitude is over- or under-covered.

**The $\phi_i$ spacing.** Successive points are rotated azimuthally by $2\pi/\varphi \approx 222.49°$, which is equivalent modulo a full turn to the more familiar golden angle of $360^\circ(1 - 1/\varphi) \approx 137.51^\circ$. This irrational rotation is chosen to prevent any azimuthal bunching.

## The Golden Ratio and Why It Works

Why does $\varphi$ produce uniquely good azimuthal coverage?

Every irrational number $x$ has a continued fraction expansion

$$
x = a_0 + \cfrac{1}{a_1 + \cfrac{1}{a_2 + \cfrac{1}{a_3 + \cdots}}},
$$

and the speed at which the rational convergents $p_k/q_k$ approximate $x$ is governed by how large the partial quotients $a_k$ are. A large partial quotient at position $k$ means the number is "almost rational" near that convergent; this forces the sequence $\{kx \bmod 1\}$ to cluster badly around a finite set of values for multiples $k \approx q_k$.

The golden ratio has the continued fraction $\varphi = [1; 1, 1, 1, \ldots]$,  all partial quotients equal to 1, the smallest possible value. This makes $\varphi$ the **most irrational** real number in a precise sense: its rational approximants converge as slowly as the theory allows. The consequence for the sequence $\{i/\varphi \bmod 1\}$ is that no rational approximation ever causes systematic clustering.

This is formalised by the **three-distance theorem** (Steinhaus, 1958): for any $N$ terms of the sequence $\{k \alpha \bmod 1\}$ with irrational $\alpha$, the $N$ points partition $[0,1)$ into gaps of exactly two or three distinct lengths. For $\alpha = 1/\varphi$ the gaps remain as small and as equal as the theory allows at every $N$, no other irrational does better (I think, send me a postcard if I am incorrect).

The consequence is the spiral you see in sunflower seed heads, pine cones, and daisy florets. Interestingly, billions of years of selection pressure may have found the same optimum.

## Implementation

Here is a clean Python implementation of both approaches:

```python
import numpy as np


def random_sphere_points(n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Uniform random points on the unit sphere.
    Uses the arccos(2u - 1) transformation to avoid pole clustering.
    """
    cos_theta = np.random.uniform(-1.0, 1.0, n)
    theta     = np.arccos(cos_theta)
    phi       = np.random.uniform(0.0, 2.0 * np.pi, n)
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return x, y, z


def fibonacci_sphere_points(n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Fibonacci lattice points on the unit sphere.
    Deterministic, maximally uniform coverage.
    """
    golden_ratio = (1.0 + np.sqrt(5.0)) / 2.0
    i         = np.arange(n)
    cos_theta = 1.0 - 2.0 * (i + 0.5) / n
    theta     = np.arccos(cos_theta)
    phi       = 2.0 * np.pi * i / golden_ratio
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return x, y, z
```

This is also slightly faster than the random version in practice since it avoids random number generation overhead entirely.

## Visualising the Difference

The 2D longitude-latitude projection proved to be revealing. Plot $\phi$ on the horizontal axis and $\theta$ on the vertical.

```python
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

N = 500

fig = plt.figure(figsize=(16, 10))
fig.patch.set_facecolor('#0d1117')

methods = {
    'Random (correct arccos method)': random_sphere_points(N),
    'Fibonacci Lattice':               fibonacci_sphere_points(N),
}

for col, (title, (x, y, z)) in enumerate(methods.items()):
    theta = np.arccos(np.clip(z, -1, 1))
    phi   = np.arctan2(y, x)

    # 3D scatter
    ax3 = fig.add_subplot(2, 2, col + 1, projection='3d')
    ax3.set_facecolor('#0d1117')
    ax3.scatter(x, y, z, s=8, c=theta, cmap='plasma', alpha=0.85, linewidths=0)
    ax3.set_title(title, color='white', fontsize=12, pad=8)
    ax3.set_box_aspect([1, 1, 1])
    for pane in [ax3.xaxis.pane, ax3.yaxis.pane, ax3.zaxis.pane]:
        pane.set_visible(False)
    ax3.grid(False)
    for axis in [ax3.xaxis, ax3.yaxis, ax3.zaxis]:
        axis.set_ticklabels([])

    # 2D projection
    ax2 = fig.add_subplot(2, 2, col + 3)
    ax2.set_facecolor('#0d1117')
    ax2.scatter(np.degrees(phi), np.degrees(theta),
                s=5, c=theta, cmap='plasma', alpha=0.75, linewidths=0)
    ax2.set_xlabel('Azimuthal angle φ (°)', color='#8b949e', fontsize=10)
    ax2.set_ylabel('Polar angle θ (°)',     color='#8b949e', fontsize=10)
    ax2.set_xlim(-180, 180)
    ax2.set_ylim(0, 180)
    ax2.tick_params(colors='#8b949e')
    for spine in ax2.spines.values():
        spine.set_edgecolor('#30363d')

plt.tight_layout()
plt.savefig('sphere_comparison.png', dpi=150, bbox_inches='tight',
            facecolor='#0d1117')
plt.show()
```

The 2D projection of the Fibonacci lattice immediately reveals the spiral: points lie on two families of diagonal lines at gradients set by the Fibonacci numbers, exactly as in a sunflower head. The random distribution, even with the correct marginal, shows obvious dense regions and conspicuous gaps.

## Measuring Uniformity

Lets quantify this better. A natural measure of coverage quality is the **nearest-neighbour distance distribution**: for a perfectly uniform point set all nearest-neighbour distances would be identical.

```python
from scipy.spatial import cKDTree

def uniformity_stats(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> dict:
    pts  = np.column_stack([x, y, z])
    tree = cKDTree(pts)
    # k=2 because k=1 is the point itself
    dists, _ = tree.query(pts, k=2)
    nn = dists[:, 1]
    return {
        'mean':      nn.mean(),
        'std':       nn.std(),
        'min':       nn.min(),
        'max':       nn.max(),
        'max/min':   nn.max() / nn.min(),
        'CV':        nn.std() / nn.mean(),
    }
```

For $N = 500$ a representative comparison looks like this:

| Metric | Random | Fibonacci |
|:-------|-------:|----------:|
| Mean NN distance | 0.113 | 0.113 |
| Std NN distance  | 0.032 | 0.007 |
| Max/Min ratio    | ~5.1  | ~1.33  |
| Coeff. of variation | ~28% | ~6% |

Both methods have the same mean nearest-neighbour distance (controlled entirely by $N$). But the Fibonacci lattice has a standard deviation four times smaller and a max/min ratio nearly four times smaller. Crucially, this advantage does not degrade as $N$ grows: the max/min ratio for Fibonacci converges to approximately 1.33 regardless of $N$ while for random sampling it grows without bound as voids inevitably appear somewhere on the sphere.

## The Interactive Visualiser

The visualiser below lets you explore this directly. Toggle between the 3D sphere and the 2D $(\phi, \theta)$ map. Each point is coloured by its nearest-neighbour distance: **blue = clustered**, **red = isolated (in a void)**. For a perfectly uniform distribution all points would be the same colour.

Click "New Random Sample" to draw a fresh set of random points and watch the voids move to different locations.

<iframe
  src="{{ '/assets/interactive/sphere_sampling_visualiser.html' | relative_url }}"
  title="Sphere sampling visualiser"
  loading="lazy"
  style="width: 100%; height: 860px; border: 1px solid #dbe4ee; border-radius: 12px; background: #ffffff;">
</iframe>

If you want to open it on its own page, use [this standalone version]({{ '/assets/interactive/sphere_sampling_visualiser.html' | relative_url }}).

A few things to look for:

- At $N = 100$, the random voids are large and obvious in both the 3D and 2D views
- At $N = 1000$, the random sample looks good to the eye in 3D but the 2D map still shows clear structure in the clustering
- The Fibonacci 2D map shows the characteristic diagonal-stripe pattern at every $N$, this is the golden-angle spiral projected flat, and it is the geometric reason coverage is good
- The max/min ratio in the stats panel stays near 1.33 for Fibonacci and grows towards 5 or more for random as $N$ increases

## Application: ESP Surface Reconstruction

To make this concrete, consider the ESP reconstruction problem directly. The electrostatic potential on a sphere of radius $r$ around a molecule can be expanded in real spherical harmonics

$$
V(r, \theta, \phi) = \sum_{l=0}^{l_{\max}} \sum_{m=-l}^{l} c_{lm} \, Y_l^m(\theta, \phi).
$$

With $N$ sampled values $V_i = V(\theta_i, \phi_i)$ this becomes a linear system

$$
\mathbf{Y} \, \mathbf{c} = \mathbf{V},
$$

where $Y_{ij} = Y_{l(j)}^{m(j)}(\theta_i, \phi_i)$. The numerical stability of the least-squares solution depends on the condition number of $\mathbf{Y}^\top \mathbf{Y}$, which in turn depends on how orthogonal the spherical harmonics are when evaluated at the sample points.

For a perfectly uniform grid, the spherical harmonics are near-orthogonal up to $l_{\max}^2 \approx N$, and the condition number stays close to 1. For a random distribution, clusters and voids break this approximate orthogonality, and the condition number grows rapidly with $l_{\max}$.

The practical outcome:

- **Random sampling**: RMSD between true and reconstructed ESP surface degrades past $l \approx 20$ regardless of $N$, because there are always voids that misrepresent the surface
- **Fibonacci sampling**: reconstruction stays stable to much larger $l$, limited by the true information content of the surface rather than by sampling artefacts
- The improvement in RMSD at high $l$ is roughly an order of magnitude

## Roberts' Offset

Something I only just found out (whilst writing this post!) is that Roberts (2018) showed that applying a small offset $\varepsilon$ to the index improves the pole behaviour:

$$
\cos\theta_i = 1 - \frac{2(i + \varepsilon)}{N + 2\varepsilon - 1}, \quad \varepsilon = 0.36.
$$

The default formula places the two polar points at $\cos\theta = \pm(1 - 1/N)$, slightly inside the poles. Roberts' offset pushes them outward fractionally, reducing the maximum nearest-neighbour distance at the poles and improving the overall max/min ratio by around 10–15%. For most applications the basic formula is sufficient, but the change is a single line:

```python
def fibonacci_sphere_improved(n: int, eps: float = 0.36) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fibonacci lattice with Roberts' (2018) pole correction."""
    golden_ratio = (1.0 + np.sqrt(5.0)) / 2.0
    i         = np.arange(n)
    cos_theta = 1.0 - 2.0 * (i + eps) / (n + 2.0 * eps - 1.0)
    theta     = np.arccos(np.clip(cos_theta, -1.0, 1.0))
    phi       = 2.0 * np.pi * i / golden_ratio
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return x, y, z
```

## Conclusions

The key take home messages:

1. **Uniform in angle is wrong**: THis was the main mistake I made, always use the $\arccos(2u-1)$ transformation. Sampling $\theta$ uniformly concentrates four times too many points near the poles.

2. **Correct random sampling still has voids**: The Poisson point process guarantees that roughly 37% of equal-area cells are completely empty at any $N$. This is not fixable by using more points; it is a fundamental property of random sampling.

3. **The Fibonacci lattice eliminates voids**: By using the golden angle, which is the most irrational rotation possible (I think), successive points spread across the sphere with a max/min nearest-neighbour ratio of approximately 1.33 regardless of $N$.

4. **The consequences are real**: For spherical harmonic reconstruction, random sampling causes numerical instability past $l \approx 20$. Fibonacci sampling maintains accuracy to much larger $l$.

5. **The mathematics is beautiful**: The three-distance theorem, the golden ratio, the spiral geometry of sunflower seeds: these are not decorative analogies. The sunflower found the same answer that optimal sphere sampling requires.

## References

[1] W. H. Press, S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery, *Numerical Recipes in C: The Art of Scientific Computing*, 2nd ed., Cambridge University Press, 1992.

[2] M. Roberts, "How to evenly distribute points on a sphere more effectively than the canonical Fibonacci lattice," extremelearning.com.au, 2018.

[3] R. Swinbank and R. J. Purser, "Fibonacci grids: A novel approach to global modelling," *Quarterly Journal of the Royal Meteorological Society* **132**, 1769–1793, 2006.

[4] A. Álvaro González, "Measurement of areas on a sphere using Fibonacci and latitude–longitude lattices," *Mathematical Geosciences* **42**, 49, 2010.

[5] H. Steinhaus, "Sur la division des corps matériels en parties," *Bull. Acad. Polon. Sci.* **4**, 801–804, 1958.
