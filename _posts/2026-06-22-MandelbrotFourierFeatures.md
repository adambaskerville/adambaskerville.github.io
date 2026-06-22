---
layout: post
title: "T>T: Can a Neural Network Predict the Mandelbrot set?"
date: 2025-09-27
excerpt: "watch me stumble around attempting to train a neural network to reconstruct the Mandelbrot set."
tags: [science, mathematics, programming, python, fractal, mandelbrot, fourier, neural-networks, visualisation]
comments: false
math: true
---

# Try the code yourself!

Click the following button to launch an ipython notebook on Google Colab which implements the code developed in this post:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/adambaskerville/adambaskerville.github.io/blob/master/_posts/MandelbrotFourierFeaturesCode/MandelbrotFourierFeatures.ipynb)

A few years ago I remember watching a YouTube video where it was tested whether a neural network could be trained to predict the structure of the Mandelbrot set. 

I cannot find this video and only have vague memories of it (as Bob Mortimer once said, my memories are like fingerprints on an abandoned handrail) but I thought it would be fun to tackle this problem blind and see how far I can get. Come along for the ride if you want as I know I will be learning something along the way.

My first thought was:

> Surely I can just feed $(x, y)$ coordinates into a neural network and ask it whether the point belongs to the Mandelbrot set.

I thought this was a reasonable starting point.

---

## Where I Started

The object we want to learn is the Mandelbrot set. For each complex number

$$
    c = x + iy
$$

we start with $z_0 = 0$ and iterate

$$
    z_{n+1} = z_n^2 + c.
$$

If the orbit stays bounded, the point is in the set. If it escapes, it is not.

My first instinct was the following machine-learning formulation:

$$
    f(x, y) =
    \begin{cases}
    1 & c \in \mathcal{M} \\
    0 & c \notin \mathcal{M}.
    \end{cases}
$$

So the plan in my head was simple; given coordinates, predict set membership:

1. Generate lots of coordinates.
2. Label them as inside or outside.
3. Train an MLP.
4. Admire the fractal.

---

## First Pitfall: Binary Labels Are a Terrible

The first thing I realised is that the only interesting part of this problem is the boundary, and the boundary is exactly where the target function is discontinuous...

Far away from the set, everything is clearly `0`. Deep inside the large cardioid, everything is clearly `1`. The difficult part, is the infinitely intricate boundary where those two values meet.

That is not an easy target for a smooth optimiser. If I use hard labels, the network gets rewarded for learning the broad interior/exterior split and receives almost no nuanced guidance about how to shape the edge. It just learns "there is some blob over here" and the boundary gets rounded off.

Rather than train on a hard in/out label, I wanted something continuous which still respects the geometry of the problem.

The usual smooth escape-time quantity is

$$
    \mu(c) = n + 1 - \frac{\log(\log |z_n|)}{\log 2},
$$

which is a softened version of "how fast did this point escape?". Here $n$ is the first iteration where the orbit exceeds radius 2. If a point escapes immediately, $\mu(c)$ is small. If it hangs around for many iterations before escaping, $\mu(c)$ is larger. THis is useful because points just outside the Mandelbrot set often take a long time to escape, while points far away escape quickly. So $\mu(c)$ carries geometric information about the boundary that the binary label discard.


Instead of only asking "inside or outside?", we ask "how close to the boundary does this point behave?" where $n$ is the first iteration at which $\|z_n\| > 2$. From that I defined a target

$$
    u(c) =
    \begin{cases}
    1 & \text{if the orbit never escapes within } N_{\max} \text{ iterations} \\
    \exp\!\left(-\mu(c)/\tau\right) & \text{otherwise.}
    \end{cases}
$$

- points inside the set stay near `1`
- points outside decay smoothly away from the boundary
- the optimiser now sees a gradient instead of a sharp cliff edge

That turned the problem from classification into regression, which felt more sensible. Here is the code I ended up putting together:

```python
import numpy as np


def smooth_escape_target(
    coords: np.ndarray,
    max_iter: int = 120,
    tau: float = 5.0,
) -> np.ndarray:
    """
    Continuous Mandelbrot target:
    1 inside the set (up to max_iter), exp(-mu / tau) outside.
    """
    c = coords[:, 0] + 1j * coords[:, 1]
    z = np.zeros_like(c)
    active = np.ones(c.shape, dtype=bool)
    smooth = np.zeros(c.shape, dtype=np.float64)

    for n in range(max_iter):
        z[active] = z[active] * z[active] + c[active]
        escaped = active & (np.abs(z) > 2.0)

        if np.any(escaped):
            smooth[escaped] = (
                n + 1 - np.log(np.log(np.abs(z[escaped]))) / np.log(2.0)
            )

        active &= ~escaped
        if not np.any(active):
            break

    target = np.exp(-smooth / tau)
    target[active] = 1.0
    return target.astype(np.float32)
```

---

## Second Pitfall: Uniform Sampling Is Not Helpful

The next choice I made was to sample points uniformly over the usual viewing window. That also sounds reasonable but there is an issue.

If you sample uniformly over a large rectangle, most points are not interesting as they are not near the Mandelbrot boundary. They are just safely outside the set and escape quickly. So the dataset is dominated by points fwhere there is not much going on. I had accidentally built a dataset that strongly encouraged the network to get very good at the least interesting part of the task.

I did not want to throw away the global sampling entirely because the model still needs to know where the full set sits in the plane. But I also did not want the intricate details to be afterthoughts (more like after-memories?)

So I switched to a mixed sampling strategy:

- a broad global sample for the full geometry
- extra samples from a few boundary-heavy regions

```python
def sample_training_coordinates(
    n_global: int,
    n_focus: int,
    seed: int = 0,
) -> np.ndarray:
    rng = np.random.default_rng(seed)

    # global context
    global_x = rng.uniform(-2.0, 1.0, n_global)
    global_y = rng.uniform(-1.25, 1.25, n_global)

    # deliberately oversample difficult boundary regions
    boxes = [
        (-0.95, -0.55,  0.00, 0.35),
        (-1.80, -1.20, -0.10, 0.10),
        ( 0.20,  0.45,  0.45, 0.75),
    ]

    per_box = n_focus // len(boxes)
    focus_chunks = []

    for x0, x1, y0, y1 in boxes:
        xs = rng.uniform(x0, x1, per_box)
        ys = rng.uniform(y0, y1, per_box)
        focus_chunks.append(np.column_stack([xs, ys]))

    coords = np.vstack([
        np.column_stack([global_x, global_y]),
        *focus_chunks,
    ])

    rng.shuffle(coords)
    return coords.astype(np.float32)
```

This is not profound, it is just me realising that if the hard part of the problem is localised, the data should reflect that.

Turning those coordinates into a PyTorch dataset is then straightforward:

```python
import torch
from torch.utils.data import TensorDataset


def make_dataset(
    n_global: int = 180_000,
    n_focus: int = 60_000,
    max_iter: int = 120,
    tau: float = 5.0,
    seed: int = 0,
) -> TensorDataset:
    coords = sample_training_coordinates(n_global, n_focus, seed=seed)
    target = smooth_escape_target(coords, max_iter=max_iter, tau=tau)

    x = torch.from_numpy(coords)
    y = torch.from_numpy(target[:, None])
    return TensorDataset(x, y)
```

At this stage I felt pleased, the target was smoother, the data was more focused, and everything looked ready for a vanilla MLP.

---

## Third Pitfall: Assuming "Bigger Network" Means "Sharper Fractal"

My instinct at this point was still basically: if the network is underfitting, make the network larger. That is a common reflex of mine, and sometimes it is the right one. Here it was not really the core issue.

If I feed raw coordinates straight into an MLP,

$$
    (x, y) \longrightarrow \text{MLP}(x, y),
$$

the model has to build every bulb, notch, filament, and mini-brot from smooth nonlinear transformations of just two numbers. In principle, yes, it can approximate the target. In practice, what it learns first is the low-frequency envelope.

This is the phenomenon referred to as spectral bias. Standard networks trained with gradient descent tend to learn smooth large-scale structure faster than sharp high-frequency detail:

- the main cardioid is low-frequency
- the large bulbs are still fairly low-frequency
- the filamentary boundary is high-frequency
- the recursively repeated detail is higher-frequency still

So the network was not failing randomly, it was doing exactly what its inductive bias told it to do. It would learn the bean-shaped bulk first and then smear the edge into a blur.

Instead I should have been asking the more useful question:

> Can I change the coordinates so the network does not have to invent all the fine structure from scratch?

That turned out to be the right question (finally).

---

## Fourier Features

Instead of feeding the network raw coordinates, I pass them through a random Fourier map first:

$$
    \gamma(\mathbf{v}) =
    \left[
    \sin(2\pi \mathbf{B}\mathbf{v}),
    \cos(2\pi \mathbf{B}\mathbf{v})
    \right],
    \qquad
    \mathbf{v} = \begin{bmatrix} x \\ y \end{bmatrix}.
$$

Here $\mathbf{B}$ is a random Gaussian matrix, sampled once and then frozen.

This was the shift that made the whole thing click for me. I am no longer asking the network to manufacture high-frequency spatial basis functions out of thin air. I am handing it a coordinate system that already contains oscillations at many scales.

If we feed the network just $(x, y)$, then every pattern in the final image has to be built indirectly by stacking smooth layers on top of those two numbers. The network starts with only a very plain description of position: (horizontal coordinate, vertical coordinate). From that, it has to create tiny repeating bulbs, thin filaments, sharp notches, and all the other fine boundary detail. 

Fourier features change this. Instead of giving the network only $(x, y)$, we first map $(x, y)$ into many sine and cosine functions. Each of those channels already oscillates across the plane. Some oscillate slowly, some rapidly, some along one direction, some along another. The Fourier map gives the network access to functions that are already good at representing high-frequency variation. Since the Mandelbrot boundary contains structure at many scales, that is a much better match than raw coordinates alone.


The scale of $\mathbf{B}$ matters:

- if it is too small, the encoding is still too smooth
- if it is too large, optimisation gets noisy and awkward
- somewhere in the middle, the network suddenly starts to "see" the boundary

Using PyTorch:

```python
import numpy as np
import torch
import torch.nn as nn


class FourierFeatures(nn.Module):
    def __init__(
        self,
        input_dim: int = 2,
        mapping_size: int = 256,
        scale: float = 8.0,
    ) -> None:
        super().__init__()
        B = torch.randn(input_dim, mapping_size) * scale
        self.register_buffer("B", B)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x_proj = 2.0 * np.pi * (x @ self.B)
        return torch.cat([torch.sin(x_proj), torch.cos(x_proj)], dim=-1)
```

There are two details worth noticing:

- `B` is a buffer, not a trainable parameter
- the output doubles in size because we keep both sine and cosine channels

That first point matters. I do not want the model to slowly drift away from the frequency basis I deliberately chose.

---

## Building the Two Models

Once I decided to use Fourier features, I wanted the comparison with the vanilla model to be as clean as possible. So I kept everything else the same and changed only the input encoding.

```python
class FractalNet(nn.Module):
    def __init__(
        self,
        use_fourier: bool = True,
        mapping_size: int = 256,
        fourier_scale: float = 8.0,
        hidden_dim: int = 256,
        depth: int = 4,
    ) -> None:
        super().__init__()

        if use_fourier:
            self.mapping = FourierFeatures(
                input_dim=2,
                mapping_size=mapping_size,
                scale=fourier_scale,
            )
            in_dim = 2 * mapping_size
        else:
            self.mapping = nn.Identity()
            in_dim = 2

        layers = [nn.Linear(in_dim, hidden_dim), nn.GELU()]

        for _ in range(depth - 1):
            layers.extend([nn.Linear(hidden_dim, hidden_dim), nn.GELU()])

        layers.extend([nn.Linear(hidden_dim, 1), nn.Sigmoid()])
        self.mlp = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.mlp(self.mapping(x))
```

That way, if one model does better, I can honestly attribute the improvement to the encoding.

And yes, turning a 2-vector into a 512-vector before the MLP sees it feels a little excessive at first. Those extra channels are carrying specific spatial frequencies, which is exactly what the raw-coordinate model lacks.

---

## Training

The loss function ended up being quite ordinary along with the optimiser and training in general.

```python
from itertools import cycle
from torch.utils.data import DataLoader


def train_model(
    model: nn.Module,
    dataset: TensorDataset,
    steps: int = 12_000,
    batch_size: int = 4096,
    lr: float = 2e-4,
    device: str = "cpu",
) -> list[float]:
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=True, drop_last=True)
    iterator = cycle(loader)

    model.to(device)
    optimiser = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-6)
    loss_fn = nn.MSELoss()
    history = []

    for _ in range(steps):
        x, y = next(iterator)
        x = x.to(device)
        y = y.to(device)

        optimiser.zero_grad(set_to_none=True)
        pred = model(x)
        loss = loss_fn(pred, y)
        loss.backward()
        optimiser.step()

        history.append(loss.item())

    return history
```

The key take home lesson for me from this project was that the bottleneck was not some optimisation trick, it was in the representation.

To render the learned field on a grid:

```python
@torch.no_grad()
def evaluate_grid(
    model: nn.Module,
    width: int = 800,
    height: int = 600,
    device: str = "cpu",
) -> np.ndarray:
    xs = torch.linspace(-2.0, 1.0, width)
    ys = torch.linspace(-1.25, 1.25, height)
    yy, xx = torch.meshgrid(ys, xs, indexing="ij")

    coords = torch.stack([xx.reshape(-1), yy.reshape(-1)], dim=-1).to(device)
    pred = model(coords).reshape(height, width)
    return pred.cpu().numpy()
```

The vanilla model still gets the broad silhouette first, but the Fourier-feature model starts sharpening structures much earlier and much more convincingly, fantastic!

---

## Interactive Visualiser

I wanted a way to show that behaviour interactively, and my first implementation used a spectral surrogate because it was lightweight and easy to make visually clean. But there is a more honest approach, which is the one below: just run the models offline, save their prediction snapshots during training, and let the browser iterate through those saved frames.

So the visualiser below is built from **actual recorded training runs**. The browser is not training anything itself, it is just stepping through snapshots from two real models:

- a vanilla coordinate MLP
- a Fourier-feature MLP

Both were trained against the same continuous Mandelbrot target, and their outputs were saved every few hundred optimisation steps. The true target stays fixed in the background and the current model prediction is overlaid on top with a reveal slider, which makes the comparison much more compact than showing two separate panels.

<style>
  .mandelbrot-visualiser-frame {
    width: 100%;
    height: 1600px;
    border: 1px solid #dbe4ee;
    border-radius: 12px;
    background: #f8fbff;
  }

  @media (max-width: 900px) {
    .mandelbrot-visualiser-frame {
      height: 1180px;
    }
  }
</style>

<iframe
  class="mandelbrot-visualiser-frame"
  src="{{ '/assets/interactive/mandelbrot_fourier_visualiser.html' | relative_url }}"
  title="Mandelbrot Fourier Feature Convergence Visualiser"
  loading="lazy">
</iframe>


---

## What I Learned

I think I had three bad instincts:

- assuming the binary target was "the clean formulation"
- assuming uniformly sampled data was automatically the fairest data
- assuming a blurry result meant "I need a bigger network"

The real fix was more subtle. I needed:

- a target that gave the optimiser usable gradients
- a dataset that did not bury the hard regions
- an input encoding that did not force the model to synthesise every high-frequency structure from raw coordinates alone

---

## Conclusions

I enjoyed this exercise as it scratched an itch I have had since seeing that video, but we also answered the initial question. Yes a neural network can predict the general shape of the Mandebrot set fairly accurately by including Fourier features. Of course the Mandelbrot set is infinitely recursive which our model cannot predict, but it is nice to see a simple ML model can predict a very complex geometric object

