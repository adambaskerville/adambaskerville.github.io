#!/usr/bin/env python3
"""Generate image assets for the Mandelbrot convergence visualiser."""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "assets" / "interactive" / "mandelbrot_convergence"
WIDTH = 256
HEIGHT = 192
FRAMES = 25
MAX_ITER = 120
TAU = 5.0
TRAINING_STEPS = 4000
SNAPSHOT_EVERY = TRAINING_STEPS // FRAMES
SEED = 0

PALETTE = np.array(
    [
        [7, 11, 33],
        [31, 65, 114],
        [120, 159, 181],
        [240, 203, 104],
        [252, 248, 232],
    ],
    dtype=float,
)


def smooth_escape_target(coords: np.ndarray) -> np.ndarray:
    """Continuous Mandelbrot target used by both the post and visualiser."""
    c = coords[:, 0] + 1j * coords[:, 1]
    z = np.zeros_like(c)
    active = np.ones(c.shape, dtype=bool)
    smooth = np.zeros(c.shape, dtype=np.float64)

    for n in range(MAX_ITER):
        z[active] = z[active] * z[active] + c[active]
        escaped = active & (np.abs(z) > 2.0)
        if np.any(escaped):
            smooth[escaped] = n + 1 - np.log(np.log(np.abs(z[escaped]))) / np.log(2.0)
        active &= ~escaped
        if not np.any(active):
            break

    target = np.exp(-smooth / TAU)
    target[active] = 1.0
    return target.astype(np.float32)


def make_grid(width: int = WIDTH, height: int = HEIGHT) -> tuple[np.ndarray, np.ndarray]:
    xs = np.linspace(-2.0, 1.0, width, dtype=np.float32)
    ys = np.linspace(-1.25, 1.25, height, dtype=np.float32)
    xx, yy = np.meshgrid(xs, ys)
    coords = np.column_stack([xx.ravel(), yy.ravel()]).astype(np.float32)
    return coords, smooth_escape_target(coords).reshape(height, width)


def xavier_uniform(n_in: int, n_out: int, rng: np.random.Generator) -> np.ndarray:
    limit = np.sqrt(6.0 / (n_in + n_out))
    return rng.uniform(-limit, limit, (n_in, n_out)).astype(np.float32)


def sigmoid(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-x))


def fourier_features(x: np.ndarray, B: np.ndarray) -> np.ndarray:
    proj = 2.0 * np.pi * (x @ B)
    return np.concatenate([np.sin(proj), np.cos(proj)], axis=1).astype(np.float32)


def init_mlp(input_dim: int, hidden_dim: int, seed: int) -> list[np.ndarray]:
    rng = np.random.default_rng(seed)
    return [
        xavier_uniform(input_dim, hidden_dim, rng),
        np.zeros((1, hidden_dim), dtype=np.float32),
        xavier_uniform(hidden_dim, hidden_dim, rng),
        np.zeros((1, hidden_dim), dtype=np.float32),
        xavier_uniform(hidden_dim, 1, rng),
        np.zeros((1, 1), dtype=np.float32),
    ]


def forward(params: list[np.ndarray], x: np.ndarray) -> tuple[np.ndarray, tuple[np.ndarray, ...]]:
    W1, b1, W2, b2, W3, b3 = params
    z1 = x @ W1 + b1
    a1 = np.tanh(z1)
    z2 = a1 @ W2 + b2
    a2 = np.tanh(z2)
    z3 = a2 @ W3 + b3
    y = sigmoid(z3)
    return y, (x, a1, a2, y)


def backward(params: list[np.ndarray], cache: tuple[np.ndarray, ...], y_true: np.ndarray) -> list[np.ndarray]:
    W1, b1, W2, b2, W3, b3 = params
    x, a1, a2, y = cache

    dloss = 2.0 * (y - y_true) / len(x)
    dz3 = dloss * y * (1.0 - y)
    gW3 = a2.T @ dz3
    gb3 = dz3.sum(axis=0, keepdims=True)

    da2 = dz3 @ W3.T
    dz2 = da2 * (1.0 - a2 * a2)
    gW2 = a1.T @ dz2
    gb2 = dz2.sum(axis=0, keepdims=True)

    da1 = dz2 @ W2.T
    dz1 = da1 * (1.0 - a1 * a1)
    gW1 = x.T @ dz1
    gb1 = dz1.sum(axis=0, keepdims=True)

    return [gW1, gb1, gW2, gb2, gW3, gb3]


def adam_step(
    params: list[np.ndarray],
    grads: list[np.ndarray],
    moment1: list[np.ndarray],
    moment2: list[np.ndarray],
    step: int,
    lr: float,
) -> None:
    beta1 = 0.9
    beta2 = 0.999
    eps = 1e-8

    for i, (param, grad) in enumerate(zip(params, grads)):
        moment1[i] = beta1 * moment1[i] + (1.0 - beta1) * grad
        moment2[i] = beta2 * moment2[i] + (1.0 - beta2) * (grad * grad)
        mhat = moment1[i] / (1.0 - beta1**step)
        vhat = moment2[i] / (1.0 - beta2**step)
        param -= lr * mhat / (np.sqrt(vhat) + eps)


def interpolate_palette(values: np.ndarray) -> np.ndarray:
    values = np.clip(values, 0.0, 1.0)
    pos = values * (len(PALETTE) - 1)
    lo = np.floor(pos).astype(int)
    hi = np.clip(lo + 1, 0, len(PALETTE) - 1)
    t = (pos - lo)[..., None]
    rgb = PALETTE[lo] * (1.0 - t) + PALETTE[hi] * t
    return rgb.astype(np.uint8)


def write_ppm(path: Path, values: np.ndarray) -> None:
    rgb = interpolate_palette(values)
    with path.open("wb") as handle:
        handle.write(f"P6\n{rgb.shape[1]} {rgb.shape[0]}\n255\n".encode())
        handle.write(rgb.tobytes())


def ppm_to_png(ppm_path: Path, png_path: Path) -> None:
    subprocess.run(["magick", str(ppm_path), str(png_path)], check=True)
    ppm_path.unlink()


def write_frame(path: Path, values: np.ndarray) -> None:
    ppm_path = path.with_suffix(".ppm")
    write_ppm(ppm_path, values)
    ppm_to_png(ppm_path, path)


def train_recorded_model(
    model_name: str,
    train_x: np.ndarray,
    train_y: np.ndarray,
    eval_x: np.ndarray,
    eval_target: np.ndarray,
    hidden_dim: int,
    lr: float,
    seed: int,
) -> list[dict[str, float | str]]:
    params = init_mlp(train_x.shape[1], hidden_dim, seed)
    moment1 = [np.zeros_like(p) for p in params]
    moment2 = [np.zeros_like(p) for p in params]
    rng = np.random.default_rng(SEED + seed)
    records: list[dict[str, float | str]] = []

    for step in range(1, TRAINING_STEPS + 1):
        idx = rng.integers(0, len(train_x), 1024)
        batch_x = train_x[idx]
        batch_y = train_y[idx]
        pred, cache = forward(params, batch_x)
        grads = backward(params, cache, batch_y)
        adam_step(params, grads, moment1, moment2, step, lr)

        if step % SNAPSHOT_EVERY == 0:
            eval_pred, _ = forward(params, eval_x)
            frame = np.clip(eval_pred.reshape(HEIGHT, WIDTH), 0.0, 1.0)
            mae = float(np.mean(np.abs(eval_pred - eval_target)))
            train_mse = float(np.mean((pred - batch_y) ** 2))
            filename = f"{model_name}_{len(records):02d}.png"
            write_frame(OUTPUT_DIR / filename, frame)
            records.append(
                {
                    "src": filename,
                    "step": step,
                    "mae": mae,
                    "loss": train_mse,
                }
            )

    return records


def write_frame_data(vanilla_records: list[dict[str, float | str]], fourier_records: list[dict[str, float | str]]) -> None:
    lines = [
        "window.MANDELBROT_FRAME_DATA = {",
        f"  FRAME_COUNT: {len(vanilla_records)},",
        "  vanillaStats: [],",
        "  fourierStats: [],",
        "};",
        "",
    ]

    for model_name, records in (("vanilla", vanilla_records), ("fourier", fourier_records)):
        for rec in records:
            lines.append(
                "window.MANDELBROT_FRAME_DATA."
                f"{model_name}Stats.push({{src: \"{rec['src']}\", step: {rec['step']}, mae: {rec['mae']:.6f}, loss: {rec['loss']:.6f}}});"
            )

    (OUTPUT_DIR / "frame_data.js").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    eval_coords, target_image = make_grid()
    eval_target = target_image.reshape(-1, 1).astype(np.float32)
    write_frame(OUTPUT_DIR / "target.png", target_image)

    fourier_B = np.random.default_rng(SEED + 7).standard_normal((2, 128), dtype=np.float32) * 8.0

    vanilla_records = train_recorded_model(
        "vanilla",
        eval_coords,
        eval_target,
        eval_coords,
        eval_target,
        hidden_dim=128,
        lr=2e-3,
        seed=11,
    )
    fourier_records = train_recorded_model(
        "fourier",
        fourier_features(eval_coords, fourier_B),
        eval_target,
        fourier_features(eval_coords, fourier_B),
        eval_target,
        hidden_dim=128,
        lr=1e-3,
        seed=12,
    )

    write_frame_data(vanilla_records, fourier_records)


if __name__ == "__main__":
    main()
