"""
Interactive comparison of Hartree-Fock (HF) vs Fully Correlated (FC) conditional
electron density for helium, as a function of the first electron's position r₁.

WHAT THIS SCRIPT DOES
─────────────────────
Given two pre-optimised wavefunction files (.sv format), it evaluates and
visualises ρ(r₂ | r₁): the probability density of finding electron 2 at
position r₂, given that electron 1 is fixed at r₁.

In HF, the two-electron wavefunction factorises — Ψ_HF(r₁,r₂) = ψ(r₁)ψ(r₂) —
so the r₁ dependence cancels identically. ρ_HF is the same regardless of where
electron 1 sits. The FC wavefunction contains explicit r₁₂ dependence, so
ρ_FC develops a depletion (correlation hole) centred on electron 1.

WHY VISUALISE THIS
──────────────────
The correlation hole is the entire reason correlation energy exists. Seeing it
as a physical dip in density — moving as electron 1 moves — makes the concept
concrete in a way equations alone cannot.

OUTPUT
──────
Three self-contained HTML files in assets/interactive/:
  orbital_density_2d.html        2D heatmap with HF/FC toggle + r₁ slider
  orbital_density_3d.html        3D surface with HF/FC toggle + r₁ slider
  orbital_density_visualiser.html  wrapper iframe with 2D/3D tab switcher

Usage:
    cd _posts/OrbitalElectronDensityCode
    python density_visualizer.py
"""

import re
import time
import numpy as np
import plotly.graph_objects as go


# ── File paths ────────────────────────────────────────────────────────────────
# The .sv files are Maple save files produced by optimisation runs.
# HF file: restricted HF orbital for He, 25 Laguerre basis functions.
# FC file: fully correlated wavefunction, 5456 three-index Laguerre terms,
#          expanded in perimetric coordinates (u, v, w).

HF_FILE = "hf_wavefunctions/e_e_2_25_RHF_OPT_SM1.sv"
FC_FILE = "fc_wavefunctions/e_e_infinity_2856_A_syevd_BOBYQA_Z2_R1_D.sv"


# ── Grid parameters ───────────────────────────────────────────────────────────
# GRID_RANGE: physical extent of the visualisation plane in bohr.
#   2 bohr captures the bulk of the helium density — the peak sits near 0.5 bohr.
#   Going further wastes resolution on near-zero density.
#
# GRID_N: points per axis → 100×100 = 10 000 evaluation points per density.
#   Fine enough to show the correlation hole clearly; coarse enough to compute
#   all r₁ positions in a few seconds on a laptop.
#
# R1_POSITIONS: the discrete r₁ values pre-computed and stored as separate
#   Plotly traces. The slider switches between them via visibility toggling,
#   so no runtime evaluation happens in the browser.
#   Spacing is denser near the nucleus where the density changes fastest.

GRID_RANGE   = 2
GRID_N       = 100
R1_POSITIONS = [0.0, 0.3, 0.6, 1.0, 1.5, 2.0, 2.5, 3.0]   # bohr


# ── Parsers ───────────────────────────────────────────────────────────────────
# Both parsers read plain text out of Maple .sv (save) files.
# Regex is used rather than a full Maple parser because the relevant data
# appears in a predictable line format — full parsing would be overkill.


def parse_hf(path):
    """
    Parse HF .sv file. Returns (A, coeffs) where:
      A        : the orbital exponent (controls radial decay rate)
      coeffs   : array of length max_n+1; coeffs[k] multiplies L_k(A·r)

    The HF orbital has the form:
        ψ(r) = exp(-A·r/2) · Σₖ cₖ · L_k(A·r)

    where L_k are Laguerre polynomials. The exponential prefactor is extracted
    separately at evaluation time (see hf_density_2d) so the coefficients here
    are purely the polynomial weights.

    File structure: the orbital is stored in a block headed 'F1:=' and
    terminated by '+0:', with individual terms like '(c)*L(n,u)'.
    """
    with open(path) as fh:
        content = fh.read()

    # Exponent A appears as a single line: "A:=5.236076:"
    A = float(re.search(r"\nA:=([\d.eE+\-]+):", content).group(1))

    # Extract the F1 block between 'F1:=' and '+0:' (the terminator Maple appends)
    f1_str = re.search(r"F1:=(.*?)\+0:", content, re.DOTALL).group(1)

    # Each term is written as '(coefficient)*L(order,u)'
    terms = re.findall(r"\(([-+]?[\d.eE+\-]+)\)\*L\((\d+),u\)", f1_str)

    max_n = max(int(n) for _, n in terms)
    coeffs = np.zeros(max_n + 1)
    for c, n in terms:
        coeffs[int(n)] = float(c)

    return A, coeffs


def parse_fc(path):
    """
    Parse FC .sv file. Returns (A, B, C, coeffs, iidx, jidx, kidx) where:
      A, B, C  : exponents for the three perimetric coordinates u, v, w
      coeffs   : array of term weights, shape (N_terms,)
      iidx, jidx, kidx : Laguerre orders for each term, shape (N_terms,)

    Each term in the wavefunction is:
        cᵢⱼₖ · L_i(u) · L_j(v) · L_k(w)

    and the full wavefunction is:
        Ψ = exp(-A(r₁+r₂)) · Σᵢⱼₖ cᵢⱼₖ · L_i(u) · L_j(v) · L_k(w)

    Storing the three index arrays separately (iidx, jidx, kidx) rather than
    as a list of tuples lets us index the precomputed Laguerre tables with
    numpy integer arrays in the inner loop — avoiding Python-level element access.

    Note: the FC file terminator is '\\n0:\\n' (a bare '0:' on its own line),
    not '+0:' as in the HF file — a quirk of how Maple saves the two-variable
    expansion differently from the one-variable one.
    """
    with open(path) as fh:
        content = fh.read()

    A = float(re.search(r"\nA:=([\d.eE+\-]+):", content).group(1))
    B = float(re.search(r"\nB:=([\d.eE+\-]+):", content).group(1))
    C = float(re.search(r"\nC:=([\d.eE+\-]+):", content).group(1))

    # The F block holds all 5456 terms between 'F:=\n' and '\n0:\n'
    f_str = re.search(r"F:=\n(.*?)\n0:\n", content, re.DOTALL).group(1)
    pattern = r"\(([-+]?[\d.eE+\-]+)\)\*L\((\d+),u\)\*L\((\d+),v\)\*L\((\d+),w\)"
    raw = re.findall(pattern, f_str)

    # Separate arrays per index type: lets numpy index T_u[iidx[t]] in one shot
    coeffs = np.array([float(c) for c, _, _, _ in raw], dtype=np.float64)
    iidx   = np.array([int(i)   for _, i, _, _ in raw], dtype=np.int32)
    jidx   = np.array([int(j)   for _, _, j, _ in raw], dtype=np.int32)
    kidx   = np.array([int(k)   for _, _, _, k in raw], dtype=np.int32)

    return A, B, C, coeffs, iidx, jidx, kidx


# ── Laguerre table (vectorised recurrence) ────────────────────────────────────


def laguerre_table(max_n, x):
    """
    Build table T where T[n, p] = L_n(x[p]) for n = 0 … max_n.

    WHY A TABLE
    ───────────
    The FC wavefunction sums over ~5000 terms, each needing L_i(u), L_j(v),
    L_k(w) evaluated at every grid point. Naively that would recompute the same
    Laguerre values many times for repeated indices. Building the full table once
    per coordinate array costs O(max_n × N) and then every term lookup is O(N)
    with no recomputation.

    WHY RECURRENCE
    ──────────────
    The three-term Laguerre recurrence:
        L_0(x) = 1
        L_1(x) = 1 - x
        L_n(x) = ((2n-1-x)·L_{n-1} - (n-1)·L_{n-2}) / n

    is numerically stable and exactly equivalent to the closed-form definition.
    It avoids computing factorials or summing series, and maps cleanly to numpy
    operations on the full 1-D array x at once — no Python loop over grid points.

    WHY ROWS ARE POLYNOMIAL ORDER
    ──────────────────────────────
    Storing T[n, :] (order in rows, grid in columns) means T[n] is already the
    contiguous array needed for the inner loop. Cache-friendly for the access
    pattern T_u[iidx[t]], T_v[jidx[t]], T_w[kidx[t]].
    """
    N = len(x)
    T = np.empty((max_n + 1, N), dtype=np.float64)
    T[0] = 1.0
    if max_n >= 1:
        T[1] = 1.0 - x
    for n in range(2, max_n + 1):
        T[n] = ((2 * n - 1 - x) * T[n - 1] - (n - 1) * T[n - 2]) / n
    return T


# ── Density evaluators ────────────────────────────────────────────────────────


def hf_density_2d(x_grid, z_grid, A_hf, coeffs_hf):
    """
    Evaluate the HF conditional density on the 2D (x₂, z₂) slice.

    PHYSICS
    ───────
    Ψ_HF(r₁,r₂) = ψ(r₁)ψ(r₂).  The conditional density is:
        ρ_HF(r₂ | r₁) = |Ψ|² / ∫|Ψ|²dr₂ = |ψ(r₂)|²

    r₁ cancels exactly — the density is purely a function of r₂.
    This is the mathematical expression of HF's lack of correlation:
    each electron is independent of the other.

    RADIAL VOLUME FACTOR
    ────────────────────
    We display 4π·r₂²·|ψ(r₂)|² rather than |ψ(r₂)|² alone. This is the
    radial density — the probability of finding the electron in a shell
    of radius r₂, integrated over all angles. It makes the peak location
    visually obvious (near 0.5 bohr for He) and matches what chemists
    call the radial distribution function.
    """
    r2 = np.sqrt(x_grid**2 + z_grid**2).ravel()
    u  = A_hf * r2

    max_n = len(coeffs_hf) - 1
    T_u   = laguerre_table(max_n, u)
    psi   = coeffs_hf @ T_u          # dot product over basis functions
    psi  *= np.exp(-0.5 * u)         # attach exponential prefactor

    return (4 * np.pi * r2**2 * psi**2).reshape(x_grid.shape)


def fc_density_2d(r1_val, x_grid, z_grid, A, B, C, coeffs, iidx, jidx, kidx):
    """
    Evaluate the FC conditional density on the 2D (x₂, z₂) slice with
    electron 1 fixed on the z-axis at distance r1_val from the nucleus.

    PERIMETRIC COORDINATES
    ──────────────────────
    Rather than Cartesian or spherical coordinates in 6D, the FC wavefunction
    is expressed in perimetric coordinates (u, v, w):
        u = A·(r₁₂ + r₂ - r₁)
        v = B·(r₁₂ + r₁ - r₂)
        w = C·(r₁  + r₂ - r₁₂)

    where r₁₂ = |r₁ - r₂|. By the triangle inequality all three are
    non-negative for any physical configuration, so there is no domain
    boundary issue and Laguerre polynomials (defined on [0,∞)) are a
    natural basis.

    For helium's ground state A = B (symmetric under electron interchange)
    and C = 2A.

    VALIDITY MASK
    ─────────────
    Floating-point arithmetic can produce tiny negative values near the
    triangular boundary (e.g. r₁₂ ≈ r₁ + r₂ up to rounding). We allow a
    tolerance of 1e-9 and clip to zero rather than discarding the point —
    the physical density there is genuinely zero anyway and this avoids
    gaps in the heatmap near the nucleus.

    INNER LOOP
    ──────────
    The sum over 5456 terms is a Python loop — the bottleneck. Each
    iteration uses prebuilt Laguerre tables indexed by the integer arrays
    iidx/jidx/kidx, so the per-term work is three O(N) numpy multiplies.
    Moving this to Cython or numba would give ~100× speedup; at N=10 000
    and 8 r₁ positions the current ~0.8s total is fast enough.
    """
    shape = x_grid.shape
    r2  = np.sqrt(x_grid**2 + z_grid**2).ravel()
    # electron 1 is on the z-axis: x₁=0, z₁=r1_val
    r12 = np.sqrt(x_grid**2 + (z_grid - r1_val)**2).ravel()

    u_arr = A * (r12 + r2  - r1_val)
    v_arr = B * (r12 + r1_val - r2)
    w_arr = C * (r1_val + r2  - r12)

    # Allow 1e-9 tolerance for floating-point near the triangular boundary
    valid  = (u_arr >= -1e-9) & (v_arr >= -1e-9) & (w_arr >= -1e-9)
    u_safe = np.where(valid, np.clip(u_arr, 0.0, None), 0.0)
    v_safe = np.where(valid, np.clip(v_arr, 0.0, None), 0.0)
    w_safe = np.where(valid, np.clip(w_arr, 0.0, None), 0.0)

    # Build Laguerre tables once per coordinate — reused by all terms
    T_u = laguerre_table(int(iidx.max()), u_safe)
    T_v = laguerre_table(int(jidx.max()), v_safe)
    T_w = laguerre_table(int(kidx.max()), w_safe)

    psi = np.zeros(len(r2), dtype=np.float64)
    for t in range(len(coeffs)):
        psi += coeffs[t] * T_u[iidx[t]] * T_v[jidx[t]] * T_w[kidx[t]]

    psi *= np.exp(-A * (r1_val + r2))
    psi  = np.where(valid, psi, 0.0)   # zero out invalid (unphysical) points

    return (4 * np.pi * r2**2 * psi**2).reshape(shape)


# ── Normalisation ─────────────────────────────────────────────────────────────


def normalise(Z, eps=1e-30):
    """
    Scale density to [0, 1] for display.

    WHY NORMALISE
    ─────────────
    HF and FC densities have different absolute magnitudes and different peak
    heights at each r₁ position. Normalising each independently to [0, 1]
    lets the colorscale always span the full range — showing the shape and
    relative structure clearly — without one plot washing out the other.

    The eps guard prevents division by zero for the r₁=0 HF case where the
    density at the origin can be numerically very small.
    """
    zmax = Z.max()
    return Z / (zmax if zmax > eps else 1.0)


# ── Custom HTML writers ───────────────────────────────────────────────────────
# WHY CUSTOM HTML RATHER THAN PLOTLY UPDATEMENUS
# ───────────────────────────────────────────────
# Plotly's built-in updatemenus/animate system has a known reliability issue:
# when toggling between traces via restyle inside updatemenus, the visibility
# list is not applied consistently — both buttons end up showing the same data.
# The root cause is that Plotly's internal animation state machine conflicts
# with updatemenus restyle when the figure contains many pre-baked traces.
#
# The fix: bypass updatemenus entirely. We write plain HTML <button> elements
# and drive Plotly.restyle() directly from JavaScript onclick handlers. This
# calls the same underlying API but outside Plotly's animation scheduler, so
# the trace visibility is applied immediately and deterministically.
#
# TRACE LAYOUT STRATEGY
# ─────────────────────
# All density arrays are pre-computed at each r₁ position and stored as
# separate Plotly traces (all hidden except the active one). The slider and
# buttons manipulate only the `visible` property via restyle — no data moves
# between browser and Python at interaction time. This makes all transitions
# instant regardless of grid size.
#
# Alternative considered: update a single trace's `z` data on each click.
# Rejected because embedding z arrays in button/slider args bloats the HTML
# proportionally to NR × grid² — manageable but wasteful compared to the
# visibility-only approach which uses the same data but stored once per trace.


def write_2d_html(path, x_vals, z_vals, hf_density, fc_densities, r1_positions):
    """
    Write a self-contained 2D heatmap HTML file.

    TRACE LAYOUT (total: 1 + 2·NR traces)
    ───────────────────────────────────────
      0          : HF heatmap                 — visible on load
      1 … NR     : FC heatmaps, one per r₁   — all hidden on load
      NR+1 … 2NR : e₁ scatter markers        — all hidden on load

    Having one marker trace per r₁ position means the slider only needs
    to toggle visibility, not update coordinates. The marker traces are
    tiny (one point each) so the overhead is negligible.

    WHY SHOW e₁ MARKER IN HF MODE
    ───────────────────────────────
    In HF the density is literally independent of r₁, so it makes no
    difference where the marker sits — but showing it anyway makes the
    independence visible: user drags slider, marker moves, density does
    not. That contrast is the whole pedagogical point.
    """
    import json
    import os
    NR         = len(r1_positions)
    hf_norm    = normalise(hf_density)
    colorscale = "Plasma"
    tick_vals  = np.linspace(-GRID_RANGE, GRID_RANGE, 5)

    traces = []

    # Trace 0: HF heatmap — the single fixed density, visible at start
    traces.append(go.Heatmap(
        x=x_vals.tolist(), y=z_vals.tolist(), z=hf_norm.tolist(),
        colorscale=colorscale, zmin=0, zmax=1, visible=True,
        showscale=True,
        colorbar=dict(x=1.01, len=0.9, title=dict(text="ρ (norm.)", side="right")),
    ))

    # Traces 1…NR: one FC heatmap per r₁ position, all hidden
    for ri in range(NR):
        traces.append(go.Heatmap(
            x=x_vals.tolist(), y=z_vals.tolist(),
            z=normalise(fc_densities[ri]).tolist(),
            colorscale=colorscale, zmin=0, zmax=1, visible=False,
            showscale=True,
            colorbar=dict(x=1.01, len=0.9, title=dict(text="ρ (norm.)", side="right")),
        ))

    # Traces NR+1…2NR: one e₁ marker per r₁ position, all hidden
    # Separate traces (rather than one trace with updated coords) so visibility
    # toggling via restyle is the only JavaScript operation needed.
    for ri, r1 in enumerate(r1_positions):
        traces.append(go.Scatter(
            x=[0.0], y=[float(r1)], mode="markers", visible=False,
            marker=dict(symbol="x", size=14, color="white",
                        line=dict(width=2, color="black")),
            showlegend=False,
        ))

    axis_style = dict(
        range=[-GRID_RANGE, GRID_RANGE],
        tickvals=tick_vals.tolist(),
        ticktext=[f"{v:.1f}" for v in tick_vals],
        showgrid=False, zeroline=False,
    )

    # Nucleus: a gold circle shape anchored to the data axes so it scales
    # correctly if the user zooms. Drawn as a layout shape rather than a trace
    # because it never needs to be hidden or updated.
    _nb = dict(type="circle", x0=-0.08, y0=-0.08, x1=0.08, y1=0.08,
               fillcolor="gold", line_color="black", line_width=2, opacity=0.9)

    fig = go.Figure(data=traces)
    fig.update_layout(
        title=dict(
            text=("<b>ρ(r₂ | r₁) — conditional electron density (helium)</b><br>"
                  "<sup>HF never changes · FC hole tracks electron 1 (white ✕) · "
                  "gold circle = nucleus</sup>"),
            x=0.5, xanchor="center", font=dict(size=14),
        ),
        height=520,
        paper_bgcolor="#0f0f0f", plot_bgcolor="#0f0f0f",
        font=dict(color="white"),
        shapes=[dict(**_nb, xref="x", yref="y")],
        legend=dict(visible=False),
        margin=dict(t=90, b=60, l=60, r=80),
        # scaleanchor ensures x and z axes have equal physical spacing —
        # critical so the correlation hole appears as a circle, not an ellipse
        xaxis=dict(title="x₂ (bohr)", scaleanchor="y", scaleratio=1, **axis_style),
        yaxis=dict(title="z₂ (bohr)", **axis_style),
    )

    # Serialise the entire figure to JSON once. The HTML page loads this as a
    # JS variable and passes it to Plotly.newPlot — no server required.
    fig_json      = fig.to_json()
    slider_labels = json.dumps([f"{r:.1f}" for r in r1_positions])
    NR_js         = NR

    # ── HTML template ────────────────────────────────────────────────────────
    # Controls are plain HTML above the plot div. The slider is a native
    # <input type="range"> — simpler and more reliable than Plotly's built-in
    # slider widget, which has its own state-machine issues.
    #
    # JS VISIBILITY LOGIC
    # ───────────────────
    # visForHF(ri) and visForFC(ri) build a boolean array over all traces.
    # applyVis() calls Plotly.restyle with that array, toggling exactly the
    # right traces in one API call. State is two variables: `mode` and `ri`.
    # No animation frames, no updatemenus, no Plotly internal scheduler.
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>
    * {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{ background: #0f0f0f; font-family: -apple-system, sans-serif; color: white; }}
    #controls {{
      display: flex; align-items: center; gap: 16px;
      padding: 10px 20px; background: #161616;
      border-bottom: 1px solid #2a2a2a; flex-wrap: wrap;
    }}
    .mode-btn {{
      background: #252525; color: #ccc;
      border: 1px solid #555; border-radius: 5px;
      padding: 6px 20px; font-size: 13px; cursor: pointer;
    }}
    .mode-btn.active {{
      background: #4a4fc4; border-color: #7b7ff0; color: #fff; font-weight: 600;
    }}
    #slider-wrap {{ display: flex; align-items: center; gap: 10px; flex: 1; min-width: 200px; }}
    #slider-wrap label {{ font-size: 13px; white-space: nowrap; }}
    #r1-slider {{ flex: 1; accent-color: #7b7ff0; }}
    #plot {{ width: 100%; }}
  </style>
</head>
<body>
  <div id="controls">
    <button class="mode-btn active" id="btn-hf" onclick="setMode('hf')">Hartree-Fock</button>
    <button class="mode-btn"       id="btn-fc" onclick="setMode('fc')">Fully correlated</button>
    <div id="slider-wrap">
      <label id="slider-label">r₁ = 0.0 bohr</label>
      <input type="range" id="r1-slider" min="0" max="{NR_js-1}" step="1" value="0"
             oninput="onSlider(this.value)">
    </div>
  </div>
  <div id="plot"></div>
  <script>
    var figData = {fig_json};
    var NR = {NR_js};
    var sliderLabels = {slider_labels};
    var mode = 'hf';   // current wavefunction mode
    var ri   = 0;      // current r1 slider index

    Plotly.newPlot('plot', figData.data, figData.layout, {{responsive: true}});

    // Build a visibility array for all (1 + 2*NR) traces.
    // HF mode: show trace 0 (HF heatmap) + marker at current ri.
    // The marker is shown even in HF mode so the user can see that moving
    // electron 1 has no effect on the HF density — the key pedagogical point.
    function visForHF(ri) {{
      var v = new Array(1 + 2*NR).fill(false);
      v[0]          = true;   // HF heatmap
      v[1 + NR + ri] = true;  // e₁ marker
      return v;
    }}

    // FC mode: show FC heatmap at ri + matching e₁ marker.
    function visForFC(ri) {{
      var v = new Array(1 + 2*NR).fill(false);
      v[1 + ri]      = true;  // FC heatmap for this r₁
      v[1 + NR + ri] = true;  // e₁ marker at same r₁
      return v;
    }}

    function applyVis() {{
      var v = (mode === 'hf') ? visForHF(ri) : visForFC(ri);
      Plotly.restyle('plot', {{'visible': v}});
    }}

    function setMode(m) {{
      mode = m;
      document.getElementById('btn-hf').classList.toggle('active', m === 'hf');
      document.getElementById('btn-fc').classList.toggle('active', m === 'fc');
      applyVis();
    }}

    function onSlider(val) {{
      ri = parseInt(val);
      document.getElementById('slider-label').textContent =
        'r₁ = ' + sliderLabels[ri] + ' bohr';
      applyVis();  // update both modes — marker moves in HF too
    }}
  </script>
</body>
</html>"""

    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(html)


def write_3d_html(path, x_vals, z_vals, hf_density, fc_densities, r1_positions):
    """
    Write a self-contained 3D surface HTML file.

    Same design philosophy as write_2d_html. Differences:
    - go.Surface instead of go.Heatmap; z-axis encodes density height.
    - A separate always-visible nucleus Scatter3d trace (trace 1) that
      cannot be a layout shape in a 3D scene.
    - Correlation hole is a literal dip in the surface, not just a colour
      change — makes the concept immediately visceral.

    TRACE LAYOUT (total: 2 + 2·NR traces)
    ───────────────────────────────────────
      0          : HF surface               — visible on load
      1          : nucleus marker           — always visible
      2 … NR+1   : FC surfaces, one per r₁ — all hidden on load
      NR+2 … 2NR+1 : e₁ markers            — all hidden on load

    ASPECT RATIO
    ────────────
    z:0.55 compresses the density axis relative to x/y. Without this the
    surface is very tall and the correlation hole hard to see. The flatter
    aspect ratio also improves the default camera angle.
    """
    import json
    import os
    NR         = len(r1_positions)
    hf_norm    = normalise(hf_density)
    colorscale = "Plasma"

    def mk_surface(z_data, visible):
        return go.Surface(
            x=x_vals.tolist(), y=z_vals.tolist(), z=z_data.tolist(),
            colorscale=colorscale, cmin=0, cmax=1, visible=visible,
            showscale=True,
            colorbar=dict(x=1.01, len=0.80, thickness=16,
                          title=dict(text="ρ (norm.)", side="right")),
            # Lighting tuned to make the correlation hole visible as a shadow
            # even when viewed from shallow angles
            lighting=dict(ambient=0.65, diffuse=0.85, specular=0.2, roughness=0.6),
            lightposition=dict(x=1, y=-1, z=3),
        )

    # Trace 0: HF surface (visible on load); trace 1: nucleus (always visible)
    traces = [
        mk_surface(hf_norm, True),
        go.Scatter3d(x=[0.0], y=[0.0], z=[0.02], mode="markers", visible=True,
                     marker=dict(size=7, color="gold", symbol="circle",
                                 line=dict(width=1, color="black")),
                     showlegend=False),
    ]

    # Traces 2…NR+1: FC surfaces, one per r₁
    for ri in range(NR):
        traces.append(mk_surface(normalise(fc_densities[ri]), False))

    # Traces NR+2…2NR+1: e₁ markers elevated to z=1.12 (just above surface peak)
    # so they float visibly above the surface regardless of viewing angle
    for ri, r1 in enumerate(r1_positions):
        traces.append(go.Scatter3d(
            x=[0.0], y=[float(r1)], z=[1.12],
            mode="markers+text", visible=False,
            marker=dict(size=7, color="white", symbol="x",
                        line=dict(width=2, color="black")),
            text=["e₁"], textposition="top center",
            textfont=dict(size=11, color="white"),
            showlegend=False,
        ))

    fig = go.Figure(data=traces)
    fig.update_layout(
        title=dict(
            text=("<b>ρ(r₂ | r₁) — 3D conditional electron density (helium)</b><br>"
                  "<sup>Drag to rotate · FC hole is a literal dip · "
                  "white ✕ = electron 1 · gold = nucleus</sup>"),
            x=0.5, xanchor="center", font=dict(size=14),
        ),
        height=540,
        paper_bgcolor="#0f0f0f",
        font=dict(color="white"),
        scene=dict(
            xaxis=dict(title="x₂ (bohr)", range=[-GRID_RANGE, GRID_RANGE],
                       backgroundcolor="#1a1a1a", gridcolor="#333",
                       showgrid=True, zeroline=False),
            yaxis=dict(title="z₂ (bohr)", range=[-GRID_RANGE, GRID_RANGE],
                       backgroundcolor="#1a1a1a", gridcolor="#333",
                       showgrid=True, zeroline=False),
            zaxis=dict(title="ρ (norm.)", range=[0, 1.2],
                       backgroundcolor="#1a1a1a", gridcolor="#333",
                       showgrid=True, zeroline=False),
            bgcolor="#111111",
            camera=dict(eye=dict(x=1.6, y=-1.9, z=1.2), up=dict(x=0, y=0, z=1)),
            aspectmode="manual", aspectratio=dict(x=1, y=1, z=0.55),
        ),
        legend=dict(visible=False),
        margin=dict(t=90, b=20, l=10, r=90),
    )

    fig_json      = fig.to_json()
    slider_labels = json.dumps([f"{r:.1f}" for r in r1_positions])
    NR_js         = NR

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>
    * {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{ background: #0f0f0f; font-family: -apple-system, sans-serif; color: white; }}
    #controls {{
      display: flex; align-items: center; gap: 16px;
      padding: 10px 20px; background: #161616;
      border-bottom: 1px solid #2a2a2a; flex-wrap: wrap;
    }}
    .mode-btn {{
      background: #252525; color: #ccc;
      border: 1px solid #555; border-radius: 5px;
      padding: 6px 20px; font-size: 13px; cursor: pointer;
    }}
    .mode-btn.active {{
      background: #4a4fc4; border-color: #7b7ff0; color: #fff; font-weight: 600;
    }}
    #slider-wrap {{ display: flex; align-items: center; gap: 10px; flex: 1; min-width: 200px; }}
    #slider-wrap label {{ font-size: 13px; white-space: nowrap; }}
    #r1-slider {{ flex: 1; accent-color: #7b7ff0; }}
    #plot {{ width: 100%; }}
  </style>
</head>
<body>
  <div id="controls">
    <button class="mode-btn active" id="btn-hf" onclick="setMode('hf')">Hartree-Fock</button>
    <button class="mode-btn"       id="btn-fc" onclick="setMode('fc')">Fully correlated</button>
    <div id="slider-wrap">
      <label id="slider-label">r₁ = 0.0 bohr</label>
      <input type="range" id="r1-slider" min="0" max="{NR_js-1}" step="1" value="0"
             oninput="onSlider(this.value)">
    </div>
  </div>
  <div id="plot"></div>
  <script>
    var figData = {fig_json};
    var NR = {NR_js};
    var sliderLabels = {slider_labels};
    var mode = 'hf';
    var ri   = 0;

    Plotly.newPlot('plot', figData.data, figData.layout, {{responsive: true}});

    // Trace indices:
    //   0        = HF surface
    //   1        = nucleus (always on)
    //   2..NR+1  = FC surfaces
    //   NR+2..2NR+1 = e₁ markers
    function visForHF(ri) {{
      var v = new Array(2 + 2*NR).fill(false);
      v[0] = true; v[1] = true;   // HF surface + nucleus
      v[2 + NR + ri] = true;      // e₁ marker moves even in HF mode
      return v;
    }}
    function visForFC(ri) {{
      var v = new Array(2 + 2*NR).fill(false);
      v[1] = true;              // nucleus always
      v[2 + ri] = true;         // FC surface for this r₁
      v[2 + NR + ri] = true;    // matching e₁ marker
      return v;
    }}

    function applyVis() {{
      var v = (mode === 'hf') ? visForHF(ri) : visForFC(ri);
      Plotly.restyle('plot', {{'visible': v}});
    }}

    function setMode(m) {{
      mode = m;
      document.getElementById('btn-hf').classList.toggle('active', m === 'hf');
      document.getElementById('btn-fc').classList.toggle('active', m === 'fc');
      applyVis();
    }}

    function onSlider(val) {{
      ri = parseInt(val);
      document.getElementById('slider-label').textContent =
        'r₁ = ' + sliderLabels[ri] + ' bohr';
      applyVis();
    }}
  </script>
</body>
</html>"""

    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(html)


# ── Combined HTML output ──────────────────────────────────────────────────────


def write_all_html(x_vals, z_vals, hf_density, fc_densities, r1_positions, output_dir):
    """
    Write all three HTML files to output_dir.

    WHY THREE FILES + WRAPPER
    ─────────────────────────
    The 2D and 3D files are self-contained and large (each embeds all density
    arrays as JSON). Keeping them separate means the browser only loads one
    at a time. The wrapper uses iframes and lazy-loads the 3D iframe — its
    WebGL canvas only initialises when the div is actually visible, avoiding
    a zero-size canvas bug that prevents 3D surfaces from rendering.
    """
    import os
    os.makedirs(output_dir, exist_ok=True)

    path_2d   = os.path.join(output_dir, "orbital_density_2d.html")
    path_3d   = os.path.join(output_dir, "orbital_density_3d.html")
    path_main = os.path.join(output_dir, "orbital_density_visualiser.html")

    write_2d_html(path_2d, x_vals, z_vals, hf_density, fc_densities, r1_positions)
    print(f"  2D → {path_2d}")
    write_3d_html(path_3d, x_vals, z_vals, hf_density, fc_densities, r1_positions)
    print(f"  3D → {path_3d}")

    # Wrapper: two iframes, 3D src set to '' until first click so WebGL
    # initialises in a visible container rather than a display:none one
    wrapper = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>HF vs FC Electron Density — Helium</title>
  <style>
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body { background: #0f0f0f; font-family: -apple-system, sans-serif; }
    #toggle-bar {
      display: flex; align-items: center; justify-content: center;
      gap: 8px; padding: 10px 0; background: #161616;
      border-bottom: 1px solid #2a2a2a;
    }
    #toggle-bar span { font-size: 13px; color: #888; margin-right: 4px; }
    .view-btn {
      background: #252525; color: #ccc;
      border: 1px solid #444; border-radius: 5px;
      padding: 6px 18px; font-size: 13px; cursor: pointer;
    }
    .view-btn.active {
      background: #4a4fc4; border-color: #7b7ff0; color: #fff; font-weight: 600;
    }
    iframe { width: 100%; height: 640px; border: none; display: block; }
    #frame-3d { display: none; }
  </style>
</head>
<body>
  <div id="toggle-bar">
    <span>View:</span>
    <button class="view-btn active" id="btn-2d" onclick="showView('2d')">2D Heatmap</button>
    <button class="view-btn"        id="btn-3d" onclick="showView('3d')">3D Surface</button>
  </div>
  <iframe id="frame-2d" src="orbital_density_2d.html"></iframe>
  <iframe id="frame-3d" src=""></iframe>
  <script>
    function showView(mode) {
      var f2d = document.getElementById('frame-2d');
      var f3d = document.getElementById('frame-3d');
      f2d.style.display = mode === '2d' ? 'block' : 'none';
      f3d.style.display = mode === '3d' ? 'block' : 'none';
      document.getElementById('btn-2d').classList.toggle('active', mode === '2d');
      document.getElementById('btn-3d').classList.toggle('active', mode === '3d');
      if (mode === '3d' && !f3d.src) f3d.src = 'orbital_density_3d.html';
    }
  </script>
</body>
</html>"""

    with open(path_main, "w", encoding="utf-8") as fh:
        fh.write(wrapper)
    print(f"Wrapper → {path_main}")


# ── Main ──────────────────────────────────────────────────────────────────────


def main():
    """
    Pipeline:
      1. Parse wavefunction files → exponents + coefficient arrays
      2. Build (x, z) grid
      3. Evaluate HF density once (r₁-independent, computed once)
      4. Evaluate FC density for each r₁ position in R1_POSITIONS
      5. Write HTML output files
    """
    print("Parsing HF wavefunction …", flush=True)
    A_hf, coeffs_hf = parse_hf(HF_FILE)
    print(f"  A = {A_hf:.6f},  {len(coeffs_hf)} basis functions")

    print("Parsing FC wavefunction …", flush=True)
    t0 = time.perf_counter()
    A_fc, B_fc, C_fc, coeffs_fc, iidx, jidx, kidx = parse_fc(FC_FILE)
    print(f"  A = {A_fc:.6f},  B = {B_fc:.6f},  C = {C_fc:.6f}")
    print(f"  {len(coeffs_fc)} terms  ({time.perf_counter()-t0:.1f}s)")

    # 2D grid over the (x₂, z₂) plane. meshgrid produces (GRID_N, GRID_N) arrays
    # that are passed directly to the vectorised density evaluators.
    x_vals = np.linspace(-GRID_RANGE, GRID_RANGE, GRID_N)
    z_vals = np.linspace(-GRID_RANGE, GRID_RANGE, GRID_N)
    X, Z = np.meshgrid(x_vals, z_vals)

    # HF: computed once — identical for all r₁
    print("\nComputing HF density …", flush=True)
    t0 = time.perf_counter()
    hf_dens = hf_density_2d(X, Z, A_hf, coeffs_hf)
    print(f"  Done in {time.perf_counter()-t0:.2f}s")

    # FC: computed once per r₁ position and cached — the bottleneck is the
    # 5456-term inner loop (~0.1s per position on a modern laptop)
    print(f"\nComputing FC densities ({len(R1_POSITIONS)} r₁ positions) …", flush=True)
    fc_dens_list = []
    for r1_val in R1_POSITIONS:
        t0 = time.perf_counter()
        dens = fc_density_2d(r1_val, X, Z, A_fc, B_fc, C_fc,
                             coeffs_fc, iidx, jidx, kidx)
        print(f"  r₁={r1_val:.1f} → {time.perf_counter()-t0:.1f}s", flush=True)
        fc_dens_list.append(dens)

    out_dir = "../../assets/interactive"
    print("\nWriting HTML files …", flush=True)
    write_all_html(x_vals, z_vals, hf_dens, fc_dens_list, R1_POSITIONS, out_dir)


if __name__ == "__main__":
    main()
