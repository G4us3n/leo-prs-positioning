# Roadmap

This file maps the planned analyses onto the modular pipeline so each new
task has a clear home. Tasks are listed in a sensible implementation order.

## 1. Multi-numerology overlays (Section 12 / RMSE plots)

Overlay a second, analog numerology on the existing RMSE-vs-SNR plots.

- **Where:** `stage12_rmse_sweep.m`, driven by a list of numerologies in
  `config.m`.
- **Plan:** turn the single PRS configuration into a small array of configs
  (e.g. the current FR1 30 kHz plus an analog one), run the sweep for each, and
  plot them on the same axes.
- **Style convention to follow:** blue = current numerology, orange/yellow =
  analog; solid = Doppler compensation, dashed = CRLB, dotted = no Doppler
  compensation.

## 2. SNR-after-FFT verification (Section 6 / 7)

Compute SNR after the FFT, matching the reference paper's formula, and compare
with the current per-sample method.

- **Where:** a new helper in `src/lib/` (e.g. `snr_after_fft.m`) called from
  the PRS / ranging stages.
- **Plan:** implement the post-FFT SNR definition, compute it alongside the
  existing per-sample SNR, and print/plot the two for consistency.

## 3. Multi-satellite selection (Section 8)

Use all satellites in visibility instead of only the best-SNR link.

- **Where:** `stage08_epochs.m` (keep all visible satellites per city, not just
  the max-SNR one) and `stage10_montecarlo.m` (range every visible satellite).
- **Plan:** replace the `[~, bestSat] = max(...)` selection with the full list
  of visible satellites and produce one ranging measurement per satellite.

## 4. Least-squares positioning over the trajectory

Estimate the user position at each time step (e.g. −200 to +100 s) by solving a
simple least-squares problem from all ~12 ranging measurements.

- **Where:** a new stage, e.g. `stage15_least_squares.m`, plus a solver helper
  in `src/lib/` (`ls_position.m`).
- **Plan:** at each epoch, build the geometry/design matrix from the visible
  satellites, run a least-squares (or weighted least-squares using the UERE
  covariance) update, and record the estimated position.

## 5. Positioning-accuracy evolution plot

Graph the positioning accuracy metric over the user trajectory.

- **Where:** same new stage as task 4, with a dedicated plotting helper.
- **Plan:** plot horizontal/3D error (or R95) versus time to see how accuracy
  evolves as the geometry changes during the flyover.

## 6. Alternative constellation geometries

Test geometries beyond the baseline 4×3 "train".

- **Where:** `config.m` only — no code changes needed.
- **Plan:** vary `plane_raan_offsets` (wider transversal/left-right coverage),
  `train_nu_offsets` (larger along-track separation), `nPlanes` and
  `nSatsPerPlane`, and compare DOP and positioning accuracy across at least two
  configurations.

## Conventions

- Keep new tunables in `config.m`, never hard-coded inside stages.
- New stages follow the `ctx = stageNN_name(cfg, ctx)` signature and are added
  to `main.m` behind a `cfg.run.*` flag.
- Reusable maths (solvers, SNR formulas) go in `src/lib/`.
- Generated artefacts go in `results/` (git-ignored).
