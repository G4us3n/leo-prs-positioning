# LEO-PRS Positioning Simulation

MATLAB simulation of LEO satellite-to-ground positioning using the 5G-NR
**Positioning Reference Signal (PRS)**. The project models a 12-satellite LEO
constellation over Milan, builds a realistic Ka-band link budget, generates a
3GPP-compliant PRS waveform, and estimates Time-of-Arrival (ToA) ranging
accuracy through Monte Carlo cross-correlation. Results are compared against
the Cramér–Rao Lower Bound (CRLB) and analysed through Dilution-of-Precision
(DOP) geometry.

The error and accuracy framework follows Edjekouane et al., *"User Equivalent
Range Error and Positioning Accuracy Analysis for ToA-Based Techniques Using
PRS and SSB in 5G/6G NTN"*, IEEE OJ-COMS, 2025.

## Requirements

MATLAB R2023b or later with the following toolboxes:

- Satellite Communications Toolbox (orbits, link budget, access)
- 5G Toolbox (PRS waveform generation)
- Phased Array System Toolbox (antenna modelling)
- Mapping Toolbox (geographic visualisation)
- Signal Processing Toolbox (cross-correlation)
- Navigation Toolbox (optional DOP helpers)

## Quick start

From the repository root, in MATLAB:

```matlab
run_simulation
```

This adds `src/` to the path and runs the full pipeline with the default
configuration. The final context struct is returned as `ctx` in the base
workspace.

To customise a run:

```matlab
addpath(genpath('src'));
cfg = config();
cfg.nPlanes   = 6;        % try a different geometry
cfg.run.gif   = true;     % also render the coverage GIF
ctx = main(cfg);
```

## Repository layout

```
leo-prs-positioning/
├── run_simulation.m        Root launcher (sets path, calls main)
├── src/
│   ├── config.m            Single source of truth for all parameters
│   ├── main.m              Pipeline driver: chains the stages in order
│   ├── stages/             One function per simulation section
│   │   ├── stage01_constellation.m   Scenario, orbits, ground stations
│   │   ├── stage02_antenna.m         Phased-array transmit pattern
│   │   ├── stage03_noise.m           Noise floor + atmospheric model
│   │   ├── stage04_heatmaps.m        Geographic coverage snapshots
│   │   ├── stage05_flyover.m         Time-series link evolution
│   │   ├── stage06_dashboard.m       Ground-track map + Milan dashboard
│   │   ├── stage07_prs_waveform.m    5G-NR PRS generation (TS 38.211)
│   │   ├── stage08_epochs.m          Reference-epoch selection
│   │   ├── stage09_gif.m             Optional SNR-coverage animation
│   │   ├── stage10_montecarlo.m      MC ranging with Doppler compensation
│   │   ├── stage11_xcorr.m           Cross-correlation visualisation
│   │   ├── stage12_rmse_sweep.m      RMSE-vs-SNR vs CRLB curve
│   │   ├── stage13_dop.m             Geometric DOP analysis
│   │   └── stage14_viewer.m          Optional interactive 3D viewer
│   └── lib/                Shared helper functions (added as needed)
├── results/                Generated figures, GIFs and data (git-ignored)
├── scripts/                One-off experiment scripts
└── docs/                   Notes, roadmap and references
```

## How the pipeline works

Each stage is a function with the signature

```matlab
ctx = stageNN_name(cfg, ctx)
```

It reads what it needs from the shared context struct `ctx` and the parameter
struct `cfg`, runs one logical section, and returns `ctx` with new fields
added. `main.m` calls the stages in order, and `cfg.run.*` flags let you turn
individual stages on or off (the expensive GIF and the interactive 3D viewer
are off by default).

This makes the parameters you change most often — constellation geometry, PRS
numerology, Monte Carlo settings — all live in one place (`config.m`) instead
of being scattered through a long script.

## Configuration highlights

All of the following are set in `src/config.m`:

- **Constellation geometry** — `nPlanes`, `nSatsPerPlane`, `plane_raan_offsets`,
  `train_nu_offsets`, `altitude_m`, `inclination_deg`.
- **PRS numerology** — `prs.subcarrierSpacing`, `prs.nSizeGrid`,
  `prs.combSize`, `prs.numPRSSymbols`.
- **Monte Carlo** — `mc.nTrials`, `mc.snrSweep_dB`, `mc.accTarget_m`.
- **Stage toggles** — `run.heatmaps`, `run.gif`, `run.monteCarlo`, etc.

## Roadmap

See [`docs/ROADMAP.md`](docs/ROADMAP.md) for the planned analyses and where
each one slots into this structure.

## License

MIT — see [`LICENSE`](LICENSE).
