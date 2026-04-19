# Spindle_dynamics_toolbox

MATLAB toolbox accompanying **Chen et al., *PNAS* 2025** on the temporal dynamics of human sleep spindles. It lets you detect spindles as time-frequency peaks, then fit a point-process generalized-linear-model (GLM) that jointly captures multiple factors influencing moment-to-moment spindle occurrence:

- **Sleep stage** (discrete N2 / N3)
- **Slow-oscillation (SO) power** — a continuous measure of sleep depth
- **Slow-oscillation (SO) phase** — cortical up/down state
- **Spindle history** — the timing of previous spindles (short-term and infraslow)

The key finding of the paper: **short-term timing patterns are the dominant determinant of spindle timing**, more so than sleep depth or cortical up/down state. Each subject has a near-fingerprint-like refractory + excitatory structure (typical peaks: ~1.8 s refractory period, ~2.9 s excitatory period, peak at ~3.5 s post-spindle), with variability increasing over age.

## Please cite the paper when using this toolbox

> Shuqiang Chen, Mingjian He, Ritchie E. Brown, Uri T. Eden, Michael J. Prerau.
> "Individualized temporal patterns drive human sleep spindle timing."
> *Proceedings of the National Academy of Sciences*, 122(2): e2405276121 (2025).
> doi: [10.1073/pnas.2405276121](https://doi.org/10.1073/pnas.2405276121)

A machine-readable citation is in [`CITATION.cff`](CITATION.cff) — GitHub's "Cite this repository" button uses it.

## What's in the toolbox

| Component | What it does |
|---|---|
| `Example_Script.m` | End-to-end example: loads example data, detects spindles, fits GLM, plots the history effect and phase preferences |
| `quick_start.m` | Minimal runnable template for applying the pipeline to your own data |
| `helper_function/TFsigma_peak_detector/` | Time-frequency peak detection in the sigma band (12–16 Hz) — spindle detector |
| `helper_function/multitaper/` | Multitaper spectrogram (shipped, MEX-accelerated) |
| `helper_function/major_function/` | GLM deviance / model-comparison / plotting routines (KS plots, history-curve plot, phase-preference plots) |
| `helper_function/utils/` | Design-matrix builders, spectral power extractors, colorbar / scale helpers |

## Pipeline overview

```
       EEG (polysomnography)
                │
                ▼
  ┌──────────────────────────────┐
  │  Multitaper spectrogram      │  multitaper_spectrogram (MEX)
  └──────────────┬───────────────┘
                 ▼
  ┌──────────────────────────────┐
  │  TF-sigma peak detection     │  extract_*.m, find_frequency_peaks
  │  12–16 Hz peaks → spindles   │
  └──────────────┬───────────────┘
                 ▼
  ┌──────────────────────────────┐
  │  Covariates                  │
  │    • sleep stage             │  (discrete)
  │    • SO power                │  (continuous; compute_mtspect_power)
  │    • SO phase                │  (0 – 2π; cortical up/down state)
  │    • spindle history         │  (past spindle times; build_design_mt)
  └──────────────┬───────────────┘
                 ▼
  ┌──────────────────────────────┐
  │  Point-process GLM           │  conditional-intensity regression
  │  λ(t | covariates, history)  │  (MATLAB glmfit / fitglm)
  └──────────────┬───────────────┘
                 ▼
  ┌──────────────────────────────┐
  │  Analysis & visualization    │
  │    • KS plot (goodness-of-fit)│  KSplot
  │    • history curve            │  plot_hist_curve
  │    • SO-phase preference      │  plot_sop_prefphase
  │    • stage × phase preference │  plot_stage_prefphase
  │    • deviance explained       │  compute_dev_exp
  └──────────────────────────────┘
```

## Install

```bash
git clone https://github.com/preraulab/Spindle_dynamics_toolbox.git
cd Spindle_dynamics_toolbox
```

In MATLAB:

```matlab
addpath(genpath('/path/to/Spindle_dynamics_toolbox'));
```

## Required MATLAB toolboxes

- **Signal Processing Toolbox**
- **Statistics and Machine Learning Toolbox** (for `glmfit` / `fitglm`)
- **Parallel Computing Toolbox** (optional, for speed in the TF-peak detector)

## Running the demo

```matlab
Example_Script
```

Runs the full pipeline on example data included under `example_data/`. Produces the figures from the paper (spindle train, conditional-intensity estimates per model variant, history curve, SO-phase preference plots).

For your own data: start from `quick_start.m` and substitute your EEG + scoring.

## Security warning on macOS

The toolbox ships pre-built MEX files for the multitaper spectrogram. On modern macOS, Apple's Gatekeeper may block execution:

> *"multitaper_spectrogram_coder_mex.... not opened ..."*

Fix:

1. **System Settings → Privacy & Security**
2. Scroll to the blocked-MEX notice
3. Click **Allow Anyway**
4. Restart MATLAB

See the upstream [multitaper_toolbox](https://github.com/preraulab/multitaper_toolbox) for details on the multitaper implementation.

## Key model details

### Point-process GLM framework

Spindles are treated as a marked point process. The conditional intensity function λ(t | History, Covariates) gives the instantaneous spindle rate given past spindle times and covariate values. For each participant, we fit:

```
log λ(t) = β₀ + β_stage · stage(t) + β_SOP · SOP(t) + β_SOφ · f(SO_phase(t)) + β_hist · h(history(t))
```

where `f(SO_phase)` and `h(history)` are basis expansions (spline / indicator functions).

### History effect interpretation

The history function `h(·)` reveals each participant's characteristic timing structure:

- **Refractory period** — regions where `exp(β_hist · h(τ)) < 1` (suppressed spindle rate at lag τ relative to baseline)
- **Excitatory period** — regions where `exp(β_hist · h(τ)) > 1` (elevated spindle rate at lag τ)
- **Infraslow clustering** — extending the history window to ~90 s reveals a secondary bump around ~55 s in some participants

### Goodness-of-fit

Use the **Kolmogorov-Smirnov plot** (`KSplot.m`) — after time-rescaling under the fitted model, inter-spindle intervals should be uniformly distributed. Deviations from the KS confidence bands indicate model misspecification.

### Model comparison

Use `compute_dev_exp.m` to compare nested models (e.g., stage-only vs stage+SOP vs stage+SOP+SO_phase vs stage+SOP+SO_phase+history). The paper demonstrates that adding the history term produces by far the largest deviance-explained improvement.

## Datasets used in the paper

1. **Wamsley et al.** — 17 healthy controls (ages 26–45), two-night recordings
2. **MESA** (Multi-Ethnic Study of Atherosclerosis) — 1,008 participants, single-night recordings
   - Middle-aged group: 433 participants (ages 54–65)
   - Older group: 575 participants (ages 66–94)

Both datasets are external — the `example_data/` folder provides a small single-subject demo for getting started.

## Documentation

Full API reference: **https://preraulab.github.io/Spindle_dynamics_toolbox/**

## License

BSD 3-Clause. See [`LICENSE`](LICENSE).
