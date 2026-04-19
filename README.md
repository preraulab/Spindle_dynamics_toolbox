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

---

## Initial Clone

To clone the `Spindle_dynamics_toolbox.git` repository to your own computer or remote workstation, start by creating a directory for it. Then, within that directory in your terminal, run the clone command.

#### 1. create an empty folder:
```bash
mkdir Spindle_dynamics_toolbox
```

#### 2. clone the repository to the folder you created:
```bash
git clone git@github.com:preraulab/Spindle_dynamics_toolbox.git Spindle_dynamics_toolbox
```

## Pre-Run Setup and Notes

#### 1. Required MATLAB toolboxes
To run this repository, the following MATLAB toolboxes must be installed:

- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox
- Parallel Computing Toolbox

To install the toolboxes above, in MATLAB: **HOME → Add-Ons → Get Add-Ons**. Detailed instructions [here](https://www.mathworks.com/help/matlab/matlab_env/get-add-ons.html).

#### 2. Possible Security & Privacy Issues
This repository includes `.mex` files to enable fast execution of the multitaper spectrogram. macOS may block execution due to user's security settings. If you see a security warning window pop up like:

> *"multitaper_spectrogram_coder_mex.... not opened ..."*

Follow these steps to allow execution:

1. Go to **System Settings → Privacy & Security**.
2. Scroll down and look for a message indicating that the MEX file was blocked.
3. Click **"Allow Anyway"**
4. **Restart** MATLAB.

Details of the multitaper spectrogram can be found [here](https://github.com/preraulab/multitaper_toolbox/tree/master?tab=readme-ov-file#matlab-implementation).

---

## Table of Contents
* [Background and Toolbox Overview](#background-and-toolbox-overview)
* [Pipeline Overview](#pipeline-overview)
* [Quick Start](#quick-start)
* [Raw Data to Model Fitting](#raw-data-to-model-fitting)
* [Model Results And Visualizations](#model-results-and-visualizations)
* [Key Model Details](#key-model-details)
* [Repository Structure](#repository-structure)
* [Datasets Used in the Paper](#datasets-used-in-the-paper)
* [Citations](#citations)
* [Status](#status)

<br/>

## Background and Toolbox Overview

Sleep spindles are cortical electrical waveforms observed during sleep, considered critical for memory consolidation and sleep stability. Abnormalities in sleep spindles have been found in neuropsychiatric disorders and aging and are suggested to contribute to functional deficits. Numerous studies have demonstrated that spindle activity dynamically and continuously evolves over time and is mediated by a variety of intrinsic and extrinsic factors including sleep stage, slow oscillation (SO) activity (0.5–1.5 Hz), and infraslow activity. Despite these known dynamics, the relative influences on the moment-to-moment likelihood of a spindle event occurring at a specific time are not well-characterized. Moreover, standard analyses almost universally report average spindle rate (known as spindle density) over fixed stages or time periods — thus ignoring timing patterns completely. Without a systematic characterization of spindle dynamics, our ability to identify biomarkers for aging and disordered conditions remains critically limited.

Using a rigorous statistical framework, **Chen et al. 2025** demonstrate that short-term timing patterns are the dominant determinant of spindle timing, whereas sleep depth, cortical up/down-state, and long-term (infraslow) patterns — features thought to be primary drivers of spindle occurrence — are less important. They also show that these short-term timing patterns are fingerprint-like and show increased variability over age. This study provides a new lens on spindle production mechanisms, which will allow studies of the role of spindle timing patterns in memory consolidation, aging, and disease.

<h3>Toolbox Overview</h3>
<p>
    This code toolbox implements an integrated framework for characterizing <strong>instantaneous spindle temporal dynamics</strong> influenced by multiple factors, including sleep stage/depth (SO Power), cortical up/down states (SO Phase), and the past history of spindles. It consists of <strong>two major components:</strong>
</p>

<h4>1. Quick Start GUI</h4>
<ul>
    <li>
        Provides an interactive interface for users to explore data and create basic visualizations:
        <ul>
            <li>Load example data or upload their own data.</li>
            <li>Specify model factors and add interactions between them.</li>
            <li>Customize settings to generate an overview figure.</li>
        </ul>
    </li>
</ul>

<h4>2. Step-by-Step Example Script</h4>
<ul>
    <li>
        Offers researchers a detailed workflow to:
        <ul>
            <li>Dive into the toolbox's functionalities.</li>
            <li>Address specific scientific questions using example data or their own data.</li>
        </ul>
    </li>
</ul>

<p>Designed as a user-friendly toolbox, we aim to:</p>
<ul>
    <li>Streamline research workflows.</li>
    <li>Ensure reproducibility.</li>
    <li>Inspire researchers to apply and adapt it to diverse datasets for broader applications.</li>
</ul>

## Pipeline Overview

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

## Quick Start
After installing the package, execute the `quick_start` GUI function in the MATLAB command line:

```matlab
 > quick_start;
```

The following GUI should be generated (shown below, left panel):

<p align="center">
<img src="https://github.com/preraulab/Spindle_dynamics_toolbox/blob/master/image_folder/GUI2.png" width="700" />
</p>
<p align="center">
  <b> Quick Start GUI </b>
</p>

Follow the 4 steps to make your choices. Here is an example of the settings (shown above, right panel), which specifies the model primarily used in the paper:

<ol>
    <li>Load data: Load <a href="https://github.com/preraulab/Spindle_dynamics_toolbox/tree/master/example_data" target="_blank">Example Data</a>, which is the same <a href="https://www.sleepdata.org/datasets/mesa" target="_blank">MESA</a> subject as in Figure 6 from the paper.</li>
    <li>Select factors: Click on "Stage", "SOphase", "History".</li>
    <li>Choose history: Click on "Long-term (90 sec)"</li>
    <li>Select interactions: Choose "stage:SOphase interaction"</li>
</ol>

Once you click "Run the Model", an overview figure will be generated for exploration:

<p id="overview-figure" align="center">
<img src="https://github.com/preraulab/Spindle_dynamics_toolbox/blob/master/image_folder/overview_fig.png" width="900" />
</p>
<p align="center">
  <b> Overview of the spindle dynamics </b>
</p>

Use the `scrollzoompan` interface to slide around the figure. To view the exact same epoch as the paper, set **Zoom** to **90** and **Pan** to **30145**.

#### To load your own data
Click on "Load User Data" to browse and load your own data. Ensure the file meets the following requirements:

1. **File format**:
   - Must be a `.mat` file.

2. **Required variables** (variable names are case-insensitive):
   - **`EEG`**: (double, 1D vector) raw EEG data. Other accepted names: `eeg_data`, `raw_EEG`, `data`
   - **`Fs`**: (double, scalar) Sampling frequency of the EEG data in Hz. Other accepted names: `sampling_rate`
   - **`stage_val`**: (double, 1D vector) Sleep stages (`1: N3, 2: N2, 3: N1, 4: REM, 5: Wake`). Other accepted names: `stage_vals`, `stages`, `stage`, `stageval`, `stagevals`
   - **`stage_time`**: (double, 1D vector) Corresponding times for the sleep stages (seconds). Other accepted names: `stage_t`, `time_stages`, `stage_times`

## Raw Data to Model Fitting
Walk-through of the data loading, preprocessing, model specification, and model fitting steps in the example script. To run everything and generate all figures with a single command:

```matlab
 > Example_Script;
```

### Model Specification
`specify_mdl` consolidates all model specifications — model factors (`BinarySelect`), interactions (`InteractSelect`), and other options — into a single struct (`ModelSpec`).

```matlab
[ModelSpec] = specify_mdl(BinarySelect, InteractSelect, 'hist_choice', hist_choice, 'control_pt', control_pt);
```

```matlab
% Input:
%       <Required inputs>
%       - BinarySelect: (1x4 vector, double), indicates which factors are selected
%         Factors in fixed order: SOphase, stage, SOpower, history
%         e.g., [1,1,0,1] means select SOphase, stage, and history
%       - InteractSelect: (1xn cell), each entry is an interaction term in the form of A:B
%         Case- and order-insensitive; accepts separators ':', '&', '-'
%         e.g., {'stage:SOphase', 'stage:history'}
%
%       <Optional inputs>
%       - hist_choice: (string) 'short' (default, up to 15 s) or 'long' (up to 90 s)
%       - control_pt:  (1xk vector, double), spline control point locations
%                      default: [0:15:90 120 150]
%       - binsize:     (double), point-process bin size in seconds (default: 0.1)
%       - hard_cutoffs: (1x2 vector, double), frequency cutoff in Hz (default: [12 16])
%
% Output:
%       ModelSpec: struct containing all model specifications
```

### Data Preprocessing
`preprocessToDesignMatrix` extracts spindle event information (`res_table`), saves preprocessed binned data (`BinData`), and builds the design matrix `X`:

```matlab
[X, BinData, res_table] = preprocessToDesignMatrix(EEG, Fs, stage_val, stage_time, ModelSpec);
```

```matlab
% Input (all required):
%       - EEG:        (double, 1D) raw EEG data
%       - Fs:         (double, scalar) sampling frequency in Hz
%       - stage_val:  (double, 1D) stage values where 1–5 = N3, N2, N1, REM, Wake
%       - stage_time: (double, 1D) corresponding time of the stages (s)
%       - ModelSpec:  struct (from specify_mdl)
%
% Output:
%       - X:         (double, matrix) design matrix
%       - BinData:   (struct) data saved per binsize
%       - res_table: (table) event info and signals. Key fields include:
%           peak_ctimes : detected event central times (s)
%           peak_freqs  : detected event frequencies (Hz)
%           SOpow       : slow oscillation power
%           SOphase     : slow oscillation phase
```

### Model Fitting
MATLAB's `glmfit` fits the point-process GLM:

```matlab
[b, dev, stats] = glmfit(X, y, 'poisson');
```

- `X` — design matrix
- `y` — point-process binary train (response)
- `b` — fitted parameters
- `dev` — deviance
- `stats` — MATLAB struct with coefficient covariance, p-values, residuals, etc.

## Model Results And Visualizations

### 1. History Modulation Curve
The history modulation curve estimates a multiplicative modulation of the spindle event rate due to a prior event at any given time lag, answering: *How much more likely is a spindle now, given an event was observed X seconds ago?* `plot_hist_curve.m` computes and plots it, along with history features (refractory period, excitatory period, peak time, peak height, infraslow multiplier).

```matlab
[xlag, yhat, yu, yl, hist_features] = plot_hist_curve(stats, ModelSpec, BinData);
```

Output fields:
- `xlag`, `yhat` — history lag (s) and modulation value
- `yu`, `yl` — 95% CI bounds
- `hist_features.ref_period` — refractory period (s)
- `hist_features.exc_period` — excited period (s)
- `hist_features.p_time` — peak time (s)
- `hist_features.p_height` — peak height
- `hist_features.AUC_is` — infraslow multiplier (area under 40–70 s / 30)

<p align="center">
<img src="https://github.com/preraulab/Spindle_dynamics_toolbox/blob/master/image_folder/hist_fig1.png" width="800" />
</p>
<p align="center">
  <b>Figure 1: History Modulation Curve </b>
</p>

We observe that history modulation begins with a refractory period, ramps up to a peak within an excitatory period, and gradually decreases back to 1 (no modulation for events long ago).

### 2. Spindle Preferential SO Phase Shifts With Sleep Depth
Sleep spindles have been widely reported to preferentially occur in the cortical up-state. `plot_stage_prefphase.m` visualizes how the preferred phase shifts with **discrete sleep stage**:

```matlab
[pp, pp_CI] = plot_stage_prefphase(b, stats, ModelSpec);
```

- `pp` — preferred phase (rad) per stage
- `pp_CI` — 95% CI per stage

`plot_sop_prefphase.m` visualizes how the preferred phase shifts with **continuous SO power**:

```matlab
[phi0, sop0, sop_pp_mat] = plot_sop_prefphase(b, stats, ModelSpec);
```

- `phi0` — phase eval points (500 evenly from −π to π)
- `sop0` — SOP eval points (510 evenly from −5 to 30)
- `sop_pp_mat` — spindle-rate heatmap over (phi0, sop)

<p align="center">
<img src="https://github.com/preraulab/Spindle_dynamics_toolbox/blob/master/image_folder/ppshift_fig2.png" width="900" />
</p>
<p align="center">
  <b>Figure 2: Phase Shift: Stage vs. SO Power </b>
</p>

Spindles tend to occur at the SO peak (phase ~0) during light sleep, and the preferred phase shifts earlier to the SO rising phase in deeper sleep for this participant.

### 3. Model With History Greatly Improves Model Performance
If the model is correct, the time-rescaling theorem maps event times into a homogeneous Poisson process. After rescaling, Kolmogorov-Smirnov (KS) plots compare the inter-spindle-interval distribution to the model prediction. A well-fit model's KS plot follows the 45° line within its significance bounds (black). Curves outside these bounds (red) indicate lack-of-fit.

```matlab
[ks, ksT] = KSplot(CIF, y, ploton);
```

- `CIF` — conditional intensity function
- `y` — binary event train
- `ploton` — if 1, plot
- `ks` — KS statistic
- `ksT` — 0 = passes KS test, 1 = fails

<p align="center">
<img src="https://github.com/preraulab/Spindle_dynamics_toolbox/blob/master/image_folder/ks_fig3.png" width="900" />
</p>
<p align="center">
  <b>Figure 3: KS plot for models with different components </b>
</p>

The stage-only or phase-only models fail the KS test. A single-history-component model passes it, and the full all-factors model achieves an even smaller KS statistic.

### 4. Short-Term History Contributes the Most to Statistical Deviance
The framework allows quantitative comparison of relative contributions via deviance analysis — the point-process equivalent of ANOVA in linear regression. `compute_dev_exp.m` computes proportional deviance explained by each factor:

```matlab
[dev_exp_sta, dev_exp_sop] = compute_dev_exp(BinData);
```

- `dev_exp_sta` — deviance explained from `[stage, phase, history]`
- `dev_exp_sop` — deviance explained from `[SOP,   phase, history]`

<p align="center">
<img src="https://github.com/preraulab/Spindle_dynamics_toolbox/blob/master/image_folder/dev_exp_fig4.png" width="400" />
</p>
<p align="center">
  <b>Figure 4: Proportional Deviance Explained By Different Factors </b>
</p>

The history component explains the most deviance, surpassing sleep stage and SO phase.

## Key Model Details

### Point-process GLM framework

Spindles are treated as a marked point process. The conditional intensity function λ(t | History, Covariates) gives the instantaneous spindle rate given past spindle times and covariate values. For each participant:

```
log λ(t) = β₀ + β_stage · stage(t) + β_SOP · SOP(t) + β_SOφ · f(SO_phase(t)) + β_hist · h(history(t))
```

where `f(SO_phase)` and `h(history)` are basis expansions (spline / indicator functions).

### History effect interpretation

The history function `h(·)` reveals each participant's characteristic timing structure:

- **Refractory period** — `exp(β_hist · h(τ)) < 1` at lag τ (suppressed rate)
- **Excitatory period** — `exp(β_hist · h(τ)) > 1` at lag τ (elevated rate)
- **Infraslow clustering** — extending the history window to ~90 s reveals a secondary bump around ~55 s in some participants

### Goodness-of-fit

`KSplot.m` produces the time-rescaled KS diagnostic. After time-rescaling under the fitted model, inter-spindle intervals should be uniformly distributed. Deviations from the KS confidence bands indicate model misspecification.

### Model comparison

`compute_dev_exp.m` compares nested models (e.g., stage-only vs stage+SOP vs stage+SOP+SO_phase vs stage+SOP+SO_phase+history). The paper shows that adding history produces by far the largest deviance-explained improvement.

## Repository Structure
The toolbox is organized as follows, with key functions:

```
.REPO_ROOT
├── quick_start.m          Quick-start GUI for overview figure exploration
├── Example_Script.m       Main script, from raw data to model results and visualizations
│                          Uses example data in example_data/
├── example_data/          Example data
├── image_folder/          Saved figures
└── helper_function/       Helper functions and vendored toolboxes
    ├── major_function/
    │       specify_mdl.m                Model specifications
    │       preprocessToDesignMatrix.m   Data preprocessing
    │       plot_hist_curve.m            History modulation curve
    │       plot_stage_prefphase.m       SO preferred phase as a function of sleep stage
    │       plot_sop_prefphase.m         SO preferred phase as a function of SO power
    │       KSplot.m                     KS statistics, plot, and test
    │       compute_dev_exp.m            Proportional deviance explained
    │
    ├── TFsigma_peak_detector/           Spindle detection algorithm (TF-peak)
    ├── multitaper/                      Multitaper spectrogram toolbox (vendored)
    └── utils/                           Various utility functions
```

## Datasets Used in the Paper

1. **Wamsley et al.** — 17 healthy controls (ages 26–45), two-night recordings
2. **MESA** (Multi-Ethnic Study of Atherosclerosis) — 1,008 participants, single-night recordings
   - Middle-aged group: 433 participants (ages 54–65)
   - Older group: 575 participants (ages 66–94)

Both datasets are external. The `example_data/` folder provides a small single-subject demo for getting started.

## Documentation

Full API reference: **https://preraulab.github.io/Spindle_dynamics_toolbox/**

## Citations
The code in this repository is companion to the paper:

> Shuqiang Chen, Mingjian He, Ritchie E. Brown, Uri T. Eden, Michael J. Prerau.
> "Individualized temporal patterns drive human sleep spindle timing."
> *Proceedings of the National Academy of Sciences*, 122(2): e2405276121 (2025).
> doi: [10.1073/pnas.2405276121](https://doi.org/10.1073/pnas.2405276121)

which should be cited for all use.

<br/>

## License

BSD 3-Clause. See [`LICENSE`](LICENSE).

## Status

All implementations are functional but subject to refinement. Next optimization: deal with perfect predictors (e.g., Wake and REM) to save runtime.
<br/>
Last updated by Shuqiang Chen, 02/11/2025
<br/>
For questions or suggestions please email Shuqiang Chen at shuqiang@bu.edu
