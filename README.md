# CMJ bilateral analysis for ACL research

A local Python / Streamlit app for reviewing normalized countermovement jump exports, averaging two or three trials, comparing limbs, identifying reviewed movement phases, and exporting multi-subject metrics and selected publication-quality graphs.

## Run in Visual Studio Code

Install **Python 3.11 or newer** and the VS Code Python extension. Open this project folder in VS Code and use **Terminal → New Terminal**.

Windows:

```powershell
py -m venv .venv
.\.venv\Scripts\python -m pip install -r requirements.txt
.\.venv\Scripts\python -m streamlit run app.py
```

macOS / Linux:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m streamlit run app.py
```

Streamlit opens the app in your browser. File selection and downloads happen through the browser. The Python process runs locally; no external analysis API is used. Stop the app with Ctrl+C. In VS Code, select `.venv` with **Python: Select Interpreter**. You do not need to edit file paths in the source code.

## Workflow

1. Select one or more subject `.txt` / `.tsv` files. Give each file a unique subject ID and optionally record the operated limb.
2. Review the editable column mapping. For your confirmed export structure, the app selects exactly the first three complete 15-column blocks (columns 1–45) and preconfirms this mapping. Other layouts default to the first three occurrences of each limb/signal/axis as trials 1–3 and require review. Extra occurrences are excluded. Select two or three valid trials, or additional trials when available. A trial ID must refer to the same jump in every channel. Excluding one trial should exclude it across all channels. Do not average unknown EMG processing variants together.
3. Confirm mapping, source units, synchronization and coordinate conventions. Your knee axes are labeled X = flexion, Y = varus, Z = rotation. Limb synchronization is preconfirmed for your exports. Anatomical sign consistency should still be checked. An optional right-knee sign inversion is available per axis and is recorded in the manifest; do not flip signs merely to reduce observed asymmetry.
4. Inspect plots with individual trials, mean ± sample SD and Right − Left curves. Bilateral comparison requires a matching channel on both sides.
5. Review phase events separately for each trial. Enter onset, bottom/deepest position, takeoff, landing and end of landing absorption, in percent. Confirm a row only when all events are present and ordered. Partial recordings can still produce whole-cycle metrics without phase metrics. Average-curve markers show mean reviewed event positions.
6. Optionally enable exploratory subject-level parametric SPM and/or group paired SPM. See the statistical limitations below.
7. Select the graphs and PNG/SVG/PDF formats you want, click **Prepare download**, then download the ZIP. CSV tables open in Excel, R, Python or other statistical tools.

All 101 samples are retained at 0–100%. No time warping or resampling is performed. Trial averages give equal weight to each selected trial. Group inference gives equal weight to each subject's average, regardless of the number of trials.

## Input format and the provided examples

The app reads the exact two-row tab-separated format in the supplied examples:

```text
        Left Knee Angles    Right Knee Angles    ... repeated trial channels ...
ITEM    X                   X                   ...
1       4.031282            8.258310             ...
...
101     ...                 ...                 ...
```

Actual separators must be tabs. Header 1 identifies the channel; header 2 identifies the axis. Names beginning `Left ` / `Right ` or `L_` / `R_` are recognized. There must be exactly 101 finite numeric samples, indexed 1–101 or 0–100. Missing values, irregular row widths and partial curves are rejected instead of silently interpolated. UTF-8 text is supported.

Inspection of the supplied files found:

- `SUBID1 ACL.txt` and `SUBID2 ACL.txt` are byte-for-byte identical. They cannot be used as evidence from two independent subjects. The GUI flags duplicates and excludes subsequent copies from group inference; descriptive exports remain possible.
- Each file has 98 signal columns: three complete 15-channel blocks, followed by a partial block and additional EMG columns with no trial or processing labels. You confirmed the first 45 columns as the three CMJ trials; this recognized layout is selected automatically. All 98 columns remain available in the mapping editor.
- Bilateral knee angles X/Y/Z, hip angle X, vertical GRF Z, knee power X and gluteus medius EMG are available in those first three blocks. Only right gluteus maximus is present; the app does not invent its left counterpart.
- The selected GRF traces start at zero. The algorithm therefore refuses to claim a complete force-derived CMJ sequence from these samples. Confirm the export's actual event boundaries and whether it contains the full movement before entering phase events.
- You confirmed GRF in newtons, EMG acquired in millivolts then normalized to MVC, synchronized left/right jump events, and knee X/Y/Z = flexion/varus/rotation. These are now saved defaults. Normalized EMG is dimensionless (MVC ratio or %MVC); the scale remains unspecified until you select it, with no numeric rescaling. Knee-power units, angle units and actual trial durations remain unconfirmed. Force magnitudes in the samples are unusually small for whole-body GRF in newtons; verify export scaling before interpreting absolute forces. No conversion is applied automatically.

The private subject files are not bundled with the source package.

## Phase definitions and limits

Automatic candidates use the **sum of left and right GRF Z** and a configurable threshold relative to its peak. A candidate needs initial sustained contact, at least three consecutive below-threshold samples, and sustained contact after flight. Takeoff is the first below-threshold sample and landing the first return above threshold. These are sampled threshold estimates, not interpolated physical event times. Negative or incorrectly oriented force requires correction at source.

When a candidate flight exists, the maximum mean left/right knee X angle before takeoff is proposed as a **bottom-position proxy**, assuming positive flexion. It is not center-of-mass minimum or zero COM velocity. Onset and the end of landing absorption require manual input. Every phase remains unconfirmed until reviewed.

The exported phases are descent, ascent, flight and landing absorption. Force-derived unweighting/braking/propulsion require calibrated force, body weight and actual time or COM velocity; they are not inferred here. Percent-normalized curves also do not provide valid jump height, seconds-based duration, impulse, mechanical work or integrated EMG in seconds. If normalization uses a nonlinear joint-angle coordinate rather than linear elapsed-time normalization, even normalized areas describe that coordinate, not a time average. Both limbs must use the same coordinate before comparison.

## Exported results

| File | Contents |
|---|---|
| `all_subjects_trial_metrics.csv` | Per-trial, per-channel whole-cycle and reviewed-phase metrics |
| `all_subjects_summary_metrics.csv` | Mean, sample SD and contributing count of each trial-level metric |
| `all_subjects_average_waveforms.csv` | 101-point mean and sample SD across selected trials |
| `all_subjects_limb_differences.csv` | Pointwise Right − Left, signed symmetry and waveform RMSE |
| `inference_results.csv`, `spm_clusters.csv` | Test summaries, errors, and parametric cluster endpoints/p-values |
| `analysis_manifest.json` | Source hashes, selected trials, operated limb, units, phase events, sign changes, inference settings and SPM cluster results |
| Subject `column_mapping.csv`, `phase_review.csv` | Auditable channel mapping and reviewed event entries |
| Subject `selected_trial_waveforms.csv`, `trial_metrics.csv`, `summary_metrics.csv` | Selected individual curves and separate metrics for each subject; metric rows include operated-limb metadata |
| Subject/group `spm_*_field.csv` | SPM t field and threshold when parametric inference succeeds |
| Subject/group `exact_*_node_p.csv` | Max-t corrected node p-values for exact comparisons |
| Group `knee_*_input.csv` | Exact subject-average curves used for group tests |
| Selected images | Mean/SD comparisons and enabled SPM plots in selected formats |

Metrics include mean, minimum, maximum, range of motion (max − min), RMS, peak absolute magnitude, its location in percent, and area in source units × percentage points. Mean and RMS are trapezoid-weighted over the normalized coordinate; phase endpoints are linearly interpolated when fractional. “ROM” applied to non-angle channels means their numeric range. Summary peak values are **the mean of trial peaks**, not the peak of the averaged waveform. Phase summaries may have fewer contributing trials; use the count columns. SD with a single contributing trial is undefined.

Signed symmetry is `200 × (Right − Left) / (|Right| + |Left|)`. Zero denominators produce missing values. This measure is unstable near zero and is usually less meaningful for signed angles or powers; inspect absolute differences and RMSE. It is not a clinical limb-symmetry clearance threshold. EMG outputs describe the supplied waveforms; rectification, filtering and MVC normalization must be established upstream before physiological interpretation.

## SPM design and changes from the previous script

The previous script repeats independent two-sample t tests (`ttest2`) across many gait variables with hard-coded file paths and gait labels. This implementation uses reusable paired analysis, movement-percent labels and explicit source mappings.

**Subject level:** the observations are matched jump trials. Right and left from the same trial are paired. The 101 nodes are the continuous domain, not 101 independent observations. Averaging trials creates a representative curve but does not create replicates for a t test. The app preserves the trials for inference.

The optional parametric test is `spm1d.stats.ttest_paired(right, left)`, two-sided, with random-field cluster inference. It requires at least three pairs in this app and nonzero difference variance at every node. It assumes independent trial differences across trials, a sufficiently smooth field and appropriate Gaussian residual behavior. Those assumptions cannot be established with three trials. Results are exploratory and are not return-to-sport decisions. Two trials receive descriptive and exact comparisons only.

The default exact paired sign-flip comparison enumerates all `2^n` sign combinations (n = 2–12) and uses maximum absolute t across the 101 nodes to correct the node p-values. This is an exact permutation field comparison, **distinct from parametric cluster-based SPM**. It assumes sign exchangeability and independent trial differences. At n = 2, the smallest two-sided p is 0.5; at n = 3 it is 0.25. Failure to reject is not equivalence. More truly independent repeated trials or an appropriate group study is needed for stronger evidence.

**Group level:** each unique, aligned subject supplies one trial-averaged curve per limb. The test remains paired Right versus Left. This is not operated versus unoperated when surgical sides differ. A minimum of three unique subjects is enforced, but three remains a very small study. Duplicate file hashes are excluded; you must also ensure IDs represent genuinely independent people and that different files are not repeated sessions from the same person.

**Multiplicity:** default family alpha = 0.05 is divided across X, Y and Z. Random-field inference handles the continuous field within each plane; Bonferroni handles the three plane tests. Optionally correct over `(number of uploaded subjects + one group comparison) × 3`, conservatively counting even excluded files. If both exact and parametric tests are examined, designate the primary method before interpreting results rather than choosing whichever is significant. Additional hypotheses require a planned multiplicity strategy. Group eligibility, errors and settings are recorded; undefined fields are not silently patched with noise.

## Development validation

```bash
python -m pip install pytest
python -m pytest -q
```

Tests cover parsing, pairing, phase rejection, trial-level metric aggregation, fractional-phase integration, small-sample permutation limits, parametric SPM and figure output, plus the Streamlit upload/review/export workflow. One integration test uses the original private files only when present in the onboarding workspace and otherwise skips.

`requirements.txt` specifies compatible ranges. `requirements-tested.txt` records the exact direct dependency versions exercised during development.
