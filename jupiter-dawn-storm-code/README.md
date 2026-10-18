# Jupiter dawn-storm analysis code

Local release candidate **v1.1.0**, prepared 2026-10-18. This candidate has
not been uploaded or assigned a new DOI. No software license has been selected.
See [中文说明](README_中文.md) and [license status](LICENSE_STATUS.md).

Analysis associated with *Triggering Jupiter’s Auroral Dawn Storms via
Magnetotail Cross-Field Current Instability*. Only two **2017** intervals are
included: May 16 (07:52–10:03 UTC, 7861 one-second samples) and the May 10 quiet
crossing (08:00–10:00 UTC, 7201 samples). The current calibrated PSD and its
display range are distinguished from the original CWT amplitude displays.

## What is included

| Path | Contents |
| --- | --- |
| `run_analysis.py` | Common entry point; default task validates the package |
| `code/psd_calibration.py` | Component decomposition, calibrated PSD, peaks, Welch comparisons, synthetic validation |
| `code/plot_psd_figures.py`, `plot_quiet_psd.py` | Current Fig4 and S3; displayed PSD range 3.5–50 mHz |
| `code/plot_fig3.py` | Replot the archived Fig3 CWT amplitudes |
| `code/recompute_cwt.py` | Optional separate CWT recomputation at explicit precision 12 |
| `code/plot_bandpass.py` | Recompute the S1 2–5 mHz Butterworth signal |
| `data/` | Processed field windows, complete spectra and original display arrays |
| `original/` | Original wavelet helper and two polynomial-root scripts |
| `reference_figures/` | Current reference PNGs, added after validation |
| `docs/` | Methods, data dictionary, figure coverage, validation, publishing guide |
| `CITATION.cff`, `SOURCE_MANIFEST.json`, `SHA256SUMS.txt` | Citation metadata, provenance, file integrity |

## Installation

Python 3.11 or later is recommended. Create an environment and install the
ordinary third-party dependencies; no author-specific directory is required.

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe run_analysis.py check
```

On Linux/macOS use `.venv/bin/python` in place of
`.\.venv\Scripts\python.exe`. `numpy>=2` supplies `numpy.trapezoid`;
`PyWavelets>=1.9` must expose `cwt(..., precision=...)`. The exact library
versions and CWT implementation hash are recorded by a numerical rerun.

## Quick figure reproduction

```powershell
.\.venv\Scripts\python.exe run_analysis.py --list
.\.venv\Scripts\python.exe run_analysis.py figures
```

This uses the included, real processed 2017 observation windows and archived
CWT/PSD arrays. It plots Fig3, Fig4, S1 and S3 under `output/figures/`, in PNG,
PDF and SVG formats. It does not fetch raw data. Figure rendering can vary
with fonts and Matplotlib; inspect the source tables for numerical values.

## Recompute the calibrated PSD

```powershell
.\.venv\Scripts\python.exe run_analysis.py psd
.\.venv\Scripts\python.exe run_analysis.py figures
```

The first command recomputes both event spectra from the bundled **three
field-component** columns, not from the precomputed PSD. It also runs impulse
response, precision convergence, white-noise and sinusoidal response checks.
This may take several minutes. The second command uses the newly computed
spectra under the same output directory. Reference inputs under `data/` are
not overwritten. A custom output directory is supported by `--output`.

Optional prepared daily files can be supplied with `--mag-dir /path/to/mag`.
If a daily file is absent, that event falls back to the bundled field window.
See [DATA_SOURCES.md](DATA_SOURCES.md) for names, units and coordinate format.
Use `--save-coefficients` only if the large complex CWT coefficient arrays
are needed; they are excluded from the release archive.

Other tasks are `cwt`, `bandpass`, `theory`, and `all`. `all` includes the full
PSD/CWT numerical reruns, all supported plots and both original root scripts.
The `cwt` task writes `*_display_recomputed.npz` separately; it never silently
replaces the archived arrays used by Fig3 and S3(a–c).

## Interpretation and coverage

- The perpendicular PSD has a local maximum at 4.211917906742496 mHz.
- The complete main-event parallel PSD retains its local maximum near
  3.359238216634793 mHz. That point is outside the current **display** range.
- Restricting the display does not establish absence of parallel peaks in
  the full source grid or perpendicular power dominance.
- CWT amplitude maps use `abs(W)/sqrt(scale)` in nT. They are not PSD maps.
- The package starts from prepared JSS/RTP magnetic-field components. The
  original PDS-product-to-JSS/RTP preprocessing and HST reduction/map-making
  workflows are not supplied. HST observations are GO-14634 only.
- Original root scripts are provided for inspection and execution; exact
  correspondence of their plotted branches to final Fig1 is not certified.

See [methods](docs/METHODS.md), [data dictionary](docs/DATA_DICTIONARY.md),
[figure coverage](docs/FIGURE_MAP.md), and [validation](docs/VALIDATION.md).

## Public release and citation

Use the [Chinese GitHub/Zenodo guide](docs/公开步骤_中文.md). The previously
supplied citation was `fymhku (2026), v1.0.1`, DOI
`10.5281/zenodo.23184052`. It is a historical identifier, not an assigned DOI
for this v1.1.0 candidate. After publishing, use the **new version's actual
DOI** and finalize `CITATION.cff`, the record's creator list and manuscript
reference. No ORCID, paper DOI or new software DOI has been invented.
