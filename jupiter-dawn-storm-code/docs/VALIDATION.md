# Validation — 2026-10-05

The supported `python run_analysis.py all` invocation completed with exit code
0 using the authors' locally prepared May 10 and May 16 MAG input files.
Generated outputs were checked locally and are excluded from the release ZIP.

## Execution results

- All eight supported tasks completed: field, CWT, wavelet power, band-pass,
  additional quiet plots, quiet five-panel assembly, and both theory scripts.
- Twelve PNG files were generated. The quiet five-panel image is 1800 x 4161 px.
- May 16 perpendicular spectrum: highest local peak 4.212 mHz; the active-period
  ion-frequency 25th-75th percentile interval was approximately 3.825-6.413 mHz.
- May 16 active shading: 09:17:40-09:46:44 UTC. Envelope detection separately
  reported 09:27:38-09:36:52 UTC; it did not alter the configured shading.
- May 10 window: 7201 unique, finite samples at regular 1-second cadence.
- May 16 window: 7861 unique, finite samples at regular 1-second cadence.

## Structural checks

- All 58 Python files parsed successfully.
- The executable ASTs of all seven `juno_ulf` modules matched the supplied
  originals after excluding docstrings; the numerical algorithms were retained.
- All 61 recorded source files retained their original source SHA-256 hashes.
  Historical copies in `legacy/` are byte-identical to those sources.
- Relative Markdown links resolved. Text files had no unexpected control
  characters. Scans found no matches for common GitHub token, AWS access-key,
  private-key, or credential-assignment patterns.

## Actual runtime used for execution checks

| Component | Imported version |
| --- | --- |
| Python | 3.14.7 |
| NumPy | 2.5.3 |
| SciPy | 1.18.1 |
| pandas | 3.0.6 |
| Matplotlib | 3.11.2 |
| PyWavelets | 1.8.0 |
| Pillow | 12.3.0 |

Existing local scientific libraries were used; a fresh dependency installation
and every version allowed by `requirements.txt` were not tested.

The theory scripts emitted harmless font-fallback and non-interactive `show()`
warnings; their PNG files were saved successfully. Historical exploratory
entry points were syntax-checked but not run as a complete workflow. Execution
success does not independently validate the physical interpretation, PSD
normalization or exact correspondence of every final manuscript figure.
