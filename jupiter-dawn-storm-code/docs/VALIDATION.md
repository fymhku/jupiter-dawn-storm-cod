# Validation of the v1.1.0 local candidate

Executed 2026-10-18 using the supplied processed event windows, without
reading the authors' original daily-file path. All generated numerical
products were written to a separate output directory.

- `run_analysis.py all`: exit code 0. Both PSD calculations, CWT amplitude
  recomputation, band-pass calculation, four plots and both original
  polynomial root scripts completed.
- Main and quiet calibrated PSD arrays match the current distributed
  reference spectra **exactly**, for both components (all 200 grid points).
- Main CWT amplitude differences are at floating-point rounding level:
  max Bperp 1.67e-16 nT, max Bparallel 4.44e-16 nT. Quiet amplitudes match exactly.
- Fig3, Fig4 and S3 PNGs are pixel-identical to their current Word-embedded
  reference images in the tested environment.
- S1's numerical filtered signal matches the reference exactly. Its new
  rendering is 3779×1417 pixels; the manuscript-embedded reference is
  1800×680. S1 is not claimed to be a pixel-identical layout reconstruction.
- Discrete impulse reconstruction/alignment passed; precision18-to20
  convergence passed the supplied 0.1% bound for both events; the
  128-replicate white-noise and phase-averaged sinusoid checks passed.
- The root scripts produced image outputs; their execution is not a
  certification of correspondence to final Fig1's model or selected branch.
- Expected original-theory warnings about font fallback and noninteractive
  `plt.show()` do not affect the successful image export.

Full numerical checks are supplied under `validation/`; complete input
parameters, peaks, environment and CWT source hash are in
`data/numerical_results.json`.

## Actual tested environment

```json
{
  "python": "3.12.14",
  "numpy": "2.5.4",
  "pandas": "3.0.6",
  "scipy": "1.18.1",
  "PyWavelets_distribution": "1.10.0",
  "pywt_module_string": "1.8.0"
}
```

The installed PyWavelets distribution metadata and module version string
are inconsistent in the authors' local environment; both are recorded
instead of silently treating them as identical. The public precision
keyword is present and no library implementation is monkeypatched.
Dependencies are ordinary installed libraries outside the release, not
vendored copies. A clean internet installation on another operating system
was not tested. The package has no embedded author-specific dependency path.

## Scope

These checks establish reproducibility from the supplied prepared field
windows and figure arrays in the tested environment. They do not validate
raw instrument calibration, native-PDS-to-JSS/RTP conversion, HST reduction,
statistical significance of local peaks, the theoretical mechanism's unique
identification, or complete raw-data-to-every-figure reproduction.
