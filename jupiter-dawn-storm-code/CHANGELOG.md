# Packaging changes — 2026-10-05

- Collected supplied Juno scripts, utility modules and polynomial solvers.
- Added the May 10 quiet-period plotting scripts from the manuscript folder.
- Kept a byte-identical historical source snapshot under `legacy/`.
- Made supported workflows use relative paths or environment variables.
- Set the main default event to May 16, 2-5 mHz band-pass, and aligned the
  selected interval to 09:17:40-09:46:44 UTC in the current main text and SM.
- Preserved the existing CWT grids, decomposition, filtering, wavelet-power
  calculation and polynomial coefficient arrays.
- Clarified metadata and stale quiet-panel comments; made percentile-band
  labels explicit; added a fallback font to working theory copies.
- Added a common runner, dependency lists, archive links, figure mapping,
  implementation notes and a source manifest.
- Excluded `.claude`, caches, generated figures, vendored libraries, wheels,
  THEMIS data, reference-paper text extracts and manuscripts.
- No source files in the original analysis or manuscript folders changed.
- No repository, author metadata, software license or DOI was invented.
