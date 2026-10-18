# Data sources and provenance

## Reused public observation archives

- Juno magnetometer Jupiter archive, NASA PDS Planetary Plasma Interactions
  Node: https://doi.org/10.17189/1519711
  (mission index: https://pds-ppi.igpp.ucla.edu/mission/juno).
- Hubble Space Telescope ultraviolet observations, **GO-14634**,
  Mikulski Archive for Space Telescopes:
  https://doi.org/10.17909/t9-69wa-4q90.

Archive identifiers are retained from the manuscript's supplied data-source
information. The package does not download or redistribute complete raw
daily archive files, HST exposures, manuscript files or literature PDFs.

## Included processed observation windows

`main_timeseries.csv` and `quiet_timeseries.csv` contain the actual prepared
one-second magnetic-field windows and derived quantities from the authors'
analysis, not synthetic demonstrations. Only their Br, Btheta and Bphi
columns are needed for a fresh decomposition and PSD computation.

| Event | UTC interval (inclusive) | Samples | Prepared daily filename |
| --- | --- | --- | --- |
| main | 2017-05-16 07:52:00–10:03:00 | 7861 | JUNO_MAG_JSS_1s_2017_05_16.txt |
| quiet | 2017-05-10 08:00:00–10:00:00 | 7201 | JUNO_MAG_JSS_1s_2017_05_10.txt |

SHA-256 of the locally prepared source daily files recorded at calculation:

- main: `29da5d8ce74ff3f1def57010125e927464df5b0805fed0acad649f7a14f41ba8`
- quiet: `d9b9f266cabfa10e19ca158967b5e64c7e1b81c6115d47409a1730a079d9ee92`

These are prepared-input hashes, not asserted PDS-native product identifiers.
The original native-product selections, calibration and coordinate conversion
are not implemented in the supplied scripts. Thus full reproduction from
raw PDS products remains outside this package's demonstrated scope.

## Optional daily-file contract

Each daily file has one header line, then four whitespace-delimited columns:
UTC (`YYYY-MM-DD/HH:MM:SS`), Br, Btheta, Bphi. Fields are in nT in the prepared
spherical JSS/RTP frame. Sampling must be uninterrupted at 1 Hz; the selected
windows are inclusive. A PDS file cannot be converted by simply renaming it.
Place prepared files in `data/mag/` or pass `--mag-dir` to `run_analysis.py`.

## Derived arrays and tables

- Full PSD tables: the main grid is 3–50 mHz and the quiet grid 1–50 mHz,
  each with 200 original samples. Current figures display 3.5–50 mHz.
- `*_display_arrays.npz`: archived CWT amplitudes and frequency axes for
  exact data-level replotting of current Fig3 and S3(a–c).
- `*_original_raw_spectrum.csv`: historical raw squared-CWT power for
  transparent comparison, not the current calibrated PSD.
- `*_Welch.csv`: independently computed spectral-density reference tables.
- `FigS1_reference.csv`: the archived band-pass signal and source component.

See `docs/DATA_DICTIONARY.md`. File-specific origin and hashes are in
`SOURCE_MANIFEST.json`; complete distributed-file hashes are in
`SHA256SUMS.txt`. Upstream observation archives should be cited separately
from this software; their terms are not replaced by this package.
