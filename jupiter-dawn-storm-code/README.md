# Jupiter dawn-storm analysis code

Analysis scripts associated with the Jupiter dawn-storm manuscript and its
Supplemental Material. The supported workflows cover magnetic-field
decomposition, continuous wavelet transforms, wavelet power spectra,
band-pass filtering, the 2017-05-10 quiet-period panels, and polynomial-root
calculations. Historical exploratory scripts are retained in `legacy/`.

## Contents

| Location | Purpose |
| --- | --- |
| `run_analysis.py` | Common entry point for the supported workflows |
| `juno_ulf_analysis/` | Modular MAG analysis; defaults to 2017-05-16 |
| `supplementary/` | 2017-05-10 quiet-period figures and panel-label helper |
| `theory/` | Three- and eighth-degree polynomial-root scripts |
| `legacy/` | Original standalone analysis scripts, helpers and old configuration |
| `DATA_SOURCES.md` | Public archives and the required local input format |
| `docs/FIGURE_MAP.md` | Relation between scripts and manuscript figures |
| `docs/METHODS.md` | Implemented definitions, averaging ranges and normalization |
| `docs/VALIDATION.md` | Packaging checks and actual execution results |
| `SOURCE_MANIFEST.json` | Source-file hashes and changes made in this package |

## Install

Use Python 3.10 or later in an isolated environment. From this directory:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

On Linux or macOS, use `.venv/bin/python` instead of
`.\.venv\Scripts\python.exe` in the commands below.

## Data

Prepare the two daily MAG files in `data/mag/`:

```text
JUNO_MAG_JSS_1s_2017_05_10.txt
JUNO_MAG_JSS_1s_2017_05_16.txt
```

These are locally prepared, whitespace-delimited JSS/RTP input files, with
one header line and four columns: UTC time, Br, Btheta and Bphi. This package
does not include observation files or a PDS-native-to-JSS/RTP conversion
pipeline. See [DATA_SOURCES.md](DATA_SOURCES.md) for the exact input contract.
Raw observation archives and already prepared input files should not be
confused: renaming a PDS file alone does not convert its format or frame.

If the files are in another directory, set the path for this session:

```powershell
$env:JUNO_MAG_DIR = 'E:\your-data-folder'
```

Linux/macOS: `export JUNO_MAG_DIR=/path/to/your-data-folder`.

## Run

```powershell
.\.venv\Scripts\python.exe run_analysis.py --list
.\.venv\Scripts\python.exe run_analysis.py all
```

Alternatively run one task: `field`, `cwt`, `psd`, `bandpass`, `quiet`,
`quiet-panels`, `theory`, or `theory-8`. The theory tasks need no observation
files. Outputs are written under `output/`; override with `JUNO_ULF_OUT`.
The main event can be edited in `juno_ulf_analysis/config.py`.

## Settings retained or aligned during packaging

| Workflow | UTC window | Selected interval | Frequency range |
| --- | --- | --- | --- |
| Main MAG analysis | 2017-05-16 07:52:00-10:03:00 | 09:17:40-09:46:44 | Original modular CWT/power grid: 3-50 mHz; filter: 2-5 mHz |
| Quiet-period analysis | 2017-05-10 08:00:00-10:00:00 | 08:44:00-08:56:00 in the separate band-pass plot | CWT/power: 1-50 mHz; filter: 2-5 mHz |

The main active interval follows the current text: it contains the visually
largest band-passed perturbations. The default main workflow analyzes only
2017-05-16. The original March/May two-storm configuration is preserved in
`legacy/juno_ulf_analysis_config.py`.

The numerical definitions of the existing modular functions are retained.
The PSD-named routine averages **raw squared CWT coefficients over the full
analysis window**; the 25th-75th percentile ion-frequency band alone uses
the selected active interval. See [docs/METHODS.md](docs/METHODS.md).
Running these workflows is not a claim that every manuscript figure's final
arrangement or numerical value has been independently reproduced.

## Historical scripts

`legacy/` preserves the supplied standalone scripts byte for byte. These
have old events, absolute machine-specific paths and additional dependency
requirements. They are provided for provenance and inspection rather than
as the default entry point. See [legacy/README.md](legacy/README.md).

## Citation and license

Repository URL, release DOI, publication details and author metadata can be
added after the repository and manuscript are finalized. Cite the data
archives separately from the analysis code. No software license has been
selected in this packaging step; add the authors' chosen license on GitHub.
