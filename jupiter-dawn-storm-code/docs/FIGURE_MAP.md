# Figure-to-code map

| Figure | Entry / source | Demonstrated scope |
| --- | --- | --- |
| Fig1 | `run_analysis.py theory`, `original/theory/solve*.py` | Original polynomial root scripts execute; exact final-branch/figure correspondence not certified |
| Fig2 (HST) | Archive route in `DATA_SOURCES.md` | HST reduction and auroral mapping not supplied |
| Fig3 | `run_analysis.py figures`, `code/plot_fig3.py`, `data/main_display_arrays.npz` | Replot archived current CWT amplitudes; optional numerical rerun is separate `cwt` task |
| Fig4 | `run_analysis.py psd` then `figures` | Full PSD recomputation and current figure plotting |
| S1 | `run_analysis.py bandpass` | Field decomposition, 2–5 mHz filtering, plot and numerical series |
| S2 (HST) | Archive route in `DATA_SOURCES.md` | HST reduction and auroral mapping not supplied |
| S3(a–c) | `code/plot_quiet_psd.py`, `data/quiet_display_arrays.npz` | Original archived field/CWT panels retained |
| S3(d,e) | `run_analysis.py psd` then `figures` | Quiet-event PSD recomputed and plotted in current display range |

The archived figure PNGs and numerical source arrays serve different
purposes. Plotting from a saved array demonstrates plotting reproducibility;
the separate PSD rerun demonstrates the supplied field-window-to-PSD step.
This is not a claim of raw-instrument-to-every-figure reproduction.
