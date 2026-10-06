# Scripts and manuscript figures

| Figure / analysis | Relevant scripts | Reproduction scope |
| --- | --- | --- |
| Main theory figure | `theory/solve.py`, `theory/solve_8.py` | Supplied polynomial-root calculations; correspondence to the final model/figure needs verification |
| Main HST auroral image | None in the supplied analysis package | Public observations linked in `DATA_SOURCES.md`; HST reduction/map code is not included |
| Main May 16 field/CWT figure | `juno_ulf_analysis/run_field_plot.py`, `run_cwt.py` | Constituent plots; the supplied modular CWT grid is 3-50 mHz; final combined layout is not assembled here |
| Main May 16 spectral figure | `juno_ulf_analysis/run_psd.py` | Full-window raw-coefficient wavelet power plus active-period ion-frequency quartiles |
| SM Fig. S1 band-pass | `juno_ulf_analysis/run_bandpass.py` | 2-5 mHz signal; current 09:17:40-09:46:44 shaded interval |
| SM Fig. S2 quiet auroral image | None | HST data source only |
| SM Fig. S3 quiet five-panel plot | `supplementary/stack_0510.py`, `textlabel.py` | Existing five-panel assembly, 1-50 mHz, imported from the paper workspace |
| Additional quiet plots | `supplementary/make_quiet_figs.py` | Separate 2-5 mHz band-pass and two wavelet-power plots |

This mapping makes the available code's coverage explicit; it is not a
certification of exact reproduction of every final manuscript figure.
