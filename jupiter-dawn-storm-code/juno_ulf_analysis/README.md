# Modular MAG workflow

Use the repository-level `README.md`, `requirements.txt` and `run_analysis.py`.
This copy uses repository-relative paths or the `JUNO_MAG_DIR` and
`JUNO_ULF_OUT` environment variables. The default event is 2017-05-16.

Scripts can also be run directly from the repository root:

```text
python juno_ulf_analysis/run_field_plot.py
python juno_ulf_analysis/run_cwt.py
python juno_ulf_analysis/run_psd.py
python juno_ulf_analysis/run_bandpass.py
```

The original configuration and README are archived in `../legacy/`.
