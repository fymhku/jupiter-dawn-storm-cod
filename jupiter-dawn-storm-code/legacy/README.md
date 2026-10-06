# Historical source snapshot

Standalone Python scripts and `tool/` are preserved byte for byte from the
supplied analysis folder, excluding `_pdf2txt.py`, which extracts text from
unrelated reference PDFs. Event lists and the original modular README,
configuration and dependency list are also retained.

These scripts include exploratory analyses and events absent from the
current manuscript. They contain absolute local input/output paths; edit
those paths in a personal working copy before running. They are not invoked
by the repository's supported `run_analysis.py`.

Some legacy scripts use pandas APIs removed in pandas 3, MVA routines from
`pyrfu`, `xarray`, and Excel export support. `requirements-legacy.txt` records
these additional dependencies; it is not a validated complete environment.
`main.py` imports `tool.find_period`, whose source file was not present in
the supplied folder. Retaining the source does not imply that every legacy
entry point is self-contained or executable.

The original configuration uses March 17 and May 16, and a slightly
different May 16 selected interval. The current portable configuration is
in `../juno_ulf_analysis/config.py`.
