# -*- coding: utf-8 -*-
"""Portable defaults for the 2017-05-16 analysis.

The original configuration is archived in legacy/juno_ulf_analysis_config.py.
All times are UTC. The active interval follows the current manuscript and SM.
"""
import os
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
MAG_DIR = os.environ.get("JUNO_MAG_DIR", str(_ROOT / "data" / "mag"))
OUT_DIR = os.environ.get("JUNO_ULF_OUT", str(_ROOT / "output"))
EVENTS = [dict(
    event_id="2017-05-16",
    date="2017_05_16",
    win_start="2017-05-16 07:52:00",
    win_end="2017-05-16 10:03:00",
    active=("2017-05-16 09:17:40", "2017-05-16 09:46:44"),
    bands=[(2e-3, 5e-3)],
)]
WINDOW_SEC = 20 * 60
# Retain the original modular workflow's frequency grid.
CWT = dict(wavelet="cmor1.5-1.0", n_scales=200,
           freq_range=(3e-3, 50e-3), ion_mass_amu=16.0)
PSD = dict(f_min=3e-3, f_max=50e-3, l1_normalize=False,
           percentiles=(25, 75))
BANDPASS = dict(order=4, block_min=20,
                env_smooth_sec=60, env_threshold=0.30)
