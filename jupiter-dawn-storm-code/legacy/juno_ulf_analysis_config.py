# -*- coding: utf-8 -*-
"""
Central configuration for the Juno ULF-wave analysis.

Data location
-------------
The Juno MAG JSS 1-s product (spherical RTP components, 1-s cadence) is read
from ``MAG_DIR``.  Point it at the folder that holds files named

    JUNO_MAG_JSS_1s_<YYYY_MM_DD>.txt

either by editing the default below or by setting the environment variable
``JUNO_MAG_DIR``.  The files are public PDS products.
"""

import os

# ----------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------
MAG_DIR = os.environ.get(
    "JUNO_MAG_DIR",
    r"D:\Juno_Data\Juno_data_txt\mag_tvar\MAG_JSS_1s\juno_B1s_JSSRTP",
)

_HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.environ.get("JUNO_ULF_OUT", os.path.join(_HERE, "output"))

# ----------------------------------------------------------------------
# Events used in the paper
# ----------------------------------------------------------------------
#   date      : file-name date (YYYY_MM_DD)
#   win_start/win_end : 2-h analysis window, UTC "YYYY-MM-DD HH:MM:SS"
#   active    : the "active period" (UTC) used for the ion-cyclotron band.
#               For 2017-05-16 the wave packet 09:27:38-09:36:52 was widened
#               by +/-10 min so that the standard [0.25, 0.75] ion-cyclotron
#               interval contains the B_perp spectral peak.
#   bands     : band-pass bands (Hz) plotted for this event
EVENTS = [
    dict(
        event_id="2017-03-17",
        date="2017_03_17",
        win_start="2017-03-17 08:56:30",
        win_end="2017-03-17 10:56:30",
        active=("2017-03-17 10:13:29", "2017-03-17 10:28:13"),
        bands=[(6e-3, 9e-3)],
    ),
    dict(
        event_id="2017-05-16",
        date="2017_05_16",
        win_start="2017-05-16 07:52:00",
        win_end="2017-05-16 10:03:00",
        active=("2017-05-16 09:17:38", "2017-05-16 09:46:52"),
        bands=[(6e-3, 9e-3), (2e-3, 5e-3)],
    ),
]

# ----------------------------------------------------------------------
# Analysis parameters
# ----------------------------------------------------------------------
WINDOW_SEC = 20 * 60          # background-averaging window for B_par / B_perp

CWT = dict(
    wavelet="cmor1.5-1.0",    # complex Morlet (bandwidth 1.5, centre freq 1.0)
    n_scales=200,             # number of log-spaced scales
    freq_range=(3e-3, 50e-3),  # Hz
    ion_mass_amu=16.0,        # ion mass for f_ci (16 amu = O+)
)

PSD = dict(
    f_min=3e-3,               # Hz
    f_max=50e-3,              # Hz
    l1_normalize=False,
    percentiles=(25, 75),     # symmetric inter-quantile interval for f_ci
)

BANDPASS = dict(
    order=4,                  # Butterworth order (zero-phase, sosfiltfilt)
    block_min=20,             # background block length (min) for B_perp
    env_smooth_sec=60,        # envelope smoothing (s)
    env_threshold=0.30,       # fraction of peak envelope -> packet edge
)
