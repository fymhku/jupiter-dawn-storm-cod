# -*- coding: utf-8 -*-
"""Juno ULF-wave analysis package.

Modules
-------
data        : read the Juno MAG JSS 1-s (RTP) text product
components  : B_parallel / B_perp decomposition from a non-overlapping block mean
cwt         : continuous wavelet transform (complex Morlet) + spectrogram
psd         : time-averaged wavelet power spectrum (PSD) and peak finding
plotting    : shared figure helpers (UTC time axes)
"""

__all__ = ["data", "components", "cwt", "psd", "plotting"]
__version__ = "1.0.0"
