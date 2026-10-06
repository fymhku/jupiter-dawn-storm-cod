# -*- coding: utf-8 -*-
"""Zero-phase band-pass filtering and wave-packet (active period) detection."""

from __future__ import annotations

import numpy as np
from scipy.signal import butter, hilbert, sosfiltfilt


def bandpass(x, lo: float, hi: float, fs: float = 1.0, order: int = 4):
    """Zero-phase Butterworth band-pass (forward-backward, no phase shift).

    Parameters
    ----------
    x : array-like
        Input series.
    lo, hi : float
        Band edges in Hz.
    fs : float
        Sampling rate in Hz.
    order : int
        Filter order.

    Returns
    -------
    numpy.ndarray
        Filtered series, same length as ``x``.
    """
    sos = butter(order, [lo, hi], btype="bandpass", fs=fs, output="sos")
    return sosfiltfilt(sos, np.asarray(x, dtype=float))


def envelope(x, smooth_sec: float = 60.0, fs: float = 1.0):
    """Amplitude envelope (Hilbert magnitude) with a boxcar smoothing.

    The envelope is used only to locate wave packets, so a modest smoothing
    length is enough to suppress cycle-to-cycle wiggles.
    """
    x = np.asarray(x, dtype=float)
    env = np.abs(hilbert(x))
    n = max(1, int(round(smooth_sec * fs)))
    if n > 1:
        env = np.convolve(env, np.ones(n) / n, mode="same")
    return env


def detect_last_packet(env, threshold_frac: float = 0.30):
    """Bounds of the last contiguous run above ``threshold_frac * max(env)``.

    Parameters
    ----------
    env : array-like
        Amplitude envelope.
    threshold_frac : float
        Edge level as a fraction of the peak envelope.

    Returns
    -------
    (i0, i1) : tuple of int
        Inclusive sample indices of the last wave packet.
    """
    env = np.asarray(env, dtype=float)
    thr = threshold_frac * env.max()
    above = env > thr

    runs = []
    i, n = 0, len(above)
    while i < n:
        if above[i]:
            j = i
            while j < n and above[j]:
                j += 1
            runs.append((i, j - 1))
            i = j
        else:
            i += 1

    if not runs:
        return 0, n - 1
    return runs[-1]
