# -*- coding: utf-8 -*-
"""Continuous wavelet transform (CWT) analysis of Juno magnetic-field data.

A complex Morlet wavelet (``cmor1.5-1.0``) is used on 1-s sampled data.  The
coefficients are normalised by ``1/sqrt(scale)`` before the amplitude is taken,
so that amplitude is comparable across frequencies.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pywt
from matplotlib.ticker import FuncFormatter, LogLocator

from .plotting import add_utc_axis

# e / (2*pi*m_p) expressed for B in nT  ->  0.0152 Hz per nT (proton)
_PROTON_HZ_PER_NT = 0.0152


def ion_cyclotron_frequency(b_total_nT, mass_amu: float = 16.0) -> np.ndarray:
    """Ion cyclotron frequency ``f_ci = e B / (2 pi m_i)``.

    Parameters
    ----------
    b_total_nT : array-like
        Total field magnitude in nT.
    mass_amu : float
        Ion mass in atomic mass units.  ``1`` gives the proton gyrofrequency;
        the default ``16`` corresponds to singly charged oxygen (O+), which is
        the relevant ion for the Jovian inner magnetosphere.

    Returns
    -------
    numpy.ndarray
        ``f_ci`` in Hz.
    """
    return _PROTON_HZ_PER_NT * np.asarray(b_total_nT, dtype=float) / mass_amu


def compute_cwt(data_series: pd.Series, b_total, freq_range=(3e-3, 50e-3),
                wavelet: str = "cmor1.5-1.0", n_scales: int = 200,
                ion_mass_amu: float = 16.0, sampling_rate: float = 1.0):
    """Continuous wavelet transform of one field component.

    Parameters
    ----------
    data_series : pandas.Series
        Component time series (nT) on a UTC DatetimeIndex.
    b_total : array-like
        Total field magnitude ``|B|`` (nT) at the same samples; used for f_ci.
    freq_range : (float, float)
        Analysed frequency band in Hz.
    wavelet : str
        PyWavelets wavelet name.
    n_scales : int
        Number of logarithmically spaced scales.
    ion_mass_amu : float
        Ion mass (amu) used for the cyclotron frequency.
    sampling_rate : float
        Sampling rate in Hz (the JSS 1-s product is 1 Hz).

    Returns
    -------
    coef : ndarray, shape (n_scales, n_time)
        Complex wavelet coefficients.
    amplitude : ndarray, shape (n_scales, n_time)
        ``|coef| / sqrt(scale)``.
    freqs : ndarray
        Log-spaced frequencies (Hz).
    time_sec : ndarray
        Seconds relative to the first sample.
    time_utc : pandas.DatetimeIndex
        UTC timestamps.
    f_ion : ndarray
        Ion cyclotron frequency (Hz) at each sample.
    """
    data = data_series.to_numpy(dtype=float)
    time_utc = data_series.index
    time_sec = np.asarray(
        [(t - time_utc[0]).total_seconds() for t in time_utc], dtype=float
    )

    fmin, fmax = freq_range
    freqs = np.logspace(np.log10(fmin), np.log10(fmax), n_scales)

    centre_freq = pywt.central_frequency(wavelet)
    scales = centre_freq * sampling_rate / freqs

    coef, _ = pywt.cwt(data, scales, wavelet, sampling_period=1.0 / sampling_rate)
    amplitude = np.abs(coef) / np.sqrt(scales)[:, np.newaxis]

    f_ion = ion_cyclotron_frequency(b_total, mass_amu=ion_mass_amu)
    return coef, amplitude, freqs, time_sec, time_utc, f_ion


def _sci_formatter(y, pos):
    """Format a log-frequency axis as powers of ten."""
    if y <= 0:
        return "0"
    exponent = int(np.floor(np.log10(y)))
    if -4 <= exponent < 0:
        return rf"$10^{{{exponent}}}$"
    coeff = y / (10 ** exponent)
    return rf"${coeff:.0f} \times 10^{{{exponent}}}$"


def plot_spectrogram(amplitude, freqs, time_sec, time_utc, f_ion,
                     component_name: str, title: str = "", cmap: str = "jet",
                     figsize=(14, 7)):
    """CWT amplitude spectrogram with the ion cyclotron frequency overlaid.

    Returns the ``matplotlib.figure.Figure`` (not shown, not saved).
    """
    fig, ax = plt.subplots(figsize=figsize)

    vmin = np.percentile(amplitude, 5)
    vmax = np.percentile(amplitude, 95)
    img = ax.pcolormesh(time_sec, freqs, amplitude, shading="auto",
                        cmap=cmap, vmin=vmin, vmax=vmax)

    ax.set_yscale("log")
    ax.set_ylim(freqs[0], freqs[-1])
    ax.yaxis.set_major_locator(LogLocator(numticks=6))
    ax.yaxis.set_major_formatter(FuncFormatter(_sci_formatter))
    ax.yaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_ylabel("Frequency (Hz)")
    ax.set_xlabel("Time (s)")

    ax.plot(time_sec, f_ion, color="w", lw=1.2, alpha=0.9,
            label=r"$f_{ci}$ (O$^+$)")

    cbar = fig.colorbar(img, ax=ax, pad=0.02)
    cbar.set_label("Amplitude (nT)")

    ax.set_title(title or f"CWT spectrogram - {component_name}", fontsize=12)
    ax.legend(loc="upper right", facecolor="#333333", edgecolor="white",
              labelcolor="white", fontsize=9)

    add_utc_axis(ax, time_utc, label="UTC")
    fig.tight_layout()
    return fig
