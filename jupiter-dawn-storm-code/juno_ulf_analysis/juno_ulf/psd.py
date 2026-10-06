# -*- coding: utf-8 -*-
"""Time-averaged wavelet power spectrum (referred to as PSD in the paper).

The spectrum is the time average of the squared wavelet coefficients,

    E(f) = < |W(f, t)|^2 >_t

evaluated on the CWT frequency grid.  ``find_peaks`` locates local maxima of
that curve (plateaus are handled explicitly), and the highest-energy maximum
within the analysis band is reported as the spectral peak.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np


def frequency_energy_spectrum(freqs, coefs, f_min: float, f_max: float,
                              l1_normalize: bool = False,
                              check_shape=None):
    """Time-averaged energy spectrum inside ``[f_min, f_max]``.

    Parameters
    ----------
    freqs : ndarray
        CWT frequency grid (Hz).
    coefs : ndarray, shape (n_freq, n_time)
        Complex wavelet coefficients.
    f_min, f_max : float
        Frequency bounds of interest (Hz).
    l1_normalize : bool
        If True, normalise so that the sum of the returned energy is 1.
    check_shape : tuple, optional
        ``(n_freq, n_time)`` used to validate / transpose ``coefs``.

    Returns
    -------
    target_freqs : ndarray
        Frequencies inside the band.
    energy : ndarray
        Corresponding time-averaged energy density (nT^2).
    peak_freq : float
        Frequency of the largest value of ``energy``.
    """
    coefs = np.asarray(coefs)
    if check_shape is not None and coefs.shape != tuple(check_shape):
        if check_shape is not None and coefs.shape == tuple(check_shape)[::-1]:
            coefs = coefs.T
        else:
            raise ValueError(
                f"coefs shape {coefs.shape} does not match {tuple(check_shape)}"
            )

    energy_density = np.mean(np.abs(coefs) ** 2, axis=1)

    mask = (freqs >= f_min) & (freqs <= f_max)
    target_freqs = freqs[mask]
    energy = energy_density[mask]

    if l1_normalize:
        total = energy.sum()
        if total > 0:
            energy = energy / total

    peak_freq = float(target_freqs[int(np.argmax(energy))]) if energy.size else np.nan
    return target_freqs, energy, peak_freq


def find_peaks(freq, energy):
    """Local maxima of ``energy(freq)``.

    A point is accepted when it is not smaller than both neighbours and the
    surrounding plateau is strictly greater than the values just outside it.
    Returns a list of peak frequencies (Hz).
    """
    energy = np.asarray(energy, dtype=float)
    n = len(energy)
    if n < 3:
        return []

    candidate = [False] * n
    for i in range(1, n - 1):
        if energy[i] >= energy[i - 1] and energy[i] >= energy[i + 1]:
            candidate[i] = True

    intervals = []
    i = 1
    while i < n - 1:
        if candidate[i]:
            start = i
            while i < n - 1 and candidate[i]:
                i += 1
            intervals.append((start, i - 1))
        else:
            i += 1

    peaks = []
    for start, end in intervals:
        if energy[start] > energy[start - 1] and energy[end] > energy[end + 1]:
            peaks.append(float(freq[(start + end) // 2]))
    return peaks


def highest_peak(freqs, energy):
    """Highest-energy local maximum, or ``None`` if the spectrum is monotonic."""
    peaks = find_peaks(freqs, energy)
    if not peaks:
        return None, []
    energies = [energy[int(np.argmin(np.abs(freqs - f)))] for f in peaks]
    return float(peaks[int(np.argmax(energies))]), peaks


def ion_cyclotron_band(f_ion, time_sec, active_start_s: float, active_end_s: float,
                       percentiles=(25, 75)):
    """Inter-quantile band of ``f_ion`` inside an active time window.

    Parameters
    ----------
    f_ion : array-like
        Ion cyclotron frequency (Hz) per sample.
    time_sec : array-like
        Sample times in seconds relative to the window start.
    active_start_s, active_end_s : float
        Active-period bounds in the same seconds.
    percentiles : (float, float)
        Symmetric quantiles, e.g. ``(25, 75)``.

    Returns
    -------
    (lo, hi) : tuple of float
        Band edges in Hz.
    """
    f_ion = np.asarray(f_ion, dtype=float)
    time_sec = np.asarray(time_sec, dtype=float)
    i0 = int(np.argmin(np.abs(time_sec - active_start_s)))
    i1 = int(np.argmin(np.abs(time_sec - active_end_s)))
    if i1 < i0:
        i0, i1 = i1, i0
    window = f_ion[i0:i1 + 1]
    lo, hi = np.percentile(window, percentiles)
    return float(lo), float(hi), window


def plot_energy_curve(freqs, energy, title: str = "", data_units: str = "nT",
                      l1_normalize: bool = False, figsize=(9, 5.5)):
    """Semi-log plot of the energy spectrum.  Returns the figure."""
    fig, ax = plt.subplots(figsize=figsize)
    ax.semilogx(freqs, energy, "b-", lw=1.8, label="Energy density")
    ax.set_xlabel("Frequency (Hz)")
    if l1_normalize:
        ax.set_ylabel("Normalised magnetic energy (1/Hz)")
    else:
        ax.set_ylabel(f"Magnetic energy ({data_units}$^2$/Hz)")
    ax.set_title(title or "Frequency-energy spectrum")
    ax.grid(True, which="both", ls="--", alpha=0.6)
    ax.legend(loc="upper right")
    fig.tight_layout()
    return fig
