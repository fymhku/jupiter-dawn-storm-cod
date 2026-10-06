# -*- coding: utf-8 -*-
"""Field-aligned decomposition of the magnetic field.

The background (parallel) direction is taken as the mean field vector of a
non-overlapping block of length ``window_sec``.  For every sample inside that block

    B_par(t)  = B(t) . b0                    with b0 = <B> / |<B>|
    B_perp(t) = sqrt(|B(t)|^2 - B_par(t)^2)

``B_par`` is the compressional (field-aligned) component and ``B_perp`` the
transverse component magnitude.  Using a *block* mean keeps the definition
local in time, which matters when the background field rotates appreciably
over the 2-h window.
"""

from __future__ import annotations

from datetime import timedelta

import numpy as np
import pandas as pd


def parallel_perp_by_window(Br: pd.Series, Bth: pd.Series, Bph: pd.Series,
                            window_sec: float = 1200.0):
    """Split the field into parallel and perpendicular components.

    Parameters
    ----------
    Br, Bth, Bph : pandas.Series
        Spherical field components (nT) on a common UTC DatetimeIndex.
    window_sec : float
        Length of the running background-averaging block, in seconds
        (default 1200 s = 20 min).

    Returns
    -------
    (B_par, B_perp) : tuple of pandas.Series
        Parallel component (signed, nT) and perpendicular magnitude
        (non-negative, nT) on the same index.
    """
    idx = Br.index
    if not isinstance(idx, pd.DatetimeIndex):
        raise TypeError("Br must be indexed by a pandas DatetimeIndex")

    B_par = pd.Series(index=idx, dtype=float)
    B_perp = pd.Series(index=idx, dtype=float)

    current = idx.min()
    t_end = idx.max()
    while current <= t_end:
        block_end = current + timedelta(seconds=window_sec)
        mask = (idx >= current) & (idx < block_end)
        block_idx = idx[mask]
        if len(block_idx) == 0:
            current = block_end
            continue

        br = Br.loc[block_idx].to_numpy(dtype=float)
        bt = Bth.loc[block_idx].to_numpy(dtype=float)
        bp = Bph.loc[block_idx].to_numpy(dtype=float)

        mean_vec = np.array([br.mean(), bt.mean(), bp.mean()])
        norm = np.linalg.norm(mean_vec)
        b0 = mean_vec / norm if norm > 1e-12 else np.array([0.0, 0.0, 1.0])

        b_par = br * b0[0] + bt * b0[1] + bp * b0[2]
        b_tot = np.sqrt(br ** 2 + bt ** 2 + bp ** 2)
        b_perp = np.sqrt(np.maximum(b_tot ** 2 - b_par ** 2, 0.0))

        B_par.loc[block_idx] = b_par
        B_perp.loc[block_idx] = b_perp

        current = block_end

    return B_par, B_perp
