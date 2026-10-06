# -*- coding: utf-8 -*-
"""Read locally prepared Juno MAG JSS/RTP 1-s ASCII files.

Locally prepared ASCII input format (one file per day)
--------------------------------------------------
    UT_Time                        Br(nT)          Btheta(nT)          Bphi(nT)
    2017-03-17/00:00:00       0.40889433       -1.6055669       -2.5356380
    2017-03-17/00:00:01       0.37830267       -1.6117249       -2.4714024
    ...

* cadence  : 1 s
* UT_Time  : ``YYYY-MM-DD/HH:MM:SS``
* Br/Btheta/Bphi : magnetic field in nT, spherical (r, theta, phi) components
* file name: ``JUNO_MAG_JSS_1s_<YYYY_MM_DD>.txt``
"""

from __future__ import annotations

import os

import pandas as pd

_FILENAME = "JUNO_MAG_JSS_1s_{date}.txt"
_COLUMNS = ["Br", "Btheta", "Bphi"]


def mag_file(date_str: str, mag_dir: str) -> str:
    """Full path of the daily MAG file for ``date_str`` (``YYYY_MM_DD``)."""
    return os.path.join(mag_dir, _FILENAME.format(date=date_str))


def load_mag_day(date_str: str, mag_dir: str) -> pd.DataFrame:
    """Load one full day of Juno MAG JSS 1-s data.

    Parameters
    ----------
    date_str : str
        Date in ``YYYY_MM_DD`` form (as used in the file name).
    mag_dir : str
        Folder containing the daily files.

    Returns
    -------
    pandas.DataFrame
        Indexed by UTC ``DatetimeIndex`` with columns
        ``Br``, ``Btheta``, ``Bphi`` (nT), sorted in time.
    """
    path = mag_file(date_str, mag_dir)
    if not os.path.exists(path):
        raise FileNotFoundError(f"MAG file not found: {path}")

    df = pd.read_csv(
        path,
        sep=r"\s+",          # whitespace-delimited
        header=None,
        skiprows=1,          # drop the header line
        names=["time_str", *_COLUMNS],
        dtype={"time_str": str},
    )
    df["time"] = pd.to_datetime(df["time_str"], format="%Y-%m-%d/%H:%M:%S",
                                errors="coerce")
    df = df.dropna(subset=["time"]).set_index("time").drop(columns=["time_str"])
    return df.sort_index()


def load_mag_window(date_str: str, win_start: str, win_end: str, mag_dir: str):
    """Load a UTC time window from one daily MAG file.

    Parameters
    ----------
    date_str : str
        ``YYYY_MM_DD``.
    win_start, win_end : str
        Inclusive UTC bounds, ``"YYYY-MM-DD HH:MM:SS"``.
    mag_dir : str
        Folder containing the daily files.

    Returns
    -------
    (Br, Btheta, Bphi) : tuple of pandas.Series
        Magnetic field components (nT) inside the window.
    """
    day = load_mag_day(date_str, mag_dir)
    win = day.loc[pd.Timestamp(win_start):pd.Timestamp(win_end)]
    if win.empty:
        raise ValueError(
            f"No data in {date_str} between {win_start} and {win_end}"
        )
    return win["Br"], win["Btheta"], win["Bphi"]
