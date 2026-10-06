# -*- coding: utf-8 -*-
"""Shared figure helpers."""

from __future__ import annotations

import matplotlib.dates as mdates
import matplotlib.pyplot as plt


def add_utc_axis(ax, time_utc, label: str = "UTC", top: bool = True):
    """Add a secondary time axis labelled in UTC (``HH:MM``) to ``ax``.

    The primary axis is assumed to be in seconds from the start of the window;
    the secondary axis maps those seconds onto the real UTC timestamps.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis whose x data are seconds relative to ``time_utc[0]``.
    time_utc : pandas.DatetimeIndex
        UTC timestamps of the data.
    label : str
        Label for the new axis.
    top : bool
        Put the axis on top (True) or bottom (False).

    Returns
    -------
    matplotlib.axes.Axes
        The newly created twin axis.
    """
    twin = ax.twiny() if top else ax.twiny()
    twin.set_xlim(time_utc[0], time_utc[-1])
    twin.xaxis.set_major_locator(mdates.AutoDateLocator(minticks=5, maxticks=10))
    twin.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    twin.set_xlabel(label, fontsize=10)
    if not top:
        twin.xaxis.set_ticks_position("bottom")
        twin.xaxis.set_label_position("bottom")
    return twin


def format_hhmm(ax, label: str = "UTC"):
    """Format an axis whose x data are matplotlib date numbers as ``HH:MM``."""
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.set_xlabel(label)
    return ax
