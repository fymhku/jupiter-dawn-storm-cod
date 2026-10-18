"""Clip a semilog-x polyline without changing its visible geometry or data."""
import numpy as np

DISPLAY_RANGE_HZ = (0.0035, 0.05)


def visible_curve(frequency, power):
    """Return in-view source samples plus interpolated boundary intersections.

    Matplotlib/PDF clipping can leave off-axes strokes in the PDF's geometry.
    Clip the rendering path explicitly. Boundary intersections are linear in
    log(f), exactly as for an existing semilog-x line segment. They are not new
    PSD estimates and are never used for peak finding or numerical conclusions.
    All source samples and complete CSV files remain unchanged.
    """
    frequency = np.asarray(frequency)
    power = np.asarray(power)
    assert len(frequency) == len(power) and np.all(np.diff(frequency) > 0)
    low, high = DISPLAY_RANGE_HZ
    inside = (frequency >= low) & (frequency <= high)
    view_f = frequency[inside].copy()
    view_p = power[inside].copy()
    assert len(view_f) >= 2
    if not np.isclose(view_f[0], low, rtol=1e-14, atol=0):
        view_f = np.r_[low, view_f]
        view_p = np.r_[np.interp(np.log(low), np.log(frequency), power), view_p]
    if not np.isclose(view_f[-1], high, rtol=1e-14, atol=0):
        view_f = np.r_[view_f, high]
        view_p = np.r_[view_p, np.interp(np.log(high), np.log(frequency), power)]
    locations = np.searchsorted(view_f, frequency[inside])
    assert np.array_equal(view_f[locations], frequency[inside])
    assert np.array_equal(view_p[locations], power[inside])
    return view_f, view_p
