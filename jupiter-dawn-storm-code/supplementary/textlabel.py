# -*- coding: utf-8 -*-
"""
Text labels rendered with matplotlib's own font machinery.

Target size: the parenthesis in the manuscript's existing figures is 62 px tall
for an 1800 px wide image.  The font size in points is calibrated by measuring
the rendered bitmap, so it does not depend on DPI or figure size.

Style: DejaVu Sans BOLD (upright -- no oblique/italic).

Nothing is cropped from any image: the label is drawn as text.
"""
import io
import os
import sys

WS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(WS)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

FONT = "DejaVu Sans"
WEIGHT = "bold"
STYLE = "normal"           # upright, not oblique

TARGET_PAREN_PX = 62       # parenthesis height at 1800 px image width


def _render(text, fontsize, dpi):
    """
    Render `text` on a transparent background and return the alpha array.
    (buffer_rgba() is composited, so its alpha is uniformly 255; saving with
    transparent=True is what actually preserves the coverage.)
    """
    fig = Figure(figsize=(4, 2), dpi=dpi)
    FigureCanvasAgg(fig)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")
    ax.text(0.05, 0.50, text, fontsize=fontsize, fontfamily=FONT,
            fontweight=WEIGHT, fontstyle=STYLE, color="black",
            ha="left", va="center")
    buf = io.BytesIO()
    fig.savefig(buf, format="png", transparent=True, dpi=dpi)
    buf.seek(0)
    img = plt.imread(buf)
    return img[:, :, 3] if img.ndim == 3 else img


def _bbox(alpha):
    ys, xs = np.where(alpha > 0.08)
    if len(ys) == 0:
        return 0, 0
    return ys.max() - ys.min() + 1, xs.max() - xs.min() + 1


def calibrate(dpi=400):
    """Font size whose parenthesis is TARGET_PAREN_PX tall."""
    lo, hi = 2.0, 200.0
    for _ in range(45):
        mid = 0.5 * (lo + hi)
        ph, _ = _bbox(_render("(", mid, dpi))
        if ph < TARGET_PAREN_PX:
            lo = mid
        else:
            hi = mid
    fs = 0.5 * (lo + hi)
    ph, _ = _bbox(_render("(", fs, dpi))
    ah, _ = _bbox(_render("a", fs, dpi))
    bh, _ = _bbox(_render("b", fs, dpi))
    return fs, ph, ah, bh


def label_bitmap(text, fontsize, dpi=400, pad=6):
    """RGBA bitmap of `text`, cropped tight with a small transparent margin."""
    alpha = _render(text, fontsize, dpi)
    ys, xs = np.where(alpha > 0.08)
    if len(ys) == 0:
        raise ValueError("empty label %r" % text)
    crop = alpha[ys.min():ys.max() + 1, xs.min():xs.max() + 1]
    h, w = crop.shape
    out = np.zeros((h + 2 * pad, w + 2 * pad, 4), dtype=np.uint8)
    out[pad:pad + h, pad:pad + w, 3] = np.clip(crop * 255, 0, 255).astype(np.uint8)
    return out


if __name__ == "__main__":
    fs, ph, ah, bh = calibrate()
    print("calibrated font size: %.2f pt at dpi=400" % fs)
    print("  parenthesis      : %d px  (target %d)" % (ph, TARGET_PAREN_PX))
    print("  letter 'a' height: %d px" % ah)
    print("  letter 'b' height: %d px" % bh)
    from matplotlib import font_manager
    f = font_manager.findfont(font_manager.FontProperties(
        family=FONT, weight=WEIGHT, style=STYLE))
    print("font file:\n  %s" % f)
