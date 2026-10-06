#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""CWT spectrograms of the parallel and perpendicular field components.

For every event in ``config.EVENTS`` this script

  1. loads the Juno MAG JSS 1-s window,
  2. splits the field into B_parallel and B_perp (20-min non-overlapping background blocks),
  3. runs the complex-Morlet CWT and
  4. saves one spectrogram per component.

Output: ``<OUT_DIR>/cwt/<date>_<component>_cwt.png``
"""

import os

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

import config
from juno_ulf import components, cwt, data


def main():
    out_dir = os.path.join(config.OUT_DIR, "cwt")
    os.makedirs(out_dir, exist_ok=True)

    for ev in config.EVENTS:
        print(f"\n=== {ev['event_id']} : {ev['win_start']} - {ev['win_end']} ===")
        Br, Bth, Bph = data.load_mag_window(
            ev["date"], ev["win_start"], ev["win_end"], config.MAG_DIR
        )
        B_total = np.sqrt(Br ** 2 + Bth ** 2 + Bph ** 2)
        B_par, B_perp = components.parallel_perp_by_window(
            Br, Bth, Bph, config.WINDOW_SEC
        )

        for series, label, tag in ((B_par, r"$B_{||}$", "B_parallel"),
                                   (B_perp, r"$B_\perp$", "B_perp")):
            coef, amp, freqs, t_sec, t_utc, f_ion = cwt.compute_cwt(
                series, B_total, **config.CWT
            )
            title = (f"Juno MAG JSS-1s CWT - {ev['event_id']}\n"
                     f"{ev['win_start'][11:]}-{ev['win_end'][11:]} UTC   {label}")
            fig = cwt.plot_spectrogram(amp, freqs, t_sec, t_utc, f_ion,
                                       label, title=title)
            out = os.path.join(out_dir, f"{ev['date']}_{tag}_cwt.png")
            fig.savefig(out, dpi=150, bbox_inches="tight")
            plt.close(fig)
            print(f"  saved {out}")

    print("\nCWT done.")


if __name__ == "__main__":
    main()
