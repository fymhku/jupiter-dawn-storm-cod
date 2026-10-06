#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Magnetic-field (RTP spherical) time series with the active period shaded.

For every event in ``config.EVENTS`` the three JSS-RTP components
(B_r, B_theta, B_phi) are plotted over the analysis window and the event's
active period is marked with a light-yellow band.

Output: ``<OUT_DIR>/field/<date>_B_RTP_field.png``
"""

import os

import matplotlib
matplotlib.use("Agg")

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd

import config
from juno_ulf import data


def main():
    out_dir = os.path.join(config.OUT_DIR, "field")
    os.makedirs(out_dir, exist_ok=True)

    for ev in config.EVENTS:
        print(f"\n=== {ev['event_id']} : {ev['win_start']} - {ev['win_end']} ===")
        Br, Bth, Bph = data.load_mag_window(
            ev["date"], ev["win_start"], ev["win_end"], config.MAG_DIR
        )

        act0 = pd.Timestamp(ev["active"][0])
        act1 = pd.Timestamp(ev["active"][1])

        fig, ax = plt.subplots(figsize=(12, 6))
        ax.plot(Br.index, Br.values, lw=0.8, label=r"$B_r$")
        ax.plot(Bth.index, Bth.values, lw=0.8, label=r"$B_\theta$")
        ax.plot(Bph.index, Bph.values, lw=0.8, label=r"$B_\phi$")
        ax.axvspan(act0, act1, color="yellow", alpha=0.25,
                   label="Active period")

        ax.set_ylabel("B (nT)")
        ax.set_xlabel("UTC")
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper right", fontsize=9)
        ax.set_title(f"B field (RTP) - {ev['event_id']}   "
                     f"{ev['win_start'][:10]}\n"
                     f"{ev['win_start'][11:]}-{ev['win_end'][11:]} UTC"
                     .replace("  ", " "), fontsize=12)
        fig.tight_layout()

        out = os.path.join(out_dir, f"{ev['date']}_B_RTP_field.png")
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  active period {act0:%H:%M:%S} - {act1:%H:%M:%S}; saved {out}")

    print("\nField plot done.")


if __name__ == "__main__":
    main()
