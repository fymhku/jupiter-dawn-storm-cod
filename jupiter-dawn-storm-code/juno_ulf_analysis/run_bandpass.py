#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Band-passed perpendicular field and the active period.

For every event and every band in ``config.EVENTS[...]['bands']`` this script

  1. loads the window and builds the perpendicular magnitude |B_perp|,
  2. band-pass filters it (zero-phase Butterworth) and
  3. plots the filtered series with the event's active period shaded.

The shaded active interval is supplied by ``config.EVENTS``. For
2017-05-16 it marks the interval containing the visually largest perturbations
in the band-passed series. Envelope detection is printed as a separate
diagnostic and does not replace the configured interval.

Output: ``<OUT_DIR>/bandpass/<date>_<band>_bandpass.png``
"""

import os

import matplotlib
matplotlib.use("Agg")

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd

import config
from juno_ulf import bandpass as bp
from juno_ulf import components, data


def _band_tag(lo, hi):
    return f"{lo*1e3:g}-{hi*1e3:g}mHz"


def main():
    out_dir = os.path.join(config.OUT_DIR, "bandpass")
    os.makedirs(out_dir, exist_ok=True)

    for ev in config.EVENTS:
        print(f"\n=== {ev['event_id']} : {ev['win_start']} - {ev['win_end']} ===")
        Br, Bth, Bph = data.load_mag_window(
            ev["date"], ev["win_start"], ev["win_end"], config.MAG_DIR
        )
        _, B_perp = components.parallel_perp_by_window(
            Br, Bth, Bph, config.WINDOW_SEC
        )
        tw = B_perp.index

        act0 = pd.Timestamp(ev["active"][0])
        act1 = pd.Timestamp(ev["active"][1])

        for lo, hi in ev["bands"]:
            sig = bp.bandpass(B_perp.to_numpy(), lo, hi,
                              fs=1.0, order=config.BANDPASS["order"])
            env = bp.envelope(sig, smooth_sec=config.BANDPASS["env_smooth_sec"])
            i0, i1 = bp.detect_last_packet(env, config.BANDPASS["env_threshold"])
            print(f"  [{_band_tag(lo, hi)}] detected last packet "
                  f"{tw[i0]:%H:%M:%S} - {tw[i1]:%H:%M:%S} (peak env "
                  f"{env[i0:i1+1].max():.3f} nT); shading config active period "
                  f"{act0:%H:%M:%S} - {act1:%H:%M:%S}")

            fig, ax = plt.subplots(figsize=(12, 5.5))
            ax.plot(tw, sig, color="#1f77b4", lw=1.4, label="band pass")
            ax.axvspan(act0, act1, color="#fff59d", alpha=0.45, lw=0,
                       zorder=0, label="active period")
            ax.set_ylabel(r"$|B_\perp|$ (nT)")
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
            ax.set_xlabel("UTC")
            ax.grid(True, alpha=0.3)
            ax.legend(loc="upper right", fontsize=9)
            ax.set_title(rf"$B_\perp$ band pass {lo*1e3:g}-{hi*1e3:g} mHz, "
                         rf"{ev['win_start'][11:]}-{ev['win_end'][11:]} UTC",
                         fontsize=12)
            fig.tight_layout()

            out = os.path.join(out_dir,
                               f"{ev['date']}_{_band_tag(lo, hi)}_bandpass.png")
            fig.savefig(out, dpi=150)
            plt.close(fig)
            print(f"  saved {out}")

    print("\nBand-pass done.")


if __name__ == "__main__":
    main()
