#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Time-averaged wavelet PSD with the ion-cyclotron frequency band.

For every event in ``config.EVENTS`` this script

  1. loads the window and splits the field into B_parallel / B_perp,
  2. computes the CWT and the time-averaged energy spectrum (3-50 mHz),
  3. finds the highest-energy local maximum of B_perp,
  4. derives the ion-cyclotron band from |B| inside the event's active period
     (symmetric inter-quantile interval, default [0.25, 0.75]) and
  5. plots the spectrum with that band shaded.

The band is drawn for B_perp only; B_parallel is plotted without it because its
spectrum has no local maximum in 3-50 mHz.

Output: ``<OUT_DIR>/psd/<date>_<component>_psd.png``
"""

import os

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import config
from juno_ulf import components, cwt, data, psd


def _active_seconds(ev, t_utc):
    """Active-period bounds expressed in seconds from the window start."""
    a0 = (pd.Timestamp(ev["active"][0]) - t_utc[0]).total_seconds()
    a1 = (pd.Timestamp(ev["active"][1]) - t_utc[0]).total_seconds()
    return a0, a1


def main():
    out_dir = os.path.join(config.OUT_DIR, "psd")
    os.makedirs(out_dir, exist_ok=True)

    f_min, f_max = config.PSD["f_min"], config.PSD["f_max"]
    p_lo, p_hi = config.PSD["percentiles"]

    for ev in config.EVENTS:
        print(f"\n=== {ev['event_id']} : {ev['win_start']} - {ev['win_end']} ===")
        Br, Bth, Bph = data.load_mag_window(
            ev["date"], ev["win_start"], ev["win_end"], config.MAG_DIR
        )
        B_total = np.sqrt(Br ** 2 + Bth ** 2 + Bph ** 2)
        B_par, B_perp = components.parallel_perp_by_window(
            Br, Bth, Bph, config.WINDOW_SEC
        )

        for series, label, tag, draw_band in (
            (B_par, r"$B_{||}$", "B_parallel", False),
            (B_perp, r"$B_\perp$", "B_perp", True),
        ):
            coef, amp, freqs, t_sec, t_utc, f_ion = cwt.compute_cwt(
                series, B_total, **config.CWT
            )

            target_freqs, energy, _ = psd.frequency_energy_spectrum(
                freqs, coef, f_min, f_max, l1_normalize=config.PSD["l1_normalize"]
            )
            f_peak, peak_list = psd.highest_peak(target_freqs, energy)

            a0, a1 = _active_seconds(ev, t_utc)
            fci_lo, fci_hi, window = psd.ion_cyclotron_band(
                f_ion, t_sec, a0, a1, percentiles=(p_lo, p_hi)
            )

            print(f"  {tag}: peak={'-' if f_peak is None else f'{f_peak*1e3:.3f} mHz'} "
                  f"| f_ci[{p_lo:g}-{p_hi:g}]=[{fci_lo*1e3:.3f}, {fci_hi*1e3:.3f}] mHz "
                  f"| inside={None if f_peak is None else (fci_lo <= f_peak <= fci_hi)}")

            fig = psd.plot_energy_curve(
                target_freqs, energy,
                title=f"PSD - {ev['event_id']}\n"
                      f"{ev['win_start'][11:]}-{ev['win_end'][11:]} UTC   {label}",
                data_units="nT", l1_normalize=config.PSD["l1_normalize"],
            )
            ax = fig.axes[0]
            ax.set_xlim(f_min, f_max)

            if draw_band:
                ax.axvspan(fci_lo, fci_hi, color="yellow", alpha=0.30,
                           label=rf"$f_{{ci}}$: {p_lo:g}-{p_hi:g}th percentiles")
                ax.axvline(fci_lo, color="orange", lw=1.0, alpha=0.9)
                ax.axvline(fci_hi, color="orange", lw=1.0, alpha=0.9)

            if f_peak is not None:
                y_peak = energy[int(np.argmin(np.abs(target_freqs - f_peak)))]
                ax.axvline(f_peak, color="red", ls="--", lw=1.5, alpha=0.9)
                ax.plot(f_peak, y_peak, "o", color="red", ms=7,
                        markeredgecolor="black")
                ax.text(f_peak * 1.15, y_peak, f"{f_peak*1e3:.2f} mHz",
                        fontsize=9, color="red", va="center",
                        bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                                  alpha=0.8))
            ax.legend(loc="upper right")

            out = os.path.join(out_dir, f"{ev['date']}_{tag}_psd.png")
            fig.savefig(out, dpi=150, bbox_inches="tight")
            plt.close(fig)
            print(f"  saved {out}")

    print("\nPSD done.")


if __name__ == "__main__":
    main()
