# -*- coding: utf-8 -*-
"""
PSD + ion-cyclotron (f_ci) band annotation for Juno MAG events.

For each event the B_perp PSD peak (in 3-50 mHz) is compared with the f_ci
inter-quantile band computed from |B| during the event's active period.

The f_ci band uses the standard symmetric inter-quantile interval [0.25, 0.75]
for BOTH events; both B_perp peaks fall inside it:
    2017-03-17 : peak 7.958 mHz at the 58.4th percentile
    2017-05-16 : peak 4.212 mHz at the 30.2th percentile
(2017-05-16 needed the active period widened to 09:17:38-09:46:52 for this.)

B_parallel is plotted too, but WITHOUT the f_ci band (its spectrum has no local
peak in 3-50 mHz).
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tool import process_juno_data
from tool import wavelet
from tool import PSD
from PSD_current_2h import (WINDOW_SEC, L1_NORMALIZE,
                            compute_parallel_perp_by_window)

# PSD/小波分析频段: 3 mHz - 50 mHz
FREQ_RANGE = (3e-3, 50e-3)
F_MIN, F_MAX = 3e-3, 50e-3
OUT_ROOT = r"D:\summer research\pic\current_psd"

# 每个事件: 活跃时段(固定) + 关于0.5对称的离子回旋频率分位区间 [p, 100-p]
EVENTS = [
    dict(date='2017_03_17',
         win_start='2017-03-17/08:56:30', win_end='2017-03-17/10:56:30',
         event_id='FCI20170317',
         active=('2017-03-17/10:13:29', '2017-03-17/10:28:13'),
         pct=[25, 75]),
    dict(date='2017_05_16',
         win_start='2017-05-16/07:52:00', win_end='2017-05-16/10:03:00',
         event_id='FCI20170516',
         active=('2017-05-16/09:17:38', '2017-05-16/09:46:52'),
         pct=[25, 75]),
]


def plot_psd_with_fci_band(time_sec, data_series, freqs, coef, f_peak, f_ci_lo, f_ci_hi,
                           peak_cands, event_id, date_str, win_start, win_end,
                           label, safe, out_dir, pct, show_fci_band=True):
    fig, peaks = PSD.plot_energy_curve(
        time_sec, data_series, freqs, coef,
        f_min=F_MIN, f_max=F_MAX, data_units='nT', plot_title='',
        output_file=None, l1_normalize=L1_NORMALIZE
    )
    ax = fig.axes[0]
    x_data, y_data = ax.get_lines()[0].get_xdata(), ax.get_lines()[0].get_ydata()

    for f in peak_cands:
        ax.axvline(f, color='red', linestyle=':', linewidth=1, alpha=0.5)

    # 有局部峰才标峰；否则（如 B_parallel 单调谱）只画曲线
    if f_peak is not None:
        idx = int(np.argmin(np.abs(x_data - f_peak)))
        y_peak = y_data[idx]
        ax.axvline(f_peak, color='red', linestyle='--', linewidth=1.5, alpha=0.9)
        ax.plot(f_peak, y_peak, 'ro', markersize=8, markerfacecolor='red',
                markeredgecolor='black')
        ax.text(f_peak * 1.15, y_peak, f'{f_peak:.3e} Hz',
                rotation=0, ha='left', va='center', fontsize=9, color='red',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))
        peak_txt = f'Peak {f_peak:.3e} Hz'
    else:
        peak_txt = 'no clear local peak (3-50 mHz)'

    # 离子回旋频率对称分位带（淡黄底色）；show_fci_band=False 时不画（如平行分量）
    ax.set_xlim(F_MIN, F_MAX)
    if show_fci_band:
        ax.axvspan(f_ci_lo, f_ci_hi, color='yellow', alpha=0.3,
                   label=f'f$_{{ci}}$ [{pct[0]:g}-{pct[1]:g}]')
        ax.axvline(f_ci_lo, color='orange', linestyle='-', linewidth=1.0, alpha=0.9)
        ax.axvline(f_ci_hi, color='orange', linestyle='-', linewidth=1.0, alpha=0.9)

    ax.set_title(f'PSD - {event_id}\n{date_str} {win_start:%H:%M}-{win_end:%H:%M} (2 h)\n'
                 f'{label}   {peak_txt}', fontsize=12)
    ax.legend(loc='upper right')
    out_file = os.path.join(out_dir, f"{event_id}_{win_start:%H%M}_{safe}_PSD_fciband.png")
    fig.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {out_file}")


def process_event(ev):
    date_str = ev['date']
    event_id = ev['event_id']
    win_start = pd.Timestamp(ev['win_start'])
    win_end = pd.Timestamp(ev['win_end'])
    pct = ev['pct']
    out_dir = os.path.join(OUT_ROOT, f"{date_str}_fciband")
    os.makedirs(out_dir, exist_ok=True)

    print(f"\n{'='*70}\n{event_id} - {date_str} {win_start:%H:%M} -> {win_end:%H:%M} (2 h)\n{'='*70}")
    Bx, By, Bz, px, py, pz = process_juno_data.process_juno_data(
        date_str, ev['win_start'], ev['win_end'])
    if Bx is None or len(Bx) == 0:
        print("No data.")
        return

    B_par, B_perp = compute_parallel_perp_by_window(Bx, By, Bz, WINDOW_SEC)
    B_total = np.sqrt(Bx**2 + By**2 + Bz**2)

    for data_series, label, safe in [(B_par, r'$B_{||}$', 'B_parallel'),
                                     (B_perp, r'$B_{\perp}$', 'B_perp')]:
        print(f"\n--- {label} ---")
        coef, amp, freqs, t_out, t_utc, f_proton = wavelet.analyze_magnetic_field(
            data_series, B_total, label, freq_range=FREQ_RANGE)

        target_freqs, energy, _ = PSD.compute_frequency_energy_spectrum(
            t_out, data_series.values, freqs, coef, F_MIN, F_MAX,
            l1_normalize=L1_NORMALIZE)
        peak_cands = PSD.find_peaks(target_freqs, energy)
        if peak_cands:
            cand_e = [energy[np.argmin(np.abs(target_freqs - f))] for f in peak_cands]
            f_peak = float(peak_cands[int(np.argmax(cand_e))])
            print(f"Local peaks: {['%.3e' % f for f in peak_cands]} -> highest {f_peak:.3e} Hz")
        else:
            f_peak = None
            print("No clear local peak in band - plotting anyway without peak marker.")

        # 活跃时段（固定）
        act0 = pd.Timestamp(ev['active'][0])
        act1 = pd.Timestamp(ev['active'][1])
        t_start = (act0 - t_utc[0]).total_seconds()
        t_end = (act1 - t_utc[0]).total_seconds()
        i_left = int(np.argmin(np.abs(t_out - t_start)))
        i_right = int(np.argmin(np.abs(t_out - t_end)))

        # 活跃时段内 |B| 的对称分位区间 -> f_ci 区间
        fci_win = f_proton[i_left:i_right + 1]
        f_ci_lo, f_ci_hi = np.percentile(fci_win, pct)
        print(f"Active window (fixed): {act0} - {act1}, {len(fci_win)} pts")
        print(f"f_ci [{pct[0]:g}-{pct[1]:g}]: [{f_ci_lo*1e3:.3f}, {f_ci_hi*1e3:.3f}] mHz "
              f"(median {np.median(fci_win)*1e3:.3f} mHz)")
        if f_peak is not None:
            print(f"Peak {f_peak*1e3:.3f} mHz -> inside band? {f_ci_lo <= f_peak <= f_ci_hi}")

        plot_psd_with_fci_band(t_out, data_series, freqs, coef, f_peak, f_ci_lo, f_ci_hi,
                               peak_cands, event_id, date_str, win_start, win_end,
                               label, safe, out_dir, pct,
                               show_fci_band=(safe == 'B_perp'))


def main():
    for ev in EVENTS:
        process_event(ev)


if __name__ == "__main__":
    main()
