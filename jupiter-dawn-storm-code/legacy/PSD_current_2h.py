import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.dates import DateFormatter, AutoDateLocator
from datetime import timedelta
from tool import process_juno_data
from tool import wavelet
from tool import PSD

# 配置
L1_NORMALIZE = False
FREQ_RANGE = [0.0001, 0.05]
WINDOW_SEC = 20 * 60
F_MIN, F_MAX = 1e-3, 0.1
POWER_THRESHOLD_FRAC = 0.8
HALF_WINDOW = timedelta(hours=1)          # 前后各 1 小时 -> 总长 2 小时
EVENT_FILE = r"D:\summer research\pic\current.txt"
OUT_ROOT = r"D:\summer research\pic\current_psd"


def parse_events(path):
    """读取事件文件，每行: <date> <start_time> <date> <end_time>（日期重复出现，可能非零填充、空格不规整）。"""
    events = []
    with open(path, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            date_raw = parts[0]
            t_start = parts[1]
            # 格式 A: date start date end (4 tokens)；格式 B: date start end (3 tokens)
            t_end = parts[3] if len(parts) >= 4 else parts[2]
            y, m, d = date_raw.split('_')
            date_padded = f"{int(y):04d}_{int(m):02d}_{int(d):02d}"
            date_dash = date_padded.replace('_', '-')
            start_dt = pd.Timestamp(f"{date_dash} {t_start}")
            end_dt = pd.Timestamp(f"{date_dash} {t_end}")
            events.append({'date': date_padded, 'start': start_dt, 'end': end_dt})
    return events


def center_window(start_dt, end_dt):
    """以事件中点对称拓宽：两侧各 HALF_WINDOW，总长 2 小时。"""
    mid = start_dt + (end_dt - start_dt) / 2
    return mid - HALF_WINDOW, mid + HALF_WINDOW


def compute_parallel_perp_by_window(Bx, By, Bz, window_sec):
    idx = Bx.index
    B_par = pd.Series(index=idx, dtype=float)
    B_perp = pd.Series(index=idx, dtype=float)
    t_start, t_end = idx.min(), idx.max()
    current = t_start
    while current <= t_end:
        window_end = current + timedelta(seconds=window_sec)
        mask = (idx >= current) & (idx < window_end)
        win_idx = idx[mask]
        if len(win_idx) == 0:
            current = window_end
            continue
        Bx_win = Bx.loc[win_idx].mean()
        By_win = By.loc[win_idx].mean()
        Bz_win = Bz.loc[win_idx].mean()
        B_norm = np.sqrt(Bx_win**2 + By_win**2 + Bz_win**2)
        if B_norm < 1e-12:
            b0 = np.array([0, 0, 1])
        else:
            b0 = np.array([Bx_win, By_win, Bz_win]) / B_norm
        Bx_vals = Bx.loc[win_idx].values
        By_vals = By.loc[win_idx].values
        Bz_vals = Bz.loc[win_idx].values
        B_par_vals = Bx_vals*b0[0] + By_vals*b0[1] + Bz_vals*b0[2]
        B_tot_vals = np.sqrt(Bx_vals**2 + By_vals**2 + Bz_vals**2)
        B_perp_vals = np.sqrt(np.maximum(B_tot_vals**2 - B_par_vals**2, 0))
        B_par.loc[win_idx] = B_par_vals
        B_perp.loc[win_idx] = B_perp_vals
        current = window_end
    return B_par, B_perp


def active_time_window(power_t, time_sec, frac=POWER_THRESHOLD_FRAC):
    """返回功率 >= frac*max 的连续区间（从全局峰值向两侧扩展），单位为秒。"""
    pmax = power_t.max()
    if pmax <= 0:
        return time_sec[0], time_sec[-1], 0, len(time_sec) - 1
    thr = frac * pmax
    peak_idx = int(np.argmax(power_t))
    i = peak_idx
    while i > 0 and power_t[i - 1] >= thr:
        i -= 1
    left = i
    i = peak_idx
    while i < len(power_t) - 1 and power_t[i + 1] >= thr:
        i += 1
    right = i
    return time_sec[left], time_sec[right], left, right


def plot_psd_annotated(time_sec, data_series, freqs, coef, f_peak, f_ci_lo, f_ci_hi,
                       peak_cands, event_id, date_str, win_start, win_end,
                       label, safe, out_dir):
    fig, peaks = PSD.plot_energy_curve(
        time_sec, data_series, freqs, coef,
        f_min=F_MIN, f_max=F_MAX, data_units='nT', plot_title='',
        output_file=None, l1_normalize=L1_NORMALIZE
    )
    ax = fig.axes[0]
    lines = ax.get_lines()
    x_data, y_data = lines[0].get_xdata(), lines[0].get_ydata()

    # 0) 所有局部极大峰：淡虚线标出，便于看谱结构
    for f in peak_cands:
        ax.axvline(f, color='red', linestyle=':', linewidth=1, alpha=0.5)

    # 1) 最高峰（能量最大）标注
    idx = int(np.argmin(np.abs(x_data - f_peak)))
    y_peak = y_data[idx]
    ax.axvline(f_peak, color='red', linestyle='--', linewidth=1.5, alpha=0.9)
    ax.plot(f_peak, y_peak, 'ro', markersize=8, markerfacecolor='red',
            markeredgecolor='black')
    ax.text(f_peak * 1.15, y_peak, f'{f_peak:.3e} Hz',
            rotation=0, ha='left', va='center', fontsize=9, color='red',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))

    # 2) 离子回旋频率带（活跃时段内 |B| 的 25%-75% 分位对应区间），淡黄底色
    x_min, x_max = ax.get_xlim()
    x_min = min(x_min, f_ci_lo * 0.7)
    x_max = max(x_max, f_ci_hi * 1.3)
    ax.set_xlim(x_min, x_max)
    ax.axvspan(f_ci_lo, f_ci_hi, color='yellow', alpha=0.3,
               label='ion cyclotron frequency')
    ax.axvline(f_ci_lo, color='orange', linestyle='-', linewidth=1.0, alpha=0.9)
    ax.axvline(f_ci_hi, color='orange', linestyle='-', linewidth=1.0, alpha=0.9)

    norm_tag = ' (L1 normalized)' if L1_NORMALIZE else ''
    ax.set_title(f'PSD - {event_id}\n{date_str} {win_start:%H:%M}-{win_end:%H:%M} (2 h window)\n'
                 f'{label}{norm_tag}\nPeak {f_peak:.3e} Hz', fontsize=12)
    ax.legend(loc='upper right')
    out_file = os.path.join(out_dir,
                            f"{event_id}_{date_str}_{win_start:%H%M}_{safe}_PSD_peak.png")
    fig.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {out_file}")


def plot_wavelet_power(time_sec, time_utc, freqs, power, f_proton, f_peak, f_peak_idx,
                       t_start, t_end, thr, event_id, date_str, win_start, win_end,
                       label, safe, out_dir):
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(12, 9), sharex=True,
        gridspec_kw={'height_ratios': [3, 1]}, dpi=150)
    fig.subplots_adjust(hspace=0.08)

    pos_power = power[f_peak_idx, :]
    vmin = max(np.percentile(power[power > 0], 1), np.finfo(float).tiny)
    img = ax_top.pcolormesh(time_sec, freqs, power, shading='auto', cmap='jet',
                            norm=LogNorm(vmin=vmin, vmax=power.max()))
    ax_top.set_yscale('log')
    ax_top.set_ylim(freqs[0], freqs[-1])
    ax_top.set_ylabel('Frequency (Hz)', fontsize=11)
    ax_top.axhline(f_peak, color='red', linestyle='--', linewidth=1.5,
                   label=f'Peak {f_peak:.3e} Hz')
    ax_top.axvspan(t_start, t_end, alpha=0.25, color='yellow',
                   label=f'Active window ({t_start:.0f}-{t_end:.0f} s)')
    ax_top.plot(time_sec, f_proton, color='lime', linewidth=1.8, alpha=0.95,
                label='Ion cyclotron freq.')
    ax_top.legend(loc='upper right', fontsize=8)
    ax_top.set_title(f'Wavelet Power Spectrum - {event_id}\n'
                     f'{date_str} {win_start:%H:%M}-{win_end:%H:%M} (2 h window)\n'
                     f'{label} (peak {f_peak:.3e} Hz)', fontsize=12)
    cbar = fig.colorbar(img, ax=ax_top, fraction=0.025, pad=0.02)
    cbar.set_label('Power (nT$^2$)', fontsize=10)

    ax_bot.plot(time_sec, pos_power, 'b-', linewidth=1)
    ax_bot.axhline(thr, color='gray', linestyle=':', linewidth=1.2,
                   label=f'Threshold {POWER_THRESHOLD_FRAC}*max')
    ax_bot.axvspan(t_start, t_end, alpha=0.25, color='yellow')
    ax_bot.set_xlabel('Time (s)', fontsize=11)
    ax_bot.set_ylabel(f'Power at\n{f_peak:.3e} Hz\n(nT$^2$)', fontsize=9)
    ax_bot.legend(loc='best', fontsize=8)
    ax_bot.grid(True, alpha=0.3)

    utc_ax = ax_top.twiny()
    utc_ax.set_xlim(time_utc[0], time_utc[-1])
    utc_ax.xaxis.set_major_locator(AutoDateLocator(minticks=5, maxticks=10))
    utc_ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
    utc_ax.set_xlabel('UTC Time', fontsize=10)

    out_file = os.path.join(out_dir,
                            f"{event_id}_{date_str}_{win_start:%H%M}_{safe}_wavelet_power_window.png")
    fig.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {out_file}")


def plot_b_window(Bx, By, Bz, B_total, time_utc, f_proton, t_start_utc, t_end_utc,
                  means, event_id, date_str, win_start, win_end, label, safe, out_dir):
    fig, (ax, ax_ci) = plt.subplots(
        2, 1, figsize=(12, 8), sharex=True,
        gridspec_kw={'height_ratios': [2, 1]}, dpi=150)
    fig.subplots_adjust(hspace=0.08)

    ax.plot(time_utc, Bx.values, label='Bx', linewidth=0.8)
    ax.plot(time_utc, By.values, label='By', linewidth=0.8)
    ax.plot(time_utc, Bz.values, label='Bz', linewidth=0.8)
    ax.plot(time_utc, B_total.values, 'k-', label='|B|', linewidth=1)
    ax.axvspan(t_start_utc, t_end_utc, alpha=0.25, color='yellow',
               label=f'Active window\n{t_start_utc:%H:%M:%S}-{t_end_utc:%H:%M:%S}')
    ax.set_ylabel('B (nT)', fontsize=11)
    ax.legend(loc='best', fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
    txt = (f'Window mean:\n'
           f'Bx={means["Bx"]:.3f}, By={means["By"]:.3f}, Bz={means["Bz"]:.3f} nT\n'
           f'|B|={means["Btot"]:.3f}, B_par={means["Bpar"]:.3f}, B_perp={means["Bperp"]:.3f} nT')
    ax.text(0.02, 0.98, txt, transform=ax.transAxes, va='top', fontsize=9,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85))
    ax.set_title(f'B field & ion cyclotron freq - {event_id}\n'
                 f'{date_str} {win_start:%H:%M}-{win_end:%H:%M} (2 h window)\n{label}', fontsize=12)

    ax_ci.plot(time_utc, f_proton, color='green', linewidth=1.2,
               label='Ion cyclotron freq.')
    ax_ci.axvspan(t_start_utc, t_end_utc, alpha=0.25, color='yellow')
    ax_ci.set_ylabel('f$_{ci}$ (Hz)', fontsize=11)
    ax_ci.legend(loc='best', fontsize=8)
    ax_ci.grid(True, alpha=0.3)
    ax_ci.xaxis.set_major_formatter(DateFormatter('%H:%M'))
    ax_ci.set_xlabel('UTC Time', fontsize=11)

    fig.tight_layout()
    out_file = os.path.join(out_dir,
                            f"{event_id}_{date_str}_{win_start:%H%M}_{safe}_B_window_mean.png")
    fig.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {out_file}")


def main():
    os.makedirs(OUT_ROOT, exist_ok=True)
    events = parse_events(EVENT_FILE)
    print(f"Read {len(events)} events.")
    summary_rows = []

    for i, ev in enumerate(events, start=1):
        event_id = f"EV{i:02d}"
        date_str = ev['date']
        win_start, win_end = center_window(ev['start'], ev['end'])
        out_dir = os.path.join(OUT_ROOT, f"{date_str}_ev{i:02d}")
        os.makedirs(out_dir, exist_ok=True)

        print(f"\n{'='*70}\n{event_id} - {date_str}, window "
              f"{win_start} -> {win_end} (2 h)\n{'='*70}")

        win_start_str = win_start.strftime('%Y-%m-%d/%H:%M:%S')
        win_end_str = win_end.strftime('%Y-%m-%d/%H:%M:%S')

        try:
            Bx, By, Bz, px, py, pz = process_juno_data.process_juno_data(date_str, win_start_str, win_end_str)
            if Bx is None or len(Bx) == 0:
                print("No data. Skipping.")
                continue
            B_par, B_perp = compute_parallel_perp_by_window(Bx, By, Bz, WINDOW_SEC)
            B_total = np.sqrt(Bx**2 + By**2 + Bz**2)
            time_sec = (Bx.index - Bx.index.min()).total_seconds().values
            time_utc = Bx.index
            components = [(B_par, r'$B_{||}$', 'B_parallel'),
                          (B_perp, r'$B_{\perp}$', 'B_perp')]
        except Exception as e:
            print(f"Error: {e}")
            continue

        for data_series, label, safe in components:
            print(f"\n--- {label} ---")
            try:
                coef, amp, freqs, t_out, t_utc, f_proton = wavelet.analyze_magnetic_field(
                    data_series, B_total, label, freq_range=FREQ_RANGE
                )
            except Exception as e:
                print(f"Wavelet failed: {e}")
                continue

            # 1) 全时段 PSD：用 find_peaks 找出所有局部极大峰，
            #    再取其中能量最高的作为谱峰（避免直接取全局最大值在低频边缘）
            target_freqs, energy, _ = PSD.compute_frequency_energy_spectrum(
                t_out, data_series.values, freqs, coef, F_MIN, F_MAX,
                l1_normalize=L1_NORMALIZE
            )
            peak_cands = PSD.find_peaks(target_freqs, energy)
            if not peak_cands:
                print("No clear local peak (monotonic spectrum) - skipping this component.")
                continue
            cand_energies = [energy[np.argmin(np.abs(target_freqs - f))] for f in peak_cands]
            f_peak = float(peak_cands[int(np.argmax(cand_energies))])
            print(f"Detected local peaks (Hz): "
                  f"{['%.3e' % f for f in peak_cands]} -> highest peak {f_peak:.3e} Hz")

            # 2) 小波功率谱峰值频率处的功率随时间 -> 能量主要发生区域
            power = np.abs(coef)**2
            f_peak_idx = int(np.argmin(np.abs(freqs - f_peak)))
            pos_power = power[f_peak_idx, :]
            t_start, t_end, i_left, i_right = active_time_window(pos_power, t_out)
            thr = POWER_THRESHOLD_FRAC * pos_power.max()
            fci_win = f_proton[i_left:i_right + 1]
            f_ci_pct = np.percentile(fci_win, [25, 75])
            f_ci_lo, f_ci_hi = float(f_ci_pct[0]), float(f_ci_pct[1])
            f_ci_mean = float(np.mean(fci_win))
            t_start_utc = t_utc[0] + timedelta(seconds=float(t_start))
            t_end_utc = t_utc[0] + timedelta(seconds=float(t_end))

            # 3) 该时间窗内磁场平均
            mask = (t_out >= t_start) & (t_out <= t_end)
            Bx_m, By_m, Bz_m = Bx.values[mask].mean(), By.values[mask].mean(), Bz.values[mask].mean()
            Btot_m = B_total.values[mask].mean()
            Bpar_m = data_series.values[mask].mean() if safe == 'B_parallel' else B_par.values[mask].mean()
            Bperp_m = B_perp.values[mask].mean() if safe == 'B_perp' else data_series.values[mask].mean()
            means = {'Bx': Bx_m, 'By': By_m, 'Bz': Bz_m, 'Btot': Btot_m,
                     'Bpar': Bpar_m, 'Bperp': Bperp_m}

            print(f"Peak freq (highest power): {f_peak:.6e} Hz")
            print(f"Active window: {t_start:.0f}-{t_end:.0f} s "
                  f"({t_start_utc} - {t_end_utc}), duration {t_end-t_start:.0f} s")
            print(f"f_ci in window: P25={f_ci_lo:.6e}, P75={f_ci_hi:.6e}, mean={f_ci_mean:.6e} Hz")
            print(f"Window mean B: Bx={Bx_m:.3f}, By={By_m:.3f}, Bz={Bz_m:.3f}, |B|={Btot_m:.3f} nT")

            # 4) 图
            plot_psd_annotated(t_out, data_series, freqs, coef, f_peak, f_ci_lo, f_ci_hi,
                               peak_cands, event_id, date_str, win_start, win_end,
                               label, safe, out_dir)
            plot_wavelet_power(t_out, t_utc, freqs, power, f_proton, f_peak, f_peak_idx,
                               t_start, t_end, thr, event_id, date_str, win_start, win_end,
                               label, safe, out_dir)
            plot_b_window(Bx, By, Bz, B_total, t_utc, f_proton, t_start_utc, t_end_utc,
                          means, event_id, date_str, win_start, win_end, label, safe, out_dir)

            summary_rows.append({
                'Event': event_id, 'Date': date_str,
                'Window_Start': win_start, 'Window_End': win_end,
                'Component': safe,
                'Peak_Freq_Hz': f_peak,
                'Active_Start': t_start_utc, 'Active_End': t_end_utc,
                'Duration_s': t_end - t_start,
                'Mean_f_ci_Hz': f_ci_mean,
                'f_ci_P25_Hz': f_ci_lo, 'f_ci_P75_Hz': f_ci_hi,
                'Mean_Bx_nT': Bx_m, 'Mean_By_nT': By_m, 'Mean_Bz_nT': Bz_m,
                'Mean_Btot_nT': Btot_m, 'Mean_Bpar_nT': Bpar_m, 'Mean_Bperp_nT': Bperp_m,
            })

    if summary_rows:
        df = pd.DataFrame(summary_rows)
        excel_out = os.path.join(OUT_ROOT, 'current_2h_summary.xlsx')
        df.to_excel(excel_out, index=False)
        print(f"\nExcel saved: {excel_out}")

    print("\nAll done.")


if __name__ == "__main__":
    main()
