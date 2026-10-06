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
POWER_THRESHOLD_FRAC = 0.8   # 功率 ≥ 0.8*max 的连续区间视为能量主要发生区域

# 要处理的事件: (event_id, date_str, start_time, end_time, output_dir)
EVENTS = [
    ('PJ03', '2016_12_04', '10:00:00', '12:00:00', r"D:\summer research\pic\16.12.4_noFCI"),
    ('PJ03', '2016_12_05', '15:30:00', '17:30:00', r"D:\summer research\pic\16.12.5_noFCI"),
]


def compute_parallel_perp_by_window(Bx, By, Bz, window_sec):
    """同前，返回平行/垂直 Series"""
    idx = Bx.index
    if not isinstance(idx, pd.DatetimeIndex):
        raise TypeError("Index must be DatetimeIndex")
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
    """返回功率 ≥ frac*max 的连续区间（从全局峰值向两侧扩展），单位为秒。"""
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


def plot_wavelet_power(time_sec, time_utc, freqs, power, f_proton, f_peak, f_peak_idx,
                       t_start, t_end, thr, event_id, date_str, start_time, end_time,
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
                   label=f'Active window\n({t_start:.0f}-{t_end:.0f} s)')
    ax_top.plot(time_sec, f_proton, 'w-', linewidth=1, alpha=0.9,
                label='Ion cyclotron freq.')
    ax_top.legend(loc='upper right', fontsize=8)
    ax_top.set_title(f'Wavelet Power Spectrum - {event_id}\n'
                     f'{date_str} {start_time}-{end_time}\n'
                     f'{label} (peak freq {f_peak:.3e} Hz)', fontsize=12)
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

    safe_date = date_str.replace('_', '-')
    safe_start = start_time.replace(':', '-')
    safe_end = end_time.replace(':', '-')
    out_file = os.path.join(out_dir,
                            f"{event_id}_{safe_date}_{safe_start}_{safe_end}_{safe}_wavelet_power_window.png")
    fig.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {out_file}")


def plot_b_window(Bx, By, Bz, B_total, time_utc, t_start_utc, t_end_utc,
                  means, event_id, date_str, start_time, end_time, label, safe, out_dir):
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(time_utc, Bx.values, label='Bx', linewidth=0.8)
    ax.plot(time_utc, By.values, label='By', linewidth=0.8)
    ax.plot(time_utc, Bz.values, label='Bz', linewidth=0.8)
    ax.plot(time_utc, B_total.values, 'k-', label='|B|', linewidth=1)
    ax.axvspan(t_start_utc, t_end_utc, alpha=0.25, color='yellow',
               label=f'Active window\n{t_start_utc:%H:%M:%S}-{t_end_utc:%H:%M:%S}')
    ax.set_xlabel('UTC Time', fontsize=11)
    ax.set_ylabel('B (nT)', fontsize=11)
    ax.legend(loc='best', fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
    txt = (f'Window mean:\n'
           f'Bx={means["Bx"]:.3f}, By={means["By"]:.3f}, Bz={means["Bz"]:.3f} nT\n'
           f'|B|={means["Btot"]:.3f}, B_par={means["Bpar"]:.3f}, B_perp={means["Bperp"]:.3f} nT')
    ax.text(0.02, 0.98, txt, transform=ax.transAxes, va='top', fontsize=9,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85))
    ax.set_title(f'B field in active window - {event_id}\n'
                 f'{date_str} {start_time}-{end_time}\n{label}', fontsize=12)
    fig.tight_layout()
    safe_date = date_str.replace('_', '-')
    safe_start = start_time.replace(':', '-')
    safe_end = end_time.replace(':', '-')
    out_file = os.path.join(out_dir,
                            f"{event_id}_{safe_date}_{safe_start}_{safe_end}_{safe}_B_window_mean.png")
    fig.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {out_file}")


summary_rows = []

for event_id, date_str, start_time, end_time, OUTPUT_DIR in EVENTS:
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print(f"\n{'='*70}\nProcessing {event_id} - {date_str} {start_time} to {end_time}\n{'='*70}")

    date_fmt = date_str.replace('_', '-')
    start_full = f"{date_fmt}/{start_time}"
    end_full = f"{date_fmt}/{end_time}"

    try:
        Bx, By, Bz, px, py, pz = process_juno_data.process_juno_data(date_str, start_full, end_full)
        if Bx is None or len(Bx) == 0:
            print(f"No data for {event_id}")
            continue

        B_par, B_perp = compute_parallel_perp_by_window(Bx, By, Bz, WINDOW_SEC)
        B_total = np.sqrt(Bx**2 + By**2 + Bz**2)   # Series

        t0 = Bx.index.min()
        time_sec = (Bx.index - t0).total_seconds().values
        time_utc = Bx.index

        components = [(B_par, r'$B_{||}$', 'B_parallel'),
                      (B_perp, r'$B_{\perp}$', 'B_perp')]
    except Exception as e:
        print(f"Error: {e}")
        continue

    for data_series, label, safe in components:
        print(f"\nWavelet on {label}...")
        try:
            coef, amp, freqs, t_out, t_utc, f_proton = wavelet.analyze_magnetic_field(
                data_series, B_total, label, freq_range=FREQ_RANGE
            )
        except Exception as e:
            print(f"Failed: {e}")
            continue

        # 1) 全时段 PSD，取功率能量最高的单一频率
        target_freqs, energy, f_peak = PSD.compute_frequency_energy_spectrum(
            t_out, data_series.values, freqs, coef, F_MIN, F_MAX,
            l1_normalize=L1_NORMALIZE
        )

        # 2) 小波功率谱：峰值频率处的功率随时间变化 → 定位能量主要发生的区域
        power = np.abs(coef)**2
        f_peak_idx = int(np.argmin(np.abs(freqs - f_peak)))
        pos_power = power[f_peak_idx, :]
        t_start, t_end, i_left, i_right = active_time_window(pos_power, t_out)
        f_ci_mean = float(np.mean(f_proton[i_left:i_right + 1]))
        thr = POWER_THRESHOLD_FRAC * pos_power.max()
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
        print(f"Mean ion cyclotron freq in window: {f_ci_mean:.6e} Hz")
        print(f"Window mean B: Bx={Bx_m:.3f}, By={By_m:.3f}, Bz={Bz_m:.3f}, |B|={Btot_m:.3f} nT")

        # 4) PSD 曲线图（PSD_perp.py 原图）上把该频率标注出来
        fig, peaks = PSD.plot_energy_curve(
            t_out, data_series, freqs, coef,
            f_min=F_MIN, f_max=F_MAX,
            data_units='nT', plot_title='',
            output_file=None, l1_normalize=L1_NORMALIZE
        )
        ax = fig.axes[0]
        y_min, y_max = ax.get_ylim()
        x_min, x_max = ax.get_xlim()

        lines = ax.get_lines()
        if lines:
            x_data = lines[0].get_xdata()
            y_data = lines[0].get_ydata()

            # 4a) 保留原有 find_peaks 的峰值竖线
            if not isinstance(peaks, (list, np.ndarray)):
                peaks_list = [peaks] if peaks is not None else []
            else:
                peaks_list = list(peaks)
            for f in peaks_list:
                idx = np.argmin(np.abs(x_data - f))
                energy_v = y_data[idx]
                ax.axvline(x=f, color='red', linestyle='--', linewidth=1.2, alpha=0.8)
                ax.text(f+0.0001, y_max * 0.85, f'{f:.3e} Hz',
                        rotation=0, ha='left', va='center', fontsize=8, color='red',
                        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))

            # 4b) 仅在 16.12.5 B_perp 的 PSD 图上标注活跃时间窗内的平均离子回旋频率（仅红色实线，不带数值）
            if event_id == 'PJ03' and date_str == '2016_12_05' and safe == 'B_perp':
                ax.axvline(x=f_ci_mean, color='red', linestyle='-', linewidth=1.6,
                           alpha=0.9, label='Ion cyclotron freq (window mean)')

        norm_tag = ' (L1 normalized)' if L1_NORMALIZE else ''
        ax.set_title(f'Magnetic Energy - {event_id}\n{date_str} {start_time}-{end_time}\n'
                     f'{label}{norm_tag}', fontsize=12)
        ax.legend(loc='best')

        safe_date = date_str.replace('_', '-')
        safe_start = start_time.replace(':', '-')
        safe_end = end_time.replace(':', '-')
        out_file = os.path.join(OUTPUT_DIR,
                                f"{event_id}_{safe_date}_{safe_start}_{safe_end}_{safe}_peaks.png")
        fig.savefig(out_file, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved {out_file}")

        # 5) 辅助图：小波功率谱 + 时间窗、磁场窗口平均
        plot_wavelet_power(t_out, t_utc, freqs, power, f_proton, f_peak, f_peak_idx,
                           t_start, t_end, thr, event_id, date_str, start_time,
                           end_time, label, safe, OUTPUT_DIR)
        plot_b_window(Bx, By, Bz, B_total, t_utc, t_start_utc, t_end_utc, means,
                      event_id, date_str, start_time, end_time, label, safe, OUTPUT_DIR)

        summary_rows.append({
            'Event': f"{event_id}:{date_str} {start_time} {end_time}",
            'Component': safe,
            'HighestPower_Freq_Hz': f_peak,
            'Active_Start_UTC': t_start_utc, 'Active_End_UTC': t_end_utc,
            'Duration_s': t_end - t_start,
            'Mean_Bx_nT': Bx_m, 'Mean_By_nT': By_m, 'Mean_Bz_nT': Bz_m,
            'Mean_Btot_nT': Btot_m, 'Mean_Bpar_nT': Bpar_m, 'Mean_Bperp_nT': Bperp_m,
        })

# 保存 Excel（每个事件写入对应输出目录，含两个分量）
if summary_rows:
    df = pd.DataFrame(summary_rows)
    for event_id, date_str, start_time, end_time, out_dir in EVENTS:
        key = f"{event_id}:{date_str} {start_time} {end_time}"
        event_df = df[df['Event'] == key]
        if event_df.empty:
            continue
        excel_out = os.path.join(out_dir, 'PSD_peak_window_summary.xlsx')
        event_df.to_excel(excel_out, index=False)
        print(f"\nExcel saved: {excel_out}")

print("\nAll done.")
