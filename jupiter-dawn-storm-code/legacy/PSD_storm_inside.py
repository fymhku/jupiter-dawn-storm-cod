import os
import numpy as np
import pandas as pd
from datetime import timedelta
from tool import process_juno_data
from tool import wavelet
from tool import PSD
from PSD_storm import (FREQ_RANGE, F_MIN, F_MAX, WINDOW_SEC, POWER_THRESHOLD_FRAC,
                       L1_NORMALIZE, HALF_WINDOW, EVENT_FILE,
                       parse_events, center_window,
                       compute_parallel_perp_by_window, active_time_window,
                       plot_psd_annotated, plot_wavelet_power, plot_b_window,
                       plot_cwt_image)

# 配置
OUT_ROOT = r"D:\summer research\pic\storm_psd_inside"
B_LOBE_HOURS = 10            # B_lobe: 事件开始后 N 小时内的最大 1h 平均 |B|
B_LOBE_FRAC = 0.9            # 阈值 = 0.9 * B_lobe
MIN_INSIDE_POINTS = 900      # 至少 15 分钟的连续"电流片内"数据才分析
SECOND_PEAK_DATES = {'2017_01_24'}   # 对这些日期的事件，额外分析第二大能量峰


def compute_B_lobe(Bx, By, Bz, hours=B_LOBE_HOURS):
    """在[开始, 开始+hours]内，滑动 1 小时窗求 |B| 平均，取最大值作为 B_lobe。"""
    Btot = pd.Series(np.sqrt(Bx**2 + By**2 + Bz**2), index=Bx.index)
    rolling = Btot.rolling('3600s', min_periods=1800).mean()
    lobe = rolling.max()
    if np.isnan(lobe) or lobe <= 0:
        lobe = Btot.mean()
    return float(lobe)


def longest_inside_run(inside):
    """返回最长的 True 连续段的 (start_idx, end_idx)；无则返回 None。"""
    inside = np.asarray(inside, dtype=bool)
    n = len(inside)
    best_len, best_start = 0, -1
    i = 0
    while i < n:
        if inside[i]:
            j = i
            while j < n and inside[j]:
                j += 1
            if (j - i) > best_len:
                best_len, best_start = j - i, i
            i = j
        else:
            i += 1
    if best_len == 0:
        return None
    return best_start, best_start + best_len - 1


def select_ranked_peak(peak_cands, target_freqs, energy, rank):
    """返回 energy 第 rank 大（rank=1 最大）对应的峰频率；峰不足 rank 条则返回 None。"""
    cand_e = [energy[np.argmin(np.abs(target_freqs - f))] for f in peak_cands]
    order = np.argsort(np.array(cand_e))[::-1]
    if len(order) < rank:
        return None
    return float(peak_cands[order[rank - 1]])


def process_peak(event_id, date_str, Bx, By, Bz, B_total, B_par, B_perp,
                 data_series, label, safe, coef, amp, freqs, t_out, t_utc, f_proton,
                 peak_cands, rank, tag, f_peak, use_start, use_end, out_dir,
                 B_lobe, threshold, inside_start_dt, inside_end_dt):
    """针对某一个峰频率，做活跃窗、平均磁场、出图，并返回一行汇总。"""
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

    mask = (t_out >= t_start) & (t_out <= t_end)
    Bx_m, By_m, Bz_m = Bx.values[mask].mean(), By.values[mask].mean(), Bz.values[mask].mean()
    Btot_m = B_total.values[mask].mean()
    Bpar_m = data_series.values[mask].mean() if safe == 'B_parallel' else B_par.values[mask].mean()
    Bperp_m = B_perp.values[mask].mean() if safe == 'B_perp' else data_series.values[mask].mean()
    means = {'Bx': Bx_m, 'By': By_m, 'Bz': Bz_m, 'Btot': Btot_m,
             'Bpar': Bpar_m, 'Bperp': Bperp_m}

    print(f"Peak (rank {rank}): {f_peak:.6e} Hz")
    print(f"Active window: {t_start:.0f}-{t_end:.0f} s "
          f"({t_start_utc} - {t_end_utc}), duration {t_end-t_start:.0f} s")
    print(f"f_ci in window: P25={f_ci_lo:.6e}, P75={f_ci_hi:.6e}, mean={f_ci_mean:.6e} Hz")
    print(f"Window mean B: Bx={Bx_m:.3f}, By={By_m:.3f}, Bz={Bz_m:.3f}, |B|={Btot_m:.3f} nT")

    plot_psd_annotated(t_out, data_series, freqs, coef, f_peak, f_ci_lo, f_ci_hi,
                       peak_cands, event_id, date_str, use_start, use_end,
                       label, safe, out_dir, tag)
    plot_wavelet_power(t_out, t_utc, freqs, power, f_proton, f_peak, f_peak_idx,
                       t_start, t_end, thr, event_id, date_str, use_start, use_end,
                       label, safe, out_dir, tag)
    plot_b_window(Bx, By, Bz, B_total, t_utc, f_proton, t_start_utc, t_end_utc,
                  means, event_id, date_str, use_start, use_end, label, safe, out_dir, tag)
    plot_cwt_image(amp, freqs, t_out, t_utc, f_proton, event_id, date_str,
                   use_start, use_end, label, safe, out_dir, tag)

    return {
        'Event': event_id, 'Date': date_str,
        'B_lobe_nT': B_lobe, 'Threshold_nT': threshold,
        'Inside_Start': inside_start_dt, 'Inside_End': inside_end_dt,
        'Inside_Duration_s': (inside_end_dt - inside_start_dt).total_seconds(),
        'Component': safe, 'Peak_Rank': rank, 'Peak_Freq_Hz': f_peak,
        'Active_Start': t_start_utc, 'Active_End': t_end_utc,
        'Duration_s': t_end - t_start,
        'Mean_f_ci_Hz': f_ci_mean,
        'f_ci_P25_Hz': f_ci_lo, 'f_ci_P75_Hz': f_ci_hi,
        'Mean_Bx_nT': Bx_m, 'Mean_By_nT': By_m, 'Mean_Bz_nT': Bz_m,
        'Mean_Btot_nT': Btot_m, 'Mean_Bpar_nT': Bpar_m, 'Mean_Bperp_nT': Bperp_m,
    }


def main():
    os.makedirs(OUT_ROOT, exist_ok=True)
    events = parse_events(EVENT_FILE)
    print(f"Read {len(events)} events.")
    summary_rows = []

    for i, ev in enumerate(events, start=1):
        event_id = f"EV{i:02d}"
        date_str = ev['date']
        start_dt, end_dt = ev['start'], ev['end']
        win_start, win_end = center_window(start_dt, end_dt)
        out_dir = os.path.join(OUT_ROOT, f"{date_str}_ev{i:02d}")
        os.makedirs(out_dir, exist_ok=True)

        print(f"\n{'='*70}\n{event_id} - {date_str}, "
              f"event {start_dt} -> {end_dt}, window {win_start} -> {win_end}\n{'='*70}")

        # ---- B_lobe: 事件开始后 B_LOBE_HOURS 内最大 1h 平均 |B| ----
        lobe_end = start_dt + timedelta(hours=B_LOBE_HOURS)
        try:
            Bx_l, By_l, Bz_l, _, _, _ = process_juno_data.process_juno_data(
                date_str,
                start_dt.strftime('%Y-%m-%d/%H:%M:%S'),
                lobe_end.strftime('%Y-%m-%d/%H:%M:%S'))
            if Bx_l is None or len(Bx_l) == 0:
                print("No B_lobe data. Skipping.")
                continue
            B_lobe = compute_B_lobe(Bx_l, By_l, Bz_l)
        except Exception as e:
            print(f"B_lobe error: {e}")
            continue
        threshold = B_LOBE_FRAC * B_lobe
        print(f"B_lobe = {B_lobe:.3f} nT, threshold (0.9*B_lobe) = {threshold:.3f} nT")

        # ---- 载入原始分析窗 ----
        win_start_str = win_start.strftime('%Y-%m-%d/%H:%M:%S')
        win_end_str = win_end.strftime('%Y-%m-%d/%H:%M:%S')
        try:
            Bx, By, Bz, px, py, pz = process_juno_data.process_juno_data(date_str, win_start_str, win_end_str)
            if Bx is None or len(Bx) == 0:
                print("No data. Skipping.")
                continue
        except Exception as e:
            print(f"Error: {e}")
            continue

        Btot = np.sqrt(Bx**2 + By**2 + Bz**2)
        inside = Btot.values < threshold
        run = longest_inside_run(inside)
        if run is None:
            print("No inside-current-sheet interval (|B| >= threshold throughout). Skipping.")
            continue
        i0, i1 = run
        if (i1 - i0 + 1) < MIN_INSIDE_POINTS:
            print(f"Inside-current-sheet interval too short ({i1 - i0 + 1} pts "
                  f"= {(i1 - i0 + 1) / 60:.1f} min). Skipping.")
            continue

        # 只保留最长的连续"电流片内"区间
        sub_idx = Bx.index[i0:i1 + 1]
        Bx, By, Bz = Bx.loc[sub_idx], By.loc[sub_idx], Bz.loc[sub_idx]
        inside_start_dt = sub_idx[0]
        inside_end_dt = sub_idx[-1]
        print(f"Inside current sheet: {inside_start_dt} -> {inside_end_dt}, "
              f"duration {(inside_end_dt - inside_start_dt).total_seconds() / 60:.1f} min")

        B_par, B_perp = compute_parallel_perp_by_window(Bx, By, Bz, WINDOW_SEC)
        B_total = np.sqrt(Bx**2 + By**2 + Bz**2)   # Series, 索引与 Bx 一致
        time_utc = sub_idx
        use_start, use_end = inside_start_dt, inside_end_dt
        components = [(B_par, r'$B_{||}$', 'B_parallel'),
                      (B_perp, r'$B_{\perp}$', 'B_perp')]

        for data_series, label, safe in components:
            print(f"\n--- {label} ---")
            try:
                coef, amp, freqs, t_out, t_utc, f_proton = wavelet.analyze_magnetic_field(
                    data_series, B_total, label, freq_range=FREQ_RANGE
                )
            except Exception as e:
                print(f"Wavelet failed: {e}")
                continue

            target_freqs, energy, _ = PSD.compute_frequency_energy_spectrum(
                t_out, data_series.values, freqs, coef, F_MIN, F_MAX,
                l1_normalize=L1_NORMALIZE
            )
            peak_cands = PSD.find_peaks(target_freqs, energy)
            if not peak_cands:
                print("No clear local peak (monotonic spectrum) - skipping this component.")
                continue
            print(f"Detected local peaks (Hz): {['%.3e' % f for f in peak_cands]}")

            ranks = [1] if date_str not in SECOND_PEAK_DATES else [1, 2]
            for rank in ranks:
                f_peak = select_ranked_peak(peak_cands, target_freqs, energy, rank)
                if f_peak is None:
                    print(f"No rank-{rank} peak (only {len(peak_cands)} peaks).")
                    continue
                tag = '' if rank == 1 else '_p2'
                print(f"Processing rank-{rank} peak: {f_peak:.6e} Hz")
                summary_rows.append(process_peak(
                    event_id, date_str, Bx, By, Bz, B_total, B_par, B_perp,
                    data_series, label, safe, coef, amp, freqs, t_out, t_utc, f_proton,
                    peak_cands, rank, tag, f_peak, use_start, use_end, out_dir,
                    B_lobe, threshold, inside_start_dt, inside_end_dt))

    if summary_rows:
        df = pd.DataFrame(summary_rows)
        excel_out = os.path.join(OUT_ROOT, 'storm_inside_summary.xlsx')
        df.to_excel(excel_out, index=False)
        print(f"\nExcel saved: {excel_out}")

    print("\nAll done.")


if __name__ == "__main__":
    main()
