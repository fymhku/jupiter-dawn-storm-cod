import os
import numpy as np
import pandas as pd
from datetime import timedelta, datetime
from tool import process_juno_data
from tool import wavelet
from tool import PSD
from PSD_storm import (parse_events, center_window,
                       compute_parallel_perp_by_window, active_time_window,
                       plot_psd_annotated, plot_wavelet_power, plot_b_window,
                       plot_cwt_image)

# 配置：CWT 频率范围改为 1e-3 ~ 5e-2 Hz
L1_NORMALIZE = False
WINDOW_SEC = 20 * 60
POWER_THRESHOLD_FRAC = 0.8
HALF_WINDOW = timedelta(hours=1)          # 从中心点前后各 1h -> 总长 2h
FREQ_RANGE = [1e-3, 5e-2]                  # CWT 频谱范围（按用户要求）
F_MIN, F_MAX = 1e-3, 5e-2                  # PSD 能量谱范围与 CWT 保持一致
EVENT_FILE = r"D:\summer research\计算曲率相关\py_code\storm.txt"
OUT_ROOT = r"D:\summer research\pic\storm_psd_cwt_1e-3"
TARGET_DATES = {'2017_03_17', '2017_05_16'}


def main():
    os.makedirs(OUT_ROOT, exist_ok=True)
    events = parse_events(EVENT_FILE)
    events = [ev for ev in events if ev['date'] in TARGET_DATES]
    print(f"Read {len(events)} events for dates {sorted(TARGET_DATES)}.")
    summary_rows = []

    for i, ev in enumerate(events, start=1):
        event_id = f"EV{i:02d}"
        date_str = ev['date']
        win_start, win_end = center_window(ev['start'], ev['end'])
        out_dir = os.path.join(OUT_ROOT, f"{date_str}_ev{i:02d}")
        os.makedirs(out_dir, exist_ok=True)

        print(f"\n{'='*70}\n{event_id} - {date_str}, event "
              f"{ev['start']} -> {ev['end']}, window {win_start} -> {win_end}\n{'='*70}")

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

            print(f"Peak freq (highest power): {f_peak:.6e} Hz")
            print(f"Active window: {t_start:.0f}-{t_end:.0f} s "
                  f"({t_start_utc} - {t_end_utc}), duration {t_end-t_start:.0f} s")
            print(f"f_ci in window: P25={f_ci_lo:.6e}, P75={f_ci_hi:.6e}, mean={f_ci_mean:.6e} Hz")
            print(f"Window mean B: Bx={Bx_m:.3f}, By={By_m:.3f}, Bz={Bz_m:.3f}, |B|={Btot_m:.3f} nT")

            plot_psd_annotated(t_out, data_series, freqs, coef, f_peak, f_ci_lo, f_ci_hi,
                               peak_cands, event_id, date_str, win_start, win_end,
                               label, safe, out_dir)
            plot_wavelet_power(t_out, t_utc, freqs, power, f_proton, f_peak, f_peak_idx,
                               t_start, t_end, thr, event_id, date_str, win_start, win_end,
                               label, safe, out_dir)
            plot_b_window(Bx, By, Bz, B_total, t_utc, f_proton, t_start_utc, t_end_utc,
                          means, event_id, date_str, win_start, win_end, label, safe, out_dir)
            plot_cwt_image(amp, freqs, t_out, t_utc, f_proton, event_id, date_str,
                           win_start, win_end, label, safe, out_dir)

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
        excel_out = os.path.join(OUT_ROOT, 'storm_cwt_summary.xlsx')
        df.to_excel(excel_out, index=False)
        print(f"\nExcel saved: {excel_out}")

    print("\nAll done.")


if __name__ == "__main__":
    main()
