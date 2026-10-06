import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
from datetime import timedelta
from tool import process_juno_data, wavelet, PSD
from PSD_storm import compute_parallel_perp_by_window, active_time_window
import warnings
warnings.filterwarnings('ignore')

# B field (RTP) 随时间曲线: 只画 Br / Btheta / Bphi, 不画 |B|,
# 峰值活跃阶段用 Bperp 谱峰定, 淡黄色标记
FREQ_RANGE = (3e-3, 5e-2)
WINDOW_SEC = 20 * 60
L1_NORMALIZE = False
POWER_THRESHOLD_FRAC = 0.8
OUT_ROOT = r"D:\summer research\pic\storm_psd_3e-3"

# (event_id, date_str, 窗口起, 窗口止, 固定活跃时段(起,止) 或 None, 输出日期文件夹)
# 活跃时段用垂直分量带通/PSD 得到的新结果:
#   2017_03_17: 10:13:29-10:28:13
#   2017_05_16: 09:17:38-09:46:52 (主波包前后各放宽 10 min)
EVENTS = [
    ('EV02', '2017_03_17', '2017-03-17/08:56:30', '2017-03-17/10:56:30',
     ('2017-03-17/10:13:29', '2017-03-17/10:28:13'),
     os.path.join(OUT_ROOT, '2017_03_17_ev07_star')),
    ('EV03', '2017_05_16', '2017-05-16/07:52:00', '2017-05-16/10:03:00',
     ('2017-05-16/09:17:38', '2017-05-16/09:46:52'),
     os.path.join(OUT_ROOT, '2017_05_16_ev14')),
]


def plot_rtp_field(Br, Bth, Bph, t_start_utc, t_end_utc, f_peak,
                   event_id, date_str, win_start, win_end, out_dir):
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(Br.index, Br.values, label='$B_r$', linewidth=0.8)
    ax.plot(Bth.index, Bth.values, label='$B_\\theta$', linewidth=0.8)
    ax.plot(Bph.index, Bph.values, label='$B_\\phi$', linewidth=0.8)
    if t_start_utc is not None:
        ax.axvspan(t_start_utc, t_end_utc, alpha=0.25, color='yellow',
                   label='Active period')
    ax.set_ylabel('B (nT)', fontsize=11)
    ax.set_xlabel('UTC Time', fontsize=11)
    ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
    ax.set_title(f'B field (RTP) - {event_id}  {date_str.replace("_", "-")}\n'
                 f'{win_start:%Y-%m-%d %H:%M:%S} - {win_end:%H:%M:%S}',
                 fontsize=12)
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out_file = os.path.join(out_dir,
                            f"{event_id}_{date_str}_{win_start:%H%M}_B_RTP_field.png")
    fig.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {out_file}")


def main():
    for ev in EVENTS:
        event_id, date_str, win_start_str, win_end_str, fixed_active, out_dir = ev
        os.makedirs(out_dir, exist_ok=True)
        print(f"\n=== {event_id} {date_str} {win_start_str} -> {win_end_str} ===")
        Bx, By, Bz, _, _, _ = process_juno_data.process_juno_data(date_str, win_start_str, win_end_str)
        if Bx is None or len(Bx) == 0:
            print("No data. Skipping.")
            continue

        Br, Bth, Bph = Bx, By, Bz          # JSS-RTP 产品: 已是 r/theta/phi
        t0_data, t1_data = Bx.index[0], Bx.index[-1]

        # --- 固定活跃时段: 直接用新得到的活跃区标注，跳过自动检测 ---
        if fixed_active is not None:
            t_start_utc = pd.Timestamp(fixed_active[0])
            t_end_utc = pd.Timestamp(fixed_active[1])
            print(f"Using fixed active period: "
                  f"{t_start_utc:%Y-%m-%d %H:%M:%S} - {t_end_utc:%Y-%m-%d %H:%M:%S}")
            plot_rtp_field(Br, Bth, Bph, t_start_utc, t_end_utc, None,
                           event_id, date_str, t0_data, t1_data, out_dir)
            continue

        # --- 未指定固定活跃时段: 保留原自动检测方案 ---
        B_total = np.sqrt(Bx**2 + By**2 + Bz**2)
        _, B_perp = compute_parallel_perp_by_window(Bx, By, Bz, WINDOW_SEC)

        label = r'$B_\perp$'
        coef, amp, freqs, t_out, t_utc, f_proton = wavelet.analyze_magnetic_field(
            B_perp, B_total, label, freq_range=FREQ_RANGE)

        target_freqs, energy, _ = PSD.compute_frequency_energy_spectrum(
            t_out, B_perp.values, freqs, coef, FREQ_RANGE[0], FREQ_RANGE[1],
            l1_normalize=L1_NORMALIZE)
        peak_cands = PSD.find_peaks(target_freqs, energy)
        if not peak_cands:
            print("No clear local peak in B_perp (monotonic) - plot field without marking.")
            plot_rtp_field(Br, Bth, Bph, None, None, None,
                           event_id, date_str, t_utc[0], t_utc[-1], out_dir)
            continue

        cand_energies = [energy[np.argmin(np.abs(target_freqs - f))] for f in peak_cands]
        f_peak = float(peak_cands[int(np.argmax(cand_energies))])
        power = np.abs(coef)**2
        f_peak_idx = int(np.argmin(np.abs(freqs - f_peak)))
        pos_power = power[f_peak_idx, :]
        t_start, t_end, _, _ = active_time_window(pos_power, t_out,
                                                  frac=POWER_THRESHOLD_FRAC)
        t_start_utc = t_utc[0] + timedelta(seconds=float(t_start))
        t_end_utc = t_utc[0] + timedelta(seconds=float(t_end))
        print(f"B_perp peak {f_peak:.3e} Hz, active window "
              f"{t_start_utc:%H:%M:%S} - {t_end_utc:%H:%M:%S}")

        plot_rtp_field(Br, Bth, Bph, t_start_utc, t_end_utc, f_peak,
                       event_id, date_str, t_utc[0], t_utc[-1], out_dir)


if __name__ == "__main__":
    main()
