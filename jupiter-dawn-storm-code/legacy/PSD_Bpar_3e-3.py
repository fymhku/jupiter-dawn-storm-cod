import os
import numpy as np
import matplotlib.pyplot as plt
from tool import process_juno_data, wavelet, PSD
from PSD_storm import compute_parallel_perp_by_window
import warnings
warnings.filterwarnings('ignore')

# B_parallel PSD, 3e-3 ~ 5e-2, 真实呈现（不标离子回旋、不寻峰）
FREQ_RANGE = (3e-3, 5e-2)
WINDOW_SEC = 20 * 60
OUT_ROOT = r"D:\summer research\pic\storm_psd_3e-3"

# (event_id, date_str, 窗口起, 窗口止, 输出日期文件夹)
EVENTS = [
    ('EV02', '2017_03_17', '2017-03-17/08:56:30', '2017-03-17/10:56:30',
     os.path.join(OUT_ROOT, '2017_03_17_ev07_star')),
    ('EV03', '2017_05_16', '2017-05-16/07:52:00', '2017-05-16/10:03:00',
     os.path.join(OUT_ROOT, '2017_05_16_ev14')),
]


def plot_bpar_psd(target_freqs, energy, event_id, date_str, win_start, win_end, out_dir):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.semilogx(target_freqs, energy, 'b-', linewidth=2, label=r'$B_{||}$ energy density')
    ax.set_xscale('log')
    ax.set_xlabel('Frequency (Hz)', fontsize=12)
    ax.set_ylabel('Magnetic Energy (nT$^2$/Hz)', fontsize=12)
    bpar_label = r'$B_{||}$'
    ax.set_title(f'PSD - {event_id}  {date_str.replace("_", "-")}\n'
                 f'{win_start:%Y-%m-%d %H:%M} - {win_end:%H:%M}\n'
                 f'{bpar_label}  (freq range {FREQ_RANGE[0]:.1e} - {FREQ_RANGE[1]:.1e} Hz)',
                 fontsize=12)
    ax.grid(True, which='both', linestyle='--', alpha=0.6)
    ax.legend(loc='upper right', fontsize=10)
    fig.tight_layout()

    out_file = os.path.join(
        out_dir,
        f"{event_id}_{date_str}_{win_start:%H%M}_B_parallel_PSD.png")
    fig.savefig(out_file, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {out_file}")


def main():
    for ev in EVENTS:
        event_id, date_str, win_start_str, win_end_str, out_dir = ev
        os.makedirs(out_dir, exist_ok=True)
        print(f"\n=== {event_id} {date_str} {win_start_str} -> {win_end_str} ===")
        Bx, By, Bz, _, _, _ = process_juno_data.process_juno_data(date_str, win_start_str, win_end_str)
        if Bx is None or len(Bx) == 0:
            print("No data. Skipping.")
            continue
        B_par, _ = compute_parallel_perp_by_window(Bx, By, Bz, WINDOW_SEC)
        B_total = np.sqrt(Bx**2 + By**2 + Bz**2)

        coef, amp, freqs, t_out, t_utc, f_proton = wavelet.analyze_magnetic_field(
            B_par, B_total, r'$B_{||}$', freq_range=FREQ_RANGE)
        target_freqs, energy, _ = PSD.compute_frequency_energy_spectrum(
            t_out, B_par.values, freqs, coef, FREQ_RANGE[0], FREQ_RANGE[1])

        win_start_dt = t_utc[0]
        win_end_dt = t_utc[-1]
        plot_bpar_psd(target_freqs, energy, event_id, date_str,
                      win_start_dt, win_end_dt, out_dir)


if __name__ == "__main__":
    main()
