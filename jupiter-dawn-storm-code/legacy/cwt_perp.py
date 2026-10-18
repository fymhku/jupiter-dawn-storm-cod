import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tool import process_juno_data
from tool import wavelet
from datetime import timedelta

# 配置
INPUT_FILE = r"D:\summer research\pic\aurora date_modify_0.9.txt"
OUTPUT_DIR = r"D:\summer research\pic\aurora pic_0.9"
FREQ_RANGE = [0.0001, 0.05]
WINDOW_SEC = 20 * 60   # 20分钟

os.makedirs(OUTPUT_DIR, exist_ok=True)

def compute_parallel_perp_by_window(Bx, By, Bz, window_sec):
    """
    Bx, By, Bz : pandas Series, 索引为 UTC datetime
    返回: B_parallel, B_perp (同索引的 Series)
    """
    idx = Bx.index
    if not isinstance(idx, pd.DatetimeIndex):
        raise TypeError("索引必须是 DatetimeIndex")
    
    B_par = pd.Series(index=idx, dtype=float)
    B_perp = pd.Series(index=idx, dtype=float)
    
    t_start = idx.min()
    t_end = idx.max()
    current = t_start
    while current <= t_end:
        window_end = current + timedelta(seconds=window_sec)
        mask = (idx >= current) & (idx < window_end)
        win_idx = idx[mask]
        if len(win_idx) == 0:
            current = window_end
            continue
        
        # 窗口平均磁场
        Bx_win = Bx.loc[win_idx].mean()
        By_win = By.loc[win_idx].mean()
        Bz_win = Bz.loc[win_idx].mean()
        B_norm = np.sqrt(Bx_win**2 + By_win**2 + Bz_win**2)
        if B_norm < 1e-12:
            b0 = np.array([0, 0, 1])
        else:
            b0 = np.array([Bx_win, By_win, Bz_win]) / B_norm
        
        # 窗口内数据
        Bx_vals = Bx.loc[win_idx].values
        By_vals = By.loc[win_idx].values
        Bz_vals = Bz.loc[win_idx].values
        
        B_par_vals = Bx_vals*b0[0] + By_vals*b0[1] + Bz_vals*b0[2]
        B_tot_vals = np.sqrt(Bx_vals**2 + By_vals**2 + Bz_vals**2)
        B_perp_vals = np.sqrt(B_tot_vals**2 - B_par_vals**2)
        
        B_par.loc[win_idx] = B_par_vals
        B_perp.loc[win_idx] = B_perp_vals
        
        current = window_end
    
    return B_par, B_perp

with open(INPUT_FILE, 'r') as f:
    lines = f.readlines()

for line in lines:
    if ':' not in line:
        continue
    event_part, datetime_part = line.strip().split(':', 1)
    event_id = event_part.strip()
    parts = datetime_part.split()
    if len(parts) < 3:
        continue
    date_str, start_time, end_time = parts[0], parts[1], parts[2]

    print(f"\n{'='*60}\nProcessing {event_id} - {date_str} {start_time} to {end_time}\n{'='*60}")
    event_dir = os.path.join(OUTPUT_DIR, event_id)
    os.makedirs(event_dir, exist_ok=True)

    date_fmt = date_str.replace('_', '-')
    start_full = f"{date_fmt}/{start_time}"
    end_full = f"{date_fmt}/{end_time}"

    try:
        Bx, By, Bz, px, py, pz = process_juno_data.process_juno_data(date_str, start_full, end_full)
        if Bx is None or len(Bx) == 0:
            print(f"No data for {event_id}")
            continue
        
        # 计算滑动窗口平行/垂直分量
        B_par, B_perp = compute_parallel_perp_by_window(Bx, By, Bz, WINDOW_SEC)
        B_total = np.sqrt(Bx**2 + By**2 + Bz**2)   # 仍然是 Series
        
        # 相对时间（秒）
        t0 = Bx.index.min()
        time_sec = (Bx.index - t0).total_seconds().values
        
        print(f"Data points: {len(Bx)}, duration: {time_sec[-1]-time_sec[0]:.1f} sec")
    except Exception as e:
        print(f"Error: {e}")
        continue

    # 分析两个分量：直接传入 Series，不要用 .values
    for data_series, label, safe in [(B_par, r'$B_{||}$ ', 'B_parallel'),
                                     (B_perp, r'$B_{\perp}$', 'B_perp')]:
        print(f"\nWavelet on {label}...")
        try:
            coef, amp, freqs, t_out, t_utc, f_proton = wavelet.analyze_magnetic_field(
                data_series,      # 传入 Series
                B_total,          # 传入 Series
                label,
                freq_range=FREQ_RANGE
            )
        except Exception as e:
            print(f"Failed: {e}")
            continue
        
        fig = wavelet.plot_spectrogram(amp, freqs, t_out, t_utc, f_proton, label)
        fig.axes[0].set_title(f'Juno Wavelet - {event_id}\n{date_str} {start_time}-{end_time}\n{label} ', fontsize=12)
        
        safe_date = date_str.replace('_', '-')
        safe_start = start_time.replace(':', '-')
        safe_end = end_time.replace(':', '-')
        out_file = os.path.join(event_dir, f"{event_id}_{safe_date}_{safe_start}_{safe_end}_{safe}.png")
        fig.savefig(out_file, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved {out_file}")

print("\nAll done.")