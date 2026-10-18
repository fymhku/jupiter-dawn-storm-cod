import numpy as np
import matplotlib.pyplot as plt
import pywt
import matplotlib

def compute_frequency_energy_spectrum(time, data, freqs, coefs, f_min, f_max, l1_normalize=False):
    """
    Compute frequency-energy spectrum (returns energy curve data only)

    Parameters:
        time: time array (1D)
        data: raw data array (1D)
        freqs: frequency array (1D)
        coefs: wavelet coefficients matrix (2D, shape: [num_freqs × num_time_points])
        f_min: minimum frequency of interest (Hz)
        f_max: maximum frequency of interest (Hz)
        l1_normalize: if True, normalize energy so sum over the target frequency
                      range equals 1 (L1 normalization)

    Returns:
        target_freqs: frequency array within target range
        energy: corresponding energy density array (L1-normalized if requested)
        peak_freq: frequency at energy peak
    """
    # Validate wavelet coefficients dimensions
    if coefs.shape != (len(freqs), len(time)):
        if coefs.shape == (len(time), len(freqs)):
            coefs = coefs.T
        else:
            raise ValueError(f"Wavelet coefficients dimension mismatch: expected ({len(freqs)}, {len(time)}), got {coefs.shape}")
    
    # Compute average energy density per frequency (time average of |coefficients|²)
    energy_density = np.mean(np.abs(coefs)**2, axis=1)
    
    # Select data within target frequency range
    freq_mask = (freqs >= f_min) & (freqs <= f_max)
    target_freqs = freqs[freq_mask]
    energy = energy_density[freq_mask]

    # L1 归一化：使目标频率范围内的能量总和为 1
    if l1_normalize:
        total = np.sum(energy)
        if total > 0:
            energy = energy / total

    # Find peak energy frequency
    peak_idx = np.argmax(energy)
    peak_freq = target_freqs[peak_idx]
    
    return target_freqs, energy, peak_freq

def find_peaks(freq, energy):
    n = len(energy)
    if n < 3:
        return []  # 数据点不足，无法检测峰值
    
    # 步骤1：标记候选峰值点（满足 ≥ 左右邻点）
    candidate = [False] * n
    for i in range(1, n - 1):
        if energy[i] >= energy[i - 1] and energy[i] >= energy[i + 1]:
            candidate[i] = True

    # 步骤2：合并连续候选点为平台区间
    intervals = []
    i = 1
    while i < n - 1:
        if candidate[i]:
            start = i
            while i < n - 1 and candidate[i]:
                i += 1
            end = i - 1  # 平台结束索引
            intervals.append((start, end))
        else:
            i += 1

    # 步骤3和4：验证峰值平台并记录位置
    peaks = []
    for (start, end) in intervals:
        # 检查平台是否严格大于左右外部点
        left_valid = (energy[start] > energy[start - 1])
        right_valid = (energy[end] > energy[end + 1])
        
        if left_valid and right_valid:
            # 取平台中间索引（偶数长度时靠左）
            mid_index = (start + end) // 2
            peaks.append(freq[mid_index])
            
    return peaks
def plot_peak_frequency_coefficients(coefficients, amplitude, freqs, time_sec, time_utc, component_name):
    """
    绘制峰值频率对应的小波系数随时间变化
    
    参数:
        coefficients: 小波系数 (复数数组)
        amplitude: 小波振幅 (实部数组)
        freqs: 频率数组 (Hz)
        time_sec: 时间数组 (秒)
        time_utc: UTC时间数组
        component_name: 磁场分量名称
    """
    # 1. 找到整个频谱中的峰值频率
    # 计算每个频率的平均振幅
    # mean_amplitude = np.mean(amplitude, axis=1)
    
    # # 找到峰值频率的索引
    # peak_idx = np.argmax(mean_amplitude)
    # peak_freq = freqs[peak_idx]
    
    # print(f"找到峰值频率: {peak_freq:.6f} Hz")
    
    # 2. 提取该频率对应的小波系数
    peak_coeffs = coefficients[peak_idx, :]
    
    # 3. 计算实部、虚部、振幅和相位
    real_part = np.real(peak_coeffs)
    imag_part = np.imag(peak_coeffs)
    amplitude_vals = np.abs(peak_coeffs)
    phase = np.angle(peak_coeffs)  # 相位（弧度）
    
    # 4. 创建图形
    fig, axes = plt.subplots(4, 1, figsize=(15, 16), sharex=True)
    fig.suptitle(f'Peak Frequency ({peak_freq:.5f} Hz) Coefficients for {component_name}', 
                 fontsize=16, y=0.95)
    
    # 图1: 实部
    ax0 = axes[0]
    ax0.plot(time_sec, real_part, 'b-', linewidth=1.2)
    ax0.set_ylabel('Real Part (nT)', fontsize=12)
    ax0.grid(True, alpha=0.3)
    ax0.set_title('Real Component of Wavelet Coefficients', fontsize=14)
    
    # 图2: 虚部
    ax1 = axes[1]
    ax1.plot(time_sec, imag_part, 'r-', linewidth=1.2)
    ax1.set_ylabel('Imaginary Part (nT)', fontsize=12)
    ax1.grid(True, alpha=0.3)
    ax1.set_title('Imaginary Component of Wavelet Coefficients', fontsize=14)
    
    # 图3: 振幅
    ax2 = axes[2]
    ax2.plot(time_sec, amplitude_vals, 'g-', linewidth=1.2)
    ax2.set_ylabel('Amplitude (nT)', fontsize=12)
    ax2.grid(True, alpha=0.3)
    ax2.set_title('Amplitude of Wavelet Coefficients', fontsize=14)
    
    # 图4: 相位
    ax3 = axes[3]
    ax3.plot(time_sec, phase, 'm-', linewidth=1.2)
    ax3.set_ylabel('Phase (radians)', fontsize=12)
    ax3.set_xlabel('Time (seconds)', fontsize=12)
    ax3.grid(True, alpha=0.3)
    ax3.set_title('Phase of Wavelet Coefficients', fontsize=14)
    
    # 添加UTC时间轴
    utc_ax = ax3.twiny()
    utc_ax.set_xlim(time_utc[0], time_utc[-1])
    utc_ax.set_xlabel('UTC Time', fontsize=12)
    
    # 设置UTC时间刻度
    locator = mdates.AutoDateLocator()
    utc_ax.xaxis.set_major_locator(locator)
    utc_ax.xaxis.set_major_formatter(DateFormatter('%H:%M:%S'))
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.93)  # 为总标题调整空间
    plt.show()
def plot_energy_curve(time, data, freqs, coefs, f_min, f_max,
                      data_units='nT', plot_title='', output_file=None, l1_normalize=False):
    # 计算能量谱
    target_freqs, energy, peak_freq = compute_frequency_energy_spectrum(
        time, data, freqs, coefs, f_min, f_max, l1_normalize=l1_normalize
    )

    # 创建图形和坐标轴对象
    fig, ax = plt.subplots(figsize=(10, 6))

    # 在半对数坐标上绘制能量曲线
    ax.semilogx(target_freqs, energy, 'b-', linewidth=2, label='Energy Density')

    # 设置标签和标题
    ax.set_xlabel('Frequency (Hz)', fontsize=12)
    if l1_normalize:
        ax.set_ylabel('Normalized Magnetic Energy (L1, 1/Hz)', fontsize=12)
    else:
        ax.set_ylabel(f'Magnetic Energy ({data_units}²/Hz)', fontsize=12)
    ax.set_title(plot_title or 'Frequency-Energy Spectrum', fontsize=14)
    
    # 添加网格和图例
    ax.grid(True, which='both', linestyle='--', alpha=0.6)
    ax.legend(loc='upper right', fontsize=10)
    
    # 自动调整布局
    fig.tight_layout()
    
    # 可选：保存图像
    if output_file:
        fig.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Image saved to: {output_file}")
    peaks=find_peaks(target_freqs, energy)
    # 返回图形对象和峰值频率
    return fig, peaks
