import numpy as np
import pandas as pd
import pywt
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter, MinuteLocator
import matplotlib.dates as mdates
from matplotlib.ticker import LogLocator, LogFormatter, FuncFormatter
from matplotlib import colors
import matplotlib.gridspec as gridspec

def analyze_magnetic_field(data_series, b_total, component_name, freq_range=(0.001, 0.4)):
    """
    Analyze 1Hz sampled magnetic field time series data (without anti-aliasing filter)
    
    Parameters:
        data_series: pandas Series (index as UTC time, values as magnetic field component)
        b_total: numpy array of total magnetic field strength at each time point (nT)
        component_name: string, magnetic field component name (e.g., 'Bx','By','Bz','B_total')
        freq_range: analysis frequency range (Hz), default (0.001, 0.4)
    
    Returns:
        amplitude: wavelet amplitude [frequency×time]
        freqs: frequency array (Hz) - logarithmic spacing
        time_sec: time array (seconds, relative to start)
        time_utc: original UTC time array
        ion_freq: ion cyclotron frequency array (Hz)
    """
    # 1. Prepare data
    data = data_series.values
    time_utc = data_series.index
    sampling_rate = 1.0  # fixed sampling rate of 1Hz
    
    # 2. Calculate relative time (seconds)
    start_time = time_utc[0]
    time_sec = np.array([(t - start_time).total_seconds() for t in time_utc])
    
    # 3. Set wavelet parameters
    wavelet = 'cmor1.5-1.0'  # complex Morlet wavelet
    num_scales = 200  # number of scales
    
    # 4. Generate logarithmically spaced frequency array
    fmin, fmax = freq_range
    freqs = np.logspace(np.log10(fmin), np.log10(fmax), num_scales)
    
    # 5. Calculate wavelet scales
    center_freq = pywt.central_frequency(wavelet)
    scales = center_freq * sampling_rate / freqs
    
    # 6. Perform Continuous Wavelet Transform (CWT)
    coefficients, freqs_cwt = pywt.cwt(data, scales, wavelet, sampling_period=1.0/sampling_rate)
    
    # 7. Normalize wavelet coefficients
    norm_factors = 1.0 / np.sqrt(scales)
    normalized_coeffs = coefficients * norm_factors[:, np.newaxis]
    
    # 8. Calculate amplitude (magnitude of normalized coefficients)
    amplitude = np.abs(normalized_coeffs)
    
    # 9. Calculate ion cyclotron frequency (f_ci = eB/(2πm_i))
    # Constants: e = 1.602e-19 C, m_i = 1.673e-27 kg (proton mass)
    # Factor: 0.0152 ≈ (1.602e-19)/(2π*1.673e-27)
    ion_freq = 0.0152 * b_total/16  # Hz
    
    # 10. Calculate total duration and statistics
    duration_hours = time_sec[-1] / 3600
    data_std = np.std(data)
    amp_std = np.std(amplitude)
    
    print(f"\n{'='*50}")
    print(f"Magnetic Field Component Analysis: {component_name}")
    print(f"Data range: {time_utc[0]} to {time_utc[-1]}")
    print(f"Duration: {duration_hours:.2f} hours")
    print(f"Number of samples: {len(data)}")
    print(f"Sampling rate: {sampling_rate} Hz")
    print(f"Nyquist frequency: {sampling_rate/2} Hz")
    print(f"Analysis frequency range: {fmin}-{fmax} Hz (logarithmic spacing)")
    print(f"Number of frequency points: {num_scales}")
    print(f"Original data standard deviation: {data_std:.2f} nT")
    print(f"Normalized amplitude standard deviation: {amp_std:.2f} nT")
    print(f"Amplitude range: {np.min(amplitude):.2f} - {np.max(amplitude):.2f} nT")
    print(f"Average ion cyclotron frequency: {np.mean(ion_freq):.4f} Hz")
    print('='*50)
    
    return coefficients,amplitude, freqs, time_sec, time_utc, ion_freq

def sci_notation_formatter(y, pos):
    """
    Formatter function for scientific notation in frequency axis
    Formats as 10^x with superscript exponent
    """
    if y <= 0:
        return "0"
    exponent = int(np.floor(np.log10(y)))
    coeff = y / (10**exponent)
    
    # For exponent values between -1 and -4, use simplified notation
    if -4 <= exponent < 0:
        return f"$10^{{{exponent}}}$"
    return f"${coeff:.0f} \\times 10^{{{exponent}}}$"

def plot_spectrogram(amplitude, freqs, time_sec, time_utc, ion_freq, component_name):
    """
    Plot spectrogram with linear color mapping
    """
    # Create figure with adjusted aspect ratio
    fig = plt.figure(figsize=(15, 8))
    
    # Create grid layout
    grid = plt.GridSpec(1, 2, width_ratios=[0.92, 0.05], wspace=0.05)
    ax_main = plt.subplot(grid[0])
    cbar_ax = plt.subplot(grid[1])
    
    # Calculate linear color range (5-95th percentile)
    vmin = np.percentile(amplitude, 5)
    vmax = np.percentile(amplitude, 95)
    
    # Plot with linear color mapping
    # 使用完整数据范围
    img = ax_main.pcolormesh(
        time_sec, freqs, amplitude, 
        shading='auto', cmap='jet'
    )
        
    # Set axes
    ax_main.set_yscale('log')
    ax_main.set_ylim(freqs[0], freqs[-1])
    
    # Set frequency axis with scientific notation
    ax_main.yaxis.set_major_locator(LogLocator(numticks=6))
    ax_main.yaxis.set_major_formatter(FuncFormatter(sci_notation_formatter))
    ax_main.yaxis.set_minor_formatter(plt.NullFormatter())
    
    ax_main.set_ylabel('Frequency (Hz)', fontsize=12)
    
    # Add ion cyclotron frequency as white curve
    ax_main.plot(time_sec, ion_freq, 'w-', linewidth=1, alpha=0.9, 
                 label='Ion Cyclotron Frequency')
    
    # Add time axis
    ax_main.set_xlabel('Time (seconds)', fontsize=12)
    
    # Add legend
    ax_main.legend(loc='upper right', facecolor='#333333', edgecolor='white', 
                   labelcolor='white', fontsize=10)
    
    # Add title
    ax_main.set_title(f'Magnetic Field Component {component_name} Spectrogram (Linear Scale)', fontsize=14)
    
    # Add vertical colorbar on the right
    cbar = fig.colorbar(img, cax=cbar_ax, orientation='vertical', 
                        label='Amplitude (nT)')
    
    # Add UTC time axis at the top
    utc_ax = ax_main.twiny()
    utc_ax.set_xlim(time_utc[0], time_utc[-1])
    utc_ax.set_xlabel('UTC Time', fontsize=10)
    
    # Set UTC time ticks
    locator = mdates.AutoDateLocator(minticks=5, maxticks=10)
    utc_ax.xaxis.set_major_locator(locator)
    utc_ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
    
    plt.tight_layout()
    #plt.show()
    return fig


def plot_time_series(data_series, component_name):
    """Plot time series"""
    plt.figure(figsize=(15, 5))
    plt.plot(data_series.index, data_series.values, 'b-', linewidth=0.8)
    plt.title(f'Magnetic Field Component {component_name} Time Series', fontsize=14)
    plt.xlabel('UTC Time', fontsize=12)
    plt.ylabel(f'{component_name} (nT)', fontsize=12)
    plt.gca().xaxis.set_major_formatter(DateFormatter('%H:%M'))
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

def plot_frequency_slice(amplitude, freqs, time_sec, component_name, target_freq):
    """
    Plot amplitude over time at a specific frequency
    
    Parameters:
        amplitude: wavelet amplitude
        freqs: frequency array (Hz)
        time_sec: time array (seconds)
        component_name: magnetic field component name
        target_freq: target frequency (Hz)
    """
    # Find closest frequency index
    idx = np.argmin(np.abs(freqs - target_freq))
    freq_actual = freqs[idx]
    
    plt.figure(figsize=(15, 6))
    
    # Plot amplitude over time
    plt.plot(time_sec, amplitude[idx, :], 'b-', linewidth=1.5)
    
    # Set title and labels
    plt.title(f'Magnetic Field Component {component_name} Amplitude at {freq_actual:.4f} Hz', fontsize=14)
    plt.xlabel('Time (seconds)', fontsize=12)
    plt.ylabel('Amplitude (nT)', fontsize=12)
    plt.grid(True, alpha=0.3)
    
    # Add UTC time axis
    ax_top = plt.gca().twiny()
    ax_top.set_xlim(time_utc[0], time_utc[-1])
    ax_top.set_xlabel('UTC Time', fontsize=12)
    
    # Set UTC time ticks
    locator = MinuteLocator(byminute=range(0, 60, 30))
    ax_top.xaxis.set_major_locator(locator)
    ax_top.xaxis.set_major_formatter(DateFormatter('%H:%M'))
    
    plt.tight_layout()
    plt.show()


def plot_three_spectrograms(amplitudes, freqs, time_sec, time_utc, ion_freqs, component_names):
    """
    Plot three spectrograms vertically with individual colorbars on the right,
    and shared time axis (only bottom plot shows time labels)
    """
    # Create figure with appropriate size
    fig = plt.figure(figsize=(15, 18))
    
    # Create 3 rows, 2 columns (spectrogram + colorbar for each)
    gs = gridspec.GridSpec(3, 2, width_ratios=[0.92, 0.05], height_ratios=[1, 1, 1], hspace=0.2)
    
    # Create axes for each component
    axes = []
    cbar_axes = []
    
    for i in range(3):
        # Main axis for spectrogram
        ax = plt.subplot(gs[i, 0])
        
        # Calculate individual color range (5-95th percentile)
        vmin = np.percentile(amplitudes[i], 5)
        vmax = np.percentile(amplitudes[i], 95)
        
        # Plot spectrogram
        img = ax.pcolormesh(
            time_sec, freqs, amplitudes[i], 
            shading='auto', cmap='jet', 
            vmin=vmin, vmax=vmax
        )
        
        # Set axes properties
        ax.set_yscale('log')
        ax.set_ylim(freqs[0], freqs[-1])
        # ax.set_title(f'Component: {component_names[i]}', fontsize=14, pad=12)
        
        # Set frequency axis with scientific notation
        ax.yaxis.set_major_locator(LogLocator(numticks=6))
        ax.yaxis.set_major_formatter(FuncFormatter(sci_notation_formatter))
        ax.yaxis.set_minor_formatter(plt.NullFormatter())
        ax.set_ylabel('Frequency (Hz)', fontsize=12)
        
        # Add ion cyclotron frequency as white curve
        ax.plot(time_sec, ion_freqs[i], 'w-', linewidth=2, alpha=0.9, 
                label='Ion Cyclotron Frequency')
        
        # Add legend only to the first plot
        if i == 0:
            ax.legend(loc='upper right', facecolor='#333333', edgecolor='white', 
                      labelcolor='white', fontsize=10)
        
        # Add UTC time axis at the top
        utc_ax = ax.twiny()
        utc_ax.set_xlim(time_utc[0], time_utc[-1])
        
        # Set UTC time ticks
        locator = mdates.AutoDateLocator(minticks=5, maxticks=10)
        utc_ax.xaxis.set_major_locator(locator)
        utc_ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
        
        # Only add UTC label to the first plot
        if i == 0:
            utc_ax.set_xlabel('UTC Time', fontsize=10)
        else:
            utc_ax.set_xticklabels([])
        
        # Hide x-axis labels for top two plots
        if i < 2:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel('Time (seconds)', fontsize=12)
        
        # Create colorbar axis for this component
        cax = plt.subplot(gs[i, 1])
        cbar = fig.colorbar(img, cax=cax)
        cbar.set_label('Amplitude (nT)', fontsize=10)
        cax.tick_params(labelsize=9)
        
        axes.append(ax)
        cbar_axes.append(cax)
    
    # Set overall title
    plt.suptitle('Magnetic Field Components Spectrograms (Linear Scale)', fontsize=16, y=0.95)
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.96)  # Adjust for suptitle
    plt.show()
def plot_wavelet_at_ion_freq(amplitude, freqs, time_sec, time_utc, ion_freq, component_name):
    """
    Plot wavelet amplitude at ion cyclotron frequency over time
    
    Parameters:
        amplitude: wavelet amplitude matrix (frequency×time)
        freqs: frequency array (Hz)
        time_sec: time array (seconds)
        time_utc: UTC time array
        ion_freq: ion cyclotron frequency array (Hz)
        component_name: magnetic field component name
    """
    # Create figure
    fig, ax = plt.subplots(figsize=(15, 6))
    
    # Store interpolated amplitudes
    amp_at_ion = np.zeros_like(time_sec, dtype=float)
    
    # Interpolate amplitude at ion cyclotron frequency for each time point
    for t_idx, t in enumerate(time_sec):
        # Get current ion cyclotron frequency
        current_ion_freq = ion_freq[t_idx]
        
        # Clamp frequency to analysis range
        if current_ion_freq < freqs[0]:
            current_ion_freq = freqs[0]
        elif current_ion_freq > freqs[-1]:
            current_ion_freq = freqs[-1]
        
        # Interpolate in log-frequency space
        log_freqs = np.log10(freqs)
        log_target = np.log10(current_ion_freq)
        amp_at_ion[t_idx] = np.interp(log_target, log_freqs, amplitude[:, t_idx])
    
    # Plot amplitude
    ax.plot(time_sec, amp_at_ion, 'b-', linewidth=1.5)
    
    # Set labels and title
    ax.set_xlabel('Time (seconds)', fontsize=12)
    ax.set_ylabel('Amplitude (nT)', fontsize=12)
    ax.set_title(f'{component_name} - Wavelet Amplitude at Ion Cyclotron Frequency', fontsize=14)
    ax.grid(True, alpha=0.3)
    
    # Add UTC time axis at the top
    utc_ax = ax.twiny()
    utc_ax.set_xlim(time_utc[0], time_utc[-1])
    utc_ax.set_xlabel('UTC Time', fontsize=10)
    
    # Set UTC time ticks
    locator = mdates.AutoDateLocator(minticks=5, maxticks=10)
    utc_ax.xaxis.set_major_locator(locator)
    utc_ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
    
    plt.tight_layout()
    plt.show()
# ===================== Usage Example =====================
# if __name__ == "__main__":
#     # Example: create simulated data
#     np.random.seed(42)
#     rng = pd.date_range('2023-01-01 00:00', '2023-01-01 06:00', freq='1S')
    
#     # Create simulated data with various frequency components
#     t = np.arange(len(rng))
    
#     # Base signal
#     data = 50 + 3 * np.random.randn(len(rng))
    
#     # Add low frequency variation (0.001 Hz)
#     data += 2 * np.sin(2 * np.pi * 0.001 * t)
    
#     # Add ion cyclotron frequency signal (0.15 Hz)
#     data += 1.5 * np.sin(2 * np.pi * 0.15 * t)
    
#     # Add high frequency signal (0.3 Hz)
#     data += 1.0 * np.sin(2 * np.pi * 0.3 * t)
    
#     # Create Series for magnetic field component
#     b_component_series = pd.Series(data, index=rng, name='B_component')
    
#     # Create total magnetic field strength (simulated variation)
#     b_total = 10 + 2 * np.sin(2 * np.pi * 0.0005 * t)  # slow variation
    
#     # Analyze magnetic field data
#     component = 'B_component'
#     amplitude, freqs, time_sec, time_utc, ion_freq = analyze_magnetic_field(
#         b_component_series,
#         b_total,  # total magnetic field array
#         component,
#         freq_range=(0.001, 0.4)  # 0.001-0.4 Hz range
#     )
    
#     # Visualize results
#     plot_spectrogram(amplitude, freqs, time_sec, time_utc, ion_freq, component)
#     plot_time_series(b_component_series, component)
#     plot_frequency_slice(amplitude, freqs, time_sec, component, target_freq=0.15)