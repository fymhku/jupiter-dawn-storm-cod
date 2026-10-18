from tool import process_juno_data
from tool import perform_mva
from tool import fit_harris_model
from tool import calculate_corotation_velocity
from tool import convert_time_to_seconds
from tool import plot
from tool import wavelet
from tool import PSD
from tool import distinguish as dis
from tool import find_period as fp
import pandas as pd
import numpy as np
import process_juno_data_multy
# Bx,By,Bz,X,Y,Z=process_juno_data.process_juno_data("2017_03_24", "2017-03-24/11:00:00", c)
# Bl, Bm, Bn, Vn, eigenvals=perform_mva.perform_mva(Bx, By, Bz, X, Y, Z)
start_time = "2017-01-29/06:00:00"
end_time = "2017-01-29/09:00:00"
Bx, By, Bz, x, y, z = process_juno_data_multy.process_juno_data(start_time, end_time)
index=Bx.index
plot.plot_magnetic_field_components(
    datetime=index,
    bx=Bx.values,
    by=By.values,
    bz=Bz.values,
    title_prefix='Magnetic Field',
    save_figure=True,
    save_path='D:\summer research\pic\\17.1.29\mag.png'
)
# B_t=np.sqrt(Bx.values**2+By.values**2+Bz.values**2)
# B_total_series = pd.Series(B_t, index=index, name='B_total')
# result=fp.analyze_magnetic_field(B_total_series, "2017-03-24/11:00:00","2017-03-24/15:00:00")
# print(f"   开始时间: {result['target_period_start'].strftime('%Y-%m-%d/%H:%M:%S')}")
# print(f"   结束时间: {result['target_end'].strftime('%Y-%m-%d/%H:%M:%S')}")
# # freq_range=(0.0001, 0.05)
# dis.highlight_signal_region(Bx.values)

#coefficients=[wavelet.analyze_magnetic_field(Bx,B_t, 'Bx', freq_range)[0],wavelet.analyze_magnetic_field(By,B_t, r'$B_r$', freq_range)[0],wavelet.analyze_magnetic_field(Bz,B_t, 'Bz', freq_range)[0]]
# # _, freqs, time_sec, time_utc,proton_freq = wavelet.analyze_magnetic_field(Bx,B_t, 'Bx', freq_range)
# # component_names=['Bx','By','Bz']
# # proton_freqs=[proton_freq,proton_freq,proton_freq]
# # wavelet.plot_three_spectrograms(amplitudes, freqs, time_sec, time_utc, proton_freqs, component_names)

#coefficients,amplitude, freqs, time_sec, time_utc,proton_freq = wavelet.analyze_magnetic_field(Bx,B_t, 'By', freq_range)
#wavelet.plot_wavelet_at_ion_freq(amplitude, freqs, time_sec, time_utc, proton_freq, 'Bx')
# wavelet.plot_spectrogram(amplitude, freqs, time_sec, time_utc,proton_freq, r'$B_{r}$')
# params, B_fit_T, Curv, curv0=fit_harris_model.fit_harris_model(Bl,Bn,["2017-05-16/05:00:00", "2017-05-16/09:00:00"],Vn)
# print("B0:",params['B0'],"L:",params['Lu'],"Bn:",params['bbn'])
# target_freqs, energy, peak_freq=PSD.compute_frequency_energy_spectrum(time_sec,                # 时间序列 (一维数组)
#     Bx.values,                # 原始数据 (一维数组，与time等长)
#     freqs,               # 小波分析的频率数组 (一维数组)
#     coefs=coefficients,               # 小波系数矩阵 (二维数组，shape: [len(freqs), len(time)])
#     f_min=1e-3,          # 目标频率最小值
#     f_max=0.1,           # 目标频率最大值
# )
# fig,peak=PSD.plot_energy_curve(time_sec, Bx.values, freqs, coefs=coefficients, f_min=1e-3, f_max=0.1, 
#                        data_units='nT', plot_title='', output_file=None)
# print(peak)
