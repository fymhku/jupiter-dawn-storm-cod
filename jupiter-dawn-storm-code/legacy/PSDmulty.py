import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tool import process_juno_data
from tool import wavelet
from tool import PSD
from scipy.signal import find_peaks

# 配置参数
INPUT_FILE = "D:\summer research\pic\PSD_aurora_0.9\\fix.txt"
OUTPUT_DIR = "D:\summer research\pic\PSD_aurora_0.9"
EXCEL_OUTPUT = "D:\summer research\pic\PSD_aurora_0.9\\fix.xlsx"
FREQ_RANGE = [0.0001, 0.05]

# 创建主输出目录
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 读取日期时间文件
with open(INPUT_FILE, 'r') as f:
    lines = f.readlines()

# 初始化数据存储列表
all_results = []

# 处理每个时间段
for line in lines:
    if ':' not in line:
        continue
        
    # 解析新格式
    event_part, datetime_part = line.strip().split(':', 1)
    event_id = event_part.strip()
    
    # 分割日期和时间段
    datetime_parts = datetime_part.split()
    if len(datetime_parts) < 3:
        continue
        
    date_str = datetime_parts[0]
    start_time = datetime_parts[1]
    end_time = datetime_parts[2]
    
    print(f"\n{'='*60}")
    print(f"Processing {event_id} - {date_str} from {start_time} to {end_time}")
    print(f"{'='*60}")
    
    # 为当前事件创建独立文件夹
    event_output_dir = os.path.join(OUTPUT_DIR, event_id)
    os.makedirs(event_output_dir, exist_ok=True)
    
    # 构建完整的时间字符串
    date_formatted = date_str.replace('_', '-')
    start_time_full = f"{date_formatted}/{start_time}"
    end_time_full = f"{date_formatted}/{end_time}"
    
    # 读取Juno数据
    try:
        Bx, By, Bz, pos_x, pos_y, pos_z = process_juno_data.process_juno_data(
            date_str, start_time_full, end_time_full
        )
        
        if Bx is None or len(Bx) == 0:
            print(f"Error: No valid data returned for {event_id}")
            continue
        
        # 计算总磁场
        B_total = np.sqrt(Bx**2 + By**2 + Bz**2)
        
    except Exception as e:
        print(f"Error reading data for {event_id}: {str(e)}")
        continue
    
    # 初始化当前事件的峰值存储字典
    event_peaks_dict = {
        'Event': f"{event_id}:{date_str} {start_time} {end_time}"
    }
    
    # 处理三个磁场分量
    component_index = 1
    for safe_component, component, data in zip(["B_r","B_theta","B_phi"], 
                                                [r'$B_r$', r'$B_\theta$', r'$B_\phi$'], 
                                                [Bx, By, Bz]):
        # 执行小波变换
        try:
            print(f"\nPerforming wavelet analysis on {component}...")
            coefficients, amplitude, freqs, time_sec, time_utc, proton_freq = wavelet.analyze_magnetic_field(
                data, 
                B_total, 
                component,
                freq_range=FREQ_RANGE
            )
        except Exception as e:
            print(f"Wavelet transform failed for {event_id}-{component}: {str(e)}")
            continue
        
        comp_name = component
        
        # 创建新图形
        fig, peak_freq = PSD.plot_energy_curve(time_sec, data, freqs, coefficients,  
                                               f_min=1e-3, f_max=0.1, 
                                               data_units='nT', plot_title='', 
                                               output_file=None)
        
        # 获取第一个子图（假设是频率-能量谱图）
        ax = fig.axes[0]
        
        # 获取频率-能量谱数据
        lines = ax.get_lines()
        if len(lines) > 0:
            line = lines[0]  # 假设第一条线是频率-能量曲线
            x_data = line.get_xdata()
            y_data = line.get_ydata()
            
            # 提取峰值频率和对应的能量密度
            if peak_freq is not None:
                if isinstance(peak_freq, (list, np.ndarray)):
                    # 如果peak_freq是数组，获取所有峰值
                    for i, freq in enumerate(peak_freq):
                        # 找到最接近峰值频率的点
                        idx = np.argmin(np.abs(x_data - freq))
                        energy_density = y_data[idx]
                        
                        # 在Excel中存储峰值信息
                        col_name_freq = f"Peak{component_index}_Freq_Hz"
                        col_name_energy = f"Peak{component_index}_Energy_nT2_Hz"
                        
                        # 格式化频率为科学计数法
                        freq_str = f"{freq:.3e}"
                        energy_str = f"{energy_density:.3e}"
                        
                        # 存储到字典中
                        event_peaks_dict[col_name_freq] = freq_str
                        event_peaks_dict[col_name_energy] = energy_str
                        
                        # 标记峰值
                        ax.plot(freq, energy_density, 'ro', markersize=8, 
                                markerfacecolor='red', markeredgecolor='black')
                        ax.annotate(f'{freq:.3e} Hz', 
                                   xy=(freq, energy_density),
                                   xytext=(10, 10),
                                   textcoords='offset points',
                                   fontsize=9,
                                   bbox=dict(boxstyle="round,pad=0.3", 
                                            facecolor='yellow', 
                                            alpha=0.7))
                        
                        component_index += 1
                else:
                    # 如果peak_freq是单个值
                    idx = np.argmin(np.abs(x_data - peak_freq))
                    energy_density = y_data[idx]
                    
                    col_name_freq = f"Peak{component_index}_Freq_Hz"
                    col_name_energy = f"Peak{component_index}_Energy_nT2_Hz"
                    
                    freq_str = f"{peak_freq:.3e}"
                    energy_str = f"{energy_density:.3e}"
                    
                    event_peaks_dict[col_name_freq] = freq_str
                    event_peaks_dict[col_name_energy] = energy_str
                    
                    ax.plot(peak_freq, energy_density, 'ro', markersize=10, 
                            markerfacecolor='red', markeredgecolor='black')
                    ax.annotate(f'Peak: {peak_freq:.3e} Hz', 
                               xy=(peak_freq, energy_density),
                               xytext=(15, 15),
                               textcoords='offset points',
                               fontsize=10,
                               bbox=dict(boxstyle="round,pad=0.3", 
                                        facecolor='yellow', 
                                        alpha=0.8))
                    
                    component_index += 1
        
        # 设置标题
        ax.set_title(
            f'Frequency-Energy Spectrum with Peak Detection - {event_id}\n'
            f'{date_str} {start_time}-{end_time}\n'
            f'Component: {comp_name}', 
            fontsize=12
        )
        
        # 添加图例
        ax.legend(loc='best')
        
        # 保存图像
        safe_date = date_str.replace('_', '-')
        safe_time_start = start_time.replace(':', '-')
        safe_time_end = end_time.replace(':', '-')
        filename = f"{event_id}_{safe_date}_{safe_time_start}_{safe_time_end}_{safe_component}_peaks.png"
        output_path = os.path.join(event_output_dir, filename)
        
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        
        print(f"Saved {output_path}")
    
    # 将当前事件的结果添加到总列表中
    all_results.append(event_peaks_dict)

# 创建DataFrame并保存为Excel
if all_results:
    df = pd.DataFrame(all_results)
    
    # 重新排列列的顺序，使频率和能量密度成对出现
    columns_order = ['Event']
    
    # 找出所有的峰值列
    peak_cols = [col for col in df.columns if col.startswith('Peak')]
    
    # 按峰值编号排序
    peak_numbers = sorted(set([int(col.split('_')[0].replace('Peak', '')) for col in peak_cols]))
    
    # 按顺序添加频率和能量密度列
    for peak_num in peak_numbers:
        freq_col = f"Peak{peak_num}_Freq_Hz"
        energy_col = f"Peak{peak_num}_Energy_nT2_Hz"
        
        if freq_col in df.columns:
            columns_order.append(freq_col)
        if energy_col in df.columns:
            columns_order.append(energy_col)
    
    # 重新排列DataFrame的列
    df = df[columns_order]
    
    # 保存到Excel
    df.to_excel(EXCEL_OUTPUT, index=False)
    print(f"\n{'='*60}")
    print(f"Excel文件已保存到: {EXCEL_OUTPUT}")
    print(f"总共处理了 {len(df)} 个事件")
    print(f"{'='*60}")
    
    # 打印前几行数据预览
    print("\n数据预览:")
    print(df.head())
else:
    print("\n没有找到任何有效数据，无法创建Excel文件。")

print("\nProcessing completed!")