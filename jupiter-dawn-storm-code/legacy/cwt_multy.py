import os
import numpy as np
import matplotlib.pyplot as plt
from tool import process_juno_data  # 替换为实际数据读取模块
from tool import wavelet  # 替换为实际小波分析模块

# 配置参数
INPUT_FILE = "D:\summer research\pic\\aurora date_modify_0.9.txt"  # 输入文件名
OUTPUT_DIR = "D:\summer research\pic\\aurora pic_0.9"    # 输出目录
FREQ_RANGE = [0.0005, 0.05]     # 小波分析频率范围 (Hz)

# 创建主输出目录
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 读取日期时间文件
with open(INPUT_FILE, 'r') as f:
    lines = f.readlines()

# 处理每个时间段
for line in lines:
    # 解析新格式: PJ12:2018_03_30 18:00:00 22:00:00
    if ':' not in line:
        continue
        
    # 分割事件编号和日期时间部分
    event_part, datetime_part = line.strip().split(':', 1)
    event_id = event_part.strip()
    
    # 分割日期和时间段
    datetime_parts = datetime_part.split()
    if len(datetime_parts) < 3:
        continue
        
    date_str = datetime_parts[0]      # 2018_03_30 (保持下划线格式)
    start_time = datetime_parts[1]    # 18:00:00
    end_time = datetime_parts[2]      # 22:00:00
    
    print(f"\n{'='*60}")
    print(f"Processing {event_id} - {date_str} from {start_time} to {end_time}")
    print(f"{'='*60}")
    
    # 为当前事件创建独立文件夹
    event_output_dir = os.path.join(OUTPUT_DIR, event_id)
    os.makedirs(event_output_dir, exist_ok=True)
    
    # 构建完整的时间字符串（符合process_juno_data要求的格式）
    date_formatted = date_str.replace('_', '-')  # 转换为 YYYY-MM-DD
    start_time_full = f"{date_formatted}/{start_time}"
    end_time_full = f"{date_formatted}/{end_time}"
    
    # 读取Juno数据 - 使用原始date_str构建文件名，完整时间字符串筛选数据
    try:
        print(f"Calling process_juno_data with:")
        print(f"  date_str = {date_str}")
        print(f"  start_time = {start_time_full}")
        print(f"  end_time = {end_time_full}")
        
        # 调用数据读取函数
        Bx, By, Bz, pos_x, pos_y, pos_z = process_juno_data.process_juno_data(
            date_str, start_time_full, end_time_full
        )
        
        # 检查返回数据是否有效
        if Bx is None or len(Bx) == 0:
            print(f"Error: No valid data returned for {event_id}")
            continue
        
        # 打印数据统计信息
        print(f"Data loaded successfully:")
        print(f"  Bx: {len(Bx)} points, min={np.min(Bx):.2f}, max={np.max(Bx):.2f}")
        print(f"  By: {len(By)} points, min={np.min(By):.2f}, max={np.max(By):.2f}")
        print(f"  Bz: {len(Bz)} points, min={np.min(Bz):.2f}, max={np.max(Bz):.2f}")
        
        # 计算总磁场
        B_total = np.sqrt(Bx**2 + By**2 + Bz**2)
        
    except Exception as e:
        print(f"Error reading data for {event_id}: {str(e)}")
        continue
    
    # 处理三个磁场分量
    for safe_component,component, data in zip(["B_r","B_theta","B_phi"],[r'$B_r$', r'$B_\theta$', r'$B_\phi$'], [Bx, By, Bz]):
        # 执行小波变换
        try:
            print(f"\nPerforming wavelet analysis on {component}...")
            coefficients,amplitude, freqs, time_sec, time_utc, proton_freq = wavelet.analyze_magnetic_field(
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
        fig = wavelet.plot_spectrogram(
            amplitude, freqs, time_sec, time_utc, proton_freq, comp_name
        )
        
        # 设置更详细的标题（包含事件编号）
        fig.axes[0].set_title(
            f'Juno Wavelet Analysis - {event_id}\n{date_str} {start_time}-{end_time}\nComponent: {comp_name}', 
            fontsize=12
        )
        
        # 创建安全文件名
        safe_date = date_str.replace('_', '-')
        safe_time_start = start_time.replace(':', '-')
        safe_time_end = end_time.replace(':', '-')
        filename = f"{event_id}_{safe_date}_{safe_time_start}_{safe_time_end}_{safe_component}.png"
        output_path = os.path.join(event_output_dir, filename)
        
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close(fig)  # 显式关闭图形释放内存
        
        print(f"Saved {output_path}")

print("\nProcessing completed!")