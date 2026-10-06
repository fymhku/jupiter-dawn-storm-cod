import os
import pandas as pd
import numpy as np
from datetime import datetime, timedelta

def process_juno_data(start_datetime_str, end_datetime_str):
    """
    连续读取Juno卫星在指定时间段内的磁场和位置数据
    
    参数:
        start_datetime_str (str): 开始时间戳，格式为"YYYY-MM-DD/HH:MM:SS"
        end_datetime_str (str): 结束时间戳，格式为"YYYY-MM-DD/HH:MM:SS"
    
    返回:
        tuple: 包含六个pandas Series的元组 (Bx_series, By_series, Bz_series, 
                pos_x_series, pos_y_series, pos_z_series)
    """
    # 设置文件路径
    pos_dir = r"D:\Juno_Data\Juno_data_txt\pos\1s"
    mag_dir = r"D:\Juno_Data\Juno_data_txt\mag_tvar\MAG_JSS_1s\juno_B1s_JSSRTP"
    
    # 将时间字符串转换为datetime对象
    start_dt = pd.to_datetime(start_datetime_str, format='%Y-%m-%d/%H:%M:%S')
    end_dt = pd.to_datetime(end_datetime_str, format='%Y-%m-%d/%H:%M:%S')
    
    # 计算日期范围
    date_list = []
    current_date = start_dt.replace(hour=0, minute=0, second=0)
    while current_date <= end_dt:
        date_list.append(current_date.strftime("%Y_%m_%d"))
        current_date += timedelta(days=1)
    
    print(f"处理日期范围: {start_dt.strftime('%Y-%m-%d')} 到 {end_dt.strftime('%Y-%m-%d')}")
    print(f"包含天数: {len(date_list)} 天")
    print(f"筛选时间范围: {start_dt.strftime('%Y-%m-%d/%H:%M:%S')} 至 {end_dt.strftime('%Y-%m-%d/%H:%M:%S')}")
    
    # 初始化存储所有数据的列表
    all_data = []
    
    for date_str in date_list:
        # 构建文件名
        pos_filename = f"JUNO_MAG_JSS_1s_{date_str}.txt"
        mag_filename = f"JUNO_MAG_JSS_1s_{date_str}.txt"
        
        # 构建完整文件路径
        pos_filepath = os.path.join(pos_dir, pos_filename)
        mag_filepath = os.path.join(mag_dir, mag_filename)
        
        # 检查文件是否存在
        if not os.path.exists(pos_filepath) or not os.path.exists(mag_filepath):
            if not os.path.exists(pos_filepath) and not os.path.exists(mag_filepath):
                print(f"跳过 {date_str}: 位置和磁场文件均不存在")
            elif not os.path.exists(pos_filepath):
                print(f"跳过 {date_str}: 位置文件不存在")
            else:
                print(f"跳过 {date_str}: 磁场文件不存在")
            continue
        
        print(f"\n处理日期: {date_str}")
        print(f"位置文件: {pos_filename}")
        print(f"磁场文件: {mag_filename}")
        
        try:
            # 读取磁场数据
            mag_df = pd.read_csv(
                mag_filepath, 
                delim_whitespace=True, 
                header=None, 
                names=['time_str', 'Bx', 'By', 'Bz'],
                dtype={'time_str': str, 'Bx': float, 'By': float, 'Bz': float},
                skiprows=1,
                on_bad_lines='warn'
            )
            
            # 读取位置数据
            pos_df = pd.read_csv(
                pos_filepath, 
                delim_whitespace=True, 
                header=None, 
                names=['time_str', 'pos_x', 'pos_y', 'pos_z'],
                dtype={'time_str': str, 'pos_x': float, 'pos_y': float, 'pos_z': float},
                skiprows=1,
                on_bad_lines='warn'
            )
            
            # 检查数据是否为空
            if mag_df.empty or pos_df.empty:
                if mag_df.empty:
                    print(f"警告: 磁场数据为空")
                if pos_df.empty:
                    print(f"警告: 位置数据为空")
                continue
            
            # 将时间字符串转换为datetime对象
            mag_df['datetime'] = pd.to_datetime(mag_df['time_str'], format='%Y-%m-%d/%H:%M:%S', errors='coerce')
            pos_df['datetime'] = pd.to_datetime(pos_df['time_str'], format='%Y-%m-%d/%H:%M:%S', errors='coerce')
            
            # 移除无效时间记录
            mag_df = mag_df.dropna(subset=['datetime'])
            pos_df = pos_df.dropna(subset=['datetime'])
            
            if mag_df.empty or pos_df.empty:
                print(f"警告: 转换后数据为空")
                continue
            
            # 合并磁场和位置数据
            merged = pd.merge(
                mag_df[['datetime', 'Bx', 'By', 'Bz']], 
                pos_df[['datetime', 'pos_x', 'pos_y', 'pos_z']], 
                on='datetime', 
                how='inner'
            )
            
            if merged.empty:
                print(f"警告: 合并后数据为空")
                continue
            
            # 添加到总数据
            all_data.append(merged)
            print(f"成功加载: {len(merged)} 条记录")
            
        except Exception as e:
            print(f"处理 {date_str} 时出错: {type(e).__name__} - {e}")
            continue
    
    # 检查是否有数据加载
    if not all_data:
        print("\n错误: 没有加载到任何有效数据")
        return None, None, None, None, None, None
    
    # 合并所有数据
    combined_df = pd.concat(all_data, ignore_index=False)
    
    # 按时间排序
    combined_df.sort_values('datetime', inplace=True)
    combined_df.set_index('datetime', inplace=True)
    
    print(f"\n合并后总记录数: {len(combined_df)}")
    
    # 筛选指定时间范围内的数据
    filtered_df = combined_df.loc[start_dt:end_dt]
    
    print(f"筛选后记录数: {len(filtered_df)}")
    
    # 创建各个分量的Series
    Bx_series = filtered_df['Bx']
    By_series = filtered_df['By']
    Bz_series = filtered_df['Bz']
    pos_x_series = filtered_df['pos_x']
    pos_y_series = filtered_df['pos_y']
    pos_z_series = filtered_df['pos_z']
    
    # 检查是否有数据
    if filtered_df.empty:
        print("\n警告: 筛选后数据为空!")
        return Bx_series, By_series, Bz_series, pos_x_series, pos_y_series, pos_z_series
    
    # 打印处理摘要
    print("\n数据处理完成:")
    print(f"总记录数: {len(combined_df)}")
    print(f"筛选范围内记录数: {len(filtered_df)}")
    
    # 计算时间范围
    time_range = filtered_df.index[-1] - filtered_df.index[0]
    print(f"时间范围: {filtered_df.index[0]} 到 {filtered_df.index[-1]}")
    print(f"时间跨度: {time_range}")
    
    print("\n首条记录:")
    print(f"时间: {filtered_df.index[0]}")
    print(f"磁场: Bx={Bx_series.iloc[0]:.4f}, By={By_series.iloc[0]:.4f}, Bz={Bz_series.iloc[0]:.4f}")
    print(f"位置: x={pos_x_series.iloc[0]:.4f}, y={pos_y_series.iloc[0]:.4f}, z={pos_z_series.iloc[0]:.4f}")
    
    print("\n末条记录:")
    print(f"时间: {filtered_df.index[-1]}")
    print(f"磁场: Bx={Bx_series.iloc[-1]:.4f}, By={By_series.iloc[-1]:.4f}, Bz={Bz_series.iloc[-1]:.4f}")
    print(f"位置: x={pos_x_series.iloc[-1]:.4f}, y={pos_y_series.iloc[-1]:.4f}, z={pos_z_series.iloc[-1]:.4f}")
    
    # 检查数据连续性
    time_diffs = np.diff(filtered_df.index)
    min_diff = min(time_diffs) if len(time_diffs) > 0 else timedelta(0)
    max_diff = max(time_diffs) if len(time_diffs) > 0 else timedelta(0)
    print(f"\n时间间隔统计: 最小间隔={min_diff}, 最大间隔={max_diff}")
    
    return Bx_series, By_series, Bz_series, pos_x_series, pos_y_series, pos_z_series
start_time = "2016-12-10/00:00:00"
end_time = "2016-12-15/23:59:59"
Bx, By, Bz, x, y, z = process_juno_data(start_time, end_time)