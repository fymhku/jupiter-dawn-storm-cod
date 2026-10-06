import os
import pandas as pd
from datetime import datetime

def process_juno_data(date_str, start_time, end_time):
    """
    处理Juno卫星的磁场和位置数据，并筛选指定时间范围内的数据
    
    参数:
        date_str (str): 日期字符串，格式为"YYYY_MM_DD"，例如"2016_12_10"
        start_time (str): 开始时间戳，格式为"YYYY-MM-DD/HH:MM:SS"
        end_time (str): 结束时间戳，格式为"YYYY-MM-DD/HH:MM:SS"
    
    返回:
        tuple: 包含六个pandas Series的元组 (Bx_series, By_series, Bz_series, 
                pos_x_series, pos_y_series, pos_z_series)
    """
    # 设置文件路径
    pos_dir = r"D:\Juno_Data\Juno_data_txt\pos\1s"
    mag_dir = r"D:\Juno_Data\Juno_data_txt\mag_tvar\MAG_JSS_1s\juno_B1s_JSSRTP"
    
    # 构建文件名
    pos_filename = f"JUNO_MAG_JSS_1s_{date_str}.txt"
    mag_filename = f"JUNO_MAG_JSS_1s_{date_str}.txt"
    
    # 构建完整文件路径
    pos_filepath = os.path.join(pos_dir, pos_filename)
    mag_filepath = os.path.join(mag_dir, mag_filename)
    
    print(f"位置文件路径: {pos_filepath}")
    print(f"磁场文件路径: {mag_filepath}")
    print(f"筛选时间范围: {start_time} 至 {end_time}")
    
    try:
        # 检查文件是否存在
        if not os.path.exists(pos_filepath):
            raise FileNotFoundError(f"位置文件不存在: {pos_filepath}")
        if not os.path.exists(mag_filepath):
            raise FileNotFoundError(f"磁场文件不存在: {mag_filepath}")
        
        # 使用pandas读取数据文件
        print("使用pandas读取数据文件...")
        
        # 读取磁场数据 - 跳过第一行（标题行）
        mag_df = pd.read_csv(
            mag_filepath, 
            delim_whitespace=True, 
            header=None, 
            names=['time_str', 'Bx', 'By', 'Bz'],
            dtype={'time_str': str, 'Bx': float, 'By': float, 'Bz': float},
            skiprows=1,  # 跳过标题行
            on_bad_lines='warn'
        )
        
        # 读取位置数据 - 跳过第一行（标题行）
        pos_df = pd.read_csv(
            pos_filepath, 
            delim_whitespace=True, 
            header=None, 
            names=['time_str', 'pos_x', 'pos_y', 'pos_z'],
            dtype={'time_str': str, 'pos_x': float, 'pos_y': float, 'pos_z': float},
            skiprows=1,  # 跳过标题行
            on_bad_lines='warn'
        )
        
        # 检查数据是否为空
        if mag_df.empty:
            raise ValueError("磁场数据文件为空")
        if pos_df.empty:
            raise ValueError("位置数据文件为空")
        
        # 合并两个数据集
        merged_df = pd.merge(
            mag_df, 
            pos_df, 
            on='time_str', 
            how='inner'
        )
        
        # 检查合并后的数据是否为空
        if merged_df.empty:
            raise ValueError("合并后数据为空 - 没有匹配的时间戳")
        
        # 将时间字符串转换为datetime对象
        merged_df['datetime'] = pd.to_datetime(merged_df['time_str'], format='%Y-%m-%d/%H:%M:%S', errors='coerce')
        
        # 检查是否有无效时间
        if merged_df['datetime'].isnull().any():
            invalid_count = merged_df['datetime'].isnull().sum()
            print(f"警告: 发现 {invalid_count} 条记录的时间格式无效")
            # 删除无效时间记录
            merged_df = merged_df.dropna(subset=['datetime'])
        
        # 设置datetime为索引
        merged_df.set_index('datetime', inplace=True)
        
        # 筛选指定时间范围内的数据
        start_dt = pd.to_datetime(start_time, format='%Y-%m-%d/%H:%M:%S')
        end_dt = pd.to_datetime(end_time, format='%Y-%m-%d/%H:%M:%S')
        
        print(f"合并后总记录数: {len(merged_df)}")
        filtered_df = merged_df.loc[start_dt:end_dt]
        print(f"筛选后记录数: {len(filtered_df)}")
        
        # 创建各个分量的Series
        Bx_series = filtered_df['Bx']
        By_series = filtered_df['By']
        Bz_series = filtered_df['Bz']
        pos_x_series = filtered_df['pos_x']
        pos_y_series = filtered_df['pos_y']
        pos_z_series = filtered_df['pos_z']
        
        # 打印处理摘要
        print("\n数据处理完成:")
        print(f"总记录数: {len(merged_df)}")
        print(f"筛选范围内记录数: {len(filtered_df)}")
        
        if not filtered_df.empty:
            print("\n首条记录:")
            print(f"时间: {filtered_df.index[0]}")
            print(f"磁场: Bx={Bx_series.iloc[0]:.4f}, By={By_series.iloc[0]:.4f}, Bz={Bz_series.iloc[0]:.4f}")
            print(f"位置: x={pos_x_series.iloc[0]:.4f}, y={pos_y_series.iloc[0]:.4f}, z={pos_z_series.iloc[0]:.4f}")
            
            print("\n末条记录:")
            print(f"时间: {filtered_df.index[-1]}")
            print(f"磁场: Bx={Bx_series.iloc[-1]:.4f}, By={By_series.iloc[-1]:.4f}, Bz={Bz_series.iloc[-1]:.4f}")
            print(f"位置: x={pos_x_series.iloc[-1]:.4f}, y={pos_y_series.iloc[-1]:.4f}, z={pos_z_series.iloc[-1]:.4f}")
        
        return Bx_series, By_series, Bz_series, pos_x_series, pos_y_series, pos_z_series
    
    except FileNotFoundError as e:
        print(f"错误: {e}")
        return None, None, None, None, None, None
    except pd.errors.ParserError as e:
        print(f"CSV解析错误: {e}")
        return None, None, None, None, None, None
    except ValueError as e:
        print(f"数据错误: {e}")
        return None, None, None, None, None, None
    except Exception as e:
        print(f"处理文件时出错: {type(e).__name__} - {e}")
        return None, None, None, None, None, None

