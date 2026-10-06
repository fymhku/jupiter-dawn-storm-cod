import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import os

# ================================================
# 1. 磁场处理函数（保持原有格式）
# ================================================

def calculate_total_magnetic_field(Bx_series, By_series, Bz_series):
    """计算总磁场强度 B_total = sqrt(Bx² + By² + Bz²)"""
    # 确保三个序列具有相同的时间索引
    if not (Bx_series.index.equals(By_series.index) and 
            Bx_series.index.equals(Bz_series.index)):
        # 如果索引不同，先对齐数据
        aligned_data = pd.DataFrame({
            'Bx': Bx_series,
            'By': By_series,
            'Bz': Bz_series
        }).dropna()  # 删除任何有缺失值的时间点
        
        Bx = aligned_data['Bx'].values
        By = aligned_data['By'].values
        Bz = aligned_data['Bz'].values
        index = aligned_data.index
    else:
        Bx = Bx_series.values
        By = By_series.values
        Bz = Bz_series.values
        index = Bx_series.index
    
    # 计算总磁场强度
    B_total = np.sqrt(Bx**2 + By**2 + Bz**2)
    
    # 创建新的Series
    B_total_series = pd.Series(B_total, index=index, name='B_total')
    
    return B_total_series

def parse_custom_time(time_str, reference_date=None):
    """解析自定义时间格式：YYYY-MM-DD/HH:MM:SS 或 HH:MM:SS"""
    if '/' in time_str:
        # 格式: YYYY-MM-DD/HH:MM:SS
        date_part, time_part = time_str.split('/')
        year, month, day = map(int, date_part.split('-'))
        hour, minute, second = map(int, time_part.split(':'))
        return pd.Timestamp(year, month, day, hour, minute, second)
    else:
        # 格式: HH:MM:SS
        if reference_date is None:
            raise ValueError("当时间字符串不包含日期时，必须提供reference_date")
        hour, minute, second = map(int, time_str.split(':'))
        return pd.Timestamp(reference_date.year, reference_date.month, reference_date.day,
                           hour, minute, second)

def analyze_magnetic_field(data_series, start_time_str, end_time_str, hours_window=10):
    """
    分析磁场数据，找到最大平均磁场，以及满足条件的最大时间段
    
    参数:
    data_series: pandas Series，索引为datetime，值为磁场数据
    start_time_str: 开始时间字符串，格式为"YYYY-MM-DD/HH:MM:SS"或"HH:MM:SS"
    end_time_str: 结束时间字符串，格式为"YYYY-MM-DD/HH:MM:SS"或"HH:MM:SS"
    hours_window: 时间窗口大小（小时），默认10小时
    
    返回:
    结果字典，包含所有计算值和找到的时间段
    """
    # 确保索引是datetime类型
    if not isinstance(data_series.index, pd.DatetimeIndex):
        raise ValueError("数据索引必须是DatetimeIndex类型")
    
    # 获取数据的日期部分
    first_date = data_series.index[0].date()
    
    # 解析输入的时间段
    if '/' in start_time_str:
        # 开始时间包含完整日期
        period_start = parse_custom_time(start_time_str)
        # 结束时间可能需要跨天处理
        if '/' in end_time_str:
            period_end = parse_custom_time(end_time_str)
        else:
            # 结束时间只有时间部分，使用开始时间的日期
            period_end = parse_custom_time(end_time_str, reference_date=period_start)
    else:
        # 开始时间只有时间部分，使用第一天的日期
        period_start = parse_custom_time(start_time_str, reference_date=first_date)
        # 结束时间
        if '/' in end_time_str:
            # 结束时间有完整日期，但开始时间没有，这可能是错误情况
            raise ValueError("开始时间和结束时间格式不一致")
        else:
            period_end = parse_custom_time(end_time_str, reference_date=first_date)
    
    # 如果结束时间小于开始时间，说明跨天了
    if period_end < period_start:
        period_end += timedelta(days=1)
    
    # 提取时间段开始后指定小时窗口的数据
    end_window = period_start + timedelta(hours=hours_window)
    mask_window = (data_series.index >= period_start) & (data_series.index < end_window)
    data_window = data_series[mask_window]
    
    if len(data_window) == 0:
        raise ValueError(f"在{period_start}到{end_window}之间没有找到数据")
    
    # 按小时重采样并计算每小时的均值
    hourly_means = data_window.resample('H').mean()
    
    # 计算每小时均值的绝对值，并找到最大绝对值
    hourly_abs_means = hourly_means.abs()
    max_abs_hour_mean = hourly_abs_means.max()
    
    # 记录最大绝对值对应的小时段
    max_abs_hour_time = hourly_abs_means.idxmax()
    
    # 计算阈值：最大绝对值*0.9
    threshold = max_abs_hour_mean * 0.9
    
    # 提取给定时间段内的数据
    period_mask = (data_series.index >= period_start) & (data_series.index <= period_end)
    period_data = data_series[period_mask]
    
    if len(period_data) == 0:
        raise ValueError(f"在{period_start}到{period_end}之间没有找到数据")
    
    # 找到满足条件的数据点：数据的绝对值 < threshold
    condition = (period_data.abs() < threshold)
    
    # 找到连续满足条件的最大时间段
    def find_longest_continuous_period(condition_series, data_index):
        """找到连续满足条件的最大时间段"""
        if not condition_series.any():
            return None, None, timedelta(0), 0
        
        # 初始化变量
        max_start = None
        max_end = None
        max_duration = timedelta(0)
        max_length = 0
        
        current_start = None
        current_length = 0
        
        # 遍历每个数据点
        for i, (idx, meets_condition) in enumerate(zip(data_index, condition_series)):
            if meets_condition:
                if current_start is None:
                    current_start = idx
                current_length += 1
            else:
                if current_start is not None:
                    current_end = data_index[i-1]
                    current_duration = current_end - current_start
                    
                    if current_length > max_length:
                        max_length = current_length
                        max_duration = current_duration
                        max_start = current_start
                        max_end = current_end
                    
                    current_start = None
                    current_length = 0
        
        # 检查最后一段
        if current_start is not None:
            current_end = data_index[-1]
            current_duration = current_end - current_start
            
            if current_length > max_length:
                max_length = current_length
                max_duration = current_duration
                max_start = current_start
                max_end = current_end
        
        return max_start, max_end, max_duration, max_length
    
    # 找到满足条件的最长连续时间段
    target_start, target_end, target_duration, target_length = find_longest_continuous_period(
        condition.values, period_data.index
    )
    
    # 准备结果
    result = {
        'max_abs_hour_mean': max_abs_hour_mean,
        'max_abs_hour_time': max_abs_hour_time,
        'threshold': threshold,
        'period_start': period_start,
        'period_end': period_end,
        'window_start': period_start,
        'window_end': end_window,
        'window_hours': hours_window,
        'target_period_start': target_start,
        'target_end': target_end,
        'target_period_duration': target_duration,
        'target_period_length': target_length,
        'condition_met_count': condition.sum(),
        'condition_total_count': len(condition)
    }
    
    return result

# ================================================
# 2. 文本文件处理函数
# ================================================

def parse_input_line(line):
    """解析输入文件的每一行"""
    line = line.strip()
    if not line:
        return None
    
    try:
        # 格式: PJ05:2017_03_24 11:00:00 15:00:00
        # 分割轨道编号和剩余部分
        if ':' not in line:
            raise ValueError("行中未找到冒号分隔符")
        
        orbit_id_part, rest = line.split(':', 1)
        orbit_id = orbit_id_part.strip()
        
        # 分割日期和时间部分
        rest_parts = rest.strip().split()
        
        if len(rest_parts) >= 3:
            # 日期部分
            date_part = rest_parts[0]
            start_time = rest_parts[1]
            end_time = rest_parts[2]
            
            # 将日期格式从2017_03_24转换为2017-03-24
            date_part = date_part.replace('_', '-')
            
            # 构建完整的开始和结束时间字符串
            start_time_str = f"{date_part}/{start_time}"
            end_time_str = f"{date_part}/{end_time}"
            
            return {
                'orbit_id': orbit_id,
                'date': date_part,
                'start_time': start_time,
                'end_time': end_time,
                'start_time_str': start_time_str,
                'end_time_str': end_time_str,
                'raw_line': line
            }
        else:
            raise ValueError(f"行格式不正确: {line}")
            
    except Exception as e:
        print(f"解析行失败: {line}, 错误: {e}")
        return None

def process_time_period(orbit_info):
    """处理单个时间段，找到中点前后2.5小时的时间段"""
    orbit_id = orbit_info['orbit_id']
    start_time_str = orbit_info['start_time_str']
    end_time_str = orbit_info['end_time_str']
    
    try:
        # 解析开始和结束时间
        if '/' in start_time_str:
            period_start = parse_custom_time(start_time_str)
        else:
            # 修正：使用正确的方式获取第一天的日期
            # 这里需要从数据中获取日期，但由于我们不再读取磁场数据，
            # 我们使用当前日期或从轨道信息中获取日期
            if 'date' in orbit_info:
                reference_date = pd.Timestamp(orbit_info['date'])
                period_start = parse_custom_time(start_time_str, reference_date=reference_date)
            else:
                # 如果没有日期信息，使用当前日期
                period_start = parse_custom_time(start_time_str, reference_date=pd.Timestamp.now())
            
        if '/' in end_time_str:
            period_end = parse_custom_time(end_time_str)
        else:
            period_end = parse_custom_time(end_time_str, reference_date=period_start)
        
        # 如果结束时间小于开始时间，说明跨天了
        if period_end < period_start:
            period_end += timedelta(days=1)
        
        # 计算中点时间
        midpoint = period_start + (period_end - period_start) / 2
        
        # 计算中点前后2.5小时的时间段（共5小时）
        new_start = midpoint - timedelta(hours=2.5)
        new_end = midpoint + timedelta(hours=2.5)
        
        # 格式化输出
        date_str = new_start.strftime("%Y_%m_%d")
        start_time_str_output = new_start.strftime("%H:%M:%S")
        end_time_str_output = new_end.strftime("%H:%M:%S")
        
        # 创建输出行，格式与输入相同
        output_line = f"{orbit_id}:{date_str} {start_time_str_output} {end_time_str_output}"
        
        # 打印详细信息
        print(f"  原始时间段: {period_start} 到 {period_end}")
        print(f"  中点时间: {midpoint}")
        print(f"  新时间段: {new_start} 到 {new_end}")
        
        return output_line
        
    except Exception as e:
        print(f"处理轨道 {orbit_id} 时发生错误: {str(e)}")
        return None

# ================================================
# 3. 主程序：处理txt文件中的所有时间段
# ================================================

def process_all_time_periods(input_file_path, output_file_path):
    """处理txt文件中所有的的时间段，找到中点前后2.5小时的时间段"""
    
    print("时间段中点分析程序")
    print("=" * 60)
    print(f"输入文件: {input_file_path}")
    print(f"输出文件: {output_file_path}")
    print("-" * 60)
    
    # 检查输入文件是否存在
    if not os.path.exists(input_file_path):
        print(f"错误: 输入文件 {input_file_path} 不存在")
        return False
    
    # 读取输入文件
    try:
        with open(input_file_path, 'r') as f:
            lines = f.readlines()
        print(f"读取到 {len(lines)} 行数据")
    except Exception as e:
        print(f"读取输入文件失败: {e}")
        return False
    
    # 解析输入行
    orbit_infos = []
    for line_num, line in enumerate(lines, 1):
        line = line.strip()
        if line and not line.startswith('#'):  # 跳过空行和注释
            orbit_info = parse_input_line(line)
            if orbit_info:
                orbit_infos.append(orbit_info)
            else:
                print(f"第{line_num}行无效: {line}")
    
    if not orbit_infos:
        print("错误: 输入文件中没有有效的轨道信息")
        return False
    
    print(f"找到 {len(orbit_infos)} 个有效的时间段需要处理")
    print("=" * 60)
    
    # 处理每个时间段
    results = []
    success_count = 0
    fail_count = 0
    
    for i, orbit_info in enumerate(orbit_infos, 1):
        print(f"[{i}/{len(orbit_infos)}] 处理: {orbit_info['orbit_id']}")
        
        result = process_time_period(orbit_info)
        
        if result:
            results.append(result)
            success_count += 1
            print(f"  成功生成新时间段")
        else:
            fail_count += 1
            print(f"  处理失败")
        
        # 添加空行分隔不同轨道的输出
        if i < len(orbit_infos):
            print()
    
    # 写入输出文件
    if results:
        try:
            with open(output_file_path, 'w') as f:
                for result_line in results:
                    f.write(result_line + "\n")
            
            print("=" * 60)
            print(f"处理完成!")
            print(f"成功处理: {success_count} 个时间段")
            print(f"失败处理: {fail_count} 个时间段")
            print(f"结果已保存到: {output_file_path}")
            
            # 显示前几行结果作为预览
            print("\n结果预览:")
            for i, result in enumerate(results[:5], 1):
                print(f"  {result}")
            if len(results) > 5:
                print(f"  ... 还有 {len(results)-5} 行")
            
            return True
            
        except Exception as e:
            print(f"写入输出文件失败: {e}")
            return False
    else:
        print("\n错误: 没有找到任何满足条件的时间段")
        return False

# ================================================
# 4. 运行程序
# ================================================

if __name__ == "__main__":
    # 在这里设置输入和输出文件路径
    input_file_path = r"D:\summer research\pic\aurora date_modify_0.9.txt"  # 您的输入文件路径
    output_file_path = r"D:\summer research\pic\e_Energy_time.txt"  # 输出文件路径
    
    # 处理所有时间段
    success = process_all_time_periods(input_file_path, output_file_path)
    
    if success:
        print("-" * 60)
        print("程序执行成功!")
    else:
        print("-" * 60)
        print("程序执行失败!")
    
    # 等待用户按键退出（方便查看结果）
    input("\n按Enter键退出...")