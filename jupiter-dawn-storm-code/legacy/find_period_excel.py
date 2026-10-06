import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import os
import traceback  # 用于打印详细异常

# ================================================
# 1. 磁场处理函数（完整版）
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
    
    # 计算阈值：最大绝对值*0.8
    threshold = max_abs_hour_mean * 0.8
    
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
# 2. 处理单个时间段（返回结构化字典）
# ================================================

def process_time_period(orbit_info):
    """
    处理单个时间段，返回结果字典
    
    参数:
    orbit_info: 字典，包含轨道编号、原始开始时间字符串、原始结束时间字符串
    
    返回:
    字典，包含处理结果
    """
    # 初始化结果字典
    result_dict = {
        'orbit_id': orbit_info['orbit_id'],
        'original_start': orbit_info['start_time_str'],
        'original_end': orbit_info['end_time_str'],
        'success': False,
        'target_start': None,
        'target_end': None,
        'duration': None,
        'duration_minutes': 0,
        'output_line': None,
        'target_start_str': None,
        'target_end_str': None,
        'error': None
    }
    
    try:
        # 导入process_juno_data函数
        from process_juno_data_multy import process_juno_data
        
        # 获取磁场数据
        Bx_series, By_series, Bz_series, x, y, z = process_juno_data(
            result_dict['original_start'], result_dict['original_end']
        )
        
        # 检查是否有数据
        if len(Bx_series) == 0 or len(By_series) == 0 or len(Bz_series) == 0:
            result_dict['error'] = f"轨道 {orbit_info['orbit_id']} 没有磁场数据"
            return result_dict
        
        # 检查索引是否为DatetimeIndex
        if not isinstance(Bx_series.index, pd.DatetimeIndex):
            # 尝试将索引转换为DatetimeIndex
            Bx_series.index = pd.to_datetime(Bx_series.index)
            By_series.index = pd.to_datetime(By_series.index)
            Bz_series.index = pd.to_datetime(Bz_series.index)
        
        # 计算总磁场
        B_total_series = calculate_total_magnetic_field(Bx_series, By_series, Bz_series)
        
        # 分析磁场数据
        result = analyze_magnetic_field(B_total_series, 
                                        result_dict['original_start'], 
                                        result_dict['original_end'], 
                                        hours_window=10)
        
        # 提取满足条件的时间段
        target_start = result['target_period_start']
        target_end = result['target_end']
        
        if target_start is None or target_end is None:
            result_dict['error'] = f"轨道 {orbit_info['orbit_id']} 未找到满足条件的时间段"
            return result_dict
        
        # 计算持续时间（分钟）
        duration = target_end - target_start
        duration_minutes = duration.total_seconds() / 60.0
        
        # 生成与旧版兼容的输出行（可选）
        date_str = target_start.strftime("%Y_%m_%d")
        start_time_only = target_start.strftime("%H:%M:%S")
        end_time_only = target_end.strftime("%H:%M:%S")
        output_line = f"{orbit_info['orbit_id']}:{date_str} {start_time_only} {end_time_only}"
        
        # 填充结果字典
        result_dict['success'] = True
        result_dict['target_start'] = target_start
        result_dict['target_end'] = target_end
        result_dict['duration'] = duration
        result_dict['duration_minutes'] = duration_minutes
        result_dict['output_line'] = output_line
        result_dict['target_start_str'] = target_start.strftime("%Y-%m-%d %H:%M:%S")
        result_dict['target_end_str'] = target_end.strftime("%Y-%m-%d %H:%M:%S")
        
        return result_dict
        
    except ImportError as e:
        result_dict['error'] = f"无法导入process_juno_data_multy模块: {str(e)}"
        return result_dict
    except Exception as e:
        # 打印详细堆栈信息便于调试
        traceback.print_exc()
        result_dict['error'] = f"处理轨道 {orbit_info['orbit_id']} 时发生错误: {str(e)}"
        return result_dict

# ================================================
# 3. 辅助函数：统一时间格式
# ================================================

def ensure_time_format(time_val):
    """将Excel中的时间值统一转换为'YYYY-MM-DD/HH:MM:SS'格式字符串"""
    if pd.isna(time_val):
        return None
    if isinstance(time_val, pd.Timestamp):
        return time_val.strftime("%Y-%m-%d/%H:%M:%S")
    else:
        s = str(time_val).strip()
        # 如果包含空格，替换为'/'（适配原解析函数）
        if ' ' in s:
            s = s.replace(' ', '/')
        # 确保日期部分使用连字符，时间部分用冒号
        # 这里假设原字符串已经符合要求，否则可能需要更复杂的解析
        return s

# ================================================
# 4. 主程序：从带表头的Excel读取，处理并输出Excel
# ================================================

def process_all_time_periods(input_excel_path, output_excel_path):
    """
    从Excel读取所有时间段（带表头：编号、开始时间、结束时间），
    处理后筛选持续时间>30分钟的结果，保存到新的Excel文件。
    """
    print("磁场数据分析程序 (Excel版，带表头)")
    print("=" * 60)
    print(f"输入文件: {input_excel_path}")
    print(f"输出文件: {output_excel_path}")
    print("-" * 60)

    # 检查输入文件
    if not os.path.exists(input_excel_path):
        print(f"错误: 输入文件 {input_excel_path} 不存在")
        return False

    # 读取Excel，第一行为表头
    try:
        df_input = pd.read_excel(input_excel_path, header=0)
        print(f"读取到 {len(df_input)} 行数据")
    except Exception as e:
        print(f"读取Excel文件失败: {e}")
        return False

    # 检查必备列
    required_columns = ['编号', '开始时间', '结束时间']
    if not all(col in df_input.columns for col in required_columns):
        print(f"错误: Excel文件必须包含列: {required_columns}")
        print(f"当前列名为: {list(df_input.columns)}")
        return False

    # 构建轨道信息列表
    orbit_infos = []
    for idx, row in df_input.iterrows():
        orbit_id = str(row['编号']).strip()
        start_time_raw = row['开始时间']
        end_time_raw = row['结束时间']

        # 统一时间格式
        start_time_str = ensure_time_format(start_time_raw)
        end_time_str = ensure_time_format(end_time_raw)

        if start_time_str is None or end_time_str is None:
            print(f"警告: 第{idx+2}行时间缺失，跳过")
            continue

        orbit_info = {
            'orbit_id': orbit_id,
            'start_time_str': start_time_str,
            'end_time_str': end_time_str,
            'original_row': idx
        }
        orbit_infos.append(orbit_info)

    if not orbit_infos:
        print("错误: 输入文件中没有有效的轨道信息")
        return False

    print(f"找到 {len(orbit_infos)} 个有效的时间段需要处理")
    print("=" * 60)

    # 处理每个时间段
    results = []
    success_count = 0
    fail_count = 0
    long_enough_count = 0

    for i, orbit_info in enumerate(orbit_infos, 1):
        print(f"[{i}/{len(orbit_infos)}] 处理: {orbit_info['orbit_id']}")

        result = process_time_period(orbit_info)

        # 防御性检查：如果返回None，则跳过
        if result is None:
            print(f"  警告: 轨道 {orbit_info['orbit_id']} 返回了 None，跳过")
            fail_count += 1
            continue

        if result['success']:
            success_count += 1
            duration_min = result['duration_minutes']
            print(f"  成功找到时间段，持续 {duration_min:.2f} 分钟")
            if duration_min > 30:
                long_enough_count += 1
                results.append(result)
                print(f"  持续时间超过30分钟，已记录")
            else:
                print(f"  持续时间不足30分钟，忽略")
        else:
            fail_count += 1
            print(f"  处理失败: {result.get('error', '未知错误')}")

        if i < len(orbit_infos):
            print()

    # 将结果写入Excel
    if results:
        output_rows = []
        for r in results:
            output_rows.append({
                '编号': r['orbit_id'],
                '原始开始时间': r['original_start'],
                '原始结束时间': r['original_end'],
                '筛选后开始时间': r['target_start_str'],
                '筛选后结束时间': r['target_end_str'],
                '持续时间(分钟)': round(r['duration_minutes'], 2)
            })

        df_output = pd.DataFrame(output_rows)

        try:
            df_output.to_excel(output_excel_path, index=False)

            print("=" * 60)
            print(f"处理完成!")
            print(f"成功处理: {success_count} 个时间段")
            print(f"失败处理: {fail_count} 个时间段")
            print(f"持续时间超过30分钟: {long_enough_count} 个")
            print(f"结果已保存到: {output_excel_path}")
            print("\n结果预览:")
            print(df_output.head().to_string(index=False))
            if len(df_output) > 5:
                print(f"... 共 {len(df_output)} 行")
            return True
        except Exception as e:
            print(f"写入输出Excel文件失败: {e}")
            return False
    else:
        print("\n错误: 没有找到任何满足条件（持续时间>30分钟）的时间段")
        return False

# ================================================
# 5. 运行程序
# ================================================

if __name__ == "__main__":
    # 请根据实际情况修改文件路径
    input_excel_path = "D:\summer research\pic\current_sheet\current_sheet_2022.xlsx"   # 输入Excel文件
    output_excel_path = "D:\summer research\pic\current_sheet\current_sheet_2022_modify.xlsx"  # 输出Excel文件

    success = process_all_time_periods(input_excel_path, output_excel_path)

    if success:
        print("-" * 60)
        print("程序执行成功!")
    else:
        print("-" * 60)
        print("程序执行失败!")

    input("\n按Enter键退出...")