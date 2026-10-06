import re
import os
import glob
from datetime import datetime, timedelta
import numpy as np
from pathlib import Path

def parse_time_range_file(time_range_file):
    """
    解析时间段文档
    
    格式: PJ03:2016_11_30 10:29:28 15:29:28
    
    Returns:
    list of dict: 包含事件ID、开始时间、结束时间的列表
    """
    time_ranges = []
    
    with open(time_range_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            # 解析格式: PJ03:2016_11_30 10:29:28 15:29:28
            match = re.match(r'(\w+):(\d{4}_\d{2}_\d{2})\s+(\d{2}:\d{2}:\d{2})\s+(\d{2}:\d{2}:\d{2})', line)
            if match:
                event_id = match.group(1)
                date_str = match.group(2)  # 2016_11_30
                start_time_str = match.group(3)  # 10:29:28
                end_time_str = match.group(4)  # 15:29:28
                
                # 转换日期格式
                date_part = date_str.replace('_', '-')  # 2016-11-30
                
                # 构建完整的datetime对象
                start_datetime = datetime.fromisoformat(f"{date_part}T{start_time_str}")
                end_datetime = datetime.fromisoformat(f"{date_part}T{end_time_str}")
                
                time_ranges.append({
                    'event_id': event_id,
                    'date': date_str,
                    'start_time': start_datetime,
                    'end_time': end_datetime,
                    'original_line': line
                })
            else:
                print(f"警告: 无法解析行: {line}")
    
    return time_ranges


def find_d2s_file_for_date(data_folder, date_str):
    """
    根据日期字符串查找对应的D2S文件
    
    日期格式: 2016_11_30
    文件名格式: 2016.11.30.d2s 或 2016_11_30.d2s
    """
    # 尝试不同的文件名格式
    date_formats = [
        date_str.replace('_', '.'),  # 2016.11.30
        date_str,  # 2016_11_30
        date_str.replace('_', '-'),  # 2016-11-30
    ]
    
    for date_format in date_formats:
        # 搜索匹配的文件
        patterns = [
            os.path.join(data_folder, f"{date_format}.d2s"),
            os.path.join(data_folder, f"{date_format}*.d2s"),
            os.path.join(data_folder, f"*{date_format}*.d2s"),
        ]
        
        for pattern in patterns:
            matching_files = glob.glob(pattern)
            if matching_files:
                return matching_files[0]
    
    return None


def parse_d2s_file_with_interpolation(file_path, start_time=None, end_time=None, energy_threshold=200):
    """使用线性插值处理能道边界，计算超过指定能量阈值的平均通量"""
    
    if isinstance(start_time, str):
        start_time = datetime.fromisoformat(start_time)
    if isinstance(end_time, str):
        end_time = datetime.fromisoformat(end_time)
    
    try:
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()
    except Exception as e:
        print(f"读取文件失败 {file_path}: {e}")
        return None
    
    if len(lines) < 10:
        print(f"文件 {file_path} 格式不正确，行数不足")
        return None
    
    # 解析第七行的能量信息
    ytags_line = lines[6]
    ytags_match = re.search(r'yTags="([^"]*)"', ytags_line)
    if not ytags_match:
        print(f"无法解析能量标签: {file_path}")
        return None
    
    ytags_str = ytags_match.group(1)
    energy_centers = np.array([float(x.strip()) for x in ytags_str.split(',')])
    
    # 计算每个能道的边界和权重
    channel_boundaries, channel_weights = calculate_channel_weights(energy_centers, energy_threshold)
    
    print(f"能量中心: {energy_centers} keV")
    print(f"\n能道分析 (阈值={energy_threshold} keV):")
    print(f"{'能道索引':<8} {'能量中心':<12} {'能道边界':<30} {'权重':<8} {'说明'}")
    print("-" * 70)
    
    for i, (center, boundaries, weight) in enumerate(zip(energy_centers, channel_boundaries, channel_weights)):
        if weight == 1.0:
            status = "完全包含"
        elif weight > 0:
            status = f"部分包含 ({weight*100:.1f}%)"
        else:
            status = "不包含"
        print(f"{i:<8} {center:<12.4f} ({boundaries[0]:.4f}, {boundaries[1]:.4f})     {weight:<8.4f} {status}")
    
    # 解析数据行（从第10行开始，索引为9）
    valid_fluxes = []
    timestamps = []
    invalid_count = 0
    total_lines = 0
    time_range_lines = 0
    
    # 打印前几行数据用于调试
    print(f"\n数据行示例（前3行）:")
    for i, line in enumerate(lines[9:12]):
        print(f"  行{10+i}: {line.strip()}")
    
    for line in lines[9:]:  # 从第10行开始
        if not line.strip():
            continue
        
        total_lines += 1
        
        # 解析每一行数据
        # 格式: :01:2016-11-30T00:00:49.229 +3.29E-01 +2.44E-01 ...
        parts = line.strip().split()
        
        # 第一个字段是 :01:2016-11-30T00:00:49.229
        if len(parts) < 2:
            continue
        
        # 提取时间戳 - 从第一个字段中移除 :01: 前缀
        timestamp_field = parts[0]
        # 移除开头的 :数字: 格式
        timestamp_str = re.sub(r'^:\d+:', '', timestamp_field)
        
        try:
            # 解析时间戳
            if 'T' in timestamp_str:
                # 处理带毫秒的格式: 2016-11-30T00:00:49.229
                timestamp = datetime.strptime(timestamp_str, '%Y-%m-%dT%H:%M:%S.%f')
            else:
                timestamp = datetime.fromisoformat(timestamp_str)
        except Exception as e:
            if total_lines <= 3:  # 只打印前几行的错误
                print(f"  时间解析错误 (行{9+total_lines}): {timestamp_str}, 错误: {e}")
            continue
        
        # 检查时间范围
        if start_time and timestamp < start_time:
            continue
        if end_time and timestamp > end_time:
            continue
        
        time_range_lines += 1
        
        # 提取通量数据（从第二个字段开始）
        flux_strings = parts[1:]
        
        # 检查数据完整性
        if len(flux_strings) != len(energy_centers):
            if time_range_lines <= 3:  # 只打印前几行的调试信息
                print(f"  数据长度不匹配 (行{9+total_lines}): 期望{len(energy_centers)}个通量值，实际{len(flux_strings)}个")
            continue
        
        # 检查是否包含无效数据
        has_invalid = False
        fluxes = []
        
        for flux_str in flux_strings:
            try:
                # 处理科学计数法: +3.29E-01
                flux_value = float(flux_str)
                # 检查是否为无效数据标记
                if flux_value <= -1.0E30:  # 无效数据
                    has_invalid = True
                    break
                fluxes.append(flux_value)
            except ValueError:
                has_invalid = True
                break
        
        if has_invalid:
            invalid_count += 1
            continue
        
        # 使用权重计算超过阈值的总通量
        weighted_flux = sum(flux * weight for flux, weight in zip(fluxes, channel_weights))
        
        valid_fluxes.append(weighted_flux)
        timestamps.append(timestamp)
    
    # 打印调试信息
    print(f"\n数据解析统计:")
    print(f"总数据行数: {total_lines}")
    print(f"时间范围内行数: {time_range_lines}")
    print(f"有效数据点数: {len(valid_fluxes)}")
    print(f"无效数据点数: {invalid_count}")
    
    if len(valid_fluxes) == 0:
        print(f"警告: 在指定时间范围内未找到有效数据")
        if time_range_lines > 0:
            print(f"  有{time_range_lines}行在时间范围内，但可能数据格式有问题")
        return {
            'average_flux': None,
            'std_flux': None,
            'max_flux': None,
            'min_flux': None,
            'valid_samples': 0,
            'invalid_samples': invalid_count,
            'total_lines': total_lines,
            'time_range_lines': time_range_lines,
            'time_range': (None, None),
            'error': '没有有效数据'
        }
    
    average_flux = np.mean(valid_fluxes)
    std_flux = np.std(valid_fluxes)
    max_flux = np.max(valid_fluxes)
    min_flux = np.min(valid_fluxes)
    
    print(f"\n统计结果:")
    print(f"时间范围: {timestamps[0]} 到 {timestamps[-1]}")
    print(f"超过{energy_threshold}keV平均通量: {average_flux:.4e} counts/s")
    print(f"标准差: {std_flux:.4e} counts/s")
    print(f"最大值: {max_flux:.4e} counts/s")
    print(f"最小值: {min_flux:.4e} counts/s")
    
    return {
        'average_flux': average_flux,
        'std_flux': std_flux,
        'max_flux': max_flux,
        'min_flux': min_flux,
        'valid_samples': len(valid_fluxes),
        'invalid_samples': invalid_count,
        'total_lines': total_lines,
        'time_range_lines': time_range_lines,
        'time_range': (timestamps[0], timestamps[-1]),
        'timestamps': timestamps,
        'fluxes': valid_fluxes,
        'energy_centers': energy_centers,
        'channel_weights': channel_weights
    }

def calculate_channel_weights(energy_centers, threshold):
    """
    计算每个能道对于超过阈值的贡献权重
    
    假设每个能道覆盖的能量范围是从前一个能道中心到当前能道中心的中间点，
    到当前能道中心到下一个能道中心的中间点
    
    Parameters:
    energy_centers: 能量中心数组
    threshold: 能量阈值
    
    Returns:
    channel_boundaries: 每个能道的边界
    channel_weights: 每个能道中超过阈值部分的比例 (0-1)
    """
    n_channels = len(energy_centers)
    channel_boundaries = []
    channel_weights = []
    
    for i in range(n_channels):
        # 计算能道边界
        if i == 0:
            # 第一个能道：下边界为0，上边界为与下一个能道的中间点
            lower_bound = 0
            if n_channels > 1:
                upper_bound = (energy_centers[i] + energy_centers[i+1]) / 2
            else:
                upper_bound = energy_centers[i] * 2
        elif i == n_channels - 1:
            # 最后一个能道：下边界为与前一个能道的中间点，上边界为外推
            lower_bound = (energy_centers[i-1] + energy_centers[i]) / 2
            # 假设能道宽度对称，上边界 = 能量中心 + (能量中心 - 下边界)
            upper_bound = energy_centers[i] + (energy_centers[i] - lower_bound)
        else:
            # 中间能道：边界为相邻能道中心的中间点
            lower_bound = (energy_centers[i-1] + energy_centers[i]) / 2
            upper_bound = (energy_centers[i] + energy_centers[i+1]) / 2
        
        channel_boundaries.append((lower_bound, upper_bound))
        
        # 计算该能道中超过阈值的比例
        if upper_bound <= threshold:
            # 整个能道都在阈值以下
            weight = 0.0
        elif lower_bound >= threshold:
            # 整个能道都在阈值以上
            weight = 1.0
        else:
            # 阈值在能道内部，使用线性插值计算比例
            channel_width = upper_bound - lower_bound
            if channel_width > 0:
                weight = (upper_bound - threshold) / channel_width
            else:
                weight = 0.0
        
        channel_weights.append(weight)
    
    return channel_boundaries, channel_weights


def process_time_ranges(data_folder, time_range_file, output_file, energy_threshold=200):
    """
    批量处理时间段文档中的所有时间段
    
    Parameters:
    data_folder: D2S文件所在文件夹
    time_range_file: 时间段文档路径
    output_file: 输出结果文件路径
    energy_threshold: 能量阈值 (keV)
    """
    
    # 解析时间段文档
    print(f"解析时间段文档: {time_range_file}")
    time_ranges = parse_time_range_file(time_range_file)
    print(f"找到 {len(time_ranges)} 个时间段")
    
    # 准备输出文件
    with open(output_file, 'w', encoding='utf-8') as out_f:
        # 写入文件头
        out_f.write("=" * 120 + "\n")
        out_f.write(f"电子能谱分析结果 - 能量阈值: >{energy_threshold} keV\n")
        out_f.write(f"数据文件夹: {data_folder}\n")
        out_f.write(f"处理时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        out_f.write("=" * 120 + "\n\n")
        
        # 写入列标题
        header = (f"{'事件ID':<10} {'日期':<12} {'开始时间':<20} {'结束时间':<20} "
                 f"{'平均通量(c/s)':<16} {'标准差':<14} {'最大值':<14} {'最小值':<14} "
                 f"{'有效样本':<10} {'无效样本':<10} {'状态'}")
        out_f.write(header + "\n")
        out_f.write("-" * 120 + "\n")
        
        # 处理每个时间段
        success_count = 0
        fail_count = 0
        
        for i, time_range in enumerate(time_ranges, 1):
            event_id = time_range['event_id']
            date_str = time_range['date']
            start_time = time_range['start_time']
            end_time = time_range['end_time']
            
            print(f"\n处理 [{i}/{len(time_ranges)}]: {event_id} - {date_str} {start_time} -> {end_time}")
            
            # 查找对应的D2S文件
            d2s_file = find_d2s_file_for_date(data_folder, date_str)
            
            if not d2s_file:
                print(f"  警告: 未找到日期 {date_str} 的D2S文件")
                result_line = (f"{event_id:<10} {date_str:<12} "
                             f"{start_time.strftime('%Y-%m-%d %H:%M:%S'):<20} "
                             f"{end_time.strftime('%Y-%m-%d %H:%M:%S'):<20} "
                             f"{'N/A':<16} {'N/A':<14} {'N/A':<14} {'N/A':<14} "
                             f"{'N/A':<10} {'N/A':<10} {'未找到文件'}")
                out_f.write(result_line + "\n")
                fail_count += 1
                continue
            
            print(f"  使用文件: {os.path.basename(d2s_file)}")
            
            # 解析D2S文件并计算结果
            result = parse_d2s_file_with_interpolation(
                d2s_file,
                start_time=start_time,
                end_time=end_time,
                energy_threshold=energy_threshold
            )
            
            if result is None or result.get('error'):
                error_msg = result.get('error', '解析失败') if result else '解析失败'
                print(f"  错误: {error_msg}")
                result_line = (f"{event_id:<10} {date_str:<12} "
                             f"{start_time.strftime('%Y-%m-%d %H:%M:%S'):<20} "
                             f"{end_time.strftime('%Y-%m-%d %H:%M:%S'):<20} "
                             f"{'N/A':<16} {'N/A':<14} {'N/A':<14} {'N/A':<14} "
                             f"{'N/A':<10} {'N/A':<10} {error_msg}")
                out_f.write(result_line + "\n")
                fail_count += 1
                continue
            
            # 格式化输出
            if result['average_flux'] is not None:
                flux_str = f"{result['average_flux']:.4e}"
                std_str = f"{result['std_flux']:.4e}"
                max_str = f"{result['max_flux']:.4e}"
                min_str = f"{result['min_flux']:.4e}"
                status = "成功"
            else:
                flux_str = "N/A"
                std_str = "N/A"
                max_str = "N/A"
                min_str = "N/A"
                status = "无有效数据"
            
            result_line = (f"{event_id:<10} {date_str:<12} "
                         f"{start_time.strftime('%Y-%m-%d %H:%M:%S'):<20} "
                         f"{end_time.strftime('%Y-%m-%d %H:%M:%S'):<20} "
                         f"{flux_str:<16} {std_str:<14} {max_str:<14} {min_str:<14} "
                         f"{result['valid_samples']:<10} {result['invalid_samples']:<10} {status}")
            
            out_f.write(result_line + "\n")
            
            if result['average_flux'] is not None:
                success_count += 1
                print(f"  结果: 平均通量 = {flux_str} counts/s, 有效样本 = {result['valid_samples']}")
            else:
                fail_count += 1
        
        # 写入汇总信息
        out_f.write("\n" + "=" * 120 + "\n")
        out_f.write(f"处理完成: 成功 {success_count} 个, 失败 {fail_count} 个, 总计 {len(time_ranges)} 个\n")
    
    print(f"\n处理完成！")
    print(f"成功: {success_count}, 失败: {fail_count}, 总计: {len(time_ranges)}")
    print(f"结果已保存到: {output_file}")
    
    return success_count, fail_count


def main():
    """
    主函数 - 配置参数并执行处理
    """
    # 配置参数 - 请根据实际情况修改这些路径
    data_folder = "D:\summer research\e_flux"  # D2S文件所在文件夹
    time_range_file = "D:\summer research\pic\e_Energy_time.txt"  # 时间段文档路径
    output_file = "D:\summer research\pic\e_Energy_flux.txt"  # 输出结果文件路径
    energy_threshold = 100  # 能量阈值 (keV)
    
    # 检查文件是否存在
    if not os.path.exists(data_folder):
        print(f"错误: 数据文件夹不存在: {data_folder}")
        print("请创建文件夹并将D2S文件放入其中")
        return
    
    if not os.path.exists(time_range_file):
        print(f"错误: 时间段文档不存在: {time_range_file}")
        print("请创建时间段文档，格式示例: PJ03:2016_11_30 10:29:28 15:29:28")
        return
    
    # 列出找到的D2S文件
    d2s_files = glob.glob(os.path.join(data_folder, "*.d2s"))
    print(f"在 {data_folder} 中找到 {len(d2s_files)} 个D2S文件")
    
    # 执行批量处理
    success_count, fail_count = process_time_ranges(
        data_folder=data_folder,
        time_range_file=time_range_file,
        output_file=output_file,
        energy_threshold=energy_threshold
    )
    
    print(f"\n最终统计:")
    print(f"  成功处理: {success_count} 个时间段")
    print(f"  处理失败: {fail_count} 个时间段")
    print(f"  总计: {success_count + fail_count} 个时间段")


def create_sample_time_range_file(output_path):
    """
    创建示例时间段文档（用于测试）
    """
    sample_data = """PJ01:2016_11_30 10:29:28 15:29:28
PJ02:2016_12_01 08:00:00 12:00:00
PJ03:2016_12_02 14:30:00 18:30:00
PJ04:2016_12_03 06:00:00 10:00:00
PJ05:2016_12_04 22:00:00 23:59:59
"""
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(sample_data)
    
    print(f"示例时间段文档已创建: {output_path}")


if __name__ == "__main__":
    # 如果需要创建示例时间段文档，取消下面的注释
    # create_sample_time_range_file("./sample_time_ranges.txt")
    
    # 执行主程序
    main()