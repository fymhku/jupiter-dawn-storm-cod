import os
import math
import datetime

# ==================== 用户设置区域 ====================
# 时间段输入文件路径（每行格式：轨道号:日期 开始时间 结束时间）
INPUT_FILE = r"D:\summer research\\pic\\PSD_aurora_0.9\\fix.txt"

# 计算结果输出文件路径
OUTPUT_FILE = r"D:\summer research\\pic\\PSD_aurora_0.9\\position.txt"

# 位置数据存放目录（存放 txt 文件的根目录）
DATA_DIR = r"D:\Juno_Data\Juno_data_txt\\pos\1s"
# ====================================================

def parse_time_line(time_str):
    """将时间字符串（如"2016-06-25/00:00:01"）转换为datetime对象"""
    dt_str = time_str.replace('/', ' ')
    return datetime.datetime.strptime(dt_str, "%Y-%m-%d %H:%M:%S")

def compute_avg_r(data_file, start_dt, end_dt):
    """从数据文件中计算指定时间范围内的平均R值，返回平均值，若无数据返回None"""
    total_r = 0.0
    count = 0
    try:
        with open(data_file, 'r', encoding='utf-8') as f:
            # 跳过标题行
            next(f)
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                if len(parts) < 4:
                    continue
                dt = parse_time_line(parts[0])
                if start_dt <= dt <= end_dt:
                    x = float(parts[1])
                    y = float(parts[2])
                    z = float(parts[3])
                    total_r += math.sqrt(x*x + y*y + z*z)
                    count += 1
    except Exception as e:
        print(f"处理文件 {data_file} 时出错: {e}")
        return None
    if count == 0:
        return None
    return total_r / count

def main():
    results = []
    # 读取时间段输入文件
    with open(INPUT_FILE, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # 解析格式：PJ03:2016_12_01 13:28:00 13:53:52
            if ':' not in line:
                print(f"跳过无效行：{line}")
                continue
            orbit, rest = line.split(':', 1)
            parts = rest.split()
            if len(parts) != 3:
                print(f"跳过无效行（时间格式错误）：{line}")
                continue
            date_str, start_time_str, end_time_str = parts

            # 构建数据文件路径
            file_name = f"JUNO_MAG_JSS_1s_{date_str}.txt"
            data_file = os.path.join(DATA_DIR, file_name)

            # 构建起止时间 datetime 对象
            date_for_dt = date_str.replace('_', '-')
            start_dt_str = f"{date_for_dt} {start_time_str}"
            end_dt_str = f"{date_for_dt} {end_time_str}"
            start_dt = datetime.datetime.strptime(start_dt_str, "%Y-%m-%d %H:%M:%S")
            end_dt = datetime.datetime.strptime(end_dt_str, "%Y-%m-%d %H:%M:%S")

            # 计算平均R
            avg_r = compute_avg_r(data_file, start_dt, end_dt)
            if avg_r is None:
                results.append(f"{line}\t无数据")
            else:
                results.append(f"{line}\t{avg_r:.6f}")

    # 写入输出文件
    with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
        f.write("\n".join(results))
    print(f"处理完成，结果已保存至 {OUTPUT_FILE}")

if __name__ == "__main__":
    main()