import numpy as np
import matplotlib.pyplot as plt

def highlight_signal_region(signal, absolute_tolerance=0.2, min_ratio=0.6, max_ratio=0.6):
    """
    根据信号趋势智能处理边界点，使用绝对容差
    :param signal: 输入信号序列
    :param absolute_tolerance: 绝对容差值（固定值）
    :param min_ratio: 最小值阈值比例
    :param max_ratio: 最大值阈值比例
    """
    # 计算最小值和最大值及其位置
    min_val = np.min(signal)
    max_val = np.max(signal)
    min_index = np.argmin(signal)
    max_index = np.argmax(signal)
    
    # 计算阈值
    threshold_min = min_val * min_ratio
    threshold_max = max_val * max_ratio
    
    # 创建时间轴
    t = np.arange(len(signal))
    
    # 确定信号趋势（根据min和max的相对位置）
    if min_index < max_index:
        trend = "上升趋势: 最小值在左，最大值在右"
        # 上升趋势处理
        # 1. 左边界：最小值阈值的最左边点
        low_indices = np.where(
            np.abs(signal - threshold_min) <= absolute_tolerance
        )[0]
        left_bound = low_indices[0] if low_indices.size > 0 else np.argmin(np.abs(signal - threshold_min))
        
        # 2. 右边界：最大值阈值的最右边点
        high_indices = np.where(
            np.abs(signal - threshold_max) <= absolute_tolerance
        )[0]
        right_bound = high_indices[-1] if high_indices.size > 0 else np.argmin(np.abs(signal - threshold_max))
        
    else:
        trend = "下降趋势: 最大值在左，最小值在右"
        # 下降趋势处理
        # 1. 左边界：最大值阈值的最左边点
        high_indices = np.where(
            np.abs(signal - threshold_max) <= absolute_tolerance
        )[0]
        left_bound = high_indices[0] if high_indices.size > 0 else np.argmin(np.abs(signal - threshold_max))
        
        # 2. 右边界：最小值阈值的最右边点
        low_indices = np.where(
            np.abs(signal - threshold_min) <= absolute_tolerance
        )[0]
        right_bound = low_indices[-1] if low_indices.size > 0 else np.argmin(np.abs(signal - threshold_min))
    
    # 确保边界在有效范围内
    left_bound = max(0, min(left_bound, len(signal)-1))
    right_bound = max(0, min(right_bound, len(signal)-1))
    
    # 创建绘图
    plt.figure(figsize=(12, 6))
    
    # 绘制原始信号
    plt.plot(t, signal, 'b-', label='Signal', linewidth=2)
    
    # 标记关键点
    plt.plot(t[min_index], signal[min_index], 'go', markersize=8, label=f'Min: {min_val:.2f}')
    plt.plot(t[max_index], signal[max_index], 'ro', markersize=8, label=f'Max: {max_val:.2f}')
    
    # 标记边界点
    plt.plot(t[left_bound], signal[left_bound], 'cs', markersize=8, 
             label=f'Left Bound: {signal[left_bound]:.2f}')
    plt.plot(t[right_bound], signal[right_bound], 'ms', markersize=8, 
             label=f'Right Bound: {signal[right_bound]:.2f}')
    
    # 高亮中间区域
    plt.axvspan(t[left_bound], t[right_bound], alpha=0.3, color='yellow', label='Highlighted Region')
    
    # 添加阈值线及容差范围
    plt.axhline(y=threshold_min, color='green', linestyle='--', alpha=0.7, 
                label=f'Min Threshold: {threshold_min:.2f}±{absolute_tolerance}')
    plt.axhline(y=threshold_max, color='red', linestyle='--', alpha=0.7, 
                label=f'Max Threshold: {threshold_max:.2f}±{absolute_tolerance}')
    
    # 填充容差范围
    plt.fill_between(t, threshold_min - absolute_tolerance, threshold_min + absolute_tolerance, 
                     color='green', alpha=0.1)
    plt.fill_between(t, threshold_max - absolute_tolerance, threshold_max + absolute_tolerance, 
                     color='red', alpha=0.1)
    
    # 设置图例和标题
    plt.title(f'Current Sheet')
    plt.xlabel('Time/s')
    plt.ylabel('B/nT')
    plt.legend(loc='best')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()
    # 返回结果
    return left_bound, right_bound, threshold_min, threshold_max

# ====================== 使用示例 ====================== 
# if __name__ == "__main__":
#     # 测试1：下降趋势信号（最大值在左，最小值在右）
#     print("\n测试1：下降趋势信号")
#     np.random.seed(24)
#     t_falling = np.linspace(0, 10, 500)
#     signal_falling = np.zeros_like(t_falling)
#     signal_falling[:150] = 10.0  # 初始高值
#     signal_falling[150:350] = np.linspace(10, 1, 200)  # 下降段
#     signal_falling[350:] = 1.0  # 稳定低值
#     signal_falling += np.random.normal(0, 0.4, len(t_falling))
    
#     lb_falling, rb_falling, th_min, th_max = highlight_signal_region(
#         signal_falling, absolute_tolerance=0.3
#     )
#     print(f"下降趋势：左边界={lb_falling}, 值={signal_falling[lb_falling]:.2f} (阈值={th_max:.2f})")
#     print(f"        右边界={rb_falling}, 值={signal_falling[rb_falling]:.2f} (阈值={th_min:.2f})")
#     print(f"        区间大小={rb_falling - lb_falling}点")
    
#     # 显示图表
#     plt.show()
    
#     # 测试2：上升趋势信号（最小值在左，最大值在右）
#     print("\n测试2：上升趋势信号")
#     np.random.seed(42)
#     t_rising = np.linspace(0, 10, 500)
#     signal_rising = np.zeros_like(t_rising)
#     signal_rising[:150] = 1.0  # 初始低值
#     signal_rising[150:350] = np.linspace(1, 10, 200)  # 上升段
#     signal_rising[350:] = 10.0  # 稳定高值
#     signal_rising += np.random.normal(0, 0.4, len(t_rising))
    
#     lb_rising, rb_rising, th_min, th_max = highlight_signal_region(
#         signal_rising, absolute_tolerance=0.3
#     )
#     print(f"上升趋势：左边界={lb_rising}, 值={signal_rising[lb_rising]:.2f} (阈值={th_min:.2f})")
#     print(f"        右边界={rb_rising}, 值={signal_rising[rb_rising]:.2f} (阈值={th_max:.2f})")
#     print(f"        区间大小={rb_rising - lb_rising}点")
    
#     # 显示图表
#     plt.show()
    
#     # 测试3：不同容差值的影响
#     print("\n测试3：不同容差值对下降趋势的影响")
#     tolerances = [0.1, 0.3, 0.5, 1.0]
    
#     for tol in tolerances:
#         print(f"\n容差值: {tol}")
#         lb, rb, _, _ = highlight_signal_region(signal_falling, absolute_tolerance=tol)
#         print(f"左边界={lb}, 值={signal_falling[lb]:.2f}, 右边界={rb}, 值={signal_falling[rb]:.2f}")
#         plt.show()
    
#     # 测试4：自定义阈值比例
#     print("\n测试4：自定义阈值比例")
#     lb_custom, rb_custom, th_min, th_max = highlight_signal_region(
#         signal_falling, 
#         absolute_tolerance=0.3