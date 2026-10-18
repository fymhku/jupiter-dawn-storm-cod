import numpy as np
import matplotlib.pyplot as plt
import itertools

plt.rcParams['font.sans-serif'] = ['SimHei']  # 用于中文显示
plt.rcParams['axes.unicode_minus'] = False    # 解决负号显示问题

# 定义八次方程系数函数
def get_coeffs(q):
    """
    定义八次方程的系数：a₀x⁸ + a₁x⁷ + a₂x⁶ + a₃x⁵ + a₄x⁴ + a₅x³ + a₆x² + a₇x + a₈
    返回包含9个系数的列表：[a₀, a₁, a₂, a₃, a₄, a₅, a₆, a₇, a₈]
    """
    # 示例方程的系数计算
    return [1, 0, -(1+2*q**2+q**4), 0, (q**4+2*q**2), 0, (q**2-q**4), 0, -q**4]

# 参数范围和求解设置
q_min = 0
q_max = 2
num_points = 1000
q_values = np.linspace(q_min, q_max, num_points)

# 存储选定分支的解 - 跟踪实部>0, ≤1且虚部>=0的解
selected_solution = np.full(num_points, np.nan, dtype=complex)

# 初始猜测点 - 选择实部>0, ≤1且虚部>=0的点
initial_guess = complex(0.5, 0)  # 实部在(0,1]范围内

# 求解第一个点 (q_min)
coeffs = get_coeffs(q_min)
roots = np.roots(coeffs)

# 选择满足条件的解：实部>0, ≤1且虚部>=0
valid_mask = (roots.real > 0) & (roots.real <= 1) & (roots.imag >= 0)
valid_roots = roots[valid_mask]

if len(valid_roots) > 0:
    # 选择最接近初始猜测的解
    distances = np.abs(valid_roots - initial_guess)
    selected_idx = np.argmin(distances)
    selected_solution[0] = valid_roots[selected_idx]
    print(f"在 q={q_min} 处找到了满足条件的解: {selected_solution[0]}")
else:
    # 如果没有符合条件的解，设为NaN
    selected_solution[0] = np.nan
    print(f"在 q={q_min} 处未找到满足条件的解")

# 使用跟踪方法保持解的连续性，但允许断开
prev_solution = selected_solution[0]

# 求解循环 (从第二个点开始)
for i in range(1, num_points):
    q = q_values[i]
    
    # 计算当前系数
    coeffs = get_coeffs(q)
    
    try:
        # 求解当前根
        current_roots = np.roots(coeffs)
        
        # 选择满足条件的解：实部>0, ≤1且虚部>=0
        valid_mask = (current_roots.real > 0) & (current_roots.real <= 1.0045) & (current_roots.imag >= 0)
        valid_roots = current_roots[valid_mask]
        
        if len(valid_roots) > 0:
            if not np.isnan(prev_solution):
                # 如果有前一个解，选择最接近前一个解的解
                distances = np.abs(valid_roots - prev_solution)
                selected_idx = np.argmin(distances)
                current_solution = valid_roots[selected_idx]
            else:
                # 如果没有前一个解，选择最接近初始猜测的解
                distances = np.abs(valid_roots - initial_guess)
                selected_idx = np.argmin(distances)
                current_solution = valid_roots[selected_idx]
            
            # 存储结果
            selected_solution[i] = current_solution
            prev_solution = current_solution
        else:
            # 没有满足条件的解，设为NaN
            selected_solution[i] = np.nan
            prev_solution = np.nan
        
    except Exception as e:
        print(f"在 q={q:.3f} 处求解失败: {str(e)}")
        # 设为NaN
        selected_solution[i] = np.nan
        prev_solution = np.nan

# 分离实部和虚部
real_parts = np.real(selected_solution)
imag_parts = np.imag(selected_solution)

# 在指定区间创建实部空白区域 (q=0.70 到 0.71)，虚部保持正常
blank_mask = (q_values >= 0.70) & (q_values <= 0.715)
real_parts[blank_mask] = np.nan

# 找出连续的解段 - 分别处理实部和虚部
real_segments = []
current_real_seg = []

imag_segments = []
current_imag_seg = []

for i in range(num_points):
    # 处理实部
    if not np.isnan(real_parts[i]):
        current_real_seg.append((q_values[i], real_parts[i]))
    else:
        if current_real_seg:
            real_segments.append(np.array(current_real_seg))
            current_real_seg = []
    
    # 处理虚部
    if not np.isnan(imag_parts[i]):
        current_imag_seg.append((q_values[i], imag_parts[i]))
    else:
        if current_imag_seg:
            imag_segments.append(np.array(current_imag_seg))
            current_imag_seg = []

# 添加最后一个段
if current_real_seg:
    real_segments.append(np.array(current_real_seg))
if current_imag_seg:
    imag_segments.append(np.array(current_imag_seg))

# 创建高质量图形
plt.figure(figsize=(14, 10))

# 绘制实部曲线 - 所有连续段（只添加一次标签）
first_real = True
for seg in real_segments:
    q_vals = seg[:, 0]
    reals = seg[:, 1]
    
    # 绘制实部曲线
    if first_real:
        plt.plot(q_vals, reals, 
                 color='#1f77b4', 
                 linewidth=3.0, 
                 alpha=0.9,
                 label=r'$Re(x)$')
        first_real = False
    else:
        plt.plot(q_vals, reals, 
                 color='#1f77b4', 
                 linewidth=3.0, 
                 alpha=0.9)

# 绘制虚部曲线 - 所有连续段（只添加一次标签）
first_imag = True
for seg in imag_segments:
    q_vals = seg[:, 0]
    imags = seg[:, 1]
    
    # 绘制虚部曲线
    if first_imag:
        plt.plot(q_vals, imags, 
                 color='#d62728', 
                 linewidth=2.5, 
                 alpha=0.8,
                 linestyle='--',
                 label=r'$Im(x)$')
        first_imag = False
    else:
        plt.plot(q_vals, imags, 
                 color='#d62728', 
                 linewidth=2.5, 
                 alpha=0.8,
                 linestyle='--')

# 添加网格
plt.grid(True, linestyle='--', alpha=0.6)

# 设置坐标轴范围
# 自动确定范围，但确保关键区域可见
y_min = min(np.nanmin(real_parts), np.nanmin(imag_parts)) 
y_max = max(np.nanmax(real_parts), np.nanmax(imag_parts)) 
plt.ylim(min(y_min, -0.1) if not np.isnan(y_min) else -0.1, 
         max(y_max, 1.1) if not np.isnan(y_max) else 1.1)

# 添加标题和标签
plt.title('ion cyclotron instability(ion drift)', fontsize=16, fontweight='bold')
plt.xlabel('q', fontsize=14, fontweight='bold')
plt.ylabel('x', fontsize=14, fontweight='bold')

# 添加图例
plt.legend(loc='best', fontsize=12, frameon=True, shadow=True)

# 添加方程信息
# equation_text = r"$(x^4 - qx^3 + 2x^2 + 3x + 1.5)(x^4 + qx^3 + 0.5x^2 - x + 0.8) = 0$"
# plt.figtext(0.5, 0.01, f"方程: {equation_text}", 
#             ha="center", fontsize=12, 
#             bbox={"facecolor": "white", "alpha": 0.8, "pad": 5})

# # 添加约束条件说明
# condition_text = r"约束条件: $0 < Re(x) \leq 1$ 且 $Im(x) \geq 0$"
# plt.figtext(0.5, 0.95, condition_text, 
#             ha="center", fontsize=12, 
#             bbox={"facecolor": "yellow", "alpha": 0.3, "pad": 5})

# # 添加空白区域说明
# blank_text = r"在 $q \in [0.70, 0.71]$ 区间实部设为空白"
# plt.annotate(blank_text, xy=(0.705, 0.5), xytext=(0.8, 0.7),
#              arrowprops=dict(facecolor='black', arrowstyle='->'),
#              fontsize=12, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8))

# 调整布局
plt.tight_layout(rect=[0, 0.03, 1, 0.90])

# 保存高清图像
plt.savefig('eighth_equation_constrained_branch.png', dpi=300, bbox_inches='tight')

# 显示图形
plt.show()