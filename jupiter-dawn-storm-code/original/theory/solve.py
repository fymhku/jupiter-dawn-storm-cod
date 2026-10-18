import numpy as np
import matplotlib.pyplot as plt
import itertools

plt.rcParams['font.sans-serif'] = ['SimHei']  # 用于中文显示
plt.rcParams['axes.unicode_minus'] = False    # 解决负号显示问题

# 定义三次方程 - 这里使用一个通用形式，您可以根据需要修改
def cubic_equation(x, q):
    """
    定义参数化的三次方程
    示例方程: x^3 - q*x^2 + 2*x + 1 = 0
    您可以修改此函数以适应您的实际方程
    """
    # 实际方程 - 这里使用示例方程
    return x**3 - (1-q**2)*x**2 - 2*q**2*x +1.375*q**2

# 参数范围和求解设置
q_min = 0
q_max = 2
num_points = 1000
q_values = np.linspace(q_min, q_max, num_points)

# 存储选定分支的解
selected_solution = np.full(num_points, np.nan, dtype=complex)

# 第一个点的初始猜测
initial_guess = complex(0.5, 0.5)  # 初始猜测在0-1范围内

# 求解第一个点 (q_min)
coeffs = [1, -1+q_min**2, -2*q_min**2, 1.375*q_min**2]  # 多项式系数: x^3, x^2, x, 常数项
roots = np.roots(coeffs)

# 选择虚部>=0的解
valid_roots = roots[roots.imag >= 0]
if len(valid_roots) > 0:
    # 选择最接近初始猜测的解
    distances = np.abs(valid_roots - initial_guess)
    selected_idx = np.argmin(distances)
    selected_solution[0] = valid_roots[selected_idx]
else:
    # 如果没有虚部>=0的解，使用实数解
    real_roots = roots[roots.imag == 0]
    if len(real_roots) > 0:
        selected_solution[0] = real_roots[0]

# 使用跟踪方法保持解的连续性
prev_solution = selected_solution[0]

# 求解循环 (从第二个点开始)
for i in range(1, num_points):
    q = q_values[i]
    
    # 计算当前系数
    coeffs = [1, -1+q**2, -2*q**2, 1.375*q**2]
    
    try:
        # 求解当前根
        current_roots = np.roots(coeffs)
        
        # 选择虚部>=0的解
        valid_roots = current_roots[current_roots.imag >= 0]
        
        if len(valid_roots) == 0:
            # 如果没有虚部>=0的解，使用实数解
            real_roots = current_roots[current_roots.imag == 0]
            if len(real_roots) > 0:
                # 选择最接近前一个解的解
                distances = np.abs(real_roots - prev_solution)
                selected_idx = np.argmin(distances)
                current_solution = real_roots[selected_idx]
            else:
                # 如果没有合适的解，使用前一个解
                current_solution = prev_solution
        elif len(valid_roots) == 1:
            current_solution = valid_roots[0]
        else:
            # 多个虚部>=0的解，选择最接近前一个解的解
            distances = np.abs(valid_roots - prev_solution)
            selected_idx = np.argmin(distances)
            current_solution = valid_roots[selected_idx]
        
        # 存储结果
        selected_solution[i] = current_solution
        prev_solution = current_solution
        
    except Exception as e:
        print(f"在 q={q:.3f} 处求解失败: {str(e)}")
        # 保留前一个解
        selected_solution[i] = prev_solution

# 分离实部和虚部
real_parts = np.real(selected_solution)
imag_parts = np.imag(selected_solution)

# 创建高质量图形
plt.figure(figsize=(14, 8))
#plt.style.use('seaborn-darkgrid')  # 使用美观的样式

# 绘制实部曲线
plt.plot(q_values, real_parts, 
         color='#1f77b4', 
         linewidth=3.0, 
         alpha=0.9,
         label=r'$Re(x)$')

# 绘制虚部曲线
plt.plot(q_values, imag_parts, 
         color='#d62728', 
         linewidth=2.5, 
         alpha=0.8,
         linestyle='--',
         label=r'$Im(x)$')

# 添加网格
plt.grid(True, linestyle='--', alpha=0.6)

# 添加零参考线
plt.axhline(0, color='gray', linestyle='-', alpha=0.7, linewidth=0.8)
plt.axvline(0, color='gray', linestyle='-', alpha=0.7, linewidth=0.8)

# 设置坐标轴范围
plt.ylim(-0.1, 1.1)  # 根据实部和虚部都在0-1范围内设置

# 添加标题和标签
plt.title('ion cyclotron instability(temperature anisotropy)', fontsize=16, fontweight='bold', pad=20)
plt.xlabel('q', fontsize=14, fontweight='bold')
plt.ylabel('x', fontsize=14, fontweight='bold')

# 添加图例
plt.legend(loc='best', fontsize=12, frameon=True, shadow=True)

# 添加方程信息
equation_text = r"$x^3 - qx^2 + 2x + 1 = 0$"
# plt.figtext(0.5, 0.01, f"方程: {equation_text}\n只显示虚部≥0的解分支", 
#             ha="center", fontsize=12, 
#             bbox={"facecolor": "white", "alpha": 0.8, "pad": 5})

# 调整布局
plt.tight_layout(rect=[0, 0.03, 1, 0.95])

# 保存高清图像
plt.savefig('cubic_equation_selected_branch.png', dpi=300, bbox_inches='tight')

# 显示图形
plt.show()

# 可选：在复平面上绘制解的变化轨迹
# plt.figure(figsize=(10, 8))

# # 绘制解轨迹
# plt.plot(real_parts, imag_parts, 
#          color='#2ca02c', 
#          linewidth=2.5, 
#          alpha=0.8,
#          label='解轨迹')

# # 标记起点和终点
# plt.scatter(real_parts[0], imag_parts[0], 
#             color='#1f77b4', 
#             s=100, 
#             edgecolor='black', 
#             zorder=3,
#             label=f'起点 (q={q_min})')
# plt.scatter(real_parts[-1], imag_parts[-1], 
#             color='#d62728', 
#             s=100, 
#             edgecolor='black', 
#             marker='s',
#             zorder=3,
#             label=f'终点 (q={q_max})')

# # 添加虚部=0的参考线
# plt.axhline(0, color='black', linewidth=0.8, linestyle='-', alpha=0.7)

# # 美化复平面图
# plt.title(f'复平面上的解轨迹: {equation_text}', fontsize=14)
# plt.xlabel('实部', fontsize=12)
# plt.ylabel('虚部', fontsize=12)
# plt.grid(True, linestyle='--', alpha=0.6)
# plt.axvline(0, color='black', linewidth=0.8)

# # 设置坐标轴范围
# plt.xlim(min(real_parts)-0.1, max(real_parts)+0.1)
# plt.ylim(-0.1, max(imag_parts)+0.1)

# # 添加图例
# plt.legend(fontsize=10, loc='best')

# # 保存图像
# plt.savefig('complex_plane_selected_branch.png', dpi=300, bbox_inches='tight')
# plt.show()

# # 显示解的统计数据
# print("\n解分支的统计数据:")
# print(f"实部范围: {np.min(real_parts):.4f} 到 {np.max(real_parts):.4f}")
# print(f"虚部范围: {np.min(imag_parts):.4f} 到 {np.max(imag_parts):.4f}")
# print(f"平均实部: {np.mean(real_parts):.4f}")
# print(f"平均虚部: {np.mean(imag_parts):.4f}")