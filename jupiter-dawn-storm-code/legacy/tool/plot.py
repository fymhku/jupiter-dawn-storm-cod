import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np

def plot_magnetic_field_components(datetime, bx, by, bz,
                                   title_prefix='Magnetic Field',
                                   save_figure=False,
                                   save_path=None):
    """
    将磁场三分量 (B_r, B_theta, B_phi) 绘制在同一幅图中，用不同颜色区分。

    参数:
    datetime -- 时间序列 (datetime对象数组)
    bx, by, bz -- 磁场分量数据 (nT)
    title_prefix -- 图表标题前缀 (默认为'Magnetic Field')
    save_figure -- 是否保存图表 (默认为False)
    save_path -- 图表保存完整路径 (若save_figure为True且为None，则自动生成)
    """
    # 创建DataFrame便于处理
    df = pd.DataFrame({
        'Datetime': datetime,
        'B_r': bx,
        'B_theta': by,
        'B_phi': bz
    })

    # 创建图形和坐标轴
    fig, ax = plt.subplots(figsize=(12, 6))

    # 绘制三个分量，使用不同颜色并添加标签
    ax.plot(df['Datetime'], df['B_r'], 'r-', linewidth=1.5, label=r'$B_r$')
    ax.plot(df['Datetime'], df['B_theta'], 'g-', linewidth=1.5, label=r'$B_\theta$')
    ax.plot(df['Datetime'], df['B_phi'], 'b-', linewidth=1.5, label=r'$B_\phi$')

    # 设置标题和轴标签
    ax.set_title(f'{title_prefix} - Components', fontsize=14)
    ax.set_ylabel('Magnetic Field (nT)', fontsize=12)
    ax.set_xlabel('Time (UTC)', fontsize=12)

    # 添加图例
    ax.legend(loc='best', fontsize=10)

    # 添加网格
    ax.grid(True, linestyle='--', alpha=0.6)

    # 设置时间轴格式（只显示时分秒）
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    ax.xaxis.set_major_locator(mdates.AutoDateLocator())
    ax.set_xlim(df['Datetime'].min(), df['Datetime'].max())

    # 自动旋转日期标签
    fig.autofmt_xdate()

    # 紧凑布局
    fig.tight_layout()

    # 保存图片
    if save_figure:
        if save_path is None:
            save_path = f"{title_prefix.replace(' ', '_')}_components.png"
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"图表已保存至: {save_path}")

    plt.show()
    plt.close()