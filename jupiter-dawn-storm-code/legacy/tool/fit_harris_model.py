import numpy as np
import pandas as pd
from scipy import signal
from tool import convert_time_to_seconds
def fit_harris_model(Bl, Bn, fit_tint, Vn):
    """
    拟合Harris电流片模型并计算曲率
    
    参数:
    Bl -- pandas Series, 切向磁场分量 (时间索引为datetime)
    Bn -- pandas Series, 法向磁场分量 (时间索引为datetime)
    fit_tint -- 时间区间 [start, end], 用于拟合
    Vn -- 法向速度 (km/s)
    
    返回:
    params -- 字典, 包含拟合参数
    B_fit_T -- pandas Series, 拟合的切向磁场
    Curv -- pandas Series, 曲率时间序列
    curv0 -- 中心点曲率
    """
    
    # 1. 提取拟合时间段数据
    Bl_fit = Bl
    Bn_fit = Bn
    
    t_seconds = convert_time_to_seconds.convert_time_strings_to_seconds("2017-05-01",Bl_fit.index)
    
    # 3. 计算初始参数
    B0 = (abs(Bl_fit.iloc[0]) + abs(Bl_fit.iloc[-1])) / 2
    min_idx = np.argmin(np.abs(Bl_fit.values))
    tc = t_seconds[min_idx]
    bbn = np.mean(Bn_fit.values) 
     # 使用平均值而不是中心点值
    
    t_seconds=np.array(t_seconds)
    # 4. 迭代优化alpha
    alpha = -0.99
    SS = float('inf')
    
    while SS >= 0.1:
        # 计算Lu和L
        Lu = (t_seconds[-1] - t_seconds[0]) #/ (2 * np.arctanh(alpha))
        L = Vn * Lu
        
        # Harris模型拟合
        B_fit = B0 * np.tanh((t_seconds - tc) / Lu)
        
        # 计算标准化均方根误差
        SS = np.sqrt(np.mean((B_fit - Bl_fit.values)**2)) / B0
        
        # 检查是否满足条件
        if SS >= 1.5:
            alpha += 0.005
            if alpha >= 1:
                raise ValueError('alpha超过合理范围')
    
    # 5. 创建拟合磁场的时间序列
    B_fit_T = pd.Series(B_fit, index=Bl_fit.index, name='B_fit')
    
    # 6. 计算曲率
    tanh_term = np.tanh((t_seconds - tc) / Lu)
    cosh_term = np.cosh((t_seconds - tc) / Lu)
    
    # 曲率计算公式
    numerator = np.abs(B0)
    denominator = bbn * L * cosh_term**2 * \
                (1 + (B0**2 * tanh_term**2 / bbn**2)**1.5)
    Curv = numerator / np.abs(denominator)
    
    # 中心点曲率 (t = tc)
    curv0 = np.abs(B0 / (bbn * L))
    
    print(f'中心曲率: {curv0:.4e}')
    
    # 7. 创建曲率时间序列
    Curv = pd.Series(Curv, index=Bl_fit.index, name='Curvature')
    
    # 8. 打包参数
    params = {
        'B0': B0,
        'tc': tc,
        'bbn': bbn,
        'alpha': alpha,
        'SS': SS,
        'L': L,
        'Lu': Lu
    }
    
    return params, B_fit_T, Curv, curv0