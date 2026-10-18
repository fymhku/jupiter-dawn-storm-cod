import numpy as np
import pandas as pd
from pyrfu import pyrf
from tool import calculate_corotation_velocity
import xarray as xr
def perform_mva(Bx, By, Bz, X, Y, Z):
    """
    执行最小方差分析(MVA)并计算相关参数
    
    参数:
    Bx, By, Bz -- pandas Series, 三个方向的磁场分量
    X, Y, Z -- pandas Series, 三个方向的位置坐标
    
    返回:
    Bl, Bm, Bn -- pandas Series, MVA坐标系下的磁场分量
    Vn -- 法向速度 (km/s)
    eigenvals -- MVA特征值数组
    """
    
    # 1. 组合磁场数据
    
    df = pd.DataFrame({
        'Bx': Bx,
        'By': By,
        'Bz': Bz
    })

    # 2. 删除包含 NaN 的行
    df_clean = df.dropna()
    BB = xr.DataArray(df_clean)


    
    # 2. 执行最小方差分析(MVA)
    BB_LMN, eigenvals, v = pyrf.mva(BB)
    print(f'MVA特征值: {eigenvals}')
    
    BB = np.column_stack((Bx.values, By.values, Bz.values))
    # 3. 计算MVA坐标系分量
    Bl = pd.Series(np.dot(BB, v[:,0]), index=Bx.index, name='Bl')
    Bm = pd.Series(np.dot(BB, v[:,1]), index=Bx.index, name='Bm')
    Bn = pd.Series(np.dot(BB, v[:,2]), index=Bx.index, name='Bn')
    
    # 4. 计算法向速度Vn (基于共转速度投影)
    Vcoro,theta = calculate_corotation_velocity.calculate_corotation_velocity(X, Y)
    Vn = np.abs(np.dot(Vcoro, v[:,2]))
    
    return Bl, Bm, Bn, Vn, eigenvals

