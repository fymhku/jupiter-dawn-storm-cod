import numpy as np
import pandas as pd

def calculate_corotation_velocity(X, Y):
    """
    计算木星共转速度模型（假设自转周期10小时）
    
    参数:
        X, Y: Pandas Series 包含位置坐标 (单位: Rj)
    
    返回:
        Vcoro: 共转速度向量 [Vx, Vy, Vz] (Rj/s)
        theta: 方位角 (弧度)
    """
    # 确保输入是 Pandas Series
    if not isinstance(X, pd.Series) or not isinstance(Y, pd.Series):
        raise TypeError("X and Y must be Pandas Series")
    
    # 检查索引是否对齐
    if not X.index.equals(Y.index):
        # 尝试对齐索引
        common_idx = X.index.intersection(Y.index)
        if len(common_idx) == 0:
            raise ValueError("X and Y have no common index values")
        X = X.loc[common_idx]
        Y = Y.loc[common_idx]
    
    # 计算平均位置
    X_mean = X.mean()
    Y_mean = Y.mean()
    
    # 计算平均赤道面距离 (Rj)
    R = np.sqrt(X_mean**2 + Y_mean**2)
    
    # 计算方位角（使用arctan2处理所有象限）
    theta = -np.pi - np.arcsin(Y_mean/R)
    
    # 计算角速度 (rad/s)
    omega = 2 * np.pi / (10 * 3600)  # 10小时周期
    
    # 计算共转速度向量 (Rj/s)
    Vcoro = np.array([
        -omega * 20 * np.sin(theta),  # Vx
        omega * 20 * np.cos(theta),   # Vy
        0                            # Vz
    ])
    
    return Vcoro, theta