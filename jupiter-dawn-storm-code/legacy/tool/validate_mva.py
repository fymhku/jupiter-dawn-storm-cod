def validate_mva(eigenvals):
    """
    验证 MVA（最小方差分析）特征值有效性条件：
    1. 最大特征值 > 10 * 中间特征值
    2. 中间特征值 > 3 * 最小特征值
    
    参数:
        eigenvals : list[float] 或 np.ndarray
            包含三个特征值的列表/数组（需按降序排列：[λ_max, λ_mid, λ_min]）
    
    返回:
        bool: 所有验证条件满足时返回 True，否则返回 False
    """
    # 确保特征值按降序排列（原MATLAB代码假设输入已排序）
    sorted_vals = sorted(eigenvals, reverse=True)
    λ_max, λ_mid, λ_min = sorted_vals
    
    # 计算特征值比率
    ratio_max_mid = λ_max / λ_mid
    ratio_mid_min = λ_mid / λ_min
    
    # 检查验证条件
    isValid = (ratio_max_mid > 10) and (ratio_mid_min > 3)
    
    # 如果验证失败，打印诊断信息（原MATLAB行为）
    if not isValid:
        print('MVA特征值验证失败：')
        print(f'λ1/λ2 = {ratio_max_mid:.2f} (需>10)')
        print(f'λ2/λ3 = {ratio_mid_min:.2f} (需>3)')
    
    return isValid