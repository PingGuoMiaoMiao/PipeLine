"""
Creep Strain Hardening - 各向同性应变硬化蠕变材料模型
蠕变率公式：ε̇_cr = C1 * σ^C2 * t^C3 * exp(-C4/T)
"""

import numpy as np
from typing import Optional


class CreepStrainHardeningMaterial:
    """
    各向同性应变硬化蠕变材料模型
    
    蠕变率公式：
    ε̇_cr = C1 * σ^C2 * t^C3 * exp(-C4/T)
    
    其中：
    - C1: 蠕变系数
    - C2: 应力指数
    - C3: 时间指数
    - C4: 温度相关参数（通常为Q/R，Q为激活能，R为气体常数）
    - σ: Mises等效应力
    - t: 时间
    - T: 温度（绝对温度，K）
    """
    
    def __init__(self, C1: float, C2: float, C3: float, C4: float):
        """
        初始化蠕变材料
        
        参数:
            C1: 蠕变系数
            C2: 应力指数（通常≥1）
            C3: 时间指数（通常<1）
            C4: 温度相关参数（通常为Q/R）
        """
        self.C1 = C1
        self.C2 = C2
        self.C3 = C3
        self.C4 = C4
    
    def compute_creep_strain_rate(self, stress: np.ndarray,
                                  equivalent_creep_strain: float,
                                  time: float,
                                  temperature: float) -> float:
        """
        计算蠕变应变率
        
        参数:
            stress: 应力向量（6×1，Voigt记号）
            equivalent_creep_strain: 等效蠕变应变（历史变量）
            time: 当前时间
            temperature: 温度（绝对温度，K）
        
        返回:
            等效蠕变应变率
        """
        # 计算Mises等效应力
        mises_stress = self._compute_mises_stress(stress)
        
        # 防止数值问题
        if mises_stress < 1e-6:
            return 0.0
        
        # 时间项（使用实际时间）
        # 对于Strain Hardening模型，时间项使用实际时间
        # ε̇_cr = C1 * σ^C2 * t^C3 * exp(-C4/T)
        if time < 1e-10:
            time = 1e-10  # 避免零的幂次
        
        # 蠕变应变率公式
        # ε̇_cr = C1 * σ^C2 * t^C3 * exp(-C4/T)
        # 时间项使用实际时间
        time_term = time ** self.C3 if self.C3 != 0 else 1.0
        
        # 温度项
        if temperature > 1e-6:
            temp_term = np.exp(-self.C4 / temperature)
        else:
            temp_term = 0.0
        
        # 蠕变应变率
        creep_rate = self.C1 * (mises_stress ** self.C2) * time_term * temp_term
        
        return creep_rate
    
    def _compute_mises_stress(self, stress: np.ndarray) -> float:
        """
        计算Mises等效应力
        
        参数:
            stress: 应力向量（6×1，Voigt记号: [σ11, σ22, σ33, σ12, σ13, σ23]）
        
        返回:
            Mises等效应力
        """
        # 应力偏量
        mean_stress = np.mean(stress[0:3])
        s = np.zeros(6)
        s[0:3] = stress[0:3] - mean_stress
        s[3:6] = stress[3:6]
        
        # Mises应力: σ_mises = sqrt(3/2 * s_ij * s_ij)
        mises_sq = (3.0 / 2.0) * (
            s[0]**2 + s[1]**2 + s[2]**2 +
            2.0 * (s[3]**2 + s[4]**2 + s[5]**2)
        )
        
        return np.sqrt(mises_sq)
    
    def compute_creep_strain_increment(self, stress: np.ndarray,
                                      equivalent_creep_strain: float,
                                      time: float,
                                      time_increment: float,
                                      temperature: float,
                                      method: str = 'implicit') -> float:
        """
        计算蠕变应变增量
        
        参数:
            stress: 当前应力
            equivalent_creep_strain: 当前等效蠕变应变
            time: 当前时间
            time_increment: 时间增量
            temperature: 温度（绝对温度，K）
            method: 积分方法（'explicit'或'implicit'）
        
        返回:
            等效蠕变应变增量
        """
        if method == 'explicit':
            # 显式Euler方法
            creep_rate = self.compute_creep_strain_rate(
                stress, equivalent_creep_strain, time, temperature)
            delta_creep_strain = creep_rate * time_increment
        else:
            # 隐式方法（简化的向后Euler）
            # 使用当前应力预测
            creep_rate = self.compute_creep_strain_rate(
                stress, equivalent_creep_strain, time, temperature)
            delta_creep_strain = creep_rate * time_increment
        
        return delta_creep_strain
    
    def compute_creep_strain_direction(self, stress: np.ndarray) -> np.ndarray:
        """
        计算蠕变应变方向（单位向量）
        
        参数:
            stress: 应力向量
        
        返回:
            蠕变应变方向（6×1）
        """
        # 应力偏量
        mean_stress = np.mean(stress[0:3])
        s = np.zeros(6)
        s[0:3] = stress[0:3] - mean_stress
        s[3:6] = stress[3:6]
        
        # Mises应力
        mises = self._compute_mises_stress(stress)
        
        if mises > 1e-10:
            # 归一化的应力偏量方向（3/2因子来自Mises定义）
            direction = (3.0 / 2.0) * s / mises
        else:
            direction = np.zeros(6)
        
        return direction
    
    def compute_creep_derivative(self, stress: np.ndarray,
                                equivalent_creep_strain: float,
                                time: float,
                                temperature: float,
                                E: float) -> float:
        """
        计算蠕变对刚度的导数（用于隐式积分）
        
        参数:
            stress: 当前应力
            equivalent_creep_strain: 当前等效蠕变应变
            time: 当前时间
            temperature: 温度
            E: 弹性模量
        
        返回:
            蠕变模量（用于切线刚度矩阵）
        """
        mises_stress = self._compute_mises_stress(stress)
        
        if mises_stress < 1e-6 or equivalent_creep_strain < 1e-10:
            return 0.0
        
        # 蠕变率
        creep_rate = self.compute_creep_strain_rate(
            stress, equivalent_creep_strain, time, temperature)
        
        if creep_rate < 1e-10:
            return 0.0
        
        # 蠕变模量的近似（简化处理）
        # H_cr ≈ dσ/dε_cr
        # 这里使用简化公式
        if self.C3 != 0 and equivalent_creep_strain > 1e-10:
            H_cr = mises_stress / (self.C3 * equivalent_creep_strain)
        else:
            H_cr = E * 1e-6  # 很小的值
        
        return H_cr
