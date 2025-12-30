"""
Integration Point - 积分点类
用于存储积分点的状态变量（应力、应变、塑性应变等）
"""

import numpy as np
from typing import Optional


class IntegrationPoint:
    """积分点状态变量"""
    
    def __init__(self):
        """初始化积分点"""
        # 应力（6×1，Voigt记号）
        self.stress = np.zeros(6)
        
        # 应变（6×1）
        self.strain = np.zeros(6)
        
        # 塑性应变（6×1）
        self.plastic_strain = np.zeros(6)
        
        # 等效塑性应变（标量历史变量）
        self.equivalent_plastic_strain = 0.0
        
        # 是否进入塑性
        self.is_plastic = False
        
        # 蠕变应变（6×1）
        self.creep_strain = np.zeros(6)
        
        # 等效蠕变应变（标量历史变量）
        self.equivalent_creep_strain = 0.0
        
        # 上一增量步的值（用于增量计算）
        self.stress_old = np.zeros(6)
        self.strain_old = np.zeros(6)
        self.plastic_strain_old = np.zeros(6)
        self.equivalent_plastic_strain_old = 0.0
        self.creep_strain_old = np.zeros(6)
        self.equivalent_creep_strain_old = 0.0
    
    def save_state(self):
        """保存当前状态为上一增量步的状态"""
        self.stress_old = self.stress.copy()
        self.strain_old = self.strain.copy()
        self.plastic_strain_old = self.plastic_strain.copy()
        self.equivalent_plastic_strain_old = self.equivalent_plastic_strain
        self.creep_strain_old = self.creep_strain.copy()
        self.equivalent_creep_strain_old = self.equivalent_creep_strain
    
    def update_stress(self, stress: np.ndarray, 
                     equivalent_plastic_strain: float,
                     is_plastic: bool):
        """更新应力和塑性应变"""
        self.stress = stress.copy()
        self.equivalent_plastic_strain = equivalent_plastic_strain
        self.is_plastic = is_plastic
    
    def update_strain(self, strain: np.ndarray):
        """更新应变"""
        self.strain = strain.copy()
    
    def update_plastic_strain(self, plastic_strain: np.ndarray):
        """更新塑性应变"""
        self.plastic_strain = plastic_strain.copy()
    
    def update_creep_strain(self, creep_strain: np.ndarray,
                           equivalent_creep_strain: float):
        """更新蠕变应变"""
        self.creep_strain = creep_strain.copy()
        self.equivalent_creep_strain = equivalent_creep_strain

