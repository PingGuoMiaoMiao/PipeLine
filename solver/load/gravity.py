"""
Gravity Load - 重力载荷
"""

import numpy as np


class GravityLoad:
    """重力载荷类"""
    
    @staticmethod
    def compute_distributed_load(density: float, area: float, gravity: np.ndarray) -> float:
        """
        计算单位长度的重力载荷
        
        参数:
            density: 材料密度 (kg/m³)
            area: 横截面积 (m²)
            gravity: 重力加速度向量 (m/s²)
        
        返回:
            单位长度载荷 (N/m)
        """
        gravity_magnitude = np.linalg.norm(gravity)
        return density * area * gravity_magnitude

