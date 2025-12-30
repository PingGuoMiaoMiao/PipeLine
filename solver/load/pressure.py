"""
Pressure Load - 压力载荷
"""

import math
import numpy as np


class PressureLoad:
    """压力载荷类"""
    
    @staticmethod
    def compute_internal_pressure_force(pressure: float, radius_inner: float) -> float:
        """
        计算内压产生的等效轴向力（端盖效应）
        
        参数:
            pressure: 内压 (Pa)
            radius_inner: 内半径 (m)
        
        返回:
            等效轴向力 (N)
        """
        A_internal = math.pi * radius_inner**2
        return pressure * A_internal

