"""
Shape Function Fourier - 环向傅里叶级数形函数
"""

import numpy as np
from typing import Tuple


class FourierShapeFunction:
    """环向傅里叶级数形函数"""
    
    @staticmethod
    def compute_shape_function(phi: float, n: int, mode: str = 'cos') -> float:
        """
        计算n阶傅里叶形函数值
        
        参数:
            phi: 环向角度 (rad)
            n: 阶数（2=椭圆化, 3=翘曲）
            mode: 'cos' 或 'sin'
        
        返回:
            形函数值
        """
        if mode == 'cos':
            return np.cos(n * phi)
        elif mode == 'sin':
            return np.sin(n * phi)
        else:
            raise ValueError(f"未知模式: {mode}")
    
    @staticmethod
    def compute_radial_displacement(phi: float, 
                                   oval_cos: float, oval_sin: float,
                                   warp_cos: float, warp_sin: float) -> float:
        """
        计算环向径向位移（考虑椭圆化和翘曲）
        
        参数:
            phi: 环向角度 (rad)
            oval_cos: 椭圆化cos(2φ)项系数
            oval_sin: 椭圆化sin(2φ)项系数
            warp_cos: 翘曲cos(3φ)项系数
            warp_sin: 翘曲sin(3φ)项系数
        
        返回:
            径向位移
        """
        delta_r = (oval_cos * np.cos(2 * phi) + oval_sin * np.sin(2 * phi) +
                   warp_cos * np.cos(3 * phi) + warp_sin * np.sin(3 * phi))
        return delta_r
    
    @staticmethod
    def generate_circumferential_points(num_points: int = 20) -> np.ndarray:
        """
        生成环向积分点角度
        
        参数:
            num_points: 环向划分点数
        
        返回:
            角度数组 (rad)
        """
        return np.linspace(0, 2 * np.pi, num_points, endpoint=False)
    
    @staticmethod
    def compute_section_coordinates(radius_mean: float, phi: float,
                                   oval_cos: float = 0.0, oval_sin: float = 0.0,
                                   warp_cos: float = 0.0, warp_sin: float = 0.0) -> Tuple[float, float]:
        """
        计算变形后的截面坐标
        
        参数:
            radius_mean: 平均半径
            phi: 环向角度 (rad)
            oval_cos, oval_sin: 椭圆化系数
            warp_cos, warp_sin: 翘曲系数
        
        返回:
            (x, y) 截面坐标
        """
        # 计算径向位移
        delta_r = FourierShapeFunction.compute_radial_displacement(
            phi, oval_cos, oval_sin, warp_cos, warp_sin)
        
        # 变形后的半径
        r = radius_mean + delta_r
        
        # 截面坐标
        x = r * np.cos(phi)
        y = r * np.sin(phi)
        
        return x, y
