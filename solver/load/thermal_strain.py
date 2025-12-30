"""
Thermal Strain Load - 热应变载荷
"""


class ThermalStrainLoad:
    """热应变载荷类"""
    
    @staticmethod
    def compute_thermal_force(E: float, A: float, alpha: float, 
                             temperature: float, ref_temperature: float) -> float:
        """
        计算热膨胀产生的等效轴向力
        
        参数:
            E: 弹性模量 (Pa)
            A: 横截面积 (m²)
            alpha: 热膨胀系数 (1/°C)
            temperature: 当前温度 (°C)
            ref_temperature: 参考温度 (°C)
        
        返回:
            等效热载荷 (N)
        """
        delta_T = temperature - ref_temperature
        return E * A * alpha * delta_T

