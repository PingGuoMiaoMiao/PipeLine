"""
Elastic Material - 弹性材料
"""


class ElasticMaterial:
    """弹性材料类"""
    
    def __init__(self, mat_id: int, **kwargs):
        """
        初始化弹性材料
        
        参数:
            mat_id: 材料ID
            **kwargs: 材料参数（E, nu, density, alpha, T_ref等）
        """
        self.mat_id = mat_id
        self.E = kwargs.get('E', 200e9)  # 弹性模量 (Pa)
        self.nu = kwargs.get('nu', 0.3)  # 泊松比
        self.density = kwargs.get('density', 7800.0)  # 密度 (kg/m³)
        self.alpha = kwargs.get('alpha', 1.2e-5)  # 热膨胀系数 (1/°C)
        self.T_ref = kwargs.get('T_ref', 25.0)  # 参考温度 (°C)
        self.G = self.E / (2.0 * (1.0 + self.nu))  # 剪切模量 (Pa)
    
    def __repr__(self):
        return f"ElasticMaterial({self.mat_id}, E={self.E/1e9:.1f}GPa, nu={self.nu})"
    
    def to_dict(self):
        """转换为字典"""
        return {
            'E': self.E,
            'nu': self.nu,
            'density': self.density,
            'alpha': self.alpha,
            'T_ref': self.T_ref
        }

