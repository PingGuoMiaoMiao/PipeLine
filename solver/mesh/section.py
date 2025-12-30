"""
Section - 截面类
"""


class Section:
    """截面类"""
    
    def __init__(self, section_id: int, section_type: str = "PIPE", **kwargs):
        """
        初始化截面
        
        参数:
            section_id: 截面ID
            section_type: 截面类型（如"PIPE"）
            **kwargs: 截面参数（如diameter, wall_thickness等）
        """
        self.id = section_id
        self.type = section_type
        self.parameters = kwargs
    
    def __repr__(self):
        return f"Section({self.id}, {self.type}, {self.parameters})"
    
    def get(self, key, default=None):
        """获取截面参数"""
        return self.parameters.get(key, default)
    
    def set(self, key, value):
        """设置截面参数"""
        self.parameters[key] = value

