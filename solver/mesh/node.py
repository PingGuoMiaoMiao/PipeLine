"""
Node - 节点类
"""

import numpy as np


class Node:
    """节点类"""
    
    def __init__(self, node_id: int, x: float, y: float, z: float):
        """
        初始化节点
        
        参数:
            node_id: 节点ID
            x, y, z: 节点坐标
        """
        self.id = node_id
        self.x = x
        self.y = y
        self.z = z
        self.displacement = np.zeros(6)  # [UX, UY, UZ, ROTX, ROTY, ROTZ]
    
    def __repr__(self):
        return f"Node({self.id}, ({self.x}, {self.y}, {self.z}))"
    
    def get_position(self):
        """获取节点位置"""
        return np.array([self.x, self.y, self.z])
    
    def get_displacement(self):
        """获取节点位移"""
        return self.displacement

