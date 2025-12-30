"""
Displacement Boundary - 位移边界条件
"""

from typing import Dict, List


class DisplacementBoundary:
    """位移边界条件类"""
    
    def __init__(self):
        self.boundary_conditions: Dict[int, List[int]] = {}
    
    def add(self, node_id: int, dof: int):
        """
        添加边界条件
        
        参数:
            node_id: 节点ID
            dof: 自由度（0=UX, 1=UY, 2=UZ, 3=ROTX, 4=ROTY, 5=ROTZ）
        """
        if node_id not in self.boundary_conditions:
            self.boundary_conditions[node_id] = []
        if dof not in self.boundary_conditions[node_id]:
            self.boundary_conditions[node_id].append(dof)
    
    def apply_to_matrix(self, K: 'np.ndarray', dof_start_map: Dict[int, int], 
                       penalty: float = 1e12):
        """
        应用边界条件到刚度矩阵（置大数法）
        
        参数:
            K: 全局刚度矩阵
            dof_start_map: 节点ID到DOF起始索引的映射
            penalty: 惩罚系数（默认1e12）
        """
        import numpy as np
        
        for node_id, dof_list in self.boundary_conditions.items():
            if node_id in dof_start_map:
                dof_start = dof_start_map[node_id]
                for dof in dof_list:
                    dof_idx = dof_start + dof
                    K[dof_idx, :] = 0
                    K[:, dof_idx] = 0
                    K[dof_idx, dof_idx] = penalty

