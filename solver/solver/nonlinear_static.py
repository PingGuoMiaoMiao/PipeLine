"""
Nonlinear Static Solver - 非线性静力求解器（当前为线性）
支持不同DOF数的单元（PIPE288: 6 DOF/节点，ELBOW290: 10 DOF/节点）
"""

import numpy as np
from typing import Dict, List
from solver.mesh.node import Node


class NonlinearStaticSolver:
    """非线性静力求解器（当前实现为线性求解）"""
    
    def __init__(self, nodes: Dict[int, Node], 
                 elements: Dict,
                 boundary_conditions: Dict[int, List[int]],
                 dof_per_node: int = 6):
        """
        初始化求解器
        
        参数:
            nodes: 节点字典
            elements: 单元字典（可以是PIPE288或ELBOW290）
            boundary_conditions: 边界条件字典（节点ID -> [DOF列表]）
            dof_per_node: 每个节点的DOF数（PIPE288=6，ELBOW290=10）
        """
        self.nodes = nodes
        self.elements = elements
        self.boundary_conditions = boundary_conditions
        self.dof_per_node = dof_per_node
        
        # 创建节点ID到DOF索引的映射
        sorted_node_ids = sorted(nodes.keys())
        self.node_id_to_dof_start = {
            node_id: idx * dof_per_node 
            for idx, node_id in enumerate(sorted_node_ids)
        }
        self.dof_count = len(nodes) * dof_per_node
    
    def assemble_stiffness_matrix(self) -> np.ndarray:
        """组装全局刚度矩阵"""
        K_global = np.zeros((self.dof_count, self.dof_count))
        
        for elem in self.elements.values():
            K_elem = elem.compute_stiffness_matrix()
            
            # 获取单元节点的DOF索引
            # 根据单元类型确定节点数和DOF数
            if hasattr(elem, 'node1'):
                # PIPE288（2节点，6 DOF/节点）
                if hasattr(elem, 'node2') and not hasattr(elem, 'node3'):
                    dof1_start = self.node_id_to_dof_start[elem.node1.id]
                    dof2_start = self.node_id_to_dof_start[elem.node2.id]
                    
                    # 组装（12×12）
                    for i in range(6):
                        for j in range(6):
                            K_global[dof1_start + i, dof1_start + j] += K_elem[i, j]
                            K_global[dof1_start + i, dof2_start + j] += K_elem[i, 6 + j]
                            K_global[dof2_start + i, dof1_start + j] += K_elem[6 + i, j]
                            K_global[dof2_start + i, dof2_start + j] += K_elem[6 + i, 6 + j]
                
                # ELBOW290（3节点，10 DOF/节点）
                elif hasattr(elem, 'node2') and hasattr(elem, 'node3'):
                    dof1_start = self.node_id_to_dof_start[elem.node1.id]
                    dof2_start = self.node_id_to_dof_start[elem.node2.id]
                    dof3_start = self.node_id_to_dof_start[elem.node3.id]
                    
                    # 组装（30×30）
                    for i in range(10):
                        for j in range(10):
                            K_global[dof1_start + i, dof1_start + j] += K_elem[i, j]
                            K_global[dof1_start + i, dof2_start + j] += K_elem[i, 10 + j]
                            K_global[dof1_start + i, dof3_start + j] += K_elem[i, 20 + j]
                            K_global[dof2_start + i, dof1_start + j] += K_elem[10 + i, j]
                            K_global[dof2_start + i, dof2_start + j] += K_elem[10 + i, 10 + j]
                            K_global[dof2_start + i, dof3_start + j] += K_elem[10 + i, 20 + j]
                            K_global[dof3_start + i, dof1_start + j] += K_elem[20 + i, j]
                            K_global[dof3_start + i, dof2_start + j] += K_elem[20 + i, 10 + j]
                            K_global[dof3_start + i, dof3_start + j] += K_elem[20 + i, 20 + j]
        
        return K_global
    
    def apply_boundary_conditions(self, K: np.ndarray, F: np.ndarray,
                                  penalty: float = 1e12):
        """
        应用边界条件（置大数法）
        
        注意：边界条件只应用于标准DOF（前6个），截面变形DOF不受约束
        """
        for node_id, dof_list in self.boundary_conditions.items():
            if node_id in self.node_id_to_dof_start:
                dof_start = self.node_id_to_dof_start[node_id]
                for dof in dof_list:
                    # 只约束标准DOF（0-5）
                    if dof < 6 and dof_start + dof < len(K):
                        dof_idx = dof_start + dof
                        K[dof_idx, :] = 0
                        K[:, dof_idx] = 0
                        K[dof_idx, dof_idx] = penalty
                        F[dof_idx] = 0
    
    def solve(self, load_vector: np.ndarray) -> np.ndarray:
        """求解线性方程组 K * U = F"""
        # 组装全局刚度矩阵
        K_global = self.assemble_stiffness_matrix()
        
        # 应用边界条件
        self.apply_boundary_conditions(K_global, load_vector)
        
        # 求解
        try:
            U = np.linalg.solve(K_global, load_vector)
        except np.linalg.LinAlgError:
            print("警告: 矩阵求解失败，尝试使用最小二乘")
            U = np.linalg.lstsq(K_global, load_vector, rcond=None)[0]
        
        return U
