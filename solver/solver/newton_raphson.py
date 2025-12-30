"""
Newton-Raphson Solver - Newton-Raphson非线性迭代求解器
支持弹塑性分析
"""

import numpy as np
from typing import Dict, List, Optional, Tuple
from solver.mesh.node import Node


class NewtonRaphsonSolver:
    """Newton-Raphson非线性迭代求解器"""
    
    def __init__(self, nodes: Dict[int, Node],
                 elements: Dict,
                 boundary_conditions: Dict[int, List[int]],
                 dof_per_node: int = 6,
                 max_iterations: int = 20,
                 tolerance: float = 1e-6):
        """
        初始化Newton-Raphson求解器
        
        参数:
            nodes: 节点字典
            elements: 单元字典
            boundary_conditions: 边界条件字典
            dof_per_node: 每个节点的DOF数
            max_iterations: 最大迭代次数
            tolerance: 收敛容差
        """
        self.nodes = nodes
        self.elements = elements
        self.boundary_conditions = boundary_conditions
        self.dof_per_node = dof_per_node
        self.max_iterations = max_iterations
        self.tolerance = tolerance
        
        # 创建节点ID到DOF索引的映射
        sorted_node_ids = sorted(nodes.keys())
        self.node_id_to_dof_start = {
            node_id: idx * dof_per_node 
            for idx, node_id in enumerate(sorted_node_ids)
        }
        self.dof_count = len(nodes) * dof_per_node
        
        # 当前位移和增量
        self.displacements = np.zeros(self.dof_count)
        self.displacement_increment = np.zeros(self.dof_count)
    
    def solve(self, load_vector: np.ndarray,
              incremental: bool = True,
              load_factor: float = 1.0) -> Tuple[np.ndarray, Dict]:
        """
        求解非线性方程组
        
        参数:
            load_vector: 载荷向量
            incremental: 是否使用增量加载
            load_factor: 载荷因子（增量加载时使用）
        
        返回:
            (displacements, info): 位移和求解信息
        """
        info = {
            'converged': False,
            'iterations': 0,
            'residual_norms': [],
            'convergence_history': []
        }
        
        # 缩放载荷向量
        F_ext = load_vector * load_factor
        
        # 初始位移
        U = self.displacements.copy()
        
        # Newton-Raphson迭代
        for iteration in range(self.max_iterations):
            # 组装残差向量 R = F_ext - F_int
            R = self._compute_residual(U, F_ext)
            
            # 计算残差范数
            residual_norm = np.linalg.norm(R)
            info['residual_norms'].append(residual_norm)
            
            # 检查收敛
            if residual_norm < self.tolerance:
                info['converged'] = True
                info['iterations'] = iteration + 1
                break
            
            # 组装切线刚度矩阵
            K_tangent = self._assemble_tangent_stiffness(U)
            
            # 应用边界条件
            self._apply_boundary_conditions(K_tangent, R, penalty=1e12)
            
            # 求解线性方程组: K_tan * ΔU = R
            try:
                delta_U = np.linalg.solve(K_tangent, R)
            except np.linalg.LinAlgError:
                print(f"警告: 迭代{iteration+1}，矩阵求解失败")
                delta_U = np.linalg.lstsq(K_tangent, R, rcond=None)[0]
            
            # 更新位移
            U += delta_U
            
            # 输出收敛信息
            print(f"  迭代 {iteration+1}: 残差范数 = {residual_norm:.6e}")
            info['convergence_history'].append({
                'iteration': iteration + 1,
                'residual_norm': residual_norm,
                'delta_U_norm': np.linalg.norm(delta_U)
            })
        
        if not info['converged']:
            print(f"警告: 未收敛（残差 = {info['residual_norms'][-1]:.6e}）")
        
        # 更新位移
        self.displacements = U.copy()
        if incremental:
            self.displacement_increment = U - self.displacements
        
        return U, info
    
    def _compute_residual(self, displacements: np.ndarray,
                         external_load: np.ndarray) -> np.ndarray:
        """
        计算残差向量 R = F_ext - F_int
        
        参数:
            displacements: 当前位移
            external_load: 外部载荷
        
        返回:
            残差向量
        """
        # 组装内部力向量
        F_int = self._assemble_internal_force(displacements)
        
        # 残差 = 外部力 - 内部力
        residual = external_load - F_int
        
        return residual
    
    def _assemble_internal_force(self, displacements: np.ndarray) -> np.ndarray:
        """
        组装内部力向量
        
        参数:
            displacements: 当前位移
        
        返回:
            内部力向量
        """
        F_int = np.zeros(self.dof_count)
        
        for elem in self.elements.values():
            # 获取单元节点的位移
            node_displacements = self._get_element_displacements(elem, displacements)
            
            # 计算单元内部力
            F_elem = elem.compute_internal_force(node_displacements)
            
            # 组装到全局向量
            self._assemble_element_force(elem, F_elem, F_int)
        
        return F_int
    
    def _assemble_tangent_stiffness(self, displacements: np.ndarray) -> np.ndarray:
        """
        组装切线刚度矩阵
        
        参数:
            displacements: 当前位移
        
        返回:
            切线刚度矩阵
        """
        K_tangent = np.zeros((self.dof_count, self.dof_count))
        
        for elem in self.elements.values():
            # 获取单元节点的位移
            node_displacements = self._get_element_displacements(elem, displacements)
            
            # 计算单元切线刚度矩阵
            K_elem = elem.compute_tangent_stiffness(node_displacements)
            
            # 组装到全局矩阵
            self._assemble_element_stiffness(elem, K_elem, K_tangent)
        
        return K_tangent
    
    def _get_element_displacements(self, elem, global_displacements: np.ndarray) -> np.ndarray:
        """获取单元节点的位移"""
        if hasattr(elem, 'node1'):
            # PIPE288（2节点）
            if hasattr(elem, 'node2') and not hasattr(elem, 'node3'):
                dof1_start = self.node_id_to_dof_start[elem.node1.id]
                dof2_start = self.node_id_to_dof_start[elem.node2.id]
                node_disp = np.zeros(12)
                node_disp[0:6] = global_displacements[dof1_start:dof1_start+6]
                node_disp[6:12] = global_displacements[dof2_start:dof2_start+6]
                return node_disp
            # ELBOW290（3节点）
            elif hasattr(elem, 'node2') and hasattr(elem, 'node3'):
                dof1_start = self.node_id_to_dof_start[elem.node1.id]
                dof2_start = self.node_id_to_dof_start[elem.node2.id]
                dof3_start = self.node_id_to_dof_start[elem.node3.id]
                node_disp = np.zeros(30)
                node_disp[0:10] = global_displacements[dof1_start:dof1_start+10]
                node_disp[10:20] = global_displacements[dof2_start:dof2_start+10]
                node_disp[20:30] = global_displacements[dof3_start:dof3_start+10]
                return node_disp
        
        return np.zeros(12)  # 默认
    
    def _assemble_element_force(self, elem, F_elem: np.ndarray, F_global: np.ndarray):
        """组装单元力向量到全局向量"""
        if hasattr(elem, 'node1'):
            if hasattr(elem, 'node2') and not hasattr(elem, 'node3'):
                # PIPE288
                dof1_start = self.node_id_to_dof_start[elem.node1.id]
                dof2_start = self.node_id_to_dof_start[elem.node2.id]
                F_global[dof1_start:dof1_start+6] += F_elem[0:6]
                F_global[dof2_start:dof2_start+6] += F_elem[6:12]
            elif hasattr(elem, 'node2') and hasattr(elem, 'node3'):
                # ELBOW290
                dof1_start = self.node_id_to_dof_start[elem.node1.id]
                dof2_start = self.node_id_to_dof_start[elem.node2.id]
                dof3_start = self.node_id_to_dof_start[elem.node3.id]
                F_global[dof1_start:dof1_start+10] += F_elem[0:10]
                F_global[dof2_start:dof2_start+10] += F_elem[10:20]
                F_global[dof3_start:dof3_start+10] += F_elem[20:30]
    
    def _assemble_element_stiffness(self, elem, K_elem: np.ndarray, K_global: np.ndarray):
        """组装单元刚度矩阵到全局矩阵"""
        if hasattr(elem, 'node1'):
            if hasattr(elem, 'node2') and not hasattr(elem, 'node3'):
                # PIPE288
                dof1_start = self.node_id_to_dof_start[elem.node1.id]
                dof2_start = self.node_id_to_dof_start[elem.node2.id]
                for i in range(6):
                    for j in range(6):
                        K_global[dof1_start + i, dof1_start + j] += K_elem[i, j]
                        K_global[dof1_start + i, dof2_start + j] += K_elem[i, 6 + j]
                        K_global[dof2_start + i, dof1_start + j] += K_elem[6 + i, j]
                        K_global[dof2_start + i, dof2_start + j] += K_elem[6 + i, 6 + j]
            elif hasattr(elem, 'node2') and hasattr(elem, 'node3'):
                # ELBOW290
                dof1_start = self.node_id_to_dof_start[elem.node1.id]
                dof2_start = self.node_id_to_dof_start[elem.node2.id]
                dof3_start = self.node_id_to_dof_start[elem.node3.id]
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
    
    def _apply_boundary_conditions(self, K: np.ndarray, F: np.ndarray,
                                   penalty: float = 1e12):
        """应用边界条件（置大数法）"""
        for node_id, dof_list in self.boundary_conditions.items():
            if node_id in self.node_id_to_dof_start:
                dof_start = self.node_id_to_dof_start[node_id]
                for dof in dof_list:
                    if dof < 6 and dof_start + dof < len(K):
                        dof_idx = dof_start + dof
                        K[dof_idx, :] = 0
                        K[:, dof_idx] = 0
                        K[dof_idx, dof_idx] = penalty
                        F[dof_idx] = 0
    
    def reset(self):
        """重置求解器状态"""
        self.displacements = np.zeros(self.dof_count)
        self.displacement_increment = np.zeros(self.dof_count)

