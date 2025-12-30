"""
Creep Time Integration - 蠕变时间积分求解器
支持蠕变-弹性-塑性耦合分析
"""

import numpy as np
from typing import Dict, List, Optional, Tuple
from solver.mesh.node import Node
from solver.solver.newton_raphson import NewtonRaphsonSolver


class CreepTimeIntegrationSolver:
    """蠕变时间积分求解器"""
    
    def __init__(self, nodes: Dict[int, Node],
                 elements: Dict,
                 boundary_conditions: Dict[int, List[int]],
                 dof_per_node: int = 6,
                 max_iterations: int = 20,
                 tolerance: float = 1e-6):
        """
        初始化蠕变时间积分求解器
        
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
        
        # 当前位移和时间
        self.displacements = np.zeros(self.dof_count)
        self.current_time = 0.0
        
        # 时间历史
        self.time_history = []
        self.displacement_history = []
    
    def solve_time_step(self, load_vector: np.ndarray,
                       time: float,
                       time_increment: float,
                       temperature: float,
                       newton_raphson_solver: Optional[NewtonRaphsonSolver] = None) -> Tuple[np.ndarray, Dict]:
        """
        求解一个时间步
        
        参数:
            load_vector: 载荷向量
            time: 当前时间
            time_increment: 时间增量
            temperature: 温度（绝对温度，K）
            newton_raphson_solver: Newton-Raphson求解器（可选）
        
        返回:
            (displacements, info): 位移和求解信息
        """
        # 如果没有提供求解器，创建一个
        if newton_raphson_solver is None:
            newton_raphson_solver = NewtonRaphsonSolver(
                self.nodes, self.elements, self.boundary_conditions,
                self.dof_per_node, self.max_iterations, self.tolerance)
        
        # 保存积分点状态
        self._save_integration_points_state()
        
        # 迭代求解（考虑蠕变）
        info = {
            'converged': False,
            'iterations': 0,
            'residual_norms': [],
            'time': time,
            'time_increment': time_increment
        }
        
        U = self.displacements.copy()
        
        for iteration in range(self.max_iterations):
            # 更新单元中的蠕变应变
            self._update_creep_strain(U, time, time_increment, temperature)
            
            # 计算残差
            R = self._compute_residual(U, load_vector)
            
            # 计算残差范数
            residual_norm = np.linalg.norm(R)
            info['residual_norms'].append(residual_norm)
            
            # 检查收敛
            if residual_norm < self.tolerance:
                info['converged'] = True
                info['iterations'] = iteration + 1
                break
            
            # 组装切线刚度矩阵
            K_tangent = self._assemble_tangent_stiffness(U, time, temperature)
            
            # 应用边界条件
            self._apply_boundary_conditions(K_tangent, R, penalty=1e12)
            
            # 求解线性方程组
            try:
                delta_U = np.linalg.solve(K_tangent, R)
            except np.linalg.LinAlgError:
                print(f"警告: 迭代{iteration+1}，矩阵求解失败")
                delta_U = np.linalg.lstsq(K_tangent, R, rcond=None)[0]
            
            # 更新位移
            U += delta_U
            
            # 输出收敛信息
            print(f"  迭代 {iteration+1}: 残差范数 = {residual_norm:.6e}")
        
        if not info['converged']:
            print(f"警告: 未收敛（残差 = {info['residual_norms'][-1]:.6e}）")
        
        # 更新位移
        self.displacements = U.copy()
        self.current_time = time + time_increment
        
        # 保存历史
        self.time_history.append(self.current_time)
        self.displacement_history.append(U.copy())
        
        return U, info
    
    def _save_integration_points_state(self):
        """保存积分点状态"""
        for elem in self.elements.values():
            if hasattr(elem, 'integration_points'):
                for ip in elem.integration_points:
                    ip.save_state()
    
    def _update_creep_strain(self, displacements: np.ndarray,
                            time: float,
                            time_increment: float,
                            temperature: float):
        """更新蠕变应变"""
        for elem in self.elements.values():
            if hasattr(elem, 'update_creep_strain'):
                # 新版本需要node_id_to_dof_start作为第一个参数
                elem.update_creep_strain(self.node_id_to_dof_start, displacements,
                                       time, time_increment, temperature)
    
    def _compute_residual(self, displacements: np.ndarray,
                         external_load: np.ndarray) -> np.ndarray:
        """计算残差向量"""
        # 组装内部力向量
        F_int = np.zeros(self.dof_count)
        
        for elem in self.elements.values():
            if hasattr(elem, 'compute_internal_force'):
                node_displacements = self._get_element_displacements(elem, displacements)
                F_elem = elem.compute_internal_force(node_displacements)
                self._assemble_element_force(elem, F_elem, F_int)
        
        # 残差 = 外部力 - 内部力
        residual = external_load - F_int
        
        return residual
    
    def _assemble_tangent_stiffness(self, displacements: np.ndarray,
                                   time: float,
                                   temperature: float) -> np.ndarray:
        """组装切线刚度矩阵"""
        K_tangent = np.zeros((self.dof_count, self.dof_count))
        
        for elem in self.elements.values():
            if hasattr(elem, 'compute_tangent_stiffness'):
                node_displacements = self._get_element_displacements(elem, displacements)
                K_elem = elem.compute_tangent_stiffness(node_displacements)
                self._assemble_element_stiffness(elem, K_elem, K_tangent)
        
        return K_tangent
    
    def _get_element_displacements(self, elem, global_displacements: np.ndarray) -> np.ndarray:
        """获取单元节点的位移"""
        if hasattr(elem, 'node1'):
            if hasattr(elem, 'node2') and not hasattr(elem, 'node3'):
                # PIPE288
                dof1_start = self.node_id_to_dof_start[elem.node1.id]
                dof2_start = self.node_id_to_dof_start[elem.node2.id]
                node_disp = np.zeros(12)
                node_disp[0:6] = global_displacements[dof1_start:dof1_start+6]
                node_disp[6:12] = global_displacements[dof2_start:dof2_start+6]
                return node_disp
            elif hasattr(elem, 'node2') and hasattr(elem, 'node3'):
                # ELBOW290
                dof1_start = self.node_id_to_dof_start[elem.node1.id]
                dof2_start = self.node_id_to_dof_start[elem.node2.id]
                dof3_start = self.node_id_to_dof_start[elem.node3.id]
                node_disp = np.zeros(30)
                node_disp[0:10] = global_displacements[dof1_start:dof1_start+10]
                node_disp[10:20] = global_displacements[dof2_start:dof2_start+10]
                node_disp[20:30] = global_displacements[dof3_start:dof3_start+10]
                return node_disp
        return np.zeros(12)
    
    def _assemble_element_force(self, elem, F_elem: np.ndarray, F_global: np.ndarray):
        """组装单元力向量"""
        if hasattr(elem, 'node1'):
            if hasattr(elem, 'node2') and not hasattr(elem, 'node3'):
                dof1_start = self.node_id_to_dof_start[elem.node1.id]
                dof2_start = self.node_id_to_dof_start[elem.node2.id]
                F_global[dof1_start:dof1_start+6] += F_elem[0:6]
                F_global[dof2_start:dof2_start+6] += F_elem[6:12]
            elif hasattr(elem, 'node2') and hasattr(elem, 'node3'):
                dof1_start = self.node_id_to_dof_start[elem.node1.id]
                dof2_start = self.node_id_to_dof_start[elem.node2.id]
                dof3_start = self.node_id_to_dof_start[elem.node3.id]
                F_global[dof1_start:dof1_start+10] += F_elem[0:10]
                F_global[dof2_start:dof2_start+10] += F_elem[10:20]
                F_global[dof3_start:dof3_start+10] += F_elem[20:30]
    
    def _assemble_element_stiffness(self, elem, K_elem: np.ndarray, K_global: np.ndarray):
        """组装单元刚度矩阵"""
        if hasattr(elem, 'node1'):
            if hasattr(elem, 'node2') and not hasattr(elem, 'node3'):
                dof1_start = self.node_id_to_dof_start[elem.node1.id]
                dof2_start = self.node_id_to_dof_start[elem.node2.id]
                for i in range(6):
                    for j in range(6):
                        K_global[dof1_start + i, dof1_start + j] += K_elem[i, j]
                        K_global[dof1_start + i, dof2_start + j] += K_elem[i, 6 + j]
                        K_global[dof2_start + i, dof1_start + j] += K_elem[6 + i, j]
                        K_global[dof2_start + i, dof2_start + j] += K_elem[6 + i, 6 + j]
            elif hasattr(elem, 'node2') and hasattr(elem, 'node3'):
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
        """应用边界条件"""
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
        self.current_time = 0.0
        self.time_history = []
        self.displacement_history = []

