"""
Optimized Creep Solver - 优化的蠕变求解器
整合稀疏矩阵、并行计算和自适应时间步
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from typing import Dict, List, Optional, Tuple
from solver.mesh.node import Node
from solver.solver.sparse_solver import SparseSolver
from solver.solver.adaptive_time_step import AdaptiveTimeStepController, ErrorEstimator
from solver.solver.parallel_assembly import ParallelAssembly


class OptimizedCreepSolver:
    """优化的蠕变时间积分求解器"""
    
    def __init__(self,
                 nodes: Dict[int, Node],
                 elements: Dict,
                 boundary_conditions: Dict[int, List[int]],
                 dof_per_node: int = 6,
                 max_iterations: int = 20,
                 tolerance: float = 1e-6,
                 use_sparse: bool = True,
                 use_parallel: bool = True,
                 num_processes: int = None,
                 adaptive_time_step: bool = True,
                 initial_time_step: float = 100.0):
        """
        初始化优化的蠕变求解器
        
        参数:
            nodes: 节点字典
            elements: 单元字典
            boundary_conditions: 边界条件字典
            dof_per_node: 每个节点的DOF数
            max_iterations: 最大迭代次数
            tolerance: 收敛容差
            use_sparse: 是否使用稀疏矩阵
            use_parallel: 是否使用并行计算
            num_processes: 并行进程数
            adaptive_time_step: 是否使用自适应时间步
            initial_time_step: 初始时间步长
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
        
        # 优化选项
        self.use_sparse = use_sparse
        self.use_parallel = use_parallel and len(elements) > 10
        self.adaptive_time_step = adaptive_time_step
        
        # 初始化优化组件
        if self.use_sparse:
            self.sparse_solver = SparseSolver(
                nodes, elements, boundary_conditions, dof_per_node)
        
        if self.use_parallel:
            self.parallel_assembly = ParallelAssembly(num_processes)
        
        if self.adaptive_time_step:
            self.time_step_controller = AdaptiveTimeStepController(
                initial_time_step=initial_time_step,
                min_time_step=1.0,
                max_time_step=1000.0,
                safety_factor=0.8,
                max_increase_factor=1.5,
                min_decrease_factor=0.5,
                target_iterations=5
            )
        
        # 当前状态
        self.displacements = np.zeros(self.dof_count)
        self.current_time = 0.0
        self.time_history = []
        self.displacement_history = []
        
        # 性能统计
        self.performance_stats = {
            'assembly_time': [],
            'solve_time': [],
            'total_time': [],
            'iterations': [],
            'time_steps': []
        }
    
    def solve_time_step(self,
                       load_vector: np.ndarray,
                       time: float,
                       time_increment: float,
                       temperature: float,
                       compute_stiffness_func,
                       compute_force_func,
                       update_creep_func) -> Tuple[np.ndarray, Dict]:
        """
        求解一个时间步（优化版）
        
        参数:
            load_vector: 载荷向量
            time: 当前时间
            time_increment: 时间增量
            temperature: 温度（绝对温度，K）
            compute_stiffness_func: 计算单元刚度矩阵的函数
            compute_force_func: 计算单元载荷向量的函数
            update_creep_func: 更新蠕变应变的函数
        
        返回:
            (displacements, info): 位移和求解信息
        """
        import time as time_module
        
        step_start_time = time_module.time()
        
        info = {
            'converged': False,
            'iterations': 0,
            'residual_norms': [],
            'time': time,
            'time_increment': time_increment,
            'assembly_time': 0.0,
            'solve_time': 0.0
        }
        
        # 保存积分点状态
        self._save_integration_points_state()
        
        U = self.displacements.copy()
        
        # Newton-Raphson迭代
        for iteration in range(self.max_iterations):
            # 更新蠕变应变
            update_creep_func(U, time, time_increment, temperature)
            
            # 装配刚度矩阵（优化）
            assembly_start = time_module.time()
            if self.use_sparse:
                K_tangent = self._assemble_sparse_stiffness(U, compute_stiffness_func)
            else:
                K_tangent = self._assemble_dense_stiffness(U, compute_stiffness_func)
            info['assembly_time'] += time_module.time() - assembly_start
            
            # 计算内部力
            if self.use_parallel:
                F_int = self.parallel_assembly.assemble_force_parallel(
                    self.elements, compute_force_func,
                    self.node_id_to_dof_start, self.dof_count)
            else:
                F_int = self._assemble_internal_force(U, compute_force_func)
            
            # 计算残差
            R = load_vector - F_int
            residual_norm = np.linalg.norm(R)
            info['residual_norms'].append(residual_norm)
            
            # 检查收敛
            if residual_norm < self.tolerance:
                info['converged'] = True
                info['iterations'] = iteration + 1
                break
            
            # 应用边界条件
            if self.use_sparse:
                K_tangent, R = self.sparse_solver.apply_boundary_conditions_sparse(
                    K_tangent, R)
            else:
                self._apply_boundary_conditions(K_tangent, R)
            
            # 求解线性方程组
            solve_start = time_module.time()
            if self.use_sparse:
                delta_U = self.sparse_solver.solve_sparse(K_tangent, R)
            else:
                delta_U = np.linalg.solve(K_tangent, R)
            info['solve_time'] += time_module.time() - solve_start
            
            # 更新位移
            U += delta_U
        
        if not info['converged']:
            print(f"警告: 未收敛（残差 = {info['residual_norms'][-1]:.6e}）")
        
        # 更新状态
        self.displacements = U.copy()
        self.current_time = time + time_increment
        self.time_history.append(self.current_time)
        self.displacement_history.append(U.copy())
        
        # 自适应时间步调整
        if self.adaptive_time_step:
            new_time_step = self.time_step_controller.compute_next_time_step(
                info['iterations'], info['converged'],
                info['residual_norms'][-1] if info['residual_norms'] else 1e10,
                self.tolerance)
            info['next_time_step'] = new_time_step
        
        # 性能统计
        total_time = time_module.time() - step_start_time
        info['total_time'] = total_time
        self.performance_stats['assembly_time'].append(info['assembly_time'])
        self.performance_stats['solve_time'].append(info['solve_time'])
        self.performance_stats['total_time'].append(total_time)
        self.performance_stats['iterations'].append(info['iterations'])
        
        return U, info
    
    def _assemble_sparse_stiffness(self, displacements, compute_func) -> sparse.csr_matrix:
        """装配稀疏刚度矩阵"""
        def compute_elem_stiffness(elem):
            if hasattr(elem, 'compute_tangent_stiffness'):
                node_disp = self._get_element_displacements(elem, displacements)
                return elem.compute_tangent_stiffness(node_disp)
            else:
                return elem.compute_stiffness_matrix()
        
        return self.sparse_solver.assemble_sparse_stiffness_matrix(compute_elem_stiffness)
    
    def _assemble_dense_stiffness(self, displacements, compute_func) -> np.ndarray:
        """装配dense刚度矩阵"""
        K_tangent = np.zeros((self.dof_count, self.dof_count))
        
        for elem in self.elements.values():
            if hasattr(elem, 'compute_tangent_stiffness'):
                node_disp = self._get_element_displacements(elem, displacements)
                K_elem = elem.compute_tangent_stiffness(node_disp)
            else:
                K_elem = elem.compute_stiffness_matrix()
            
            self._assemble_element_stiffness(elem, K_elem, K_tangent)
        
        return K_tangent
    
    def _assemble_internal_force(self, displacements, compute_func) -> np.ndarray:
        """装配内部力向量"""
        F_int = np.zeros(self.dof_count)
        
        for elem in self.elements.values():
            if hasattr(elem, 'compute_internal_force'):
                node_disp = self._get_element_displacements(elem, displacements)
                F_elem = elem.compute_internal_force(node_disp)
            else:
                F_elem = np.zeros(12)  # 默认
            
            self._assemble_element_force(elem, F_elem, F_int)
        
        return F_int
    
    def _get_element_displacements(self, elem, global_displacements: np.ndarray) -> np.ndarray:
        """获取单元节点的位移"""
        if hasattr(elem, 'node1'):
            if hasattr(elem, 'node2') and not hasattr(elem, 'node3'):
                dof1_start = self.node_id_to_dof_start[elem.node1.id]
                dof2_start = self.node_id_to_dof_start[elem.node2.id]
                node_disp = np.zeros(12)
                node_disp[0:6] = global_displacements[dof1_start:dof1_start+6]
                node_disp[6:12] = global_displacements[dof2_start:dof2_start+6]
                return node_disp
        return np.zeros(12)
    
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
    
    def _assemble_element_force(self, elem, F_elem: np.ndarray, F_global: np.ndarray):
        """组装单元力向量"""
        if hasattr(elem, 'node1'):
            if hasattr(elem, 'node2') and not hasattr(elem, 'node3'):
                dof1_start = self.node_id_to_dof_start[elem.node1.id]
                dof2_start = self.node_id_to_dof_start[elem.node2.id]
                F_global[dof1_start:dof1_start+6] += F_elem[0:6]
                F_global[dof2_start:dof2_start+6] += F_elem[6:12]
    
    def _apply_boundary_conditions(self, K: np.ndarray, F: np.ndarray, penalty: float = 1e12):
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
    
    def _save_integration_points_state(self):
        """保存积分点状态"""
        for elem in self.elements.values():
            if hasattr(elem, 'integration_points'):
                for ip in elem.integration_points:
                    ip.save_state()
    
    def get_performance_stats(self) -> Dict:
        """获取性能统计"""
        stats = self.performance_stats.copy()
        if stats['total_time']:
            stats['avg_assembly_time'] = np.mean(stats['assembly_time'])
            stats['avg_solve_time'] = np.mean(stats['solve_time'])
            stats['avg_total_time'] = np.mean(stats['total_time'])
            stats['total_time'] = sum(stats['total_time'])
            stats['avg_iterations'] = np.mean(stats['iterations'])
        return stats
    
    def reset(self):
        """重置求解器状态"""
        self.displacements = np.zeros(self.dof_count)
        self.current_time = 0.0
        self.time_history = []
        self.displacement_history = []
        self.performance_stats = {
            'assembly_time': [],
            'solve_time': [],
            'total_time': [],
            'iterations': [],
            'time_steps': []
        }
        if self.adaptive_time_step:
            self.time_step_controller.reset()

