"""
Optimized Creep Solver V2 - 优化的蠕变求解器（性能优化版）
包含以下优化：
1. 稀疏矩阵存储（scipy.sparse）
2. 装配过程优化（避免重复遍历）
3. 自适应时间步进
4. 保证结果精度
"""

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from typing import Dict, List, Optional, Tuple
import time
from solver.solver.newton_raphson import NewtonRaphsonSolver
from solver.mesh.node import Node


class OptimizedCreepSolverV2:
    """
    优化的蠕变求解器
    
    性能优化特性：
    1. 稀疏矩阵存储和求解
    2. 高效的矩阵装配（使用COO格式预分配）
    3. 自适应时间步进
    4. 避免重复计算
    """
    
    def __init__(self, nodes: Dict[int, Node], elements: Dict, 
                 boundary_conditions: Dict, dof_per_node: int = 6,
                 max_iterations: int = 50, tolerance: float = 1e-6,
                 use_sparse: bool = True, use_adaptive_step: bool = True):
        """
        初始化优化的蠕变求解器
        
        参数:
            nodes: 节点字典
            elements: 单元对象字典
            boundary_conditions: 边界条件字典
            dof_per_node: 每个节点的自由度
            max_iterations: 最大迭代次数
            tolerance: 收敛容差
            use_sparse: 是否使用稀疏矩阵
            use_adaptive_step: 是否使用自适应时间步
        """
        self.nodes = nodes
        self.elements = elements
        self.boundary_conditions = boundary_conditions
        self.dof_per_node = dof_per_node
        self.max_iterations = max_iterations
        self.tolerance = tolerance
        self.use_sparse = use_sparse
        self.use_adaptive_step = use_adaptive_step
        
        # 计算总自由度数
        self.dof_count = len(nodes) * dof_per_node
        
        # 创建节点ID到DOF索引的映射
        sorted_node_ids = sorted(nodes.keys())
        self.node_id_to_dof_start = {
            node_id: idx * dof_per_node 
            for idx, node_id in enumerate(sorted_node_ids)
        }
        
        # 性能统计
        self.stats = {
            'assembly_time': 0.0,
            'solve_time': 0.0,
            'total_time': 0.0,
            'num_iterations': 0,
            'num_time_steps': 0
        }
        
        # 预分配稀疏矩阵结构（COO格式）
        if use_sparse:
            self._preallocate_sparse_structure()
    
    def _preallocate_sparse_structure(self):
        """
        预分配稀疏矩阵结构（COO格式）
        避免在每次迭代中重新分配内存
        """
        # 估算非零元素数量
        # 对于管单元，每个单元贡献12x12或30x30的块
        # 估算：每个节点平均连接2-3个单元，每个单元贡献约dof_per_node^2个非零元素
        estimated_nnz = len(self.elements) * (self.dof_per_node * 2) ** 2
        
        # 存储非零元素的行、列索引和数据（COO格式）
        self._coo_rows = []
        self._coo_cols = []
        self._coo_data = []
        self._coo_capacity = estimated_nnz * 2  # 预分配额外空间
        
        # 预分配数组
        self._coo_rows = np.zeros(self._coo_capacity, dtype=np.int32)
        self._coo_cols = np.zeros(self._coo_capacity, dtype=np.int32)
        self._coo_data_array = np.zeros(self._coo_capacity, dtype=np.float64)
        self._coo_index = 0
    
    def _reset_coo_arrays(self):
        """重置COO数组（用于新的装配）"""
        self._coo_index = 0
    
    def _add_to_coo(self, row: int, col: int, value: float):
        """
        添加元素到COO数组
        
        参数:
            row: 行索引
            col: 列索引
            value: 元素值
        """
        if self._coo_index >= len(self._coo_rows):
            # 如果超出容量，扩展数组
            new_size = len(self._coo_rows) * 2
            self._coo_rows = np.resize(self._coo_rows, new_size)
            self._coo_cols = np.resize(self._coo_cols, new_size)
            self._coo_data_array = np.resize(self._coo_data_array, new_size)
        
        self._coo_rows[self._coo_index] = row
        self._coo_cols[self._coo_index] = col
        self._coo_data_array[self._coo_index] = value
        self._coo_index += 1
    
    def assemble_stiffness_matrix(self, tangent: bool = True) -> sp.csr_matrix:
        """
        装配刚度矩阵（稀疏格式）
        
        参数:
            tangent: 是否使用切线刚度矩阵
        
        返回:
            稀疏刚度矩阵（CSR格式）
        """
        start_time = time.time()
        
        # 重置COO数组
        self._reset_coo_arrays()
        
        # 遍历所有单元，装配到COO格式
        for elem_id, elem in self.elements.items():
            # 获取单元节点
            node_ids = elem.node_ids if hasattr(elem, 'node_ids') else []
            if len(node_ids) < 2:
                continue
            
            # 计算单元刚度矩阵
            if tangent and hasattr(elem, 'compute_tangent_stiffness'):
                K_elem = elem.compute_tangent_stiffness()
            elif hasattr(elem, 'compute_stiffness_matrix'):
                K_elem = elem.compute_stiffness_matrix()
            else:
                continue
            
            # 获取节点对应的DOF起始索引
            dof_starts = []
            for node_id in node_ids:
                if node_id in self.node_id_to_dof_start:
                    dof_starts.append(self.node_id_to_dof_start[node_id])
                else:
                    dof_starts.append(-1)
            
            # 组装到全局矩阵（只添加非零元素）
            elem_dof_count = K_elem.shape[0]
            for i in range(elem_dof_count):
                global_i = dof_starts[i // self.dof_per_node] + (i % self.dof_per_node)
                if global_i < 0 or global_i >= self.dof_count:
                    continue
                
                for j in range(elem_dof_count):
                    global_j = dof_starts[j // self.dof_per_node] + (j % self.dof_per_node)
                    if global_j < 0 or global_j >= self.dof_count:
                        continue
                    
                    value = K_elem[i, j]
                    if abs(value) > 1e-15:  # 只存储非零元素
                        self._add_to_coo(global_i, global_j, value)
        
        # 创建CSR格式稀疏矩阵
        if self._coo_index > 0:
            # 截取实际使用的部分
            rows = self._coo_rows[:self._coo_index]
            cols = self._coo_cols[:self._coo_index]
            data = self._coo_data_array[:self._coo_index]
            
            # 转换为CSR格式
            K_sparse = sp.coo_matrix((data, (rows, cols)), 
                                     shape=(self.dof_count, self.dof_count))
            K_csr = K_sparse.tocsr()
        else:
            K_csr = sp.csr_matrix((self.dof_count, self.dof_count))
        
        self.stats['assembly_time'] += time.time() - start_time
        return K_csr
    
    def assemble_load_vector(self, internal_pressure: Dict = None,
                            gravity: np.ndarray = None,
                            temperature: float = None,
                            ref_temperature: float = None,
                            displacements: np.ndarray = None) -> np.ndarray:
        """
        装配载荷向量（优化版）
        
        参数:
            internal_pressure: 内压字典
            gravity: 重力加速度
            temperature: 当前温度
            ref_temperature: 参考温度
            displacements: 当前位移（用于非线性）
        
        返回:
            载荷向量
        """
        F = np.zeros(self.dof_count)
        
        # 遍历所有单元，装配载荷向量
        for elem_id, elem in self.elements.items():
            # 获取单元节点
            node_ids = elem.node_ids if hasattr(elem, 'node_ids') else []
            if len(node_ids) < 2:
                continue
            
            # 计算单元载荷向量
            if hasattr(elem, 'compute_load_vector'):
                F_elem = elem.compute_load_vector(
                    internal_pressure.get(elem_id, 0.0) if internal_pressure else 0.0,
                    gravity, temperature, ref_temperature
                )
            else:
                continue
            
            # 获取节点对应的DOF起始索引
            dof_starts = []
            for node_id in node_ids:
                if node_id in self.node_id_to_dof_start:
                    dof_starts.append(self.node_id_to_dof_start[node_id])
                else:
                    dof_starts.append(-1)
            
            # 组装到全局向量
            for i in range(len(F_elem)):
                node_idx = i // self.dof_per_node
                dof_idx = i % self.dof_per_node
                global_i = dof_starts[node_idx] + dof_idx
                if global_i >= 0 and global_i < self.dof_count:
                    F[global_i] += F_elem[i]
        
        return F
    
    def apply_boundary_conditions(self, K: sp.csr_matrix, F: np.ndarray) -> Tuple[sp.csr_matrix, np.ndarray]:
        """
        应用边界条件（稀疏矩阵版）
        
        参数:
            K: 刚度矩阵（稀疏）
            F: 载荷向量
        
        返回:
            (K_modified, F_modified): 修改后的矩阵和向量
        """
        # 大数法应用边界条件
        large_number = 1e12
        
        # 转换为LIL格式以便修改
        K_lil = K.tolil()
        
        for node_id, bc in self.boundary_conditions.items():
            if node_id not in self.node_id_to_dof_start:
                continue
            
            dof_start = self.node_id_to_dof_start[node_id]
            
            # 只约束标准DOF（前6个或10个）
            num_constrained_dofs = min(6, self.dof_per_node)  # 对于ELBOW290，只约束前6个
            
            for i in range(num_constrained_dofs):
                dof_idx = dof_start + i
                if dof_idx >= self.dof_count:
                    continue
                
                if i in bc.get('displacement', {}):
                    # 固定位移
                    disp_value = bc['displacement'][i]
                    
                    # 修改对角线元素
                    K_lil[dof_idx, dof_idx] = large_number
                    F[dof_idx] = large_number * disp_value
                    
                    # 修改非对角线元素（设为0）
                    K_lil[dof_idx, :] = 0
                    K_lil[dof_idx, dof_idx] = large_number
        
        # 转换回CSR格式
        return K_lil.tocsr(), F
    
    def solve_linear_system(self, K: sp.csr_matrix, F: np.ndarray) -> np.ndarray:
        """
        求解线性方程组（稀疏矩阵版）
        
        参数:
            K: 刚度矩阵（稀疏）
            F: 载荷向量
        
        返回:
            位移向量
        """
        start_time = time.time()
        
        # 使用稀疏矩阵求解器
        if self.use_sparse:
            displacements = spsolve(K, F)
        else:
            # 转换为密集矩阵（用于对比）
            K_dense = K.toarray()
            displacements = np.linalg.solve(K_dense, F)
        
        self.stats['solve_time'] += time.time() - start_time
        return displacements
    
    def compute_adaptive_time_step(self, current_step: int, 
                                   converged_iterations: int,
                                   previous_step_size: float,
                                   min_step: float = 1e-3,
                                   max_step: float = 100.0,
                                   growth_factor: float = 1.5,
                                   reduction_factor: float = 0.5) -> float:
        """
        计算自适应时间步长
        
        参数:
            current_step: 当前时间步
            converged_iterations: 收敛时的迭代次数
            previous_step_size: 前一个时间步的步长
            min_step: 最小时间步长
            max_step: 最大时间步长
            growth_factor: 增长因子
            reduction_factor: 缩减因子
        
        返回:
            新的时间步长
        """
        if not self.use_adaptive_step:
            return previous_step_size
        
        # 根据迭代次数调整时间步
        if converged_iterations <= 3:
            # 迭代次数少，可以增大时间步
            new_step = previous_step_size * growth_factor
        elif converged_iterations >= self.max_iterations * 0.8:
            # 迭代次数多，减小时间步
            new_step = previous_step_size * reduction_factor
        else:
            # 保持不变
            new_step = previous_step_size
        
        # 限制在范围内
        new_step = max(min_step, min(max_step, new_step))
        
        return new_step
    
    def solve_creep_analysis(self, total_time: float, initial_time_step: float = 1.0,
                            internal_pressure: Dict = None,
                            gravity: np.ndarray = None,
                            temperature: float = None,
                            ref_temperature: float = None,
                            initial_displacements: Optional[np.ndarray] = None) -> Dict:
        """
        求解蠕变分析（优化版）
        
        参数:
            total_time: 总时间
            initial_time_step: 初始时间步长
            internal_pressure: 内压字典
            gravity: 重力加速度
            temperature: 当前温度
            ref_temperature: 参考温度
            initial_displacements: 初始位移
        
        返回:
            结果字典，包含位移历史、时间历史等
        """
        total_start_time = time.time()
        
        # 初始化
        if initial_displacements is None:
            displacements = np.zeros(self.dof_count)
        else:
            displacements = initial_displacements.copy()
        
        time_history = [0.0]
        displacement_history = [displacements.copy()]
        
        current_time = 0.0
        current_time_step = initial_time_step
        step_count = 0
        
        # 初始平衡（不考虑蠕变）
        print(f"初始平衡求解...")
        K = self.assemble_stiffness_matrix(tangent=False)
        F = self.assemble_load_vector(internal_pressure, gravity, temperature, ref_temperature)
        K, F = self.apply_boundary_conditions(K, F)
        displacements = self.solve_linear_system(K, F)
        
        # 时间步进循环
        while current_time < total_time:
            step_count += 1
            self.stats['num_time_steps'] += 1
            
            # 确定实际时间步长
            dt = min(current_time_step, total_time - current_time)
            new_time = current_time + dt
            
            print(f"时间步 {step_count}: t = {new_time:.2f}, dt = {dt:.2f}")
            
            # Newton-Raphson迭代求解
            converged = False
            iter_count = 0
            u_current = displacements.copy()
            
            for iteration in range(self.max_iterations):
                iter_count += 1
                self.stats['num_iterations'] += 1
                
                # 装配切线刚度矩阵
                K_t = self.assemble_stiffness_matrix(tangent=True)
                
                # 计算内部力（包括蠕变应变效应）
                F_int = self.compute_internal_force(u_current, temperature, dt, new_time)
                
                # 装配外部载荷
                F_ext = self.assemble_load_vector(internal_pressure, gravity, 
                                                 temperature, ref_temperature, u_current)
                
                # 计算残差
                residual = F_ext - F_int
                
                # 应用边界条件
                K_t, residual = self.apply_boundary_conditions(K_t, residual)
                
                # 求解增量
                delta_u = self.solve_linear_system(K_t, residual)
                
                # 更新位移
                u_current += delta_u
                
                # 检查收敛
                residual_norm = np.linalg.norm(residual)
                if residual_norm < self.tolerance:
                    converged = True
                    break
            
            if not converged:
                print(f"警告：时间步 {step_count} 未收敛")
                # 减小时间步并重试
                current_time_step *= 0.5
                continue
            
            # 更新蠕变应变
            self.update_creep_strain(u_current, temperature, dt, new_time)
            
            # 保存结果
            displacements = u_current.copy()
            time_history.append(new_time)
            displacement_history.append(displacements.copy())
            
            # 更新当前时间
            current_time = new_time
            
            # 自适应调整时间步
            current_time_step = self.compute_adaptive_time_step(
                step_count, iter_count, current_time_step
            )
            
            print(f"  收敛：{iter_count} 次迭代，残差 = {residual_norm:.2e}")
        
        self.stats['total_time'] = time.time() - total_start_time
        
        return {
            'displacements': displacements,
            'time_history': np.array(time_history),
            'displacement_history': np.array(displacement_history),
            'stats': self.stats
        }
    
    def compute_internal_force(self, displacements: np.ndarray,
                              temperature: float, dt: float, time: float) -> np.ndarray:
        """
        计算内部力（包括蠕变效应）
        
        参数:
            displacements: 当前位移
            temperature: 温度
            dt: 时间增量
            time: 当前时间
        
        返回:
            内部力向量
        """
        F_int = np.zeros(self.dof_count)
        
        # 遍历所有单元
        for elem_id, elem in self.elements.items():
            # 获取单元节点
            node_ids = elem.node_ids if hasattr(elem, 'node_ids') else []
            if len(node_ids) < 2:
                continue
            
            # 获取单元位移
            elem_disp = []
            for node_id in node_ids:
                if node_id in self.node_id_to_dof_start:
                    dof_start = self.node_id_to_dof_start[node_id]
                    node_disp = displacements[dof_start:dof_start+self.dof_per_node]
                    elem_disp.extend(node_disp)
            
            if len(elem_disp) == 0:
                continue
            
            # 计算单元内部力
            if hasattr(elem, 'compute_internal_force'):
                F_elem = elem.compute_internal_force(np.array(elem_disp))
            elif hasattr(elem, 'compute_residual'):
                F_elem = elem.compute_residual(np.array(elem_disp))
            else:
                # 简化：使用刚度矩阵
                K_elem = elem.compute_stiffness_matrix()
                F_elem = K_elem @ np.array(elem_disp)
            
            # 组装到全局向量
            dof_starts = []
            for node_id in node_ids:
                if node_id in self.node_id_to_dof_start:
                    dof_starts.append(self.node_id_to_dof_start[node_id])
                else:
                    dof_starts.append(-1)
            
            for i in range(len(F_elem)):
                node_idx = i // self.dof_per_node
                dof_idx = i % self.dof_per_node
                global_i = dof_starts[node_idx] + dof_idx
                if global_i >= 0 and global_i < self.dof_count:
                    F_int[global_i] += F_elem[i]
        
        return F_int
    
    def update_creep_strain(self, displacements: np.ndarray,
                           temperature: float, dt: float, time: float):
        """
        更新蠕变应变
        
        参数:
            displacements: 当前位移
            temperature: 温度
            dt: 时间增量
            time: 当前时间
        """
        # 遍历所有单元，更新蠕变应变
        for elem_id, elem in self.elements.items():
            if hasattr(elem, 'update_creep_strain'):
                elem.update_creep_strain(
                    self.node_id_to_dof_start, displacements,
                    time, dt, temperature
                )

