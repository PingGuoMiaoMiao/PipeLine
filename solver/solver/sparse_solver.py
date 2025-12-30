"""
Sparse Matrix Solver - 稀疏矩阵求解器优化
使用scipy.sparse提高大规模问题的计算效率
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from typing import Dict, List, Tuple, Optional
from solver.mesh.node import Node


class SparseSolver:
    """稀疏矩阵求解器（优化版）"""
    
    def __init__(self, nodes: Dict[int, Node],
                 elements: Dict,
                 boundary_conditions: Dict[int, List[int]],
                 dof_per_node: int = 6):
        """
        初始化稀疏矩阵求解器
        
        参数:
            nodes: 节点字典
            elements: 单元字典
            boundary_conditions: 边界条件字典
            dof_per_node: 每个节点的DOF数
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
        
        # 预分配稀疏矩阵结构（CSR格式）
        self._preallocate_sparsity_pattern()
    
    def _preallocate_sparsity_pattern(self):
        """预分配稀疏矩阵的稀疏结构"""
        # 估算非零元素数量
        # 每个单元最多贡献 (dof_per_node * num_nodes_per_elem)^2 个非零元素
        max_nonzeros = 0
        
        for elem in self.elements.values():
            if hasattr(elem, 'node1'):
                if hasattr(elem, 'node2') and not hasattr(elem, 'node3'):
                    # PIPE288: 2节点，每节点6 DOF = 12 DOF
                    max_nonzeros += 12 * 12
                elif hasattr(elem, 'node2') and hasattr(elem, 'node3'):
                    # ELBOW290: 3节点，每节点10 DOF = 30 DOF
                    max_nonzeros += 30 * 30
        
        # 使用估算值（通常实际非零元素会更少）
        self.estimated_nonzeros = max_nonzeros
    
    def assemble_sparse_stiffness_matrix(self, compute_func) -> sparse.csr_matrix:
        """
        组装稀疏刚度矩阵
        
        参数:
            compute_func: 计算单元刚度矩阵的函数
        
        返回:
            CSR格式的稀疏矩阵
        """
        # 使用COO格式收集数据（更高效）
        row_indices = []
        col_indices = []
        values = []
        
        for elem in self.elements.values():
            K_elem = compute_func(elem)
            
            # 获取单元节点的DOF索引
            if hasattr(elem, 'node1'):
                if hasattr(elem, 'node2') and not hasattr(elem, 'node3'):
                    # PIPE288
                    dof1_start = self.node_id_to_dof_start[elem.node1.id]
                    dof2_start = self.node_id_to_dof_start[elem.node2.id]
                    dof_indices = [dof1_start + i for i in range(6)] + \
                                 [dof2_start + i for i in range(6)]
                    
                    # 添加非零元素
                    for i, dof_i in enumerate(dof_indices):
                        for j, dof_j in enumerate(dof_indices):
                            if abs(K_elem[i, j]) > 1e-15:  # 过滤很小的值
                                row_indices.append(dof_i)
                                col_indices.append(dof_j)
                                values.append(K_elem[i, j])
                
                elif hasattr(elem, 'node2') and hasattr(elem, 'node3'):
                    # ELBOW290
                    dof1_start = self.node_id_to_dof_start[elem.node1.id]
                    dof2_start = self.node_id_to_dof_start[elem.node2.id]
                    dof3_start = self.node_id_to_dof_start[elem.node3.id]
                    dof_indices = ([dof1_start + i for i in range(10)] +
                                  [dof2_start + i for i in range(10)] +
                                  [dof3_start + i for i in range(10)])
                    
                    # 添加非零元素
                    for i, dof_i in enumerate(dof_indices):
                        for j, dof_j in enumerate(dof_indices):
                            if abs(K_elem[i, j]) > 1e-15:
                                row_indices.append(dof_i)
                                col_indices.append(dof_j)
                                values.append(K_elem[i, j])
        
        # 转换为COO格式，然后转为CSR格式（CSR格式更适合求解）
        coo_matrix = sparse.coo_matrix(
            (values, (row_indices, col_indices)),
            shape=(self.dof_count, self.dof_count)
        )
        
        # 转换为CSR格式
        csr_matrix = coo_matrix.tocsr()
        
        return csr_matrix
    
    def solve_sparse(self, K: sparse.csr_matrix, F: np.ndarray) -> np.ndarray:
        """
        求解稀疏线性方程组
        
        参数:
            K: 稀疏刚度矩阵（CSR格式）
            F: 载荷向量
        
        返回:
            位移向量
        """
        try:
            U = spsolve(K, F)
        except Exception as e:
            print(f"警告: 稀疏矩阵求解失败: {e}")
            # 回退到dense方法
            U = np.linalg.solve(K.toarray(), F)
        
        return U
    
    def apply_boundary_conditions_sparse(self, K: sparse.csr_matrix,
                                        F: np.ndarray,
                                        penalty: float = 1e12) -> Tuple[sparse.csr_matrix, np.ndarray]:
        """
        应用边界条件（稀疏矩阵版本）
        
        参数:
            K: 稀疏刚度矩阵
            F: 载荷向量
            penalty: 惩罚系数
        
        返回:
            (修改后的K, 修改后的F)
        """
        # 转换为LIL格式（更适合修改）
        K = K.tolil()
        
        for node_id, dof_list in self.boundary_conditions.items():
            if node_id in self.node_id_to_dof_start:
                dof_start = self.node_id_to_dof_start[node_id]
                for dof in dof_list:
                    if dof < 6 and dof_start + dof < self.dof_count:
                        dof_idx = dof_start + dof
                        
                        # 清零行和列
                        K[dof_idx, :] = 0
                        K[:, dof_idx] = 0
                        
                        # 设置对角元素
                        K[dof_idx, dof_idx] = penalty
                        
                        # 清零载荷
                        F[dof_idx] = 0
        
        # 转回CSR格式
        K = K.tocsr()
        
        return K, F


class OptimizedAssembly:
    """优化的装配过程"""
    
    @staticmethod
    def assemble_stiffness_vectorized(elements: Dict,
                                     node_dof_map: Dict[int, int],
                                     compute_func,
                                     dof_count: int) -> sparse.csr_matrix:
        """
        向量化的刚度矩阵装配（批量处理）
        
        参数:
            elements: 单元字典
            node_dof_map: 节点ID到DOF索引的映射
            compute_func: 计算单元刚度矩阵的函数
            dof_count: 总DOF数
        
        返回:
            稀疏刚度矩阵
        """
        # 使用列表收集数据（避免重复分配）
        row_indices = []
        col_indices = []
        values = []
        
        # 批量处理单元
        for elem in elements.values():
            K_elem = compute_func(elem)
            
            if hasattr(elem, 'node1'):
                if hasattr(elem, 'node2') and not hasattr(elem, 'node3'):
                    # PIPE288
                    dof1_start = node_dof_map[elem.node1.id]
                    dof2_start = node_dof_map[elem.node2.id]
                    
                    # 使用numpy的meshgrid加速索引计算
                    dof_indices = np.array([dof1_start + i for i in range(6)] +
                                          [dof2_start + i for i in range(6)])
                    
                    # 向量化添加非零元素
                    i_indices, j_indices = np.meshgrid(dof_indices, dof_indices, indexing='ij')
                    mask = np.abs(K_elem) > 1e-15
                    
                    row_indices.extend(i_indices[mask].tolist())
                    col_indices.extend(j_indices[mask].tolist())
                    values.extend(K_elem[mask].tolist())
        
        # 转换为稀疏矩阵
        coo_matrix = sparse.coo_matrix(
            (values, (row_indices, col_indices)),
            shape=(dof_count, dof_count)
        )
        
        return coo_matrix.tocsr()
    
    @staticmethod
    def assemble_force_vectorized(elements: Dict,
                                 node_dof_map: Dict[int, int],
                                 compute_func,
                                 dof_count: int) -> np.ndarray:
        """
        向量化的载荷向量装配
        
        参数:
            elements: 单元字典
            node_dof_map: 节点ID到DOF索引的映射
            compute_func: 计算单元载荷向量的函数
            dof_count: 总DOF数
        
        返回:
            载荷向量
        """
        F = np.zeros(dof_count)
        
        for elem in elements.values():
            F_elem = compute_func(elem)
            
            if hasattr(elem, 'node1'):
                if hasattr(elem, 'node2') and not hasattr(elem, 'node3'):
                    # PIPE288
                    dof1_start = node_dof_map[elem.node1.id]
                    dof2_start = node_dof_map[elem.node2.id]
                    
                    # 直接赋值（numpy的+=操作已经优化）
                    F[dof1_start:dof1_start+6] += F_elem[0:6]
                    F[dof2_start:dof2_start+6] += F_elem[6:12]
        
        return F

