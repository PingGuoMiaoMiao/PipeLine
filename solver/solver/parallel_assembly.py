"""
Parallel Assembly - 并行装配优化
使用多进程加速单元计算和装配过程
"""

import numpy as np
from multiprocessing import Pool, cpu_count
from functools import partial
from typing import Dict, List, Callable, Tuple
import time


class ParallelAssembly:
    """并行装配器"""
    
    def __init__(self, num_processes: int = None):
        """
        初始化并行装配器
        
        参数:
            num_processes: 进程数（None表示使用CPU核心数）
        """
        if num_processes is None:
            num_processes = max(1, cpu_count() - 1)  # 保留一个核心
        self.num_processes = num_processes
    
    def assemble_stiffness_parallel(self,
                                   elements: Dict,
                                   compute_func: Callable,
                                   node_dof_map: Dict[int, int],
                                   dof_count: int) -> np.ndarray:
        """
        并行装配刚度矩阵
        
        参数:
            elements: 单元字典
            compute_func: 计算单元刚度矩阵的函数
            node_dof_map: 节点ID到DOF索引的映射
            dof_count: 总DOF数
        
        返回:
            全局刚度矩阵（dense或sparse）
        """
        if self.num_processes <= 1 or len(elements) < 10:
            # 单元数少或单进程，使用串行
            return self._assemble_stiffness_serial(
                elements, compute_func, node_dof_map, dof_count)
        
        # 将单元分成多个任务
        element_items = list(elements.items())
        chunk_size = max(1, len(element_items) // self.num_processes)
        chunks = [element_items[i:i+chunk_size] 
                 for i in range(0, len(element_items), chunk_size)]
        
        # 并行计算
        with Pool(processes=self.num_processes) as pool:
            results = pool.starmap(
                self._compute_chunk_stiffness,
                [(chunk, compute_func, node_dof_map) for chunk in chunks]
            )
        
        # 合并结果
        K_global = np.zeros((dof_count, dof_count))
        for K_chunk in results:
            K_global += K_chunk
        
        return K_global
    
    def _compute_chunk_stiffness(self,
                                element_chunk: List,
                                compute_func: Callable,
                                node_dof_map: Dict[int, int]) -> np.ndarray:
        """计算一个单元块的刚度矩阵"""
        # 计算最大DOF数
        max_dof = max(node_dof_map.values()) + 6  # 假设每节点最多6 DOF
        K_chunk = np.zeros((max_dof, max_dof))
        
        for elem_id, elem in element_chunk:
            K_elem = compute_func(elem)
            
            if hasattr(elem, 'node1'):
                if hasattr(elem, 'node2') and not hasattr(elem, 'node3'):
                    dof1_start = node_dof_map.get(elem.node1.id, 0)
                    dof2_start = node_dof_map.get(elem.node2.id, 0)
                    
                    # 组装（简化版本）
                    K_chunk[dof1_start:dof1_start+6, dof1_start:dof1_start+6] += K_elem[0:6, 0:6]
                    K_chunk[dof1_start:dof1_start+6, dof2_start:dof2_start+6] += K_elem[0:6, 6:12]
                    K_chunk[dof2_start:dof2_start+6, dof1_start:dof1_start+6] += K_elem[6:12, 0:6]
                    K_chunk[dof2_start:dof2_start+6, dof2_start:dof2_start+6] += K_elem[6:12, 6:12]
        
        return K_chunk
    
    def _assemble_stiffness_serial(self,
                                  elements: Dict,
                                  compute_func: Callable,
                                  node_dof_map: Dict[int, int],
                                  dof_count: int) -> np.ndarray:
        """串行装配刚度矩阵"""
        K_global = np.zeros((dof_count, dof_count))
        
        for elem in elements.values():
            K_elem = compute_func(elem)
            
            if hasattr(elem, 'node1'):
                if hasattr(elem, 'node2') and not hasattr(elem, 'node3'):
                    dof1_start = node_dof_map[elem.node1.id]
                    dof2_start = node_dof_map[elem.node2.id]
                    
                    K_global[dof1_start:dof1_start+6, dof1_start:dof1_start+6] += K_elem[0:6, 0:6]
                    K_global[dof1_start:dof1_start+6, dof2_start:dof2_start+6] += K_elem[0:6, 6:12]
                    K_global[dof2_start:dof2_start+6, dof1_start:dof1_start+6] += K_elem[6:12, 0:6]
                    K_global[dof2_start:dof2_start+6, dof2_start:dof2_start+6] += K_elem[6:12, 6:12]
        
        return K_global
    
    def assemble_force_parallel(self,
                               elements: Dict,
                               compute_func: Callable,
                               node_dof_map: Dict[int, int],
                               dof_count: int) -> np.ndarray:
        """
        并行装配载荷向量
        
        参数:
            elements: 单元字典
            compute_func: 计算单元载荷向量的函数
            node_dof_map: 节点ID到DOF索引的映射
            dof_count: 总DOF数
        
        返回:
            全局载荷向量
        """
        if self.num_processes <= 1 or len(elements) < 10:
            return self._assemble_force_serial(
                elements, compute_func, node_dof_map, dof_count)
        
        # 将单元分成多个任务
        element_items = list(elements.items())
        chunk_size = max(1, len(element_items) // self.num_processes)
        chunks = [element_items[i:i+chunk_size] 
                 for i in range(0, len(element_items), chunk_size)]
        
        # 并行计算
        with Pool(processes=self.num_processes) as pool:
            results = pool.starmap(
                self._compute_chunk_force,
                [(chunk, compute_func, node_dof_map) for chunk in chunks]
            )
        
        # 合并结果
        F_global = np.zeros(dof_count)
        for F_chunk in results:
            F_global += F_chunk
        
        return F_global
    
    def _compute_chunk_force(self,
                            element_chunk: List,
                            compute_func: Callable,
                            node_dof_map: Dict[int, int]) -> np.ndarray:
        """计算一个单元块的载荷向量"""
        max_dof = max(node_dof_map.values()) + 6
        F_chunk = np.zeros(max_dof)
        
        for elem_id, elem in element_chunk:
            F_elem = compute_func(elem)
            
            if hasattr(elem, 'node1'):
                if hasattr(elem, 'node2') and not hasattr(elem, 'node3'):
                    dof1_start = node_dof_map.get(elem.node1.id, 0)
                    dof2_start = node_dof_map.get(elem.node2.id, 0)
                    
                    F_chunk[dof1_start:dof1_start+6] += F_elem[0:6]
                    F_chunk[dof2_start:dof2_start+6] += F_elem[6:12]
        
        return F_chunk
    
    def _assemble_force_serial(self,
                              elements: Dict,
                              compute_func: Callable,
                              node_dof_map: Dict[int, int],
                              dof_count: int) -> np.ndarray:
        """串行装配载荷向量"""
        F_global = np.zeros(dof_count)
        
        for elem in elements.values():
            F_elem = compute_func(elem)
            
            if hasattr(elem, 'node1'):
                if hasattr(elem, 'node2') and not hasattr(elem, 'node3'):
                    dof1_start = node_dof_map[elem.node1.id]
                    dof2_start = node_dof_map[elem.node2.id]
                    
                    F_global[dof1_start:dof1_start+6] += F_elem[0:6]
                    F_global[dof2_start:dof2_start+6] += F_elem[6:12]
        
        return F_global

