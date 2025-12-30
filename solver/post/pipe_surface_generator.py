"""
Pipe Surface Generator - 管单元表面生成器
为PIPE288直管生成圆柱体表面（用于VTK输出）
"""

import math
import numpy as np
from typing import List, Tuple, Dict
from solver.mesh.node import Node


class PipeSurfaceGenerator:
    """管单元表面生成器"""
    
    @staticmethod
    def generate_pipe288_surface_points(node1: Node, node2: Node,
                                       Ro: float, Ri: float,
                                       displacement1: np.ndarray = None,
                                       displacement2: np.ndarray = None,
                                       num_circumferential: int = 20) -> Tuple[List, List]:
        """
        为PIPE288单元生成圆柱体表面点（外表面和内表面）
        
        参数:
            node1: 节点1
            node2: 节点2
            Ro: 外半径 (m)
            Ri: 内半径 (m)
            displacement1: 节点1的位移 [ux, uy, uz] (m)
            displacement2: 节点2的位移 [ux, uy, uz] (m)
            num_circumferential: 环向划分点数
        
        返回:
            (points, cells): 点列表和单元列表（四边形面片）
        """
        if displacement1 is None:
            displacement1 = np.zeros(3)
        if displacement2 is None:
            displacement2 = np.zeros(3)
        
        # 节点位置（mm转换为m）
        p1 = np.array([node1.x, node1.y, node1.z]) / 1000.0 + displacement1
        p2 = np.array([node2.x, node2.y, node2.z]) / 1000.0 + displacement2
        
        # 单元方向向量
        dx = p2 - p1
        L = np.linalg.norm(dx)
        if L < 1e-10:
            return [], []
        
        ex = dx / L  # 轴向单位向量
        
        # 构造局部坐标系（与PIPE288单元一致）
        if abs(ex[2]) < 0.9:
            ez_global = np.array([0, 0, 1])
            ey = np.cross(ex, ez_global)
        else:
            ey_global = np.array([0, 1, 0])
            ey = np.cross(ex, ey_global)
        
        ey_norm = np.linalg.norm(ey)
        if ey_norm > 1e-10:
            ey = ey / ey_norm
        else:
            ey = np.array([0, 1, 0])
        
        ez = np.cross(ex, ey)
        ez_norm = np.linalg.norm(ez)
        if ez_norm > 1e-10:
            ez = ez / ez_norm
        else:
            ez = np.array([0, 0, 1])
        ey = np.cross(ez, ex)
        
        # 生成环向角度
        phis = np.linspace(0, 2 * math.pi, num_circumferential, endpoint=False)
        
        points = []
        cells = []
        point_idx = 0
        
        # 生成外表面和内表面的点
        for radius in [Ro, Ri]:
            for phi in phis:
                # 环向方向
                r_local = radius * (np.cos(phi) * ey + np.sin(phi) * ez)
                
                # 两个节点处的截面点
                p1_section = p1 + r_local
                p2_section = p2 + r_local
                
                # 转换回mm并添加
                points.append((p1_section[0] * 1000.0, p1_section[1] * 1000.0, p1_section[2] * 1000.0))
                points.append((p2_section[0] * 1000.0, p2_section[1] * 1000.0, p2_section[2] * 1000.0))
        
        # 生成四边形面片
        # 外表面
        for i in range(num_circumferential):
            i_next = (i + 1) % num_circumferential
            
            # 节点1处的两个点
            idx1_outer = i * 2
            idx1_outer_next = i_next * 2
            
            # 节点2处的两个点
            idx2_outer = i * 2 + 1
            idx2_outer_next = i_next * 2 + 1
            
            # 形成四边形（外表面）
            cells.append([4, idx1_outer, idx2_outer, idx2_outer_next, idx1_outer_next])
        
        # 内表面
        inner_start = num_circumferential * 2  # 外表面有 num_circumferential * 2 个点
        for i in range(num_circumferential):
            i_next = (i + 1) % num_circumferential
            
            # 节点1处的两个点
            idx1_inner = inner_start + i * 2
            idx1_inner_next = inner_start + i_next * 2
            
            # 节点2处的两个点
            idx2_inner = inner_start + i * 2 + 1
            idx2_inner_next = inner_start + i_next * 2 + 1
            
            # 形成四边形（内表面，注意方向）
            cells.append([4, idx1_inner, idx1_inner_next, idx2_inner_next, idx2_inner])
        
        # 生成端面（节点1和节点2处的圆环）
        # 节点1端面
        for i in range(num_circumferential):
            i_next = (i + 1) % num_circumferential
            idx_outer1 = i * 2
            idx_outer1_next = i_next * 2
            idx_inner1 = inner_start + i * 2
            idx_inner1_next = inner_start + i_next * 2
            # 外圆到内圆的四边形
            cells.append([4, idx_outer1, idx_outer1_next, idx_inner1_next, idx_inner1])
        
        # 节点2端面
        for i in range(num_circumferential):
            i_next = (i + 1) % num_circumferential
            idx_outer2 = i * 2 + 1
            idx_outer2_next = i_next * 2 + 1
            idx_inner2 = inner_start + i * 2 + 1
            idx_inner2_next = inner_start + i_next * 2 + 1
            # 外圆到内圆的四边形（注意方向）
            cells.append([4, idx_outer2, idx_inner2, idx_inner2_next, idx_outer2_next])
        
        return points, cells
    
    @staticmethod
    def generate_pipe288_surfaces_for_elements(pipe288_elements: Dict,
                                              nodes: Dict[int, Node],
                                              displacements: np.ndarray,
                                              node_id_to_dof_start: Dict[int, int],
                                              num_circumferential: int = 20) -> Tuple[List, List]:
        """
        为多个PIPE288单元生成表面
        
        参数:
            pipe288_elements: PIPE288单元对象字典
            nodes: 节点字典
            displacements: 位移数组
            node_id_to_dof_start: 节点ID到DOF起始索引的映射
            num_circumferential: 环向划分点数
        
        返回:
            (points, cells): 所有点的列表和所有单元的列表
        """
        all_points = []
        all_cells = []
        current_point_idx = 0
        
        for elem_id, elem in pipe288_elements.items():
            # 获取节点位移
            dof_start1 = node_id_to_dof_start.get(elem.node1.id, -1)
            dof_start2 = node_id_to_dof_start.get(elem.node2.id, -1)
            
            if dof_start1 < 0 or dof_start2 < 0:
                continue
            
            disp1 = np.zeros(3)
            disp2 = np.zeros(3)
            if dof_start1 + 3 <= len(displacements):
                disp1 = displacements[dof_start1:dof_start1+3]
            if dof_start2 + 3 <= len(displacements):
                disp2 = displacements[dof_start2:dof_start2+3]
            
            # 生成该单元的表面
            points, cells = PipeSurfaceGenerator.generate_pipe288_surface_points(
                elem.node1, elem.node2,
                elem.Ro, elem.Ri,
                disp1, disp2,
                num_circumferential
            )
            
            if len(points) > 0:
                # 添加点
                all_points.extend(points)
                
                # 添加单元（调整索引）
                for cell in cells:
                    adjusted_cell = [cell[0]] + [current_point_idx + idx for idx in cell[1:]]
                    all_cells.append(adjusted_cell)
                
                current_point_idx += len(points)
        
        return all_points, all_cells

