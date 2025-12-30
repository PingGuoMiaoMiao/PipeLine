"""
Pipe Surface Expand - 管单元曲面展开
将环向自由度扩展为三维曲面，用于VTK输出
"""

import numpy as np
from typing import List, Tuple, Dict
from solver.mesh.node import Node
from solver.element.shape_function_fourier import FourierShapeFunction


class PipeSurfaceExpander:
    """管单元曲面展开器"""
    
    @staticmethod
    def expand_elbow290_to_surface(elbow290_elements: Dict,
                                   nodes: Dict[int, Node],
                                   displacements: np.ndarray,
                                   node_dof_map: Dict[int, int],
                                   num_circumferential: int = 20) -> Tuple[List, List]:
        """
        将ELBOW290单元扩展为三维曲面（四边形面片）
        
        参数:
            elbow290_elements: ELBOW290单元对象字典
            nodes: 节点字典
            displacements: 位移数组
            node_dof_map: 节点ID到DOF起始索引的映射（ELBOW290每个节点10个DOF）
            num_circumferential: 环向划分点数
        
        返回:
            (points, cells) 点列表和单元列表
        """
        points = []
        cells = []
        
        # 遍历所有ELBOW290单元
        for elem_id, elem in elbow290_elements.items():
            # 获取节点的自由度
            node_dofs_list = []
            for node_id in [elem.node1.id, elem.node2.id, elem.node3.id]:
                if node_id in node_dof_map:
                    dof_start = node_dof_map[node_id]
                    # ELBOW290每个节点10个自由度
                    if dof_start + 10 <= len(displacements):
                        node_dofs = displacements[dof_start:dof_start + 10]
                        node_dofs_list.append(node_dofs)
                    else:
                        node_dofs_list.append(np.zeros(10))
                else:
                    node_dofs_list.append(np.zeros(10))
            
            # 为每个节点生成截面点
            node_sections = []
            node_positions = [
                np.array([elem.node1.x, elem.node1.y, elem.node1.z]) / 1000.0,
                np.array([elem.node2.x, elem.node2.y, elem.node2.z]) / 1000.0,
                np.array([elem.node3.x, elem.node3.y, elem.node3.z]) / 1000.0
            ]
            
            for node_idx in range(3):
                node_dofs = node_dofs_list[node_idx]
                node_pos = node_positions[node_idx]
                
                # 获取截面点
                section_points = elem.get_section_points(node_pos, node_dofs, num_circumferential)
                if len(section_points) == num_circumferential:
                    node_sections.append(section_points)
                else:
                    # 如果生成失败，跳过这个单元
                    node_sections = []
                    break
            
            if len(node_sections) < 2:
                continue
            
            # 连接相邻节点形成四边形面片
            for i in range(len(node_sections) - 1):
                section1 = node_sections[i]
                section2 = node_sections[i + 1]
                
                for j in range(num_circumferential):
                    j_next = (j + 1) % num_circumferential
                    
                    # 四个点形成四边形
                    p1 = section1[j]
                    p2 = section1[j_next]
                    p3 = section2[j_next]
                    p4 = section2[j]
                    
                    # 添加点到点列表
                    idx1 = len(points)
                    points.append(p1)
                    idx2 = len(points)
                    points.append(p2)
                    idx3 = len(points)
                    points.append(p3)
                    idx4 = len(points)
                    points.append(p4)
                    
                    # 添加四边形单元（4个点）
                    cells.append([4, idx1, idx2, idx3, idx4])
        
        return points, cells
    

