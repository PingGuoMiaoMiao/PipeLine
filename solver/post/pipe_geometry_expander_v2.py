"""
Pipe Geometry Expander V2 - 管单元几何展开模块（重构版）
使用Parallel Transport Frame和节点级截面生成，确保连续曲面
"""

import math
import numpy as np
from typing import List, Tuple, Dict, Optional, Set
from solver.mesh.node import Node
from solver.post.ovalization_estimator import OvalizationEstimator
from solver.post.parallel_transport_frame import ParallelTransportFrame


class PipeGeometryExpanderV2:
    """
    管单元几何展开器（V2版本）
    
    关键改进：
    1. 使用Parallel Transport Frame（平行传输标架）代替Frenet标架
    2. 以节点为单位生成截面点（每个节点只生成一次）
    3. 相邻节点之间生成四边形面片，确保跨单元连接处连续
    4. 输出单一连续曲面，无分段感
    """
    
    @staticmethod
    def generate_section_points(
        node_pos: np.ndarray,  # 节点位置 (mm)
        T: np.ndarray,  # 切向（归一化）
        N: np.ndarray,  # 法向（归一化）
        B: np.ndarray,  # 副法向（归一化）
        outer_radius: float,  # 外半径 (mm)
        inner_radius: float,  # 内半径 (mm)
        ovalization_amp: float,  # 椭圆化幅值 (mm)
        num_circumferential: int  # 环向划分点数
    ) -> Tuple[List[Tuple[float, float, float]], List[Tuple[float, float, float]]]:
        """
        为单个节点生成截面点（外表面和内表面）
        
        参数:
            node_pos: 节点位置 (mm)
            T: 切向（归一化）
            N: 法向（归一化）
            B: 副法向（归一化）
            outer_radius: 外半径 (mm)
            inner_radius: 内半径 (mm)
            ovalization_amp: 椭圆化幅值 (mm)，r(θ) = r0 + a * cos(2θ)
            num_circumferential: 环向划分点数
        
        返回:
            (outer_points, inner_points): 外表面点和内表面点的列表
        """
        # 生成环向角度
        phis = np.linspace(0, 2 * math.pi, num_circumferential, endpoint=False)
        
        outer_points = []
        inner_points = []
        
        for phi in phis:
            # 计算当前角度的半径（考虑椭圆化）
            # 外表面
            R_outer = outer_radius + ovalization_amp * np.cos(2 * phi)
            # 内表面
            R_inner = inner_radius + ovalization_amp * np.cos(2 * phi)
            
            # 在N-B平面内生成环向点
            # 使用N和B作为基向量：r = R * (cos(φ) * N + sin(φ) * B)
            r_outer_local = R_outer * (np.cos(phi) * N + np.sin(phi) * B)
            r_inner_local = R_inner * (np.cos(phi) * N + np.sin(phi) * B)
            
            # 截面点坐标
            outer_point = node_pos + r_outer_local
            inner_point = node_pos + r_inner_local
            
            outer_points.append(tuple(outer_point))
            inner_points.append(tuple(inner_point))
        
        return outer_points, inner_points
    
    @staticmethod
    def expand_pipe_centerline_to_surface(
        nodes: Dict[int, Node],
        elements: Dict[int, 'Element'],  # 类型注解，避免循环导入
        displacements: np.ndarray,  # 全局位移数组
        node_id_to_dof_start: Dict[int, int],  # 节点ID到DOF起始索引的映射
        sections: Dict[int, Dict],  # 截面参数字典 {sec_id: {'diameter': ..., 'wall_thickness': ...}}
        element_section_map: Optional[Dict[int, int]] = None,  # 单元到截面ID的映射
        num_circumferential: int = 20,
        material_properties: Optional[Dict] = None,  # 材料属性 {'E': ..., 'nu': ...}
        element_metadata: Optional[Dict[int, Dict]] = None,  # 单元元数据 {elem_id: {'R_curvature': ..., 'L': ..., ...}}
        ovalization_display_amplification: float = 1.0  # 椭圆化显示放大系数
    ) -> Tuple[List[Tuple[float, float, float]], List[List[int]]]:
        """
        将整个管系统的中心线展开为三维连续曲面（V2版本）
        
        关键特性：
        1. 以节点为单位生成截面点（每个节点只生成一次）
        2. 使用Parallel Transport Frame确保截面连续
        3. 相邻节点之间生成四边形面片
        4. 确保跨单元连接处网格拓扑连续
        
        参数:
            nodes: 节点字典
            elements: 单元字典（Element对象，包含node_ids和type）
            displacements: 全局位移数组
            node_id_to_dof_start: 节点ID到DOF起始索引的映射
            sections: 截面参数字典
            element_section_map: 单元到截面ID的映射
            num_circumferential: 环向划分点数
            material_properties: 材料属性字典
            element_metadata: 单元元数据字典（用于ELBOW290）
            ovalization_display_amplification: 椭圆化显示放大系数
        
        返回:
            (all_points, all_cells): 所有点和单元的列表
        """
        # 默认截面参数
        default_section = {'diameter': 30.0, 'wall_thickness': 1.0}
        
        # 1. 收集所有节点及其顺序（按单元连接顺序）
        node_order = PipeGeometryExpanderV2._get_node_order(nodes, elements)
        if len(node_order) < 2:
            return [], []
        
        # 2. 计算节点位置（变形后）
        node_positions = {}
        node_displacements_dict = {}
        for node_id in node_order:
            node = nodes[node_id]
            # 原始位置 (mm)
            node_pos_original = np.array([node.x, node.y, node.z])
            
            # 获取位移
            dof_start = node_id_to_dof_start.get(node_id, -1)
            node_disp = np.zeros(3)
            if dof_start >= 0 and dof_start + 3 <= len(displacements):
                node_disp = displacements[dof_start:dof_start+3] * 1000.0  # m -> mm
            
            # 变形后位置
            node_positions[node_id] = node_pos_original + node_disp
            node_displacements_dict[node_id] = node_disp
        
        # 3. 在节点之间增加插值点，使表面更光滑
        # 为每对相邻节点之间插入N个中间点
        num_segments_per_edge = 3  # 每两个节点之间插入3个中间点（共4段）
        interpolated_positions = []
        interpolated_node_map = []  # 记录每个插值点对应的原始节点索引
        
        for i in range(len(node_order) - 1):
            node1_id = node_order[i]
            node2_id = node_order[i + 1]
            p1 = node_positions[node1_id]
            p2 = node_positions[node2_id]
            
            # 添加第一个节点（如果是第一个节点对，避免重复）
            if i == 0:
                interpolated_positions.append(p1)
                interpolated_node_map.append(node1_id)
            
            # 在节点之间插入中间点
            for seg in range(1, num_segments_per_edge + 1):
                alpha = seg / (num_segments_per_edge + 1.0)
                p_interp = p1 * (1.0 - alpha) + p2 * alpha
                interpolated_positions.append(p_interp)
                interpolated_node_map.append(node1_id)  # 使用前一个节点的参数
            
            # 添加第二个节点
            interpolated_positions.append(p2)
            interpolated_node_map.append(node2_id)
        
        # 4. 计算Parallel Transport Frame标架（使用插值后的点）
        frames = ParallelTransportFrame.compute_frames_for_curve(interpolated_positions)
        
        # 5. 获取每个节点的截面参数和椭圆化幅值
        node_sections = PipeGeometryExpanderV2._get_node_sections(
            node_order, nodes, elements, sections, element_section_map,
            default_section, node_displacements_dict, node_id_to_dof_start,
            displacements, material_properties, element_metadata,
            ovalization_display_amplification
        )
        
        # 6. 为每个插值点生成截面点
        all_points = []
        section_point_indices = []  # 每个截面的 (outer_start_idx, inner_start_idx, num_points)
        current_point_idx = 0
        
        for idx, pos in enumerate(interpolated_positions):
            # 获取对应的原始节点ID，用于获取截面参数
            original_node_id = interpolated_node_map[idx]
            section_info = node_sections[original_node_id]
            outer_radius = section_info['outer_radius']
            inner_radius = section_info['inner_radius']
            ovalization_amp = section_info['ovalization_amp']
            
            # 获取标架
            T, N, B = frames[idx]
            
            # 生成截面点
            outer_points, inner_points = PipeGeometryExpanderV2.generate_section_points(
                pos, T, N, B, outer_radius, inner_radius,
                ovalization_amp, num_circumferential
            )
            
            # 记录索引范围
            outer_start = current_point_idx
            all_points.extend(outer_points)
            current_point_idx += len(outer_points)
            
            inner_start = current_point_idx
            all_points.extend(inner_points)
            current_point_idx += len(inner_points)
            
            section_point_indices.append((outer_start, inner_start, num_circumferential))
        
        # 7. 生成四边形面片（连接相邻截面）
        all_cells = []
        
        # 7.1 生成外表面和内表面的四边形（连接所有相邻截面，包括插值点）
        for i in range(len(section_point_indices) - 1):
            outer_start1, inner_start1, num_circ1 = section_point_indices[i]
            outer_start2, inner_start2, num_circ2 = section_point_indices[i + 1]
            
            # 外表面四边形
            for j in range(num_circumferential):
                j_next = (j + 1) % num_circumferential
                idx1_1 = outer_start1 + j
                idx1_2 = outer_start1 + j_next
                idx2_1 = outer_start2 + j
                idx2_2 = outer_start2 + j_next
                # 外表面：逆时针顺序（从外部看）
                all_cells.append([4, idx1_1, idx2_1, idx2_2, idx1_2])
            
            # 内表面四边形
            for j in range(num_circumferential):
                j_next = (j + 1) % num_circumferential
                idx1_1 = inner_start1 + j
                idx1_2 = inner_start1 + j_next
                idx2_1 = inner_start2 + j
                idx2_2 = inner_start2 + j_next
                # 内表面：顺时针顺序（从内部看，法向量向内）
                all_cells.append([4, idx1_1, idx1_2, idx2_2, idx2_1])
        
        # 7.2 生成端面（第一个和最后一个截面）
        if len(section_point_indices) > 0:
            # 第一个截面端面（入口）
            outer_start, inner_start, num_circ = section_point_indices[0]
            for j in range(num_circumferential):
                j_next = (j + 1) % num_circumferential
                idx_outer = outer_start + j
                idx_outer_next = outer_start + j_next
                idx_inner = inner_start + j
                idx_inner_next = inner_start + j_next
                # 端面：从外圆到内圆（逆时针顺序，法向量指向管道外部）
                all_cells.append([4, idx_outer, idx_outer_next, idx_inner_next, idx_inner])
            
            # 最后一个截面端面（出口）
            outer_start, inner_start, num_circ = section_point_indices[-1]
            for j in range(num_circumferential):
                j_next = (j + 1) % num_circumferential
                idx_outer = outer_start + j
                idx_outer_next = outer_start + j_next
                idx_inner = inner_start + j
                idx_inner_next = inner_start + j_next
                # 端面：从外圆到内圆（顺时针顺序，法向量指向管道外部）
                all_cells.append([4, idx_outer, idx_inner, idx_inner_next, idx_outer_next])
        
        return all_points, all_cells
    
    @staticmethod
    def _get_node_order(nodes: Dict[int, Node], elements: Dict[int, 'Element']) -> List[int]:
        """
        获取节点的顺序（按单元连接顺序）
        
        参数:
            nodes: 节点字典
            elements: 单元字典
        
        返回:
            node_order: 节点ID的有序列表
        """
        # 构建节点连接关系
        node_connections = {}  # node_id -> set of connected node_ids
        for node_id in nodes.keys():
            node_connections[node_id] = set()
        
        for elem in elements.values():
            if len(elem.node_ids) < 2:
                continue
            # 对于PIPE288（2节点）和ELBOW290（3节点），连接相邻节点
            for i in range(len(elem.node_ids) - 1):
                n1 = elem.node_ids[i]
                n2 = elem.node_ids[i + 1]
                if n1 in node_connections and n2 in node_connections:
                    node_connections[n1].add(n2)
                    node_connections[n2].add(n1)
        
        # 找到起始节点（只有一个连接的节点，或任意节点）
        start_node = None
        for node_id, connections in node_connections.items():
            if len(connections) == 1:
                start_node = node_id
                break
        if start_node is None:
            # 如果没有找到端点，使用第一个节点
            start_node = next(iter(nodes.keys()))
        
        # 按连接顺序遍历节点
        node_order = []
        visited = set()
        stack = [start_node]
        
        while stack:
            current = stack.pop()
            if current in visited:
                continue
            
            visited.add(current)
            node_order.append(current)
            
            # 添加未访问的相邻节点
            for neighbor in node_connections[current]:
                if neighbor not in visited:
                    stack.append(neighbor)
        
        # 如果还有未访问的节点（可能是孤立的），按ID顺序添加
        for node_id in sorted(nodes.keys()):
            if node_id not in visited:
                node_order.append(node_id)
        
        return node_order
    
    @staticmethod
    def _get_node_sections(
        node_order: List[int],
        nodes: Dict[int, Node],
        elements: Dict[int, 'Element'],
        sections: Dict[int, Dict],
        element_section_map: Optional[Dict[int, int]],
        default_section: Dict,
        node_displacements_dict: Dict[int, np.ndarray],
        node_id_to_dof_start: Dict[int, int],
        displacements: np.ndarray,
        material_properties: Optional[Dict],
        element_metadata: Optional[Dict[int, Dict]],
        ovalization_display_amplification: float
    ) -> Dict[int, Dict]:
        """
        获取每个节点的截面参数和椭圆化幅值
        
        参数:
            node_order: 节点顺序列表
            nodes: 节点字典
            elements: 单元字典
            sections: 截面参数字典
            element_section_map: 单元到截面ID的映射
            default_section: 默认截面参数
            node_displacements_dict: 节点位移字典
            node_id_to_dof_start: 节点ID到DOF起始索引的映射
            displacements: 全局位移数组
            material_properties: 材料属性字典
            element_metadata: 单元元数据字典
            ovalization_display_amplification: 椭圆化显示放大系数
        
        返回:
            node_sections: 节点截面参数字典 {node_id: {'outer_radius': ..., 'inner_radius': ..., 'ovalization_amp': ...}}
        """
        node_sections = {}
        
        # 构建节点到单元的映射
        node_to_elements = {}  # node_id -> list of (elem_id, elem)
        for elem_id, elem in elements.items():
            for node_id in elem.node_ids:
                if node_id not in node_to_elements:
                    node_to_elements[node_id] = []
                node_to_elements[node_id].append((elem_id, elem))
        
        # 为每个节点确定截面参数
        for node_id in node_order:
            # 找到包含该节点的第一个单元，使用其截面参数
            section = default_section
            ovalization_amp = 0.0
            
            if node_id in node_to_elements and len(node_to_elements[node_id]) > 0:
                elem_id, elem = node_to_elements[node_id][0]
                
                # 获取截面ID
                sec_id = 1  # 默认截面ID
                if element_section_map and elem_id in element_section_map:
                    sec_id = element_section_map[elem_id]
                
                section = sections.get(sec_id, default_section)
                
                # 计算半径
                D = section.get('diameter', 30.0)  # mm
                t = section.get('wall_thickness', 1.0)  # mm
                outer_radius = D / 2.0  # mm
                inner_radius = (D - 2 * t) / 2.0  # mm
                
                    # 如果是弯管（ELBOW290），计算椭圆化幅值
                if elem.type == 290 and material_properties and element_metadata and elem_id in element_metadata:
                    # 获取节点位移（标准DOF：前6个）
                    node_displacements_list = []
                    for nid in elem.node_ids[:3]:  # ELBOW290最多3个节点
                        dof_start = node_id_to_dof_start.get(nid, -1)
                        if dof_start >= 0 and dof_start + 6 <= len(displacements):
                            node_disp = displacements[dof_start:dof_start+6]
                            node_displacements_list.append(node_disp)
                        else:
                            node_displacements_list.append(np.zeros(6))
                    
                    metadata = element_metadata[elem_id]
                    R_curvature = metadata.get('R_curvature', 0.0)  # m
                    L_elem = metadata.get('L', 0.0)  # m
                    R_mean_m = metadata.get('R_mean', 0.0)  # m
                    
                    # 确保L_elem是标量
                    if hasattr(L_elem, '__len__') and not isinstance(L_elem, str):
                        L_elem = float(L_elem[0]) if len(L_elem) > 0 else 0.0
                    else:
                        L_elem = float(L_elem)
                    
                    if R_curvature > 1e-6 and L_elem > 1e-10:
                        t_m = (outer_radius - inner_radius) / 1000.0  # mm -> m
                        if R_mean_m <= 0:
                            R_mean_m = (outer_radius + inner_radius) / 2.0 / 1000.0  # mm -> m
                        
                        Ro_m = outer_radius / 1000.0
                        Ri_m = inner_radius / 1000.0
                        I = np.pi / 4.0 * (Ro_m**4 - Ri_m**4)
                        E = material_properties.get('E', 200e9)
                        
                        # 估计椭圆化幅值 (m)
                        # 注意：estimate_ovalization_for_element的第一个参数是node_displacements列表
                        ovalization_amp_m = OvalizationEstimator.estimate_ovalization_for_element(
                            node_displacements_list,
                            R_curvature, R_mean_m, t_m, E, I, L_elem,
                            display_amplification=ovalization_display_amplification
                        )
                        # 转换为mm
                        ovalization_amp = ovalization_amp_m * 1000.0
            else:
                # 使用默认截面
                D = default_section.get('diameter', 30.0)
                t = default_section.get('wall_thickness', 1.0)
                outer_radius = D / 2.0
                inner_radius = (D - 2 * t) / 2.0
            
            node_sections[node_id] = {
                'outer_radius': outer_radius,
                'inner_radius': inner_radius,
                'ovalization_amp': ovalization_amp
            }
        
        return node_sections

