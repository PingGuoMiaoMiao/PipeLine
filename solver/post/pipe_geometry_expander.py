"""
Pipe Geometry Expander - 管单元几何展开模块
独立的后处理模块，将管单元中心线结果展开为三维管壳曲面
完全解耦于求解器，仅基于节点位移、转角、椭圆化参数生成几何

椭圆化处理模式：显示主导、弱耦合
- 椭圆化模态不参与单元刚度矩阵
- 椭圆化幅值在后处理阶段根据弯矩估计或用户指定
"""

import math
import numpy as np
from typing import List, Tuple, Dict, Optional
from solver.mesh.node import Node
from solver.post.ovalization_estimator import OvalizationEstimator
from solver.post.frenet_frame import FrenetFrame


class PipeGeometryExpander:
    """
    管单元几何展开器
    
    功能：
    1. 基于中心线节点坐标和位移生成三维管壳曲面
    2. 支持直管和弯管
    3. 支持椭圆化变形（r(θ) = r0 + a * cos(2θ)）
    4. 输出VTK PolyData格式的四边形面片
    """
    
    @staticmethod
    def expand_pipe_segment(
        node1_pos: np.ndarray,  # 节点1位置 (mm)
        node2_pos: np.ndarray,  # 节点2位置 (mm)
        node1_disp: np.ndarray,  # 节点1位移 (mm) [ux, uy, uz]
        node2_disp: np.ndarray,  # 节点2位移 (mm) [ux, uy, uz]
        node1_rot: Optional[np.ndarray] = None,  # 节点1转角 (rad) [rotx, roty, rotz]
        node2_rot: Optional[np.ndarray] = None,  # 节点2转角 (rad) [rotx, roty, rotz]
        outer_radius: float = 15.0,  # 外半径 (mm)
        inner_radius: float = 14.0,  # 内半径 (mm)
        ovalization_amp: float = 0.0,  # 椭圆化幅值 (mm)，r(θ) = r0 + a * cos(2θ)
        num_circumferential: int = 20,  # 环向划分点数
        is_bent: bool = False,  # 是否为弯管
        node3_pos: Optional[np.ndarray] = None,  # 节点3位置（用于弯管）(mm)
        node3_disp: Optional[np.ndarray] = None  # 节点3位移（用于弯管）(mm)
    ) -> Tuple[List[Tuple[float, float, float]], List[List[int]]]:
        """
        展开单个管段为三维曲面（使用Frenet标架）
        
        参数:
            node1_pos: 节点1的原始位置 (mm)
            node2_pos: 节点2的原始位置 (mm)
            node1_disp: 节点1的位移 (mm)
            node2_disp: 节点2的位移 (mm)
            node1_rot: 节点1的转角 (rad)，可选
            node2_rot: 节点2的转角 (rad)，可选
            outer_radius: 外半径 (mm)
            inner_radius: 内半径 (mm)
            ovalization_amp: 椭圆化幅值 (mm)
            num_circumferential: 环向划分点数
            is_bent: 是否为弯管
            node3_pos: 节点3位置（用于弯管）(mm)，可选
            node3_disp: 节点3位移（用于弯管）(mm)，可选
        
        返回:
            (points, cells): 点列表和单元列表（四边形面片）
        """
        # 变形后的节点位置
        deformed_node1 = node1_pos + node1_disp
        deformed_node2 = node2_pos + node2_disp
        
        # 为两个节点分别计算Frenet标架（确保连续）
        if is_bent and node3_pos is not None:
            # 弯管：使用三点计算Frenet标架
            deformed_node3 = node3_pos + (node3_disp if node3_disp is not None else np.zeros(3))
            # 节点1的标架（两点：使用节点1和节点2）
            T1, N1, B1 = FrenetFrame.compute_frenet_frame(deformed_node1, deformed_node2)
            # 节点2的标架（三点：使用节点1、节点2、节点3）
            T2, N2, B2 = FrenetFrame.compute_frenet_frame(deformed_node1, deformed_node2, deformed_node3)
        else:
            # 直管：两个节点使用相同的标架
            T, N, B = FrenetFrame.compute_frenet_frame(deformed_node1, deformed_node2)
            T1, N1, B1 = T, N, B
            T2, N2, B2 = T, N, B
        
        # 确保N和B的方向连续（防止翻转）
        # 如果N2与N1反向，翻转N2和B2
        if np.dot(N1, N2) < 0:
            N2 = -N2
            B2 = -B2
        
        # 生成环向角度
        phis = np.linspace(0, 2 * math.pi, num_circumferential, endpoint=False)
        
        points = []
        cells = []
        
        # 为每个半径（外表面和内表面）生成点
        for radius_type, R_base in [('outer', outer_radius), ('inner', inner_radius)]:
            section_points = []
            
            # 为节点1和节点2生成截面点（使用各自的Frenet标架）
            for node_idx, (node_pos, N_frame, B_frame) in enumerate([
                (deformed_node1, N1, B1),
                (deformed_node2, N2, B2)
            ]):
                for phi in phis:
                    # 计算当前角度的半径（考虑椭圆化）
                    R = R_base + ovalization_amp * np.cos(2 * phi)
                    
                    # 在N-B平面内生成环向点
                    # 使用N和B作为基向量
                    r_local = R * (np.cos(phi) * N_frame + np.sin(phi) * B_frame)
                    
                    # 截面点坐标
                    section_point = node_pos + r_local
                    section_points.append(tuple(section_point))
            
            # 添加点到总点列表
            start_idx = len(points)
            points.extend(section_points)
            end_idx = len(points)
            
            # 生成四边形面片（连接节点1和节点2的截面）
            # 每个半径类型有 num_circumferential 个四边形
            for i in range(num_circumferential):
                i_next = (i + 1) % num_circumferential
                
                # 节点1处的两个点
                idx1_1 = start_idx + i
                idx1_2 = start_idx + i_next
                
                # 节点2处的两个点
                idx2_1 = start_idx + num_circumferential + i
                idx2_2 = start_idx + num_circumferential + i_next
                
                # 形成四边形（注意顶点顺序，确保法向量正确）
                if radius_type == 'outer':
                    # 外表面：逆时针顺序
                    cells.append([4, idx1_1, idx2_1, idx2_2, idx1_2])
                else:
                    # 内表面：顺时针顺序（法向量向内）
                    cells.append([4, idx1_1, idx1_2, idx2_2, idx2_1])
        
        # 生成端面（节点1和节点2处的圆环）
        outer_start = 0
        inner_start = num_circumferential * 2
        
        # 节点1端面
        for i in range(num_circumferential):
            i_next = (i + 1) % num_circumferential
            idx_outer1 = outer_start + i
            idx_outer1_next = outer_start + i_next
            idx_inner1 = inner_start + i
            idx_inner1_next = inner_start + i_next
            # 外圆到内圆的四边形
            cells.append([4, idx_outer1, idx_outer1_next, idx_inner1_next, idx_inner1])
        
        # 节点2端面
        for i in range(num_circumferential):
            i_next = (i + 1) % num_circumferential
            idx_outer2 = outer_start + num_circumferential + i
            idx_outer2_next = outer_start + num_circumferential + i_next
            idx_inner2 = inner_start + num_circumferential + i
            idx_inner2_next = inner_start + num_circumferential + i_next
            # 外圆到内圆的四边形（注意方向）
            cells.append([4, idx_outer2, idx_inner2, idx_inner2_next, idx_outer2_next])
        
        return points, cells
    
    @staticmethod
    def expand_pipe_centerline_to_surface(
        nodes: Dict[int, Node],
        elements: Dict[int, 'Element'],  # 类型注解，避免循环导入
        displacements: np.ndarray,  # 全局位移数组
        node_id_to_dof_start: Dict[int, int],  # 节点ID到DOF起始索引的映射
        sections: Dict[int, Dict],  # 截面参数字典 {sec_id: {'diameter': ..., 'wall_thickness': ...}}
        element_section_map: Optional[Dict[int, int]] = None,  # 单元到截面ID的映射
        ovalization_dofs: Optional[Dict[int, np.ndarray]] = None,  # 椭圆化自由度（已弃用，仅保留兼容性）
        num_circumferential: int = 20,
        material_properties: Optional[Dict] = None,  # 材料属性 {'E': ..., 'nu': ...}
        element_metadata: Optional[Dict[int, Dict]] = None,  # 单元元数据 {elem_id: {'R_curvature': ..., 'L': ..., ...}}
        ovalization_display_amplification: float = 1.0  # 椭圆化显示放大系数
    ) -> Tuple[List[Tuple[float, float, float]], List[List[int]]]:
        """
        将整个管系统的中心线展开为三维曲面
        
        参数:
            nodes: 节点字典
            elements: 单元字典（Element对象，包含node_ids和type）
            displacements: 全局位移数组
            node_id_to_dof_start: 节点ID到DOF起始索引的映射
            sections: 截面参数字典
            element_section_map: 单元到截面ID的映射，如果为None则使用默认截面
            ovalization_dofs: 椭圆化自由度字典（已弃用）
            num_circumferential: 环向划分点数
            material_properties: 材料属性字典
            element_metadata: 单元元数据字典（用于ELBOW290）
            ovalization_display_amplification: 椭圆化显示放大系数
        
        返回:
            (all_points, all_cells): 所有点和单元的列表
        """
        all_points = []
        all_cells = []
        current_point_idx = 0
        
        # 默认截面参数
        default_section = {'diameter': 30.0, 'wall_thickness': 1.0}
        
        # 遍历所有单元
        sorted_elements = sorted(elements.items())
        for elem_id, elem in sorted_elements:
            if len(elem.node_ids) < 2:
                continue
            
            # 获取节点
            node1_id = elem.node_ids[0]
            node2_id = elem.node_ids[1]
            
            if node1_id not in nodes or node2_id not in nodes:
                continue
            
            node1 = nodes[node1_id]
            node2 = nodes[node2_id]
            
            # 获取节点位置 (mm)
            node1_pos = np.array([node1.x, node1.y, node1.z])
            node2_pos = np.array([node2.x, node2.y, node2.z])
            
            # 获取节点位移
            dof_start1 = node_id_to_dof_start.get(node1_id, -1)
            dof_start2 = node_id_to_dof_start.get(node2_id, -1)
            
            if dof_start1 < 0 or dof_start2 < 0:
                continue
            
            # 提取位移（前3个DOF为平移）
            node1_disp = np.zeros(3)
            node2_disp = np.zeros(3)
            
            if dof_start1 + 3 <= len(displacements):
                node1_disp = displacements[dof_start1:dof_start1+3] * 1000.0  # m -> mm
            if dof_start2 + 3 <= len(displacements):
                node2_disp = displacements[dof_start2:dof_start2+3] * 1000.0  # m -> mm
            
            # 提取转角（可选，DOF 3-5）
            node1_rot = None
            node2_rot = None
            if dof_start1 + 6 <= len(displacements):
                node1_rot = displacements[dof_start1+3:dof_start1+6]
            if dof_start2 + 6 <= len(displacements):
                node2_rot = displacements[dof_start2+3:dof_start2+6]
            
            # 获取截面参数
            sec_id = 1  # 默认截面ID
            if element_section_map and elem_id in element_section_map:
                sec_id = element_section_map[elem_id]
            
            section = sections.get(sec_id, default_section)
            D = section.get('diameter', 30.0)  # mm
            t = section.get('wall_thickness', 1.0)  # mm
            outer_radius = D / 2.0  # mm
            inner_radius = (D - 2 * t) / 2.0  # mm
            
            # 计算单元长度（用于椭圆化估计）
            L = np.linalg.norm(node2_pos - node1_pos)  # mm
            L_m = L / 1000.0  # 转换为m
            
            # 获取椭圆化幅值（显示主导模式）
            ovalization_amp = 0.0
            is_bent = (elem.type == 290)
            
            if is_bent:
                # 弯管：根据弯矩经验公式估计椭圆化幅值
                # 椭圆化模态不参与刚度矩阵，仅用于显示
                if element_metadata and elem_id in element_metadata and material_properties:
                    metadata = element_metadata[elem_id]
                    R_curvature = metadata.get('R_curvature', 0.0)  # m
                    L_elem = metadata.get('L', L_m)  # m
                    R_mean_m = metadata.get('R_mean', 0.0)  # m
                else:
                    # 如果没有元数据，使用单元长度估算（简化）
                    R_curvature = L_m * 2.0  # 假设曲率半径为长度的2倍 (m)
                    L_elem = L_m  # m
                    R_mean_m = (outer_radius + inner_radius) / 2.0 / 1000.0  # mm -> m
                
                if R_curvature > 1e-6 and material_properties:
                    # 收集节点位移（标准DOF：前6个）
                    node_displacements = []
                    for node_id in elem.node_ids[:3]:  # ELBOW290最多3个节点
                        dof_start = node_id_to_dof_start.get(node_id, -1)
                        if dof_start >= 0 and dof_start + 6 <= len(displacements):
                            node_disp = displacements[dof_start:dof_start+6]  # [ux, uy, uz, rotx, roty, rotz] (m, rad)
                            node_displacements.append(node_disp)
                        else:
                            node_displacements.append(np.zeros(6))
                    
                    # 计算截面参数（转换为m）
                    t_m = (outer_radius - inner_radius) / 1000.0  # mm -> m
                    if R_mean_m <= 0:
                        R_mean_m = (outer_radius + inner_radius) / 2.0 / 1000.0  # mm -> m
                    
                    # 计算惯性矩 (m⁴)
                    Ro_m = outer_radius / 1000.0
                    Ri_m = inner_radius / 1000.0
                    I = np.pi / 4.0 * (Ro_m**4 - Ri_m**4)
                    E = material_properties.get('E', 200e9)
                    
                    # 估计椭圆化幅值 (m)
                    ovalization_amp_m = OvalizationEstimator.estimate_ovalization_for_element(
                        node_displacements,
                        R_curvature,
                        R_mean_m,
                        t_m,
                        E,
                        I,
                        L_elem,
                        display_amplification=ovalization_display_amplification
                    )
                    # 转换为mm（因为后续使用mm单位）
                    ovalization_amp = ovalization_amp_m * 1000.0
            
            # 获取节点3信息（用于弯管）
            node3_pos = None
            node3_disp = None
            if is_bent and len(elem.node_ids) >= 3:
                node3_id = elem.node_ids[2]
                if node3_id in nodes:
                    node3 = nodes[node3_id]
                    node3_pos = np.array([node3.x, node3.y, node3.z])
                    dof_start3 = node_id_to_dof_start.get(node3_id, -1)
                    if dof_start3 >= 0 and dof_start3 + 3 <= len(displacements):
                        node3_disp = displacements[dof_start3:dof_start3+3] * 1000.0  # m -> mm
            
            # 展开该单元（使用Frenet标架）
            points, cells = PipeGeometryExpander.expand_pipe_segment(
                node1_pos, node2_pos,
                node1_disp, node2_disp,
                node1_rot, node2_rot,
                outer_radius, inner_radius,
                ovalization_amp,
                num_circumferential,
                is_bent,
                node3_pos,
                node3_disp
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
