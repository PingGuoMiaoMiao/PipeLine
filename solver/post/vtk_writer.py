"""
VTK Writer - VTK文件写入器
"""

import math
from typing import Dict, Optional, Tuple
import numpy as np
from solver.mesh.node import Node
from solver.mesh.element import Element


class VTKWriter:
    """VTK文件写入器（PolyLine格式，支持ELBOW290曲面）"""
    
    @staticmethod
    def write_polylines(filename: str, elements: Dict[int, Element], 
                       nodes: Dict[int, Node],
                       displacements: Optional[np.ndarray] = None,
                       sections: Optional[Dict[int, Dict]] = None,
                       element_section_map: Optional[Dict[int, int]] = None,
                       ovalization_dofs: Optional[Dict[int, np.ndarray]] = None,
                       output_surface: bool = True,
                       material_properties: Optional[Dict] = None,
                       element_metadata: Optional[Dict[int, Dict]] = None,
                       ovalization_display_amplification: float = 1.0):
        """
        将单元和节点写入VTK文件
        
        参数:
            filename: 输出文件名
            elements: 单元字典
            nodes: 节点字典
            displacements: 位移数组（可选，单位：m）
            sections: 截面参数字典 {sec_id: {'diameter': ..., 'wall_thickness': ...}}
            element_section_map: 单元到截面ID的映射
            ovalization_dofs: 椭圆化自由度字典（用于ELBOW290）
            output_surface: 是否输出为表面（True）或线条（False）
        """
        if output_surface and displacements is not None:
            # 检查是否有弯管（ELBOW290）单元
            has_elbow290 = any(elem.type == 290 for elem in elements.values())
            
            if has_elbow290:
                # 弯管使用V2版本：Parallel Transport Frame + 节点级截面
                VTKWriter._write_surface_from_centerline_v2(
                    filename, elements, nodes, displacements, 
                    sections, element_section_map, ovalization_dofs,
                    material_properties,
                    element_metadata,
                    ovalization_display_amplification
                )
            else:
                # 直管使用原来的方法（基于单元的截面生成）
                VTKWriter._write_surface_from_centerline(
                    filename, elements, nodes, displacements, 
                    sections, element_section_map, ovalization_dofs,
                    material_properties,
                    element_metadata,
                    ovalization_display_amplification
                )
        else:
            # 输出为线条（向后兼容）
            VTKWriter._write_polyline_vtk_legacy(filename, elements, nodes, displacements)
    
    @staticmethod
    def _write_polyline_vtk_legacy(filename: str, elements: Dict[int, Element],
                           nodes: Dict[int, Node],
                           displacements: Optional[np.ndarray]):
        """写入PolyLine格式的VTK文件，如果提供pipe288_elements则输出为表面"""
        # 如果有PIPE288单元对象，输出为表面
        if pipe288_elements is not None and len(pipe288_elements) > 0:
            VTKWriter._write_pipe288_surface_vtk(
                filename, elements, nodes, displacements, pipe288_elements)
            return
        
        # 否则输出为线条
        with open(filename, 'w') as f:
            # VTK文件头
            f.write("# vtk DataFile Version 2.0\n")
            f.write("Pipe Elements with Displacements\n")
            f.write("ASCII\n")
            f.write("DATASET POLYDATA\n")
            f.write("\n")
            
            # 写入所有节点（变形后的坐标）
            f.write(f"POINTS {len(nodes)} double\n")
            sorted_nodes = sorted(nodes.items())
            node_id_to_index = {
                node_id: idx for idx, (node_id, _) in enumerate(sorted_nodes)
            }
            
            # 计算每个节点的DOF数（PIPE288为6，ELBOW290为10）
            dof_per_node = 6
            if displacements is not None and len(displacements) > 0:
                dof_per_node = len(displacements) // len(nodes)
            
            for idx, (node_id, node) in enumerate(sorted_nodes):
                if displacements is not None and idx * dof_per_node < len(displacements):
                    # 添加位移
                    ux = displacements[idx * dof_per_node] * 1000.0  # m -> mm
                    uy = displacements[idx * dof_per_node + 1] * 1000.0
                    uz = displacements[idx * dof_per_node + 2] * 1000.0
                    x = node.x + ux
                    y = node.y + uy
                    z = node.z + uz
                else:
                    x, y, z = node.x, node.y, node.z
                f.write(f"{x:.6e} {y:.6e} {z:.6e}\n")
            
            # 写入PolyLine
            f.write(f"\nLINES {len(elements)} {len(elements) + sum(len(e.node_ids) for e in elements.values())}\n")
            sorted_elements = sorted(elements.items())
            for elem_id, elem in sorted_elements:
                valid_node_indices = [
                    node_id_to_index[nid] 
                    for nid in elem.node_ids 
                    if nid in node_id_to_index
                ]
                if valid_node_indices:
                    f.write(f"{len(valid_node_indices)}")
                    for idx in valid_node_indices:
                        f.write(f" {idx}")
                    f.write("\n")
            
            # 写入单元类型数据
            f.write(f"\nCELL_DATA {len(elements)}\n")
            f.write("SCALARS ElementType int 1\n")
            f.write("LOOKUP_TABLE default\n")
            for elem_id, elem in sorted_elements:
                f.write(f"{elem.type}\n")
            
            # 写入节点位移数据
            if displacements is not None:
                f.write(f"\nPOINT_DATA {len(nodes)}\n")
                f.write("VECTORS Displacement double\n")
                for idx in range(len(nodes)):
                    if idx * dof_per_node < len(displacements):
                        ux = displacements[idx * dof_per_node] * 1000.0  # m -> mm
                        uy = displacements[idx * dof_per_node + 1] * 1000.0
                        uz = displacements[idx * dof_per_node + 2] * 1000.0
                    else:
                        ux = uy = uz = 0.0
                    f.write(f"{ux:.6e} {uy:.6e} {uz:.6e}\n")
                
                # 写入位移幅值
                f.write("\nSCALARS DisplacementMagnitude double 1\n")
                f.write("LOOKUP_TABLE default\n")
                for idx in range(len(nodes)):
                    if idx * dof_per_node < len(displacements):
                        ux = displacements[idx * dof_per_node] * 1000.0  # m -> mm
                        uy = displacements[idx * dof_per_node + 1] * 1000.0
                        uz = displacements[idx * dof_per_node + 2] * 1000.0
                    else:
                        ux = uy = uz = 0.0
                mag = math.sqrt(ux**2 + uy**2 + uz**2)
                f.write(f"{mag:.6e}\n")
    
    @staticmethod
    def write_polylines_with_creep(filename: str, elements: Dict[int, Element],
                                   nodes: Dict[int, Node],
                                   displacements: np.ndarray,
                                   element_objects: Dict = None):
        """
        写入包含蠕变应变信息的VTK文件
        
        参数:
            filename: 输出文件名
            elements: 单元字典
            nodes: 节点字典
            displacements: 位移数组
            element_objects: 单元对象字典（用于获取蠕变状态）
        """
        with open(filename, 'w') as f:
            # VTK文件头
            f.write("# vtk DataFile Version 2.0\n")
            f.write("Pipe Elements with Creep\n")
            f.write("ASCII\n")
            f.write("DATASET POLYDATA\n")
            f.write("\n")
            
            # 写入节点
            f.write(f"POINTS {len(nodes)} double\n")
            sorted_nodes = sorted(nodes.items())
            node_id_to_index = {
                node_id: idx for idx, (node_id, _) in enumerate(sorted_nodes)
            }
            
            dof_per_node = len(displacements) // len(nodes)
            for idx, (node_id, node) in enumerate(sorted_nodes):
                if idx * dof_per_node < len(displacements):
                    ux = displacements[idx * dof_per_node] * 1000.0
                    uy = displacements[idx * dof_per_node + 1] * 1000.0
                    uz = displacements[idx * dof_per_node + 2] * 1000.0
                    x = node.x + ux
                    y = node.y + uy
                    z = node.z + uz
                else:
                    x, y, z = node.x, node.y, node.z
                f.write(f"{x:.6e} {y:.6e} {z:.6e}\n")
            
            # 写入PolyLine
            f.write(f"\nLINES {len(elements)} {len(elements) + sum(len(e.node_ids) for e in elements.values())}\n")
            sorted_elements = sorted(elements.items())
            for elem_id, elem in sorted_elements:
                valid_node_indices = [
                    node_id_to_index[nid] 
                    for nid in elem.node_ids 
                    if nid in node_id_to_index
                ]
                if valid_node_indices:
                    f.write(f"{len(valid_node_indices)}")
                    for idx in valid_node_indices:
                        f.write(f" {idx}")
                    f.write("\n")
            
            # 写入单元类型数据
            f.write(f"\nCELL_DATA {len(elements)}\n")
            f.write("SCALARS ElementType int 1\n")
            f.write("LOOKUP_TABLE default\n")
            for elem_id, elem in sorted_elements:
                f.write(f"{elem.type}\n")
            
            # 写入蠕变应变数据
            if element_objects:
                f.write("SCALARS EquivalentCreepStrain double 1\n")
                f.write("LOOKUP_TABLE default\n")
                for elem_id, elem in sorted_elements:
                    if elem_id in element_objects:
                        creep_state = element_objects[elem_id].get_creep_state()
                        f.write(f"{creep_state['equivalent_creep_strain']:.6e}\n")
                    else:
                        f.write("0.0\n")
            
            # 写入节点位移数据
            f.write(f"\nPOINT_DATA {len(nodes)}\n")
            f.write("VECTORS Displacement double\n")
            for idx in range(len(nodes)):
                if idx * dof_per_node < len(displacements):
                    ux = displacements[idx * dof_per_node] * 1000.0
                    uy = displacements[idx * dof_per_node + 1] * 1000.0
                    uz = displacements[idx * dof_per_node + 2] * 1000.0
                else:
                    ux = uy = uz = 0.0
                f.write(f"{ux:.6e} {uy:.6e} {uz:.6e}\n")
            
            # 写入位移幅值
            f.write("\nSCALARS DisplacementMagnitude double 1\n")
            f.write("LOOKUP_TABLE default\n")
            for idx in range(len(nodes)):
                if idx * dof_per_node < len(displacements):
                    ux = displacements[idx * dof_per_node] * 1000.0
                    uy = displacements[idx * dof_per_node + 1] * 1000.0
                    uz = displacements[idx * dof_per_node + 2] * 1000.0
                else:
                    ux = uy = uz = 0.0
                mag = math.sqrt(ux**2 + uy**2 + uz**2)
                f.write(f"{mag:.6e}\n")
    
    @staticmethod
    def write_polylines_with_plasticity(filename: str, elements: Dict[int, Element],
                                       nodes: Dict[int, Node],
                                       displacements: np.ndarray,
                                       element_objects: Dict = None):
        """
        写入包含塑性状态信息的VTK文件
        
        参数:
            filename: 输出文件名
            elements: 单元字典
            nodes: 节点字典
            displacements: 位移数组
            element_objects: 单元对象字典（用于获取塑性状态）
        """
        with open(filename, 'w', encoding='utf-8') as f:
            # VTK文件头
            f.write("# vtk DataFile Version 2.0\n")
            f.write("Pipe Elements with Plasticity\n")
            f.write("ASCII\n")
            f.write("DATASET POLYDATA\n")
            f.write("\n")
            
            # 写入节点
            f.write(f"POINTS {len(nodes)} double\n")
            sorted_nodes = sorted(nodes.items())
            node_id_to_index = {
                node_id: idx for idx, (node_id, _) in enumerate(sorted_nodes)
            }
            
            dof_per_node = len(displacements) // len(nodes)
            for idx, (node_id, node) in enumerate(sorted_nodes):
                if idx * dof_per_node < len(displacements):
                    ux = displacements[idx * dof_per_node] * 1000.0
                    uy = displacements[idx * dof_per_node + 1] * 1000.0
                    uz = displacements[idx * dof_per_node + 2] * 1000.0
                    x = node.x + ux
                    y = node.y + uy
                    z = node.z + uz
                else:
                    x, y, z = node.x, node.y, node.z
                f.write(f"{x:.6e} {y:.6e} {z:.6e}\n")
            
            # 写入PolyLine
            f.write(f"\nLINES {len(elements)} {len(elements) + sum(len(e.node_ids) for e in elements.values())}\n")
            sorted_elements = sorted(elements.items())
            for elem_id, elem in sorted_elements:
                valid_node_indices = [
                    node_id_to_index[nid] 
                    for nid in elem.node_ids 
                    if nid in node_id_to_index
                ]
                if valid_node_indices:
                    f.write(f"{len(valid_node_indices)}")
                    for idx in valid_node_indices:
                        f.write(f" {idx}")
                    f.write("\n")
            
            # 写入单元类型数据
            f.write(f"\nCELL_DATA {len(elements)}\n")
            f.write("SCALARS ElementType int 1\n")
            f.write("LOOKUP_TABLE default\n")
            for elem_id, elem in sorted_elements:
                f.write(f"{elem.type}\n")
            
            # 写入塑性状态数据
            if element_objects:
                f.write("SCALARS PlasticState int 1\n")
                f.write("LOOKUP_TABLE default\n")
                for elem_id, elem in sorted_elements:
                    if elem_id in element_objects:
                        plastic_state = element_objects[elem_id].get_plastic_state()
                        f.write(f"{1 if plastic_state['is_plastic'] else 0}\n")
                    else:
                        f.write("0\n")
                
                f.write("SCALARS EquivalentPlasticStrain double 1\n")
                f.write("LOOKUP_TABLE default\n")
                for elem_id, elem in sorted_elements:
                    if elem_id in element_objects:
                        plastic_state = element_objects[elem_id].get_plastic_state()
                        f.write(f"{plastic_state['equivalent_plastic_strain']:.6e}\n")
                    else:
                        f.write("0.0\n")
            
            # 写入节点位移数据
            f.write(f"\nPOINT_DATA {len(nodes)}\n")
            f.write("VECTORS Displacement double\n")
            for idx in range(len(nodes)):
                if idx * dof_per_node < len(displacements):
                    ux = displacements[idx * dof_per_node] * 1000.0
                    uy = displacements[idx * dof_per_node + 1] * 1000.0
                    uz = displacements[idx * dof_per_node + 2] * 1000.0
                else:
                    ux = uy = uz = 0.0
                f.write(f"{ux:.6e} {uy:.6e} {uz:.6e}\n")
            
            # 写入位移幅值
            f.write("\nSCALARS DisplacementMagnitude double 1\n")
            f.write("LOOKUP_TABLE default\n")
            for idx in range(len(nodes)):
                if idx * dof_per_node < len(displacements):
                    ux = displacements[idx * dof_per_node] * 1000.0
                    uy = displacements[idx * dof_per_node + 1] * 1000.0
                    uz = displacements[idx * dof_per_node + 2] * 1000.0
                else:
                    ux = uy = uz = 0.0
                mag = math.sqrt(ux**2 + uy**2 + uz**2)
                f.write(f"{mag:.6e}\n")
    
    @staticmethod
    def _write_surface_vtk(filename: str, elements: Dict[int, Element],
                          nodes: Dict[int, Node],
                          displacements: Optional[np.ndarray],
                          elbow290_elements: Dict):
        """写入曲面格式的VTK文件（支持ELBOW290）"""
        from solver.post.pipe_surface_expand import PipeSurfaceExpander
        
        # 创建节点ID到DOF索引的映射（ELBOW290每个节点10个DOF）
        # 需要找到所有ELBOW290单元使用的节点
        elbow290_node_ids = set()
        for elem in elbow290_elements.values():
            elbow290_node_ids.add(elem.node1.id)
            elbow290_node_ids.add(elem.node2.id)
            elbow290_node_ids.add(elem.node3.id)
        
        sorted_elbow_node_ids = sorted(elbow290_node_ids)
        node_dof_map = {}
        dof_idx = 0
        for node_id in sorted_elbow_node_ids:
            node_dof_map[node_id] = dof_idx
            # ELBOW290每个节点10个DOF
            dof_idx += 10
        
        # 扩展ELBOW290单元为曲面
        points, cells = PipeSurfaceExpander.expand_elbow290_to_surface(
            elbow290_elements, nodes, displacements, node_dof_map, num_circumferential=20)
        
        with open(filename, 'w') as f:
            # VTK文件头
            f.write("# vtk DataFile Version 2.0\n")
            f.write("ELBOW290 Pipe Elements with Surface\n")
            f.write("ASCII\n")
            f.write("DATASET POLYDATA\n")
            f.write("\n")
            
            # 写入点
            f.write(f"POINTS {len(points)} double\n")
            for point in points:
                f.write(f"{point[0]:.6e} {point[1]:.6e} {point[2]:.6e}\n")
            
            # 写入多边形（四边形）
            total_cell_size = sum(1 + len(cell[1:]) for cell in cells)
            f.write(f"\nPOLYGONS {len(cells)} {total_cell_size}\n")
            for cell in cells:
                f.write(f"{cell[0]}")
                for idx in cell[1:]:
                    f.write(f" {idx}")
                f.write("\n")
            
            # 写入单元类型数据
            f.write(f"\nCELL_DATA {len(cells)}\n")
            f.write("SCALARS ElementType int 1\n")
            f.write("LOOKUP_TABLE default\n")
            for _ in cells:
                f.write("290\n")  # ELBOW290
    
    @staticmethod
    def _write_pipe288_surface_vtk(filename: str, elements: Dict[int, Element],
                                   nodes: Dict[int, Node],
                                   displacements: Optional[np.ndarray],
                                   pipe288_elements: Dict):
        """写入PIPE288单元的表面格式VTK文件"""
        from solver.post.pipe_surface_generator import PipeSurfaceGenerator
        
        # 创建节点ID到DOF索引的映射
        sorted_node_ids = sorted(nodes.keys())
        node_id_to_dof_start = {
            node_id: idx * 6 for idx, node_id in enumerate(sorted_node_ids)
        }
        
        # 生成表面
        points, cells = PipeSurfaceGenerator.generate_pipe288_surfaces_for_elements(
            pipe288_elements, nodes, displacements, node_id_to_dof_start, num_circumferential=20)
        
        with open(filename, 'w', encoding='utf-8') as f:
            # VTK文件头
            f.write("# vtk DataFile Version 2.0\n")
            f.write("PIPE288 Pipe Elements with Surface\n")
            f.write("ASCII\n")
            f.write("DATASET POLYDATA\n")
            f.write("\n")
            
            # 写入点
            f.write(f"POINTS {len(points)} double\n")
            for point in points:
                f.write(f"{point[0]:.6e} {point[1]:.6e} {point[2]:.6e}\n")
            
            # 写入多边形（四边形）
            total_cell_size = sum(1 + len(cell[1:]) for cell in cells)
            f.write(f"\nPOLYGONS {len(cells)} {total_cell_size}\n")
            for cell in cells:
                f.write(f"{cell[0]}")
                for idx in cell[1:]:
                    f.write(f" {idx}")
                f.write("\n")
            
            # 写入单元类型数据
            f.write(f"\nCELL_DATA {len(cells)}\n")
            f.write("SCALARS ElementType int 1\n")
            f.write("LOOKUP_TABLE default\n")
            for _ in cells:
                f.write("288\n")  # PIPE288
    
    @staticmethod
    def _write_surface_from_centerline(
        filename: str,
        elements: Dict[int, Element],
        nodes: Dict[int, Node],
        displacements: np.ndarray,
        sections: Optional[Dict[int, Dict]] = None,
        element_section_map: Optional[Dict[int, int]] = None,
        ovalization_dofs: Optional[Dict[int, np.ndarray]] = None,
        material_properties: Optional[Dict] = None,
        element_metadata: Optional[Dict[int, Dict]] = None,
        ovalization_display_amplification: float = 1.0
    ):
        """
        基于中心线结果生成三维表面并写入VTK文件
        
        这是完全解耦的后处理模块，不依赖求解器的内部实现
        仅基于节点坐标、位移、截面参数生成几何
        """
        from solver.post.pipe_geometry_expander import PipeGeometryExpander
        
        # 创建节点ID到DOF索引的映射
        sorted_node_ids = sorted(nodes.keys())
        
        # 确定每个节点的DOF数（PIPE288为6，ELBOW290为10）
        dof_per_node = len(displacements) // len(nodes) if len(displacements) > 0 and len(nodes) > 0 else 6
        
        node_id_to_dof_start = {
            node_id: idx * dof_per_node for idx, node_id in enumerate(sorted_node_ids)
        }
        
        # 如果没有提供截面参数，使用默认值
        if sections is None:
            sections = {1: {'diameter': 30.0, 'wall_thickness': 1.0}}
        
        # 展开中心线为三维表面
        points, cells = PipeGeometryExpander.expand_pipe_centerline_to_surface(
            nodes, elements, displacements,
            node_id_to_dof_start,
            sections, element_section_map,
            ovalization_dofs,
            num_circumferential=20,
            material_properties=material_properties,
            element_metadata=element_metadata,
            ovalization_display_amplification=ovalization_display_amplification
        )
        
        # 写入VTK文件
        with open(filename, 'w', encoding='utf-8') as f:
            # VTK文件头
            f.write("# vtk DataFile Version 2.0\n")
            f.write("Pipe Elements - Surface Geometry from Centerline\n")
            f.write("ASCII\n")
            f.write("DATASET POLYDATA\n")
            f.write("\n")
            
            # 写入点
            f.write(f"POINTS {len(points)} double\n")
            for point in points:
                f.write(f"{point[0]:.6e} {point[1]:.6e} {point[2]:.6e}\n")
            
            # 写入多边形（四边形）
            total_cell_size = sum(1 + len(cell[1:]) for cell in cells)
            f.write(f"\nPOLYGONS {len(cells)} {total_cell_size}\n")
            for cell in cells:
                f.write(f"{cell[0]}")
                for idx in cell[1:]:
                    f.write(f" {idx}")
                f.write("\n")
            
            # 写入单元类型数据
            f.write(f"\nCELL_DATA {len(cells)}\n")
            f.write("SCALARS ElementType int 1\n")
            f.write("LOOKUP_TABLE default\n")
            # 简化：所有面片都标记为同一类型
            for _ in cells:
                f.write("1\n")
    
    @staticmethod
    def _write_surface_from_centerline_v2(
        filename: str,
        elements: Dict[int, Element],
        nodes: Dict[int, Node],
        displacements: np.ndarray,
        sections: Optional[Dict[int, Dict]] = None,
        element_section_map: Optional[Dict[int, int]] = None,
        ovalization_dofs: Optional[Dict[int, np.ndarray]] = None,
        material_properties: Optional[Dict] = None,
        element_metadata: Optional[Dict[int, Dict]] = None,
        ovalization_display_amplification: float = 1.0
    ):
        """
        基于中心线结果生成三维连续表面并写入VTK文件（V2版本）
        
        使用Parallel Transport Frame和节点级截面生成，确保：
        1. 截面连续，无翻转和卷曲
        2. 跨单元连接处网格拓扑连续
        3. 管道外观连续，无分段感
        """
        from solver.post.pipe_geometry_expander_v2 import PipeGeometryExpanderV2
        
        # 创建节点ID到DOF索引的映射
        sorted_node_ids = sorted(nodes.keys())
        
        # 确定每个节点的DOF数（PIPE288为6，ELBOW290为10）
        dof_per_node = len(displacements) // len(nodes) if len(displacements) > 0 and len(nodes) > 0 else 6
        
        node_id_to_dof_start = {
            node_id: idx * dof_per_node for idx, node_id in enumerate(sorted_node_ids)
        }
        
        # 如果没有提供截面参数，使用默认值
        if sections is None:
            sections = {1: {'diameter': 30.0, 'wall_thickness': 1.0}}
        
        # 展开中心线为三维连续表面（V2版本）
        points, cells = PipeGeometryExpanderV2.expand_pipe_centerline_to_surface(
            nodes, elements, displacements,
            node_id_to_dof_start,
            sections, element_section_map,
            num_circumferential=20,
            material_properties=material_properties,
            element_metadata=element_metadata,
            ovalization_display_amplification=ovalization_display_amplification
        )
        
        # 写入VTK文件
        with open(filename, 'w', encoding='utf-8') as f:
            # VTK文件头
            f.write("# vtk DataFile Version 2.0\n")
            f.write("Pipe Elements - Continuous Surface Geometry (V2)\n")
            f.write("ASCII\n")
            f.write("DATASET POLYDATA\n")
            f.write("\n")
            
            # 写入点
            f.write(f"POINTS {len(points)} double\n")
            for point in points:
                f.write(f"{point[0]:.6e} {point[1]:.6e} {point[2]:.6e}\n")
            
            # 写入多边形（四边形）
            total_cell_size = sum(1 + len(cell[1:]) for cell in cells)
            f.write(f"\nPOLYGONS {len(cells)} {total_cell_size}\n")
            for cell in cells:
                f.write(f"{cell[0]}")
                for idx in cell[1:]:
                    f.write(f" {idx}")
                f.write("\n")
            
            # 写入单元类型数据
            f.write(f"\nCELL_DATA {len(cells)}\n")
            f.write("SCALARS ElementType int 1\n")
            f.write("LOOKUP_TABLE default\n")
            # 简化：所有面片都标记为同一类型
            for _ in cells:
                f.write("1\n")

