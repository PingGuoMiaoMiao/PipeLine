#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
主程序入口 - ELBOW290支持版本
"""

import sys
import numpy as np
from solver.mesh.cdb_parser import CDBParser
from solver.material.elastic import ElasticMaterial
from solver.element.pipe288_element import PIPE288Element
try:
    from solver.element.elbow290_element import ELBOW290Element
except ImportError:
    ELBOW290Element = None
from solver.solver.nonlinear_static import NonlinearStaticSolver
from solver.post.vtk_writer import VTKWriter


def main():
    """主函数"""
    if len(sys.argv) < 2:
        print("用法: python -m solver.main_elbow290 <input.cdb> [output.vtk]")
        print("示例: python -m solver.main_elbow290 examples/ELBOW290_PLAST.cdb output.vtk")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else input_file.replace('.cdb', '.vtk')
    
    print(f"解析CDB文件: {input_file}")
    parser = CDBParser(input_file)
    parser.parse()
    
    print(f"解析完成:")
    print(f"  - 节点数: {len(parser.nodes)}")
    print(f"  - 单元数: {len(parser.elements)}")
    print(f"  - 元素类型映射: {parser.element_type_map}")
    
    if len(parser.elements) == 0:
        print("警告: 未找到任何单元！")
        return
    
    # 检查单元类型
    pipe288_elements = {
        eid: elem for eid, elem in parser.elements.items() 
        if elem.type == 288
    }
    elbow290_elements = {
        eid: elem for eid, elem in parser.elements.items() 
        if elem.type == 290
    }
    
    if len(pipe288_elements) == 0 and len(elbow290_elements) == 0:
        print("警告: 未找到支持的单元类型，仅输出几何信息")
        VTKWriter.write_polylines(output_file, parser.elements, parser.nodes)
        return
    
    print(f"\n开始有限元分析...")
    print(f"  - PIPE288单元数: {len(pipe288_elements)}")
    print(f"  - ELBOW290单元数: {len(elbow290_elements)}")
    
    # 获取材料和截面参数
    mat_id = 1
    if mat_id not in parser.materials:
        print(f"警告: 材料ID {mat_id} 未找到，使用默认值")
        material = ElasticMaterial(mat_id, E=200e9, nu=0.3, density=7800.0, 
                                  alpha=1.2e-5, T_ref=25.0)
    else:
        mat_dict = parser.materials[mat_id]
        if 'T_ref' not in mat_dict:
            mat_dict['T_ref'] = parser.reference_temp
        material = ElasticMaterial(mat_id, **mat_dict)
    
    sec_id = 1
    if sec_id not in parser.sections:
        print(f"警告: 截面ID {sec_id} 未找到，使用默认值")
        if len(elbow290_elements) > 0:
            section = {'diameter': 90.0, 'wall_thickness': 2.0}  # ELBOW290默认值
        else:
            section = {'diameter': 30.0, 'wall_thickness': 1.0}  # PIPE288默认值
    else:
        section = parser.sections[sec_id]
    
    print(f"材料参数: E={material.E/1e9:.1f} GPa, nu={material.nu:.3f}")
    print(f"截面参数: D={section.get('diameter', 0):.1f} mm, "
          f"t={section.get('wall_thickness', 0):.1f} mm")
    print(f"载荷: 温度={parser.temperature:.1f}°C, 重力={parser.gravity}, "
          f"内压单元数={len(parser.internal_pressure)}")
    
    # 创建单元对象
    pipe288_elems = {}
    elbow290_elems = {}
    
    # 创建PIPE288单元
    for elem_id, elem in pipe288_elements.items():
        if len(elem.node_ids) >= 2:
            node1 = parser.nodes[elem.node_ids[0]]
            node2 = parser.nodes[elem.node_ids[1]]
            try:
                pipe288_elems[elem_id] = PIPE288Element(
                    elem_id, node1, node2, material, section)
            except Exception as e:
                print(f"警告: PIPE288单元 {elem_id} 创建失败: {e}")
    
    # 创建ELBOW290单元
    if ELBOW290Element is not None:
        for elem_id, elem in elbow290_elements.items():
            if len(elem.node_ids) >= 3:
                node1 = parser.nodes[elem.node_ids[0]]
                node2 = parser.nodes[elem.node_ids[1]]
                node3 = parser.nodes[elem.node_ids[2]]
                try:
                    elbow290_elems[elem_id] = ELBOW290Element(
                        elem_id, node1, node2, node3, material, section)
                except Exception as e:
                    print(f"警告: ELBOW290单元 {elem_id} 创建失败: {e}")
    else:
        print("警告: ELBOW290Element未导入，跳过ELBOW290单元")
    
    # 确定DOF数（优先使用ELBOW290的10 DOF/节点）
    if len(elbow290_elems) > 0:
        dof_per_node = 10
        all_elements = {**pipe288_elems, **elbow290_elems}
        print(f"使用ELBOW290模式（10 DOF/节点）")
    else:
        dof_per_node = 6
        all_elements = pipe288_elems
        print(f"使用PIPE288模式（6 DOF/节点）")
    
    if len(all_elements) == 0:
        print("错误: 无法创建任何单元")
        return
    
    # 创建求解器
    solver = NonlinearStaticSolver(
        parser.nodes, all_elements, parser.boundary_conditions, dof_per_node)
    
    # 组装载荷向量
    load_vector = np.zeros(solver.dof_count)
    sorted_node_ids = sorted(parser.nodes.keys())
    node_id_to_dof_start = {
        node_id: idx * dof_per_node for idx, node_id in enumerate(sorted_node_ids)
    }
    
    # PIPE288单元载荷
    for elem_id, elem in pipe288_elems.items():
        pressure = parser.internal_pressure.get(elem_id, 0.0)
        F_elem = elem.compute_load_vector(
            pressure,
            np.array(parser.gravity),
            parser.temperature,
            material.T_ref
        )
        
        dof1_start = node_id_to_dof_start[elem.node1.id]
        dof2_start = node_id_to_dof_start[elem.node2.id]
        
        load_vector[dof1_start:dof1_start+6] += F_elem[0:6]
        load_vector[dof2_start:dof2_start+6] += F_elem[6:12]
    
    # ELBOW290单元载荷
    for elem_id, elem in elbow290_elems.items():
        pressure = parser.internal_pressure.get(elem_id, 0.0)
        F_elem = elem.compute_load_vector(
            pressure,
            np.array(parser.gravity),
            parser.temperature,
            material.T_ref
        )
        
        dof1_start = node_id_to_dof_start[elem.node1.id]
        dof2_start = node_id_to_dof_start[elem.node2.id]
        dof3_start = node_id_to_dof_start[elem.node3.id]
        
        load_vector[dof1_start:dof1_start+10] += F_elem[0:10]
        load_vector[dof2_start:dof2_start+10] += F_elem[10:20]
        load_vector[dof3_start:dof3_start+10] += F_elem[20:30]
    
    print(f"求解线性方程组 (DOF数: {solver.dof_count})...")
    
    # 求解
    displacements = solver.solve(load_vector)
    
    # 更新节点位移
    for idx, node_id in enumerate(sorted_node_ids):
        if node_id in parser.nodes and node_id in node_id_to_dof_start:
            node = parser.nodes[node_id]
            dof_start = node_id_to_dof_start[node_id]
            # 只保存前6个DOF（标准DOF）
            node.displacement = np.zeros(6)
            if dof_start + 6 <= len(displacements):
                node.displacement[0:6] = displacements[dof_start:dof_start+6]
    
    # 输出位移统计
    max_disp = np.max(np.abs(displacements[0::dof_per_node]))  # UX
    max_disp_y = np.max(np.abs(displacements[1::dof_per_node]))  # UY
    max_disp_z = np.max(np.abs(displacements[2::dof_per_node]))  # UZ
    print(f"最大位移: UX={max_disp*1000:.3f} mm, "
          f"UY={max_disp_y*1000:.3f} mm, UZ={max_disp_z*1000:.3f} mm")
    
    # 输出椭圆化和翘曲统计（如果有ELBOW290）
    if len(elbow290_elems) > 0:
        # 椭圆化DOF（索引6-7）
        oval_cos = displacements[6::dof_per_node]
        oval_sin = displacements[7::dof_per_node]
        max_oval = np.max(np.sqrt(oval_cos**2 + oval_sin**2))
        # 翘曲DOF（索引8-9）
        warp_cos = displacements[8::dof_per_node]
        warp_sin = displacements[9::dof_per_node]
        max_warp = np.max(np.sqrt(warp_cos**2 + warp_sin**2))
        print(f"最大椭圆化: {max_oval*1000:.6f} mm")
        print(f"最大翘曲: {max_warp*1000:.6f} mm")
    
    # 写入VTK文件
    print(f"\n写入VTK文件: {output_file}")
    
    # 准备ELBOW290的元数据（用于椭圆化显示）
    element_metadata = {}
    material_props = {}
    if len(elbow290_elems) > 0:
        # 收集ELBOW290单元的元数据
        for elem_id, elem in elbow290_elems.items():
            element_metadata[elem_id] = {
                'R_curvature': elem.R_curvature,
                'L': elem.L,
                'R_mean': elem.R_mean
            }
        
        # 材料属性
        if len(elbow290_elems) > 0:
            first_elem = next(iter(elbow290_elems.values()))
            material_props = {
                'E': first_elem.material.E,
                'nu': first_elem.material.nu
            }
    
    VTKWriter.write_polylines(
        output_file, parser.elements, parser.nodes, displacements,
        sections=parser.sections,
        output_surface=True,
        material_properties=material_props if material_props else None,
        element_metadata=element_metadata if element_metadata else None,
        ovalization_display_amplification=1.0
    )
    
    print(f"完成！输出文件: {output_file}")
    print(f"\n可以使用ParaView打开查看: paraview {output_file}")


if __name__ == '__main__':
    main()

