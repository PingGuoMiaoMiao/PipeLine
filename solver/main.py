#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
主程序入口 - 基于管单元的高温管道蠕变高效求解
"""

import sys
import numpy as np
from solver.mesh.cdb_parser import CDBParser
from solver.material.elastic import ElasticMaterial
from solver.element.pipe288_element import PIPE288Element
from solver.solver.nonlinear_static import NonlinearStaticSolver
from solver.post.vtk_writer import VTKWriter


def main():
    """主函数"""
    if len(sys.argv) < 2:
        print("用法: python -m solver.main <input.cdb> [output.vtk]")
        print("示例: python -m solver.main examples/PIPE288_PLAST.cdb output.vtk")
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
    
    # 检查是否有PIPE288单元
    pipe288_elements = {
        eid: elem for eid, elem in parser.elements.items() 
        if elem.type == 288
    }
    if len(pipe288_elements) == 0:
        print("警告: 未找到PIPE288单元，仅输出几何信息")
        VTKWriter.write_polylines(output_file, parser.elements, parser.nodes)
        return
    
    print(f"\n开始有限元分析...")
    
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
        section = {'diameter': 30.0, 'wall_thickness': 1.0}
    else:
        section = parser.sections[sec_id]
    
    print(f"材料参数: E={material.E/1e9:.1f} GPa, nu={material.nu:.3f}")
    print(f"截面参数: D={section.get('diameter', 0):.1f} mm, "
          f"t={section.get('wall_thickness', 0):.1f} mm")
    print(f"载荷: 温度={parser.temperature:.1f}°C, 重力={parser.gravity}, "
          f"内压单元数={len(parser.internal_pressure)}")
    
    # 创建PIPE288单元对象
    pipe288_elems = {}
    for elem_id, elem in pipe288_elements.items():
        if len(elem.node_ids) >= 2:
            node1 = parser.nodes[elem.node_ids[0]]
            node2 = parser.nodes[elem.node_ids[1]]
            try:
                pipe288_elems[elem_id] = PIPE288Element(
                    elem_id, node1, node2, material, section)
            except Exception as e:
                print(f"警告: 单元 {elem_id} 创建失败: {e}")
    
    if len(pipe288_elems) == 0:
        print("错误: 无法创建PIPE288单元")
        return
    
    # 创建求解器
    solver = NonlinearStaticSolver(
        parser.nodes, pipe288_elems, parser.boundary_conditions)
    
    # 组装载荷向量
    load_vector = np.zeros(solver.dof_count)
    sorted_node_ids = sorted(parser.nodes.keys())
    node_id_to_dof_start = {
        node_id: idx * 6 for idx, node_id in enumerate(sorted_node_ids)
    }
    
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
    
    print(f"求解线性方程组 (DOF数: {solver.dof_count})...")
    
    # 求解
    displacements = solver.solve(load_vector)
    
    # 更新节点位移
    for idx, node_id in enumerate(sorted_node_ids):
        if node_id in parser.nodes and node_id in node_id_to_dof_start:
            node = parser.nodes[node_id]
            dof_start = node_id_to_dof_start[node_id]
            node.displacement = displacements[dof_start:dof_start+6]
    
    # 输出位移统计
    max_disp = np.max(np.abs(displacements[0::6]))  # UX
    max_disp_y = np.max(np.abs(displacements[1::6]))  # UY
    max_disp_z = np.max(np.abs(displacements[2::6]))  # UZ
    print(f"最大位移: UX={max_disp*1000:.3f} mm, "
          f"UY={max_disp_y*1000:.3f} mm, UZ={max_disp_z*1000:.3f} mm")
    
    # 写入VTK文件（使用新的几何展开模块）
    print(f"\n写入VTK文件: {output_file}")
    VTKWriter.write_polylines(
        output_file, parser.elements, parser.nodes, displacements,
        sections=parser.sections,
        element_section_map=None,  # 简化：所有单元使用同一截面
        ovalization_dofs=None,  # PIPE288无椭圆化自由度
        output_surface=True  # 输出为表面几何
    )
    print(f"完成！输出文件: {output_file}")
    print(f"\n可以使用ParaView打开查看: paraview {output_file}")


if __name__ == '__main__':
    main()

