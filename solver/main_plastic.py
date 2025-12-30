#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
主程序入口 - 弹塑性分析版本
支持BISO双线性等向强化模型和Newton-Raphson迭代
"""

import sys
import numpy as np
from solver.mesh.cdb_parser import CDBParser
from solver.material.elastic import ElasticMaterial
from solver.material.plastic_biso import BISOPlasticMaterial
from solver.element.pipe288_element_plastic import PIPE288ElementPlastic
from solver.solver.newton_raphson import NewtonRaphsonSolver
from solver.post.vtk_writer import VTKWriter


def main():
    """主函数"""
    if len(sys.argv) < 2:
        print("用法: python -m solver.main_plastic <input.cdb> [output.vtk] [num_steps]")
        print("示例: python -m solver.main_plastic examples/PIPE288_PLAST.cdb output_plastic.vtk 10")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else input_file.replace('.cdb', '_plastic.vtk')
    num_steps = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    
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
        print("警告: 未找到PIPE288单元")
        return
    
    print(f"\n开始弹塑性有限元分析...")
    print(f"  - 增量步数: {num_steps}")
    
    # 获取材料参数
    mat_id = 1
    if mat_id not in parser.materials:
        print(f"警告: 材料ID {mat_id} 未找到，使用默认值")
        material = ElasticMaterial(mat_id, E=200e9, nu=0.3, density=7800.0,
                                   alpha=1.2e-5, T_ref=25.0)
        # 默认弹塑性参数
        yield_stress = 25e6  # 25 MPa
        tangent_modulus = 100e9  # 100 GPa
    else:
        mat_dict = parser.materials[mat_id]
        if 'T_ref' not in mat_dict:
            mat_dict['T_ref'] = parser.reference_temp
        material = ElasticMaterial(mat_id, **mat_dict)
        
        # 从材料字典获取弹塑性参数
        yield_stress = mat_dict.get('yield_stress', 25e6)
        tangent_modulus = mat_dict.get('tangent_modulus', 100e9)
    
    # 创建弹塑性材料
    plastic_material = BISOPlasticMaterial(
        E=material.E,
        nu=material.nu,
        yield_stress=yield_stress,
        tangent_modulus=tangent_modulus
    )
    
    print(f"弹性材料参数: E={material.E/1e9:.1f} GPa, nu={material.nu:.3f}")
    print(f"弹塑性材料参数: σ_y={yield_stress/1e6:.1f} MPa, E_t={tangent_modulus/1e9:.1f} GPa")
    
    # 获取截面参数
    sec_id = 1
    if sec_id not in parser.sections:
        print(f"警告: 截面ID {sec_id} 未找到，使用默认值")
        section = {'diameter': 30.0, 'wall_thickness': 1.0}
    else:
        section = parser.sections[sec_id]
    
    print(f"截面参数: D={section.get('diameter', 0):.1f} mm, "
          f"t={section.get('wall_thickness', 0):.1f} mm")
    print(f"载荷: 温度={parser.temperature:.1f}°C, 重力={parser.gravity}, "
          f"内压单元数={len(parser.internal_pressure)}")
    
    # 创建PIPE288单元对象（支持弹塑性）
    pipe288_elems = {}
    for elem_id, elem in pipe288_elements.items():
        if len(elem.node_ids) >= 2:
            node1 = parser.nodes[elem.node_ids[0]]
            node2 = parser.nodes[elem.node_ids[1]]
            try:
                pipe288_elems[elem_id] = PIPE288ElementPlastic(
                    elem_id, node1, node2, material, section, plastic_material)
            except Exception as e:
                print(f"警告: 单元 {elem_id} 创建失败: {e}")
    
    if len(pipe288_elems) == 0:
        print("错误: 无法创建PIPE288单元")
        return
    
    # 创建Newton-Raphson求解器
    solver = NewtonRaphsonSolver(
        parser.nodes, pipe288_elems, parser.boundary_conditions,
        dof_per_node=6, max_iterations=20, tolerance=1e-6)
    
    # 组装初始载荷向量
    sorted_node_ids = sorted(parser.nodes.keys())
    node_id_to_dof_start = {
        node_id: idx * 6 for idx, node_id in enumerate(sorted_node_ids)
    }
    
    load_vector = np.zeros(solver.dof_count)
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
    
    # 增量加载
    print(f"\n开始增量加载分析...")
    total_info = []
    
    for step in range(num_steps):
        load_factor = (step + 1) / num_steps
        print(f"\n增量步 {step+1}/{num_steps} (载荷因子: {load_factor:.3f})")
        
        # Newton-Raphson迭代
        displacements, info = solver.solve(
            load_vector, incremental=True, load_factor=load_factor)
        
        total_info.append(info)
        
        if info['converged']:
            print(f"  收敛: {info['iterations']} 次迭代")
        else:
            print(f"  未收敛（残差: {info['residual_norms'][-1]:.6e}）")
    
    # 输出最终位移统计
    final_displacements = solver.displacements
    max_disp = np.max(np.abs(final_displacements[0::6]))
    max_disp_y = np.max(np.abs(final_displacements[1::6]))
    max_disp_z = np.max(np.abs(final_displacements[2::6]))
    print(f"\n最终最大位移: UX={max_disp*1000:.3f} mm, "
          f"UY={max_disp_y*1000:.3f} mm, UZ={max_disp_z*1000:.3f} mm")
    
    # 统计塑性单元
    plastic_units = []
    for elem_id, elem in pipe288_elems.items():
        plastic_state = elem.get_plastic_state()
        if plastic_state['is_plastic']:
            plastic_units.append({
                'elem_id': elem_id,
                'equivalent_plastic_strain': plastic_state['equivalent_plastic_strain']
            })
    
    print(f"\n塑性单元数: {len(plastic_units)}/{len(pipe288_elems)}")
    if len(plastic_units) > 0:
        max_plastic_strain = max(p['equivalent_plastic_strain'] for p in plastic_units)
        print(f"最大等效塑性应变: {max_plastic_strain:.6e}")
    
    # 更新节点位移
    for idx, node_id in enumerate(sorted_node_ids):
        if node_id in parser.nodes and node_id in node_id_to_dof_start:
            node = parser.nodes[node_id]
            dof_start = node_id_to_dof_start[node_id]
            node.displacement = np.zeros(6)
            if dof_start + 6 <= len(final_displacements):
                node.displacement[0:6] = final_displacements[dof_start:dof_start+6]
    
    # 写入VTK文件（包含塑性区信息）
    print(f"\n写入VTK文件: {output_file}")
    VTKWriter.write_polylines_with_plasticity(
        output_file, parser.elements, parser.nodes, final_displacements,
        pipe288_elems,
        sections=parser.sections,
        output_surface=True
    )
    print(f"完成！输出文件: {output_file}")
    print(f"\n可以使用ParaView打开查看: paraview {output_file}")


if __name__ == '__main__':
    main()

