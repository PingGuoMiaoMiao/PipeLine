#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
主程序入口 - 蠕变分析版本
支持高温蠕变-弹性-塑性耦合分析
"""

import sys
import numpy as np
from solver.mesh.cdb_parser import CDBParser
from solver.material.elastic import ElasticMaterial
from solver.material.plastic_biso import BISOPlasticMaterial
from solver.material.creep_strain_hardening import CreepStrainHardeningMaterial
from solver.element.pipe288_element_creep import PIPE288ElementCreep
from solver.solver.creep_time_integration import CreepTimeIntegrationSolver
from solver.post.vtk_writer import VTKWriter


def main():
    """主函数"""
    if len(sys.argv) < 2:
        print("用法: python -m solver.main_creep <input.cdb> [output_prefix] [total_time] [num_steps]")
        print("示例: python -m solver.main_creep examples/PIPE288_CREEP.cdb output_creep 10000 100")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_prefix = sys.argv[2] if len(sys.argv) > 2 else input_file.replace('.cdb', '_creep')
    total_time = float(sys.argv[3]) if len(sys.argv) > 3 else 10000.0
    num_steps = int(sys.argv[4]) if len(sys.argv) > 4 else 100
    
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
    
    print(f"\n开始蠕变分析...")
    print(f"  - 总时间: {total_time}")
    print(f"  - 时间步数: {num_steps}")
    print(f"  - 时间增量: {total_time/num_steps:.2f}")
    
    # 获取材料参数
    mat_id = 1
    if mat_id not in parser.materials:
        print(f"警告: 材料ID {mat_id} 未找到，使用默认值")
        material = ElasticMaterial(mat_id, E=200e9, nu=0.3, density=7800.0,
                                   alpha=1.2e-5, T_ref=25.0)
        yield_stress = 25e6
        tangent_modulus = 100e9
        # 默认蠕变参数（从算例说明中）
        C1 = 1.0e-14
        C2 = 2.0
        C3 = -1.0
        C4 = 0.0
    else:
        mat_dict = parser.materials[mat_id]
        if 'T_ref' not in mat_dict:
            mat_dict['T_ref'] = parser.reference_temp
        material = ElasticMaterial(mat_id, **mat_dict)
        
        yield_stress = mat_dict.get('yield_stress', 25e6)
        tangent_modulus = mat_dict.get('tangent_modulus', 100e9)
        
        # 蠕变参数
        C1 = mat_dict.get('C1', 1.0e-14)
        C2 = mat_dict.get('C2', 2.0)
        C3 = mat_dict.get('C3', -1.0)
        C4 = mat_dict.get('C4', 0.0)
    
    # 创建材料对象
    plastic_material = None
    if yield_stress > 0:
        plastic_material = BISOPlasticMaterial(
            E=material.E, nu=material.nu,
            yield_stress=yield_stress,
            tangent_modulus=tangent_modulus
        )
    
    creep_material = CreepStrainHardeningMaterial(C1, C2, C3, C4)
    
    print(f"弹性材料参数: E={material.E/1e9:.1f} GPa, nu={material.nu:.3f}")
    if plastic_material:
        print(f"弹塑性材料参数: σ_y={yield_stress/1e6:.1f} MPa")
    print(f"蠕变材料参数: C1={C1:.2e}, C2={C2:.1f}, C3={C3:.1f}, C4={C4:.1f}")
    
    # 获取截面参数
    sec_id = 1
    if sec_id not in parser.sections:
        print(f"警告: 截面ID {sec_id} 未找到，使用默认值")
        section = {'diameter': 30.0, 'wall_thickness': 1.0}
    else:
        section = parser.sections[sec_id]
    
    print(f"截面参数: D={section.get('diameter', 0):.1f} mm, "
          f"t={section.get('wall_thickness', 0):.1f} mm")
    
    # 温度（转换为绝对温度）
    temp_celsius = parser.temperature
    temp_kelvin = temp_celsius + 273.15
    print(f"温度: {temp_celsius:.1f}°C ({temp_kelvin:.1f} K)")
    
    # 创建单元对象
    pipe288_elems = {}
    for elem_id, elem in pipe288_elements.items():
        if len(elem.node_ids) >= 2:
            node1 = parser.nodes[elem.node_ids[0]]
            node2 = parser.nodes[elem.node_ids[1]]
            try:
                pipe288_elems[elem_id] = PIPE288ElementCreep(
                    elem_id, node1, node2, material, section,
                    plastic_material, creep_material)
            except Exception as e:
                print(f"警告: 单元 {elem_id} 创建失败: {e}")
    
    if len(pipe288_elems) == 0:
        print("错误: 无法创建单元")
        return
    
    # 创建时间积分求解器
    time_solver = CreepTimeIntegrationSolver(
        parser.nodes, pipe288_elems, parser.boundary_conditions,
        dof_per_node=6, max_iterations=20, tolerance=1e-6)
    
    # 组装载荷向量
    sorted_node_ids = sorted(parser.nodes.keys())
    node_id_to_dof_start = {
        node_id: idx * 6 for idx, node_id in enumerate(sorted_node_ids)
    }
    
    load_vector = np.zeros(time_solver.dof_count)
    for elem_id, elem in pipe288_elems.items():
        pressure = parser.internal_pressure.get(elem_id, 0.0)
        F_elem = elem.compute_load_vector(
            pressure,
            np.array(parser.gravity),
            temp_celsius,
            material.T_ref
        )
        
        dof1_start = node_id_to_dof_start[elem.node1.id]
        dof2_start = node_id_to_dof_start[elem.node2.id]
        
        load_vector[dof1_start:dof1_start+6] += F_elem[0:6]
        load_vector[dof2_start:dof2_start+6] += F_elem[6:12]
    
    # 时间步进分析
    time_increment = total_time / num_steps
    output_interval = max(1, num_steps // 10)  # 输出10个时刻
    
    print(f"\n开始时间步进分析...")
    
    for step in range(num_steps):
        current_time = step * time_increment
        next_time = current_time + time_increment
        
        print(f"\n时间步 {step+1}/{num_steps} (t={current_time:.1f} -> {next_time:.1f})")
        
        # 求解时间步
        displacements, info = time_solver.solve_time_step(
            load_vector,
            current_time,
            time_increment,
            temp_kelvin
        )
        
        if info['converged']:
            print(f"  收敛: {info['iterations']} 次迭代")
        else:
            print(f"  未收敛（残差: {info['residual_norms'][-1]:.6e}）")
        
        # 输出位移统计
        max_disp = np.max(np.abs(displacements[0::6]))
        max_disp_z = np.max(np.abs(displacements[2::6]))
        print(f"  最大位移: UX={max_disp*1000:.3f} mm, UZ={max_disp_z*1000:.3f} mm")
        
        # 统计蠕变应变
        if step % output_interval == 0 or step == num_steps - 1:
            max_creep_strain = 0.0
            for elem in pipe288_elems.values():
                creep_state = elem.get_creep_state()
                if creep_state['is_creep']:
                    max_creep_strain = max(max_creep_strain,
                                         creep_state['equivalent_creep_strain'])
            print(f"  最大等效蠕变应变: {max_creep_strain:.6e}")
            
            # 输出VTK文件
            output_file = f"{output_prefix}_t{next_time:.0f}.vtk"
            print(f"  输出VTK文件: {output_file}")
            
            # 更新节点位移
            for idx, node_id in enumerate(sorted_node_ids):
                if node_id in parser.nodes and node_id in node_id_to_dof_start:
                    node = parser.nodes[node_id]
                    dof_start = node_id_to_dof_start[node_id]
                    node.displacement = np.zeros(6)
                    if dof_start + 6 <= len(displacements):
                        node.displacement[0:6] = displacements[dof_start:dof_start+6]
            
            # 写入VTK文件（包含蠕变应变）
            VTKWriter.write_polylines_with_creep(
                output_file, parser.elements, parser.nodes, displacements,
                pipe288_elems)
    
    print(f"\n蠕变分析完成！")
    print(f"最终最大位移: UX={np.max(np.abs(displacements[0::6]))*1000:.3f} mm")
    
    # 最终蠕变应变统计
    max_creep_strain = 0.0
    for elem in pipe288_elems.values():
        creep_state = elem.get_creep_state()
        if creep_state['is_creep']:
            max_creep_strain = max(max_creep_strain,
                                 creep_state['equivalent_creep_strain'])
    print(f"最终最大等效蠕变应变: {max_creep_strain:.6e}")


if __name__ == '__main__':
    main()

