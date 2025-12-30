"""
性能测试脚本 - 对比优化前后的蠕变求解器性能
"""

import sys
import time
import numpy as np
from solver.mesh.cdb_parser import CDBParser
from solver.material.elastic import ElasticMaterial
from solver.element.pipe288_element_creep import PIPE288ElementCreep
from solver.solver.creep_time_integration import CreepTimeIntegrationSolver
from solver.solver.optimized_creep_solver_v2 import OptimizedCreepSolverV2


def run_original_solver(parser, material, elements_dict, boundary_conditions):
    """运行原始求解器"""
    print("\n" + "="*60)
    print("运行原始求解器（CreepTimeIntegrationSolver）")
    print("="*60)
    
    start_time = time.time()
    
    # 创建求解器
    solver = CreepTimeIntegrationSolver(
        parser.nodes, elements_dict, boundary_conditions, dof_per_node=6
    )
    
    # 求解
    result = solver.solve_creep_analysis(
        total_time=1000.0,
        num_steps=100,
        internal_pressure=parser.internal_pressure,
        gravity=parser.gravity,
        temperature=parser.temperature,
        ref_temperature=parser.reference_temp
    )
    
    elapsed_time = time.time() - start_time
    
    print(f"原始求解器完成")
    print(f"  总时间: {elapsed_time:.2f} 秒")
    print(f"  时间步数: {len(result.get('time_history', []))}")
    
    return result, elapsed_time


def run_optimized_solver(parser, material, elements_dict, boundary_conditions):
    """运行优化求解器"""
    print("\n" + "="*60)
    print("运行优化求解器（OptimizedCreepSolverV2）")
    print("="*60)
    
    start_time = time.time()
    
    # 创建优化求解器
    solver = OptimizedCreepSolverV2(
        parser.nodes, elements_dict, boundary_conditions, dof_per_node=6,
        use_sparse=True, use_adaptive_step=True
    )
    
    # 求解
    result = solver.solve_creep_analysis(
        total_time=1000.0,
        initial_time_step=10.0,
        internal_pressure=parser.internal_pressure,
        gravity=parser.gravity,
        temperature=parser.temperature,
        ref_temperature=parser.reference_temp
    )
    
    elapsed_time = time.time() - start_time
    
    print(f"优化求解器完成")
    print(f"  总时间: {elapsed_time:.2f} 秒")
    print(f"  时间步数: {result.get('stats', {}).get('num_time_steps', 0)}")
    print(f"  总迭代数: {result.get('stats', {}).get('num_iterations', 0)}")
    print(f"  装配时间: {result.get('stats', {}).get('assembly_time', 0):.2f} 秒")
    print(f"  求解时间: {result.get('stats', {}).get('solve_time', 0):.2f} 秒")
    
    return result, elapsed_time


def compare_results(result_original, result_optimized):
    """对比结果精度"""
    print("\n" + "="*60)
    print("结果精度对比")
    print("="*60)
    
    # 获取最终位移
    disp_original = result_original.get('displacements', np.array([]))
    disp_optimized = result_optimized.get('displacements', np.array([]))
    
    if len(disp_original) == 0 or len(disp_optimized) == 0:
        print("警告：无法对比结果（位移数组为空）")
        return
    
    # 确保长度一致
    min_len = min(len(disp_original), len(disp_optimized))
    disp_original = disp_original[:min_len]
    disp_optimized = disp_optimized[:min_len]
    
    # 计算误差
    error = np.abs(disp_original - disp_optimized)
    max_error = np.max(error)
    relative_error = max_error / (np.max(np.abs(disp_original)) + 1e-10)
    
    print(f"最大绝对误差: {max_error:.2e} m")
    print(f"最大相对误差: {relative_error*100:.4f}%")
    
    # 检查精度
    if relative_error < 0.01:  # 1%的相对误差
        print("✓ 精度检查通过（相对误差 < 1%）")
    else:
        print("⚠ 警告：相对误差较大，可能存在精度问题")


def main():
    """主函数"""
    if len(sys.argv) < 2:
        print("用法: python performance_test_creep.py <input.cdb>")
        print("示例: python performance_test_creep.py examples/PIPE288_CREEP.cdb")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    print("="*60)
    print("蠕变求解器性能对比测试")
    print("="*60)
    print(f"输入文件: {input_file}")
    
    # 解析CDB文件
    print("\n解析CDB文件...")
    parser = CDBParser(input_file)
    parser.parse()
    
    print(f"  节点数: {len(parser.nodes)}")
    print(f"  单元数: {len(parser.elements)}")
    
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
        section = {'diameter': 30.0, 'wall_thickness': 1.0}
    else:
        section = parser.sections[sec_id]
    
    # 创建单元对象
    print("\n创建单元对象...")
    elements_dict = {}
    for elem_id, elem in parser.elements.items():
        if elem.type == 288 and len(elem.node_ids) >= 2:
            node1 = parser.nodes[elem.node_ids[0]]
            node2 = parser.nodes[elem.node_ids[1]]
            try:
                elements_dict[elem_id] = PIPE288ElementCreep(
                    elem_id, node1, node2, material, section
                )
            except Exception as e:
                print(f"警告: 单元 {elem_id} 创建失败: {e}")
    
    print(f"  成功创建 {len(elements_dict)} 个单元")
    
    # 运行原始求解器
    try:
        result_original, time_original = run_original_solver(
            parser, material, elements_dict, parser.boundary_conditions
        )
    except Exception as e:
        print(f"原始求解器运行失败: {e}")
        import traceback
        traceback.print_exc()
        result_original, time_original = None, None
    
    # 运行优化求解器
    try:
        result_optimized, time_optimized = run_optimized_solver(
            parser, material, elements_dict, parser.boundary_conditions
        )
    except Exception as e:
        print(f"优化求解器运行失败: {e}")
        import traceback
        traceback.print_exc()
        result_optimized, time_optimized = None, None
    
    # 性能对比
    print("\n" + "="*60)
    print("性能对比总结")
    print("="*60)
    
    if time_original is not None and time_optimized is not None:
        speedup = time_original / time_optimized
        print(f"原始求解器时间: {time_original:.2f} 秒")
        print(f"优化求解器时间: {time_optimized:.2f} 秒")
        print(f"加速比: {speedup:.2f}x")
        
        if speedup > 1.0:
            print(f"✓ 优化成功，性能提升 {speedup:.2f} 倍")
        else:
            print(f"⚠ 优化后性能未提升（可能需要调整参数）")
    
    # 精度对比
    if result_original is not None and result_optimized is not None:
        compare_results(result_original, result_optimized)
    
    print("\n测试完成！")


if __name__ == '__main__':
    main()

