"""
PIPE288 Element with Creep - 支持蠕变的PIPE288单元
支持弹性-塑性-蠕变耦合
"""

import math
import numpy as np
from typing import Dict, Optional
from typing import Dict as DictType
from solver.mesh.node import Node
from solver.material.elastic import ElasticMaterial
from solver.material.plastic_biso import BISOPlasticMaterial
from solver.material.creep_strain_hardening import CreepStrainHardeningMaterial
from solver.integration.integration_point import IntegrationPoint


class PIPE288ElementCreep:
    """PIPE288单元类 - 支持蠕变分析"""
    
    def __init__(self, elem_id: int, node1: Node, node2: Node,
                 material: ElasticMaterial, section: Dict,
                 plastic_material: Optional[BISOPlasticMaterial] = None,
                 creep_material: Optional[CreepStrainHardeningMaterial] = None):
        """
        初始化PIPE288单元（支持蠕变）
        
        参数:
            elem_id: 单元ID
            node1: 节点1
            node2: 节点2
            material: 弹性材料对象
            section: 截面参数字典
            plastic_material: 弹塑性材料对象（可选）
            creep_material: 蠕变材料对象（可选）
        """
        self.elem_id = elem_id
        self.node1 = node1
        self.node2 = node2
        self.material = material
        self.plastic_material = plastic_material
        self.creep_material = creep_material
        self.is_plastic = plastic_material is not None
        self.is_creep = creep_material is not None
        
        # 截面参数
        D = section.get('diameter', 30.0)
        t = section.get('wall_thickness', 1.0)
        
        self.Ro = D / 2.0 / 1000.0  # 外半径 (m) = D(mm) / 2 / 1000
        self.Ri = (D - 2 * t) / 2.0 / 1000.0  # 内半径 (m) = (D-2t)(mm) / 2 / 1000
        
        # 计算截面几何参数
        self.A = math.pi * (self.Ro**2 - self.Ri**2)
        self.I = math.pi / 4.0 * (self.Ro**4 - self.Ri**4)
        self.J = math.pi / 2.0 * (self.Ro**4 - self.Ri**4)
        
        # 计算单元长度和方向
        dx = (node2.x - node1.x) / 1000.0
        dy = (node2.y - node1.y) / 1000.0
        dz = (node2.z - node1.z) / 1000.0
        self.L = math.sqrt(dx**2 + dy**2 + dz**2)
        
        if self.L > 1e-10:
            self.ex = np.array([dx, dy, dz]) / self.L
            if abs(self.ex[2]) < 0.9:
                ez_global = np.array([0, 0, 1])
                ey = np.cross(self.ex, ez_global)
            else:
                ey_global = np.array([0, 1, 0])
                ey = np.cross(self.ex, ey_global)
            ey_norm = np.linalg.norm(ey)
            if ey_norm > 1e-10:
                ey = ey / ey_norm
            else:
                ey = np.array([0, 1, 0])
            self.ez = np.cross(self.ex, ey)
            ez_norm = np.linalg.norm(self.ez)
            if ez_norm > 1e-10:
                self.ez = self.ez / ez_norm
            else:
                self.ez = np.array([0, 0, 1])
            self.ey = np.cross(self.ez, self.ex)
        else:
            raise ValueError(f"单元{elem_id}长度为零")
        
        # 积分点
        self.integration_points = []
        if self.is_creep or self.is_plastic:
            num_gauss = 2
            gauss_points = [-1.0/np.sqrt(3), 1.0/np.sqrt(3)]
            gauss_weights = [1.0, 1.0]
            for i in range(num_gauss):
                ip = IntegrationPoint()
                ip.xi = gauss_points[i]
                ip.weight = gauss_weights[i]
                self.integration_points.append(ip)
    
    def compute_stiffness_matrix(self) -> np.ndarray:
        """计算单元刚度矩阵"""
        return self._compute_elastic_stiffness_matrix()
    
    def _compute_elastic_stiffness_matrix(self) -> np.ndarray:
        """计算弹性刚度矩阵（12×12）"""
        K_local = np.zeros((12, 12))
        L = self.L
        E = self.material.E
        G = self.material.G
        
        # 轴向刚度
        k_axial = E * self.A / L
        K_local[0, 0] = k_axial
        K_local[0, 6] = -k_axial
        K_local[6, 0] = -k_axial
        K_local[6, 6] = k_axial
        
        # 扭转刚度
        k_torsion = G * self.J / L
        K_local[3, 3] = k_torsion
        K_local[3, 9] = -k_torsion
        K_local[9, 3] = -k_torsion
        K_local[9, 9] = k_torsion
        
        # 弯曲刚度
        EI = E * self.I
        k_bend = EI / L**3
        
        # Y方向弯曲
        K_local[1, 1] = 12 * k_bend
        K_local[1, 5] = 6 * L * k_bend
        K_local[1, 7] = -12 * k_bend
        K_local[1, 11] = 6 * L * k_bend
        K_local[5, 1] = 6 * L * k_bend
        K_local[5, 5] = 4 * L**2 * k_bend
        K_local[5, 7] = -6 * L * k_bend
        K_local[5, 11] = 2 * L**2 * k_bend
        K_local[7, 1] = -12 * k_bend
        K_local[7, 5] = -6 * L * k_bend
        K_local[7, 7] = 12 * k_bend
        K_local[7, 11] = -6 * L * k_bend
        K_local[11, 1] = 6 * L * k_bend
        K_local[11, 5] = 2 * L**2 * k_bend
        K_local[11, 7] = -6 * L * k_bend
        K_local[11, 11] = 4 * L**2 * k_bend
        
        # Z方向弯曲
        K_local[2, 2] = 12 * k_bend
        K_local[2, 4] = -6 * L * k_bend
        K_local[2, 8] = -12 * k_bend
        K_local[2, 10] = -6 * L * k_bend
        K_local[4, 2] = -6 * L * k_bend
        K_local[4, 4] = 4 * L**2 * k_bend
        K_local[4, 8] = 6 * L * k_bend
        K_local[4, 10] = 2 * L**2 * k_bend
        K_local[8, 2] = -12 * k_bend
        K_local[8, 4] = 6 * L * k_bend
        K_local[8, 8] = 12 * k_bend
        K_local[8, 10] = 6 * L * k_bend
        K_local[10, 2] = -6 * L * k_bend
        K_local[10, 4] = 2 * L**2 * k_bend
        K_local[10, 8] = 6 * L * k_bend
        K_local[10, 10] = 4 * L**2 * k_bend
        
        # 坐标变换
        R = np.array([
            [self.ex[0], self.ex[1], self.ex[2]],
            [self.ey[0], self.ey[1], self.ey[2]],
            [self.ez[0], self.ez[1], self.ez[2]]
        ])
        
        T = np.zeros((12, 12))
        T[0:3, 0:3] = R
        T[3:6, 3:6] = R
        T[6:9, 6:9] = R
        T[9:12, 9:12] = R
        
        K_global = T.T @ K_local @ T
        return K_global
    
    def compute_tangent_stiffness(self, node_displacements: np.ndarray) -> np.ndarray:
        """计算切线刚度矩阵"""
        K_elastic = self._compute_elastic_stiffness_matrix()
        
        # 简化：如果有蠕变或塑性，降低刚度
        if self.is_creep or self.is_plastic:
            reduction_factor = 1.0
            for ip in self.integration_points:
                # 简化处理
                if hasattr(ip, 'is_plastic') and ip.is_plastic:
                    reduction_factor = min(reduction_factor, 0.1)
            
            return K_elastic * reduction_factor
        
        return K_elastic
    
    def compute_internal_force(self, node_displacements: np.ndarray) -> np.ndarray:
        """计算单元内部力向量"""
        # 简化：使用切线刚度
        K_tangent = self.compute_tangent_stiffness(node_displacements)
        return K_tangent @ node_displacements
    
    def update_creep_strain(self, node_id_to_dof_start: DictType[int, int],
                           global_displacements: np.ndarray,
                           time: float,
                           time_increment: float,
                           temperature: float):
        """
        更新蠕变应变
        
        参数:
            node_id_to_dof_start: 节点ID到DOF起始索引的映射
            global_displacements: 全局位移向量
            time: 当前时间
            time_increment: 时间增量
            temperature: 温度（绝对温度，K）
        """
        if not self.is_creep:
            return
        
        # 从节点位移计算单元的局部应变和应力
        # 获取节点在全局位移向量中的索引
        dof_start_1 = node_id_to_dof_start.get(self.node1.id, -1)
        dof_start_2 = node_id_to_dof_start.get(self.node2.id, -1)
        
        if dof_start_1 < 0 or dof_start_2 < 0:
            return  # 节点未找到
        
        # 提取节点位移（前3个DOF是平移）
        u1 = global_displacements[dof_start_1:dof_start_1+3] if dof_start_1+3 <= len(global_displacements) else np.zeros(3)
        u2 = global_displacements[dof_start_2:dof_start_2+3] if dof_start_2+3 <= len(global_displacements) else np.zeros(3)
        
        # 计算单元轴向应变
        # 单元在局部坐标系中的轴向方向是ex
        # 局部轴向位移差
        du_local = np.dot(u2 - u1, self.ex)
        # 轴向应变
        epsilon_axial = du_local / self.L if self.L > 1e-10 else 0.0
        
        # 计算轴向应力（弹性部分）
        E = self.material.E
        axial_stress = E * epsilon_axial
        
        # 构建应力向量（Voigt记号：σ11, σ22, σ33, σ12, σ13, σ23）
        # 对于管道，主要应力是轴向应力σ11
        stress = np.zeros(6)
        stress[0] = axial_stress  # 轴向应力
        # 其他分量假设为0（简化处理）
        
        # 计算Mises等效应力（对于单轴应力，就是轴向应力的绝对值）
        mises_stress = np.abs(axial_stress)
        
        for ip in self.integration_points:
            # 如果应力太小，跳过蠕变计算
            if mises_stress < 1e-3:  # 1 kPa
                continue
            
            # 计算蠕变应变增量
            delta_equiv_creep = self.creep_material.compute_creep_strain_increment(
                stress,
                ip.equivalent_creep_strain,
                time,
                time_increment,
                temperature,
                method='explicit'
            )
            
            # 更新等效蠕变应变
            ip.equivalent_creep_strain += delta_equiv_creep
            
            # 计算蠕变应变方向
            direction = self.creep_material.compute_creep_strain_direction(stress)
            
            # 更新蠕变应变向量
            ip.creep_strain += delta_equiv_creep * direction
    
    def get_creep_state(self) -> Dict:
        """获取单元的蠕变状态"""
        if not self.is_creep:
            return {'is_creep': False, 'equivalent_creep_strain': 0.0}
        
        max_equiv_creep_strain = 0.0
        
        for ip in self.integration_points:
            max_equiv_creep_strain = max(
                max_equiv_creep_strain, ip.equivalent_creep_strain)
        
        return {
            'is_creep': True,
            'equivalent_creep_strain': max_equiv_creep_strain
        }
    
    def compute_load_vector(self, internal_pressure: float,
                           gravity: np.ndarray,
                           temperature: float,
                           ref_temperature: float) -> np.ndarray:
        """计算单元载荷向量"""
        from solver.load.pressure import PressureLoad
        from solver.load.gravity import GravityLoad
        from solver.load.thermal_strain import ThermalStrainLoad
        
        F = np.zeros(12)
        L = self.L
        
        # 内压等效轴向力
        if internal_pressure > 0:
            A_internal = math.pi * self.Ri**2
            F_axial = internal_pressure * 1e6 * A_internal  # MPa -> Pa
            F_axial_node = F_axial / 2.0
            F[0:3] += -F_axial_node * self.ex
            F[6:9] += F_axial_node * self.ex
        
        # 自重
        gravity_vec = np.array(gravity) / 1000.0
        gravity_norm = np.linalg.norm(gravity_vec)
        if gravity_norm > 1e-10:
            w = GravityLoad.compute_distributed_load(
                self.material.density, self.A, gravity)
            F_gravity = w * L / 2.0
            gravity_dir = gravity_vec / gravity_norm
            F[2] += F_gravity * gravity_dir[2]
            F[8] += F_gravity * gravity_dir[2]
        
        # 热膨胀
        if abs(temperature - ref_temperature) > 1e-6:
            F_thermal = ThermalStrainLoad.compute_thermal_force(
                self.material.E, self.A, self.material.alpha,
                temperature, ref_temperature)
            F[0:3] += -F_thermal / 2.0 * self.ex
            F[6:9] += F_thermal / 2.0 * self.ex
        
        return F

