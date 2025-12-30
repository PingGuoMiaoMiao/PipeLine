"""
PIPE288 Element with Plasticity - 支持弹塑性的PIPE288单元
"""

import math
import numpy as np
from typing import Dict, Optional
from solver.mesh.node import Node
from solver.material.elastic import ElasticMaterial
from solver.material.plastic_biso import BISOPlasticMaterial
from solver.integration.integration_point import IntegrationPoint


class PIPE288ElementPlastic:
    """PIPE288单元类 - 支持弹塑性分析"""
    
    def __init__(self, elem_id: int, node1: Node, node2: Node,
                 material: ElasticMaterial, section: Dict,
                 plastic_material: Optional[BISOPlasticMaterial] = None):
        """
        初始化PIPE288单元（支持弹塑性）
        
        参数:
            elem_id: 单元ID
            node1: 节点1
            node2: 节点2
            material: 弹性材料对象
            section: 截面参数字典（diameter, wall_thickness）
            plastic_material: 弹塑性材料对象（可选，如果提供则启用弹塑性）
        """
        self.elem_id = elem_id
        self.node1 = node1
        self.node2 = node2
        self.material = material
        self.plastic_material = plastic_material
        self.is_plastic = plastic_material is not None
        
        # 截面参数（ANSYS中为mm，转换为m）
        D = section.get('diameter', 30.0)
        t = section.get('wall_thickness', 1.0)
        
        self.Ro = D / 2.0 / 1000.0  # 外半径 (m) = D(mm) / 2 / 1000
        self.Ri = (D - 2 * t) / 2.0 / 1000.0  # 内半径 (m) = (D-2t)(mm) / 2 / 1000
        
        # 计算截面几何参数
        self.A = math.pi * (self.Ro**2 - self.Ri**2)  # 横截面积 (m²)
        self.I = math.pi / 4.0 * (self.Ro**4 - self.Ri**4)  # 惯性矩 (m⁴)
        self.J = math.pi / 2.0 * (self.Ro**4 - self.Ri**4)  # 极惯性矩 (m⁴)
        
        # 计算单元长度和方向
        dx = (node2.x - node1.x) / 1000.0
        dy = (node2.y - node1.y) / 1000.0
        dz = (node2.z - node1.z) / 1000.0
        self.L = math.sqrt(dx**2 + dy**2 + dz**2)
        
        if self.L > 1e-10:
            self.ex = np.array([dx, dy, dz]) / self.L
            # 构造局部坐标系
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
        
        # 积分点（简化：使用2个积分点）
        self.integration_points = []
        if self.is_plastic:
            # 创建积分点
            num_gauss = 2
            gauss_points = [-1.0/np.sqrt(3), 1.0/np.sqrt(3)]
            gauss_weights = [1.0, 1.0]
            for i in range(num_gauss):
                ip = IntegrationPoint()
                ip.xi = gauss_points[i]  # 局部坐标
                ip.weight = gauss_weights[i]
                self.integration_points.append(ip)
    
    def compute_stiffness_matrix(self) -> np.ndarray:
        """
        计算单元刚度矩阵（弹性或弹塑性切线刚度）
        
        注意：如果使用弹塑性，应使用compute_tangent_stiffness
        """
        if not self.is_plastic:
            return self._compute_elastic_stiffness_matrix()
        else:
            # 使用当前状态的切线刚度
            # 这里返回弹性刚度，实际应该使用compute_tangent_stiffness
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
        """
        计算切线刚度矩阵（用于Newton-Raphson迭代）
        
        参数:
            node_displacements: 节点位移（12×1）
        
        返回:
            切线刚度矩阵（12×12）
        """
        if not self.is_plastic:
            return self._compute_elastic_stiffness_matrix()
        
        # 简化：对于梁单元，主要考虑轴向和弯曲的弹塑性
        # 这里使用简化的方法：如果进入塑性，则降低刚度
        K_elastic = self._compute_elastic_stiffness_matrix()
        
        # 检查积分点是否进入塑性
        # 简化：计算单元平均应变，判断是否进入塑性
        # 实际应该在每个积分点处计算
        reduction_factor = 1.0
        for ip in self.integration_points:
            if ip.is_plastic:
                # 如果积分点进入塑性，降低刚度
                reduction_factor = min(reduction_factor, 0.1)  # 简化处理
        
        return K_elastic * reduction_factor
    
    def compute_internal_force(self, node_displacements: np.ndarray) -> np.ndarray:
        """
        计算单元内部力向量
        
        参数:
            node_displacements: 节点位移（12×1）
        
        返回:
            内部力向量（12×1）
        """
        if not self.is_plastic:
            # 弹性情况：F_int = K * U
            K = self._compute_elastic_stiffness_matrix()
            return K @ node_displacements
        
        # 弹塑性情况：需要从应力计算
        # 简化实现：使用切线刚度
        K_tangent = self.compute_tangent_stiffness(node_displacements)
        return K_tangent @ node_displacements
    
    def compute_load_vector(self, internal_pressure: float,
                           gravity: np.ndarray,
                           temperature: float,
                           ref_temperature: float) -> np.ndarray:
        """计算单元载荷向量（12×1）"""
        from solver.load.pressure import PressureLoad
        from solver.load.gravity import GravityLoad
        from solver.load.thermal_strain import ThermalStrainLoad
        
        F = np.zeros(12)
        L = self.L
        
        # 内压等效轴向力
        if internal_pressure > 0:
            F_axial = PressureLoad.compute_internal_pressure_force(
                internal_pressure * 1e6, self.Ri)  # MPa -> Pa
            F_axial_node = F_axial / 2.0
            F[0:3] += -F_axial_node * self.ex
            F[6:9] += F_axial_node * self.ex
        
        # 自重载荷
        gravity_vec = np.array(gravity) / 1000.0
        gravity_norm = np.linalg.norm(gravity_vec)
        if gravity_norm > 1e-10:
            w = GravityLoad.compute_distributed_load(
                self.material.density, self.A, gravity)
            F_gravity = w * L / 2.0
            gravity_dir = gravity_vec / gravity_norm
            F[2] += F_gravity * gravity_dir[2]
            F[8] += F_gravity * gravity_dir[2]
        
        # 热膨胀载荷
        if abs(temperature - ref_temperature) > 1e-6:
            F_thermal = ThermalStrainLoad.compute_thermal_force(
                self.material.E, self.A, self.material.alpha,
                temperature, ref_temperature)
            F[0:3] += -F_thermal / 2.0 * self.ex
            F[6:9] += F_thermal / 2.0 * self.ex
        
        return F
    
    def get_plastic_state(self) -> Dict:
        """
        获取单元的塑性状态（用于输出）
        
        返回:
            包含塑性信息的字典
        """
        if not self.is_plastic:
            return {'is_plastic': False, 'equivalent_plastic_strain': 0.0}
        
        # 获取最大等效塑性应变
        max_equiv_plastic_strain = 0.0
        any_plastic = False
        
        for ip in self.integration_points:
            if ip.is_plastic:
                any_plastic = True
                max_equiv_plastic_strain = max(
                    max_equiv_plastic_strain, ip.equivalent_plastic_strain)
        
        return {
            'is_plastic': any_plastic,
            'equivalent_plastic_strain': max_equiv_plastic_strain
        }

