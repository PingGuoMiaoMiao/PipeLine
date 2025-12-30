"""
PIPE288 Element - PIPE288单元类（线弹性分析）
"""

import math
import numpy as np
from typing import Dict
from solver.mesh.node import Node
from solver.material.elastic import ElasticMaterial


class PIPE288Element:
    """PIPE288单元类 - 线弹性分析"""
    
    def __init__(self, elem_id: int, node1: Node, node2: Node, 
                 material: ElasticMaterial, section: Dict):
        """
        初始化PIPE288单元
        
        参数:
            elem_id: 单元ID
            node1: 节点1
            node2: 节点2
            material: 材料对象
            section: 截面参数字典（diameter, wall_thickness）
        """
        self.elem_id = elem_id
        self.node1 = node1
        self.node2 = node2
        self.material = material
        
        # 截面参数（ANSYS中为mm，转换为m）
        D = section.get('diameter', 30.0)
        t = section.get('wall_thickness', 1.0)
        
        self.Ro = D / 2.0 / 1000.0  # 外半径 (m) = D(mm) / 2 / 1000
        self.Ri = (D - 2 * t) / 2.0 / 1000.0  # 内半径 (m) = (D-2t)(mm) / 2 / 1000
        
        # 计算截面几何参数
        self.A = math.pi * (self.Ro**2 - self.Ri**2)  # 横截面积 (m²)
        self.I = math.pi / 4.0 * (self.Ro**4 - self.Ri**4)  # 惯性矩 (m⁴)
        self.J = math.pi / 2.0 * (self.Ro**4 - self.Ri**4)  # 极惯性矩 (m⁴)
        
        # 计算单元长度和方向（ANSYS中坐标为mm，转换为m）
        dx = (node2.x - node1.x) / 1000.0
        dy = (node2.y - node1.y) / 1000.0
        dz = (node2.z - node1.z) / 1000.0
        self.L = math.sqrt(dx**2 + dy**2 + dz**2)  # 单元长度 (m)
        
        # 局部坐标系方向向量
        if self.L > 1e-10:
            self.ex = np.array([dx, dy, dz]) / self.L
            # 构造局部Y轴
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
            self.ey = np.cross(self.ez, self.ex)  # 重新正交化
        else:
            raise ValueError(f"单元{elem_id}长度为零")
    
    def compute_stiffness_matrix(self) -> np.ndarray:
        """计算单元刚度矩阵（全局坐标系，12×12）"""
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
        
        # Y方向弯曲
        EI = E * self.I
        k_bend_y = EI / L**3
        K_local[1, 1] = 12 * k_bend_y
        K_local[1, 5] = 6 * L * k_bend_y
        K_local[1, 7] = -12 * k_bend_y
        K_local[1, 11] = 6 * L * k_bend_y
        K_local[5, 1] = 6 * L * k_bend_y
        K_local[5, 5] = 4 * L**2 * k_bend_y
        K_local[5, 7] = -6 * L * k_bend_y
        K_local[5, 11] = 2 * L**2 * k_bend_y
        K_local[7, 1] = -12 * k_bend_y
        K_local[7, 5] = -6 * L * k_bend_y
        K_local[7, 7] = 12 * k_bend_y
        K_local[7, 11] = -6 * L * k_bend_y
        K_local[11, 1] = 6 * L * k_bend_y
        K_local[11, 5] = 2 * L**2 * k_bend_y
        K_local[11, 7] = -6 * L * k_bend_y
        K_local[11, 11] = 4 * L**2 * k_bend_y
        
        # Z方向弯曲（注意符号）
        K_local[2, 2] = 12 * k_bend_y
        K_local[2, 4] = -6 * L * k_bend_y
        K_local[2, 8] = -12 * k_bend_y
        K_local[2, 10] = -6 * L * k_bend_y
        K_local[4, 2] = -6 * L * k_bend_y
        K_local[4, 4] = 4 * L**2 * k_bend_y
        K_local[4, 8] = 6 * L * k_bend_y
        K_local[4, 10] = 2 * L**2 * k_bend_y
        K_local[8, 2] = -12 * k_bend_y
        K_local[8, 4] = 6 * L * k_bend_y
        K_local[8, 8] = 12 * k_bend_y
        K_local[8, 10] = 6 * L * k_bend_y
        K_local[10, 2] = -6 * L * k_bend_y
        K_local[10, 4] = 2 * L**2 * k_bend_y
        K_local[10, 8] = 6 * L * k_bend_y
        K_local[10, 10] = 4 * L**2 * k_bend_y
        
        # 坐标变换矩阵
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
        
        # 转换到全局坐标系
        K_global = T.T @ K_local @ T
        
        return K_global
    
    def compute_load_vector(self, internal_pressure: float, 
                          gravity: np.ndarray, 
                          temperature: float, 
                          ref_temperature: float) -> np.ndarray:
        """计算单元载荷向量（全局坐标系，12×1）"""
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
        
        # 自重载荷（ANSYS中重力单位为mm/s²，需转换为m/s²）
        gravity_vec = np.array(gravity) / 1000.0  # mm/s² -> m/s²
        gravity_norm = np.linalg.norm(gravity_vec)
        if gravity_norm > 1e-10:
            w = self.material.density * self.A * gravity_norm  # 单位长度重量 (N/m)
            F_gravity = w * L / 2.0  # 等效节点力
            gravity_dir = gravity_vec / gravity_norm
            F[2] += F_gravity * gravity_dir[2]  # 节点1 Z方向
            F[8] += F_gravity * gravity_dir[2]  # 节点2 Z方向
        
        # 热膨胀载荷
        if abs(temperature - ref_temperature) > 1e-6:
            F_thermal = ThermalStrainLoad.compute_thermal_force(
                self.material.E, self.A, self.material.alpha,
                temperature, ref_temperature)
            F[0:3] += -F_thermal / 2.0 * self.ex
            F[6:9] += F_thermal / 2.0 * self.ex
        
        return F

