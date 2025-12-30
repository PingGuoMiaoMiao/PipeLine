"""
ELBOW290 Element - ELBOW290弯管单元类
"""

import math
import numpy as np
from typing import Dict, List, Tuple
from solver.mesh.node import Node
from solver.material.elastic import ElasticMaterial
from solver.element.shape_function_fourier import FourierShapeFunction


class ELBOW290Element:
    """ELBOW290弯管单元类 - 线弹性分析（简化版）"""
    
    # 每个节点10个自由度：6个标准DOF + 4个截面变形DOF
    DOF_PER_NODE = 10
    # 标准DOF: UX, UY, UZ, ROTX, ROTY, ROTZ
    # 截面变形DOF: δ_2c (椭圆化cos), δ_2s (椭圆化sin), 
    #              δ_3c (翘曲cos), δ_3s (翘曲sin)
    
    def __init__(self, elem_id: int, node1: Node, node2: Node, node3: Node,
                 material: ElasticMaterial, section: Dict):
        """
        初始化ELBOW290单元（3节点弯管单元）
        
        参数:
            elem_id: 单元ID
            node1: 起始节点
            node2: 中间节点（用于定义弯曲）
            node3: 结束节点
            material: 材料对象
            section: 截面参数字典（diameter, wall_thickness）
        """
        self.elem_id = elem_id
        self.node1 = node1
        self.node2 = node2
        self.node3 = node3
        self.material = material
        
        # 截面参数（ANSYS中为mm，转换为m）
        D = section.get('diameter', 90.0)
        t = section.get('wall_thickness', 2.0)
        
        self.Ro = D / 2.0 / 1000.0  # 外半径 (m) = D(mm) / 2 / 1000
        self.Ri = (D - 2 * t) / 2.0 / 1000.0  # 内半径 (m) = (D-2t)(mm) / 2 / 1000
        self.R_mean = (self.Ro + self.Ri) / 2.0  # 平均半径 (m)
        
        # 计算截面几何参数
        self.A = math.pi * (self.Ro**2 - self.Ri**2)  # 横截面积 (m²)
        self.I = math.pi / 4.0 * (self.Ro**4 - self.Ri**4)  # 惯性矩 (m⁴)
        self.J = math.pi / 2.0 * (self.Ro**4 - self.Ri**4)  # 极惯性矩 (m⁴)
        
        # 计算弯管几何（曲率半径和角度）
        self._compute_bend_geometry()
        
        # 局部坐标系
        self._compute_local_coordinate_system()
    
    def _compute_bend_geometry(self):
        """计算弯管几何参数（曲率半径和角度）"""
        # 节点坐标（转换为m）
        p1 = np.array([self.node1.x, self.node1.y, self.node1.z]) / 1000.0
        p2 = np.array([self.node2.x, self.node2.y, self.node2.z]) / 1000.0
        p3 = np.array([self.node3.x, self.node3.y, self.node3.z]) / 1000.0
        
        # 计算两个向量
        v1 = p2 - p1
        v2 = p3 - p2
        
        L1 = np.linalg.norm(v1)
        L2 = np.linalg.norm(v2)
        
        if L1 < 1e-10 or L2 < 1e-10:
            raise ValueError(f"单元{self.elem_id}节点距离过小")
        
        # 归一化
        v1_norm = v1 / L1
        v2_norm = v2 / L2
        
        # 计算角度（使用点积）
        cos_angle = np.dot(v1_norm, v2_norm)
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        self.bend_angle = math.acos(cos_angle)  # 弯曲角度 (rad)
        
        # 简化的曲率半径估算（使用弧长和角度）
        # 对于弯管，假设节点在圆弧上，曲率半径 R = L / θ
        self.arc_length = L1 + L2  # 近似弧长
        if self.bend_angle > 1e-6:
            self.R_curvature = self.arc_length / self.bend_angle
        else:
            # 如果角度很小，使用大的曲率半径（接近直管）
            self.R_curvature = 1e6
        
        # 计算每个段的长度
        self.L1 = L1
        self.L2 = L2
        self.L_total = L1 + L2
        self.L = self.L_total  # 总长度，用于兼容性
    
    def _compute_local_coordinate_system(self):
        """计算局部坐标系"""
        # 节点坐标（转换为m）
        p1 = np.array([self.node1.x, self.node1.y, self.node1.z]) / 1000.0
        p2 = np.array([self.node2.x, self.node2.y, self.node2.z]) / 1000.0
        p3 = np.array([self.node3.x, self.node3.y, self.node3.z]) / 1000.0
        
        # 第一段的轴向方向
        v1 = p2 - p1
        self.ex1 = v1 / np.linalg.norm(v1)
        
        # 第二段的轴向方向
        v2 = p3 - p2
        self.ex2 = v2 / np.linalg.norm(v2)
        
        # 平均轴向方向（用于简化）
        ex_avg = (self.ex1 + self.ex2) / 2.0
        ex_avg = ex_avg / np.linalg.norm(ex_avg)
        self.ex = ex_avg
        
        # 构造局部Y轴和Z轴（在平均方向上）
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
    
    def compute_stiffness_matrix(self) -> np.ndarray:
        """
        计算单元刚度矩阵（全局坐标系，30×30）
        
        简化策略：
        1. 基于PIPE288的刚度矩阵，考虑曲率修正
        2. 添加椭圆化和翘曲的简化刚度项
        """
        # 总自由度：3节点 × 10 DOF = 30
        K = np.zeros((30, 30))
        
        E = self.material.E
        G = self.material.G
        
        # 1. 标准DOF的刚度（类似于PIPE288，但使用等效长度）
        # 简化为将弯管视为两段直管的组合
        
        # 第一段（节点1到节点2）
        K_seg1 = self._compute_segment_stiffness(E, G, self.L1)
        # 组装到总矩阵（节点1和节点2的标准DOF：0-5, 10-15）
        for i in range(6):
            for j in range(6):
                K[i, j] += K_seg1[i, j]
                K[i, j + 10] += K_seg1[i, j + 6]
                K[i + 10, j] += K_seg1[i + 6, j]
                K[i + 10, j + 10] += K_seg1[i + 6, j + 6]
        
        # 第二段（节点2到节点3）
        K_seg2 = self._compute_segment_stiffness(E, G, self.L2)
        # 组装到总矩阵（节点2和节点3的标准DOF：10-15, 20-25）
        for i in range(6):
            for j in range(6):
                K[i + 10, j + 10] += K_seg2[i, j]
                K[i + 10, j + 20] += K_seg2[i, j + 6]
                K[i + 20, j + 10] += K_seg2[i + 6, j]
                K[i + 20, j + 20] += K_seg2[i + 6, j + 6]
        
        # 2. 椭圆化刚度（显示主导模式：极弱刚度，仅保持矩阵非奇异）
        # 椭圆化模态不参与力学计算，仅用于显示
        # 设置极小的虚拟刚度，确保矩阵非奇异
        if self.R_curvature < 1e6:
            t = self.Ro - self.Ri
            # 虚拟刚度：非常小，仅为数值稳定性（比实际刚度小几个数量级）
            k_oval_virtual = E * t**3 / (self.R_curvature**2) * 1e-6  # 虚拟刚度
            k_oval_node = k_oval_virtual / 3.0
            
            # 椭圆化DOF索引：6-7（节点1），16-17（节点2），26-27（节点3）
            for node_offset in [0, 10, 20]:
                oval_idx = node_offset + 6
                K[oval_idx, oval_idx] = k_oval_node  # cos(2φ)项（虚拟）
                K[oval_idx + 1, oval_idx + 1] = k_oval_node  # sin(2φ)项（虚拟）
        
        # 3. 翘曲刚度（同样设置为虚拟刚度）
        if self.R_curvature < 1e6:
            t = self.Ro - self.Ri
            k_warp_virtual = E * t**3 / (self.R_curvature**3) * 1e-6  # 虚拟刚度
            k_warp_node = k_warp_virtual / 3.0
            
            # 翘曲DOF索引：8-9（节点1），18-19（节点2），28-29（节点3）
            for node_offset in [0, 10, 20]:
                warp_idx = node_offset + 8
                K[warp_idx, warp_idx] = k_warp_node  # cos(3φ)项（虚拟）
                K[warp_idx + 1, warp_idx + 1] = k_warp_node  # sin(3φ)项（虚拟）
        
        # 坐标变换（简化：使用平均方向）
        K_global = self._transform_to_global(K)
        
        return K_global
    
    def _compute_segment_stiffness(self, E: float, G: float, L: float) -> np.ndarray:
        """计算单段直管的刚度矩阵（12×12，类似PIPE288）"""
        K_local = np.zeros((12, 12))
        
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
        # Y方向
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
        # Z方向
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
        
        return K_local
    
    def _transform_to_global(self, K_local: np.ndarray) -> np.ndarray:
        """将刚度矩阵转换到全局坐标系"""
        # 旋转矩阵
        R = np.array([
            [self.ex[0], self.ex[1], self.ex[2]],
            [self.ey[0], self.ey[1], self.ey[2]],
            [self.ez[0], self.ez[1], self.ez[2]]
        ])
        
        # 为每个节点构建变换矩阵（10×10）
        T_node = np.zeros((10, 10))
        T_node[0:3, 0:3] = R  # 平移
        T_node[3:6, 3:6] = R  # 转动
        T_node[6:10, 6:10] = np.eye(4)  # 截面变形DOF不变换（简化）
        
        # 总变换矩阵（30×30）
        T = np.zeros((30, 30))
        for i in range(3):
            T[i*10:(i+1)*10, i*10:(i+1)*10] = T_node
        
        # 变换
        K_global = T.T @ K_local @ T
        
        return K_global
    
    def compute_load_vector(self, internal_pressure: float,
                           gravity: np.ndarray,
                           temperature: float,
                           ref_temperature: float) -> np.ndarray:
        """计算单元载荷向量（全局坐标系，30×1）"""
        F = np.zeros(30)
        
        # 标准载荷（类似PIPE288，分配到三个节点）
        # 简化：将载荷分配到三个节点
        
        # 内压等效轴向力
        if internal_pressure > 0:
            from solver.load.pressure import PressureLoad
            A_internal = math.pi * self.Ri**2
            F_axial = internal_pressure * 1e6 * A_internal
            F_axial_node = F_axial / 3.0
            # 分配到三个节点
            for node_offset in [0, 10, 20]:
                F[node_offset:node_offset+3] += -F_axial_node / 2.0 * self.ex
            for node_offset in [10, 20]:
                F[node_offset:node_offset+3] += F_axial_node / 2.0 * self.ex
        
        # 自重载荷
        gravity_vec = np.array(gravity) / 1000.0  # mm/s² -> m/s²
        gravity_norm = np.linalg.norm(gravity_vec)
        if gravity_norm > 1e-10:
            w = self.material.density * self.A * gravity_norm
            F_gravity_node = w * self.arc_length / 3.0
            gravity_dir = gravity_vec / gravity_norm
            for node_offset in [0, 10, 20]:
                F[node_offset + 2] += F_gravity_node * gravity_dir[2]
        
        # 热膨胀载荷
        if abs(temperature - ref_temperature) > 1e-6:
            from solver.load.thermal_strain import ThermalStrainLoad
            F_thermal = ThermalStrainLoad.compute_thermal_force(
                self.material.E, self.A, self.material.alpha,
                temperature, ref_temperature)
            F_thermal_node = F_thermal / 3.0
            for node_offset in [0, 10, 20]:
                F[node_offset:node_offset+3] += -F_thermal_node / 2.0 * self.ex
        
        return F
    
    def get_section_points(self, node_pos: np.ndarray, node_dofs: np.ndarray,
                          num_circumferential: int = 20) -> List[Tuple[float, float, float]]:
        """
        获取截面点坐标（用于VTK输出）
        
        参数:
            node_pos: 节点位置（m）
            node_dofs: 节点的10个自由度（前6个是标准DOF，后4个是截面变形DOF）
            num_circumferential: 环向划分点数
        
        返回:
            截面点坐标列表 [(x, y, z), ...]
        """
        
        # 位移（前3个DOF）
        displacement = node_dofs[0:3]
        
        # 截面变形DOF
        oval_cos = node_dofs[6]  # δ_2c
        oval_sin = node_dofs[7]  # δ_2s
        warp_cos = node_dofs[8]  # δ_3c
        warp_sin = node_dofs[9]  # δ_3s
        
        # 变形后的节点位置
        deformed_pos = node_pos + displacement
        
        # 生成环向角度
        phis = FourierShapeFunction.generate_circumferential_points(num_circumferential)
        
        # 计算截面点
        points = []
        for phi in phis:
            # 截面局部坐标
            x_local, y_local = FourierShapeFunction.compute_section_coordinates(
                self.R_mean, phi, oval_cos, oval_sin, warp_cos, warp_sin)
            
            # 转换到全局坐标系
            # 使用局部坐标系基向量
            point_global = (deformed_pos + 
                          x_local * self.ey + 
                          y_local * self.ez)
            
            # 转换回mm
            points.append((point_global[0] * 1000.0,
                          point_global[1] * 1000.0,
                          point_global[2] * 1000.0))
        
        return points
