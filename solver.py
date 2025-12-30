#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANSYS CDB 到 VTK 转换器 - PIPE288 线弹性分析
解析CDB文件，进行有限元分析，输出节点位移到VTK格式
"""

import sys
import re
import math
from typing import List, Dict, Tuple, Optional
import numpy as np


class Node:
    """节点类"""
    def __init__(self, node_id: int, x: float, y: float, z: float):
        self.id = node_id
        self.x = x
        self.y = y
        self.z = z
        self.displacement = np.zeros(6)  # [UX, UY, UZ, ROTX, ROTY, ROTZ]
    
    def __repr__(self):
        return f"Node({self.id}, ({self.x}, {self.y}, {self.z}))"


class Element:
    """单元类"""
    def __init__(self, elem_id: int, elem_type: int, node_ids: List[int]):
        self.id = elem_id
        self.type = elem_type  # 288 for PIPE288, 290 for ELBOW290
        self.node_ids = node_ids
    
    def __repr__(self):
        type_name = "PIPE288" if self.type == 288 else "ELBOW290" if self.type == 290 else str(self.type)
        return f"Element({self.id}, {type_name}, nodes={self.node_ids})"


class CDBParser:
    """CDB文件解析器"""
    
    def __init__(self, filename: str):
        self.filename = filename
        self.nodes: Dict[int, Node] = {}
        self.elements: Dict[int, Element] = {}
        self.element_type_map: Dict[int, int] = {}  # 元素类型ID -> 单元类型号 (288/290)
        
        # 材料参数
        self.materials: Dict[int, Dict] = {}  # 材料ID -> 材料参数
        # 截面参数
        self.sections: Dict[int, Dict] = {}  # 截面ID -> 截面参数
        # 载荷
        self.gravity: List[float] = [0.0, 0.0, 0.0]  # 重力加速度 [X, Y, Z]
        self.temperature: float = 25.0  # 温度
        self.reference_temp: float = 25.0  # 参考温度
        self.internal_pressure: Dict[int, float] = {}  # 单元ID -> 内压
        # 边界条件
        self.boundary_conditions: Dict[int, List[int]] = {}  # 节点ID -> [DOF列表]，0=UX, 1=UY, 2=UZ, 3=ROTX, 4=ROTY, 5=ROTZ
        
    def parse(self):
        """解析CDB文件"""
        with open(self.filename, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()
        
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            
            # 解析ET命令：ET, type_id, unit_type
            if line.startswith('ET,'):
                self._parse_et(line)
            
            # 解析NBLOCK（节点块）
            elif line.startswith('NBLOCK,'):
                i = self._parse_nblock(lines, i)
            
            # 解析EBLOCK（单元块）
            elif line.startswith('EBLOCK,'):
                i = self._parse_eblock(lines, i)
            
            # 解析材料参数 MPDATA
            elif line.startswith('MPDATA,'):
                self._parse_mpdata(line)
            
            # 解析参考温度 REFT
            elif line.startswith('MPTEMP,') and 'REFT' in line:
                self._parse_reft(lines, i)
            
            # 解析截面参数 SECDATA
            elif line.startswith('SECDATA,'):
                self._parse_secdata(line)
            
            # 解析重力 ACEL
            elif line.startswith('ACEL,'):
                self._parse_acel(line)
            
            # 解析温度 BFUNIF,TEMP
            elif 'BFUNIF,TEMP' in line:
                self._parse_temperature(line)
            
            # 解析参考温度 TREF
            elif line.startswith('TREF,'):
                self._parse_tref(line)
            
            # 解析边界条件 D,
            elif line.startswith('D,'):
                self._parse_boundary_condition(line)
            
            # 解析内压 SFE, PRES
            elif line.startswith('SFE,'):
                i = self._parse_sfe_pressure(lines, i)
            
            i += 1
    
    def _parse_et(self, line: str):
        """解析ET命令，建立元素类型映射"""
        # ET, type_id, unit_type
        parts = [p.strip() for p in line.split(',')]
        if len(parts) >= 3:
            try:
                type_id = int(parts[1])
                unit_type = int(parts[2])
                if unit_type in [288, 290]:  # 只关注PIPE288和ELBOW290
                    self.element_type_map[type_id] = unit_type
            except ValueError:
                pass
    
    def _parse_nblock(self, lines: List[str], start_idx: int) -> int:
        """解析NBLOCK节点块"""
        # NBLOCK格式：NBLOCK,6,SOLID,min_node_id,max_node_id
        header_line = lines[start_idx].strip()
        parts = header_line.split(',')
        
        # 下一行是格式说明：(3i9,6e21.13e3)
        format_line = lines[start_idx + 1].strip() if start_idx + 1 < len(lines) else ""
        
        # 从下一行开始读取节点数据
        i = start_idx + 2
        while i < len(lines):
            line = lines[i].strip()
            
            # 遇到结束标记或下一个命令，停止解析
            if not line or line.startswith('N,R') or line.startswith('N,') and 'LOC' in line:
                break
            
            # 解析节点数据：节点ID, 标志1, 标志2, X, Y, Z, ...
            # 格式：(3i9,6e21.13e3) 表示3个整数9位宽，6个浮点数21.13位宽
            parts = line.split()
            if len(parts) >= 3:
                try:
                    node_id = int(parts[0])
                    # flags = [int(parts[1]), int(parts[2])]
                    
                    # 读取坐标（可能有1-3个坐标）
                    x = float(parts[3]) if len(parts) > 3 else 0.0
                    y = float(parts[4]) if len(parts) > 4 else 0.0
                    z = float(parts[5]) if len(parts) > 5 else 0.0
                    
                    self.nodes[node_id] = Node(node_id, x, y, z)
                except (ValueError, IndexError):
                    pass
            
            i += 1
        
        return i - 1
    
    def _parse_eblock(self, lines: List[str], start_idx: int) -> int:
        """解析EBLOCK单元块"""
        # EBLOCK格式：EBLOCK,19,SOLID,min_elem_id,max_elem_id
        header_line = lines[start_idx].strip()
        
        # 下一行是格式说明：(19i9)
        format_line = lines[start_idx + 1].strip() if start_idx + 1 < len(lines) else ""
        
        # 从下一行开始读取单元数据
        elem_counter = 1
        i = start_idx + 2
        while i < len(lines):
            line = lines[i].strip()
            
            # 遇到结束标记（-1），停止解析
            if line == '-1' or (line.startswith('-') and line.strip() == '-1'):
                break
            
            # 解析单元数据：19个整数
            parts = line.split()
            if len(parts) >= 12:
                try:
                    # EBLOCK格式解析：
                    # 字段索引（从0开始）：
                    # 0: 材料ID
                    # 1: 单元类型ID (用于查找288/290)
                    # 2-7: 其他属性
                    # 8: 节点数
                    # 9: 单元ID（如果为0，可能需要特殊处理）
                    # 10: 可能是辅助字段或单元ID的一部分
                    # 11+: 节点ID列表（从索引11开始）
                    
                    num_nodes = int(parts[8]) if len(parts) > 8 else 2
                    elem_type_id = int(parts[1]) if len(parts) > 1 else 1
                    
                    # 获取实际的单元类型（288或290）
                    elem_type = self.element_type_map.get(elem_type_id, 288)
                    
                    # 根据ANSYS EBLOCK格式分析：
                    # 字段9（索引9）通常是单元ID，如果为0，则字段10（索引10）可能是单元ID
                    # 节点ID从字段11（索引11）开始
                    # 从实际数据看，字段9为0时，字段10是单元ID（1,2,3...），节点从字段11开始
                    
                    # 确定单元ID
                    if len(parts) > 9:
                        try:
                            field9_val = int(parts[9])
                            if field9_val > 0:
                                current_elem_id = field9_val
                            elif len(parts) > 10:
                                # 字段9为0时，使用字段10作为单元ID
                                current_elem_id = int(parts[10])
                            else:
                                current_elem_id = elem_counter
                        except (ValueError, IndexError):
                            current_elem_id = elem_counter
                    else:
                        current_elem_id = elem_counter
                    
                    # 读取节点ID（从字段11开始，即索引11）
                    node_ids = []
                    start_node_idx = 11  # 节点ID从字段12开始（索引11）
                    for j in range(start_node_idx, min(start_node_idx + num_nodes, len(parts))):
                        node_id = int(parts[j])
                        if node_id > 0:  # 有效的节点ID
                            node_ids.append(node_id)
                    
                    if node_ids:
                        self.elements[current_elem_id] = Element(current_elem_id, elem_type, node_ids)
                        elem_counter = max(elem_counter, current_elem_id) + 1
                
                except (ValueError, IndexError) as e:
                    pass
            
            i += 1
        
        return i - 1
    
    def _parse_mpdata(self, line: str):
        """解析材料参数 MPDATA"""
        # MPDATA,label,MAT, C0, LOC1, C1, LOC2, C2, ...
        # 例如: MPDATA,R5.0, 1,EX  ,       1, 1,  200000.000
        parts = [p.strip() for p in line.split(',')]
        if len(parts) >= 7:
            try:
                mat_id = int(parts[2])
                prop_name = parts[3].strip()
                # 值在parts[6]（索引6）
                value = float(parts[6])
                
                if mat_id not in self.materials:
                    self.materials[mat_id] = {}
                
                if prop_name == 'EX':
                    # ANSYS中EX单位为MPa，转换为Pa
                    self.materials[mat_id]['E'] = value * 1e6  # MPa -> Pa
                elif prop_name in ['PRXY', 'NUXY']:
                    self.materials[mat_id]['nu'] = value  # 泊松比（无量纲）
                elif prop_name == 'DENS':
                    # ANSYS中密度单位为tonne/mm³，转换为kg/m³
                    # 1 tonne/mm³ = 1000 kg / (1e-9 m³) = 1e12 kg/m³
                    self.materials[mat_id]['density'] = value * 1e12  # tonne/mm³ -> kg/m³
                elif prop_name == 'ALPX':
                    self.materials[mat_id]['alpha'] = value  # 热膨胀系数（1/°C）
            except (ValueError, IndexError) as e:
                # 调试用
                pass
    
    def _parse_reft(self, lines: List[str], start_idx: int):
        """解析参考温度"""
        # MPDATA,R5.0, 1,REFT,       1, 1,  25.0000000
        line = lines[start_idx].strip()
        parts = [p.strip() for p in line.split(',')]
        if len(parts) >= 6 and 'REFT' in parts[3]:
            try:
                mat_id = int(parts[2])
                value = float(parts[5])
                if mat_id not in self.materials:
                    self.materials[mat_id] = {}
                self.materials[mat_id]['T_ref'] = value
                self.reference_temp = value
            except (ValueError, IndexError):
                pass
    
    def _parse_secdata(self, line: str):
        """解析截面参数 SECDATA"""
        # SECDATA,  30.000    ,  1.0000    ,  20.000    ,  0.0000    ,  5.0000
        # 格式：直径, 壁厚, 圆周划分, 其他, 径向划分
        parts = [p.strip() for p in line.split(',')]
        if len(parts) >= 2:
            try:
                section_id = 1  # 默认截面ID为1
                diameter = float(parts[1])
                wall_thickness = float(parts[2]) if len(parts) > 2 else 1.0
                
                self.sections[section_id] = {
                    'diameter': diameter,
                    'wall_thickness': wall_thickness
                }
            except (ValueError, IndexError):
                pass
    
    def _parse_acel(self, line: str):
        """解析重力加速度 ACEL"""
        # ACEL,  0.00000000    ,  0.00000000    ,  9800.00000
        parts = [p.strip() for p in line.split(',')]
        if len(parts) >= 4:
            try:
                self.gravity[0] = float(parts[1])
                self.gravity[1] = float(parts[2])
                self.gravity[2] = float(parts[3])
            except (ValueError, IndexError):
                pass
    
    def _parse_temperature(self, line: str):
        """解析温度 BFUNIF,TEMP"""
        # BFUNIF,TEMP,  200.000000
        parts = [p.strip() for p in line.split(',')]
        if len(parts) >= 3:
            try:
                self.temperature = float(parts[2])
            except (ValueError, IndexError):
                pass
    
    def _parse_tref(self, line: str):
        """解析参考温度 TREF"""
        # TREF,  0.00000000
        parts = [p.strip() for p in line.split(',')]
        if len(parts) >= 2:
            try:
                # TREF通常为0，实际参考温度在MPDATA中
                pass
            except (ValueError, IndexError):
                pass
    
    def _parse_boundary_condition(self, line: str):
        """解析边界条件 D"""
        # D,      1,UX  ,  0.00000000    ,  0.00000000
        parts = [p.strip() for p in line.split(',')]
        if len(parts) >= 3:
            try:
                node_id = int(parts[1])
                dof_name = parts[2].strip()
                
                dof_map = {
                    'UX': 0, 'UY': 1, 'UZ': 2,
                    'ROTX': 3, 'ROTY': 4, 'ROTZ': 5
                }
                
                if dof_name in dof_map:
                    dof = dof_map[dof_name]
                    if node_id not in self.boundary_conditions:
                        self.boundary_conditions[node_id] = []
                    if dof not in self.boundary_conditions[node_id]:
                        self.boundary_conditions[node_id].append(dof)
            except (ValueError, IndexError):
                pass
    
    def _parse_sfe_pressure(self, lines: List[str], start_idx: int) -> int:
        """解析内压 SFE, PRES"""
        # SFE,        1,   1,PRES,1,R5.0
        #   3.00000000      0.00000000      0.00000000      0.00000000
        line = lines[start_idx].strip()
        parts = [p.strip() for p in line.split(',')]
        if len(parts) >= 2 and 'PRES' in line:
            try:
                elem_id = int(parts[1])
                # 读取下一行的压力值
                if start_idx + 1 < len(lines):
                    pressure_line = lines[start_idx + 1].strip()
                    pressure_parts = pressure_line.split()
                    if len(pressure_parts) > 0:
                        pressure = float(pressure_parts[0])
                        self.internal_pressure[elem_id] = pressure
                        return start_idx + 1  # 跳过下一行
            except (ValueError, IndexError):
                pass
        return start_idx


class PIPE288Element:
    """PIPE288单元类 - 线弹性分析"""
    
    def __init__(self, elem_id: int, node1: Node, node2: Node, 
                 material: Dict, section: Dict):
        self.elem_id = elem_id
        self.node1 = node1
        self.node2 = node2
        
        # 材料参数（所有参数使用SI单位：m, kg, Pa, s）
        self.E = material.get('E', 200e9)  # 弹性模量 (Pa)
        self.nu = material.get('nu', 0.3)  # 泊松比
        self.density = material.get('density', 7800.0)  # 密度 (kg/m³)
        self.alpha = material.get('alpha', 1.2e-5)  # 热膨胀系数 (1/°C)
        self.G = self.E / (2.0 * (1.0 + self.nu))  # 剪切模量 (Pa)
        
        # 截面参数（ANSYS中为mm，存储时保持mm，使用时转换）
        self.D = section.get('diameter', 30.0)  # 外径 (mm)
        self.t = section.get('wall_thickness', 1.0)  # 壁厚 (mm)
        
        # 转换为米（ANSYS中单位为mm）
        self.Ro = self.D / 2000.0  # 外半径 (m)
        self.Ri = (self.D - 2 * self.t) / 2000.0  # 内半径 (m)
        
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
            # 构造局部Y轴（如果X不接近Z轴，使用全局Z；否则使用全局Y）
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
        # 局部坐标系下的刚度矩阵
        K_local = np.zeros((12, 12))
        L = self.L
        
        # 轴向刚度
        k_axial = self.E * self.A / L
        K_local[0, 0] = k_axial
        K_local[0, 6] = -k_axial
        K_local[6, 0] = -k_axial
        K_local[6, 6] = k_axial
        
        # 扭转刚度
        k_torsion = self.G * self.J / L
        K_local[3, 3] = k_torsion
        K_local[3, 9] = -k_torsion
        K_local[9, 3] = -k_torsion
        K_local[9, 9] = k_torsion
        
        # Y方向弯曲（局部坐标系）
        EI = self.E * self.I
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
        F = np.zeros(12)
        L = self.L
        
        # 内压等效轴向力（ANSYS中压力单位为MPa）
        if internal_pressure > 0:
            A_internal = math.pi * self.Ri**2
            F_axial = internal_pressure * 1e6 * A_internal  # MPa -> Pa，然后乘以面积
            # 分配给两个节点（考虑端盖效应）
            F_axial_node = F_axial / 2.0
            # 沿轴向方向
            F[0:3] += -F_axial_node * self.ex
            F[6:9] += F_axial_node * self.ex
        
        # 自重载荷（均布载荷）
        w = self.density * self.A * 9.8  # 单位长度重量 (N/m)
        # 转换为等效节点力
        F_gravity = w * L / 2.0
        # 重力方向（通常是-Z方向）
        gravity_vec = np.array(gravity)
        gravity_norm = np.linalg.norm(gravity_vec)
        if gravity_norm > 1e-10:
            gravity_dir = gravity_vec / gravity_norm
            F[2] += F_gravity * gravity_dir[2]  # 节点1 Z方向
            F[8] += F_gravity * gravity_dir[2]  # 节点2 Z方向
        
        # 热膨胀载荷
        if abs(temperature - ref_temperature) > 1e-6:
            delta_T = temperature - ref_temperature
            F_thermal = self.E * self.A * self.alpha * delta_T
            F[0:3] += -F_thermal / 2.0 * self.ex
            F[6:9] += F_thermal / 2.0 * self.ex
        
        return F


class FEASolver:
    """有限元求解器"""
    
    def __init__(self, nodes: Dict[int, Node], 
                 elements: Dict[int, PIPE288Element],
                 boundary_conditions: Dict[int, List[int]]):
        self.nodes = nodes
        self.elements = elements
        self.boundary_conditions = boundary_conditions
        
        # 创建节点ID到DOF索引的映射
        sorted_node_ids = sorted(nodes.keys())
        self.node_id_to_dof_start = {node_id: idx * 6 for idx, node_id in enumerate(sorted_node_ids)}
        self.dof_count = len(nodes) * 6
    
    def solve(self, load_vector: np.ndarray) -> np.ndarray:
        """求解线性方程组 K * U = F"""
        # 组装全局刚度矩阵
        K_global = np.zeros((self.dof_count, self.dof_count))
        
        for elem in self.elements.values():
            K_elem = elem.compute_stiffness_matrix()
            # 获取单元节点的DOF索引
            dof1_start = self.node_id_to_dof_start[elem.node1.id]
            dof2_start = self.node_id_to_dof_start[elem.node2.id]
            
            # 组装到全局矩阵
            for i in range(6):
                for j in range(6):
                    K_global[dof1_start + i, dof1_start + j] += K_elem[i, j]
                    K_global[dof1_start + i, dof2_start + j] += K_elem[i, 6 + j]
                    K_global[dof2_start + i, dof1_start + j] += K_elem[6 + i, j]
                    K_global[dof2_start + i, dof2_start + j] += K_elem[6 + i, 6 + j]
        
        # 应用边界条件（置大数法）
        for node_id, dof_list in self.boundary_conditions.items():
            if node_id in self.node_id_to_dof_start:
                dof_start = self.node_id_to_dof_start[node_id]
                for dof in dof_list:
                    dof_idx = dof_start + dof
                    # 置大数
                    K_global[dof_idx, :] = 0
                    K_global[:, dof_idx] = 0
                    K_global[dof_idx, dof_idx] = 1e12
                    load_vector[dof_idx] = 0
        
        # 求解
        try:
            U = np.linalg.solve(K_global, load_vector)
        except np.linalg.LinAlgError:
            print("警告: 矩阵求解失败，尝试使用最小二乘")
            U = np.linalg.lstsq(K_global, load_vector, rcond=None)[0]
        
        return U


class VTKWriter:
    """VTK文件写入器（PolyLine格式）"""
    
    @staticmethod
    def write_polylines(filename: str, elements: Dict[int, Element], nodes: Dict[int, Node],
                       displacements: Optional[np.ndarray] = None):
        """将单元和节点写入VTK文件（PolyLine格式），可包含位移"""
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
            node_id_to_index = {node_id: idx for idx, (node_id, _) in enumerate(sorted_nodes)}
            
            for idx, (node_id, node) in enumerate(sorted_nodes):
                if displacements is not None:
                    # 添加位移
                    ux = displacements[idx * 6]
                    uy = displacements[idx * 6 + 1]
                    uz = displacements[idx * 6 + 2]
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
                valid_node_indices = [node_id_to_index[nid] for nid in elem.node_ids if nid in node_id_to_index]
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
                for idx, (node_id, node) in enumerate(sorted_nodes):
                    ux = displacements[idx * 6]
                    uy = displacements[idx * 6 + 1]
                    uz = displacements[idx * 6 + 2]
                    f.write(f"{ux:.6e} {uy:.6e} {uz:.6e}\n")
                
                # 写入位移幅值
                f.write("\nSCALARS DisplacementMagnitude double 1\n")
                f.write("LOOKUP_TABLE default\n")
                for idx in range(len(nodes)):
                    ux = displacements[idx * 6]
                    uy = displacements[idx * 6 + 1]
                    uz = displacements[idx * 6 + 2]
                    mag = math.sqrt(ux**2 + uy**2 + uz**2)
                    f.write(f"{mag:.6e}\n")


def main():
    """主函数"""
    if len(sys.argv) < 2:
        print("用法: python solver.py <input.cdb> [output.vtk]")
        print("示例: python solver.py examples/PIPE288_PLAST.cdb output.vtk")
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
    pipe288_elements = {eid: elem for eid, elem in parser.elements.items() if elem.type == 288}
    if len(pipe288_elements) == 0:
        print("警告: 未找到PIPE288单元，仅输出几何信息")
        VTKWriter.write_polylines(output_file, parser.elements, parser.nodes)
        return
    
    print(f"\n开始有限元分析...")
    
    # 获取材料和截面参数
    mat_id = 1  # 默认材料ID
    if mat_id not in parser.materials:
        print(f"警告: 材料ID {mat_id} 未找到，使用默认值")
        material = {'E': 200e9, 'nu': 0.3, 'density': 7800.0, 'alpha': 1.2e-5, 'T_ref': 25.0}
    else:
        material = parser.materials[mat_id]
        if 'T_ref' not in material:
            material['T_ref'] = parser.reference_temp
    
    sec_id = 1  # 默认截面ID
    if sec_id not in parser.sections:
        print(f"警告: 截面ID {sec_id} 未找到，使用默认值")
        section = {'diameter': 30.0, 'wall_thickness': 1.0}
    else:
        section = parser.sections[sec_id]
    
    print(f"材料参数: E={material.get('E', 0)/1e9:.1f} GPa, nu={material.get('nu', 0):.3f}")
    print(f"截面参数: D={section.get('diameter', 0):.1f} mm, t={section.get('wall_thickness', 0):.1f} mm")
    print(f"载荷: 温度={parser.temperature:.1f}°C, 重力={parser.gravity}, 内压单元数={len(parser.internal_pressure)}")
    
    # 创建PIPE288单元对象
    pipe288_elems = {}
    for elem_id, elem in pipe288_elements.items():
        if len(elem.node_ids) >= 2:
            node1 = parser.nodes[elem.node_ids[0]]
            node2 = parser.nodes[elem.node_ids[1]]
            try:
                pipe288_elems[elem_id] = PIPE288Element(elem_id, node1, node2, material, section)
            except Exception as e:
                print(f"警告: 单元 {elem_id} 创建失败: {e}")
    
    if len(pipe288_elems) == 0:
        print("错误: 无法创建PIPE288单元")
        return
    
    # 创建求解器
    solver = FEASolver(parser.nodes, pipe288_elems, parser.boundary_conditions)
    
    # 组装载荷向量
    load_vector = np.zeros(solver.dof_count)
    sorted_node_ids = sorted(parser.nodes.keys())
    node_id_to_dof_start = {node_id: idx * 6 for idx, node_id in enumerate(sorted_node_ids)}
    
    for elem_id, elem in pipe288_elems.items():
        pressure = parser.internal_pressure.get(elem_id, 0.0)
        F_elem = elem.compute_load_vector(
            pressure, 
            np.array(parser.gravity),
            parser.temperature,
            material.get('T_ref', parser.reference_temp)
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
    print(f"最大位移: UX={max_disp*1000:.3f} mm, UY={max_disp_y*1000:.3f} mm, UZ={max_disp_z*1000:.3f} mm")
    
    # 写入VTK文件
    print(f"\n写入VTK文件: {output_file}")
    VTKWriter.write_polylines(output_file, parser.elements, parser.nodes, displacements)
    print(f"完成！输出文件: {output_file}")
    print(f"\n可以使用ParaView打开查看: paraview {output_file}")


if __name__ == '__main__':
    main()
