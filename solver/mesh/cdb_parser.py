"""
CDB Parser - ANSYS CDB文件解析器
"""

from typing import List, Dict
from .node import Node
from .element import Element
from .section import Section


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
        self.boundary_conditions: Dict[int, List[int]] = {}  # 节点ID -> [DOF列表]
        
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
            
            # 解析TB BISO（弹塑性参数）
            elif line.startswith('TB,BISO,'):
                i = self._parse_tb_biso(lines, i)
            
            # 解析TB CREEP（蠕变参数）
            elif line.startswith('TB,CREE'):
                i = self._parse_tb_creep(lines, i)
            
            i += 1
    
    def _parse_et(self, line: str):
        """解析ET命令，建立元素类型映射"""
        parts = [p.strip() for p in line.split(',')]
        if len(parts) >= 3:
            try:
                type_id = int(parts[1])
                unit_type = int(parts[2])
                if unit_type in [288, 290]:
                    self.element_type_map[type_id] = unit_type
            except ValueError:
                pass
    
    def _parse_nblock(self, lines: List[str], start_idx: int) -> int:
        """解析NBLOCK节点块"""
        i = start_idx + 2
        while i < len(lines):
            line = lines[i].strip()
            
            if not line or line.startswith('N,R') or (line.startswith('N,') and 'LOC' in line):
                break
            
            parts = line.split()
            if len(parts) >= 3:
                try:
                    node_id = int(parts[0])
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
        elem_counter = 1
        i = start_idx + 2
        while i < len(lines):
            line = lines[i].strip()
            
            if line == '-1' or (line.startswith('-') and line.strip() == '-1'):
                break
            
            parts = line.split()
            if len(parts) >= 12:
                try:
                    num_nodes = int(parts[8]) if len(parts) > 8 else 2
                    elem_type_id = int(parts[1]) if len(parts) > 1 else 1
                    elem_type = self.element_type_map.get(elem_type_id, 288)
                    
                    # 确定单元ID
                    if len(parts) > 9:
                        try:
                            field9_val = int(parts[9])
                            if field9_val > 0:
                                current_elem_id = field9_val
                            elif len(parts) > 10:
                                current_elem_id = int(parts[10])
                            else:
                                current_elem_id = elem_counter
                        except (ValueError, IndexError):
                            current_elem_id = elem_counter
                    else:
                        current_elem_id = elem_counter
                    
                    # 读取节点ID
                    node_ids = []
                    start_node_idx = 11
                    for j in range(start_node_idx, min(start_node_idx + num_nodes, len(parts))):
                        node_id = int(parts[j])
                        if node_id > 0:
                            node_ids.append(node_id)
                    
                    if node_ids:
                        self.elements[current_elem_id] = Element(current_elem_id, elem_type, node_ids)
                        elem_counter = max(elem_counter, current_elem_id) + 1
                
                except (ValueError, IndexError):
                    pass
            
            i += 1
        
        return i - 1
    
    def _parse_mpdata(self, line: str):
        """解析材料参数 MPDATA"""
        parts = [p.strip() for p in line.split(',')]
        if len(parts) >= 7:
            try:
                mat_id = int(parts[2])
                prop_name = parts[3].strip()
                value = float(parts[6])
                
                if mat_id not in self.materials:
                    self.materials[mat_id] = {}
                
                if prop_name == 'EX':
                    self.materials[mat_id]['E'] = value * 1e6  # MPa -> Pa
                elif prop_name in ['PRXY', 'NUXY']:
                    self.materials[mat_id]['nu'] = value
                elif prop_name == 'DENS':
                    # ANSYS CDB中密度单位为tonne/mm³，转换为kg/m³
                    # 1 tonne/mm³ = 1000 kg / (1e-9 m³) = 1e12 kg/m³
                    self.materials[mat_id]['density'] = value * 1e12  # tonne/mm³ -> kg/m³
                elif prop_name == 'ALPX':
                    self.materials[mat_id]['alpha'] = value
            except (ValueError, IndexError):
                pass
    
    def _parse_reft(self, lines: List[str], start_idx: int):
        """解析参考温度"""
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
        parts = [p.strip() for p in line.split(',')]
        if len(parts) >= 2:
            try:
                section_id = 1
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
        parts = [p.strip() for p in line.split(',')]
        if len(parts) >= 3:
            try:
                self.temperature = float(parts[2])
            except (ValueError, IndexError):
                pass
    
    def _parse_tref(self, line: str):
        """解析参考温度 TREF"""
        pass
    
    def _parse_boundary_condition(self, line: str):
        """解析边界条件 D"""
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
        line = lines[start_idx].strip()
        parts = [p.strip() for p in line.split(',')]
        if len(parts) >= 2 and 'PRES' in line:
            try:
                elem_id = int(parts[1])
                if start_idx + 1 < len(lines):
                    pressure_line = lines[start_idx + 1].strip()
                    pressure_parts = pressure_line.split()
                    if len(pressure_parts) > 0:
                        pressure = float(pressure_parts[0])  # MPa
                        self.internal_pressure[elem_id] = pressure
                        return start_idx + 1
            except (ValueError, IndexError):
                pass
        return start_idx
    
    def _parse_tb_biso(self, lines: List[str], start_idx: int) -> int:
        """解析TB BISO弹塑性参数"""
        # TB,BISO, mat_id, ...
        line = lines[start_idx].strip()
        parts = [p.strip() for p in line.split(',')]
        if len(parts) >= 3:
            try:
                mat_id = int(parts[2])
                # 查找TBDAT行（注意是TBDAT不是TBDATA）
                i = start_idx + 1
                while i < len(lines) and i < start_idx + 10:
                    next_line = lines[i].strip()
                    if next_line.startswith('TBDAT') and not next_line.startswith('TBDATA'):
                        # TBDAT格式：TBDAT, index, yield_stress, tangent_modulus
                        tbdat_parts = [p.strip() for p in next_line.split(',')]
                        if len(tbdat_parts) >= 4:
                            yield_stress = float(tbdat_parts[2]) * 1e6  # MPa -> Pa
                            tangent_modulus = float(tbdat_parts[3]) * 1e6  # MPa -> Pa
                            if mat_id not in self.materials:
                                self.materials[mat_id] = {}
                            self.materials[mat_id]['yield_stress'] = yield_stress
                            self.materials[mat_id]['tangent_modulus'] = tangent_modulus
                            return i
                    elif next_line.startswith('TB,'):
                        # 遇到下一个TB命令，停止
                        break
                    i += 1
            except (ValueError, IndexError):
                pass
        return start_idx
    
    def _parse_tb_creep(self, lines: List[str], start_idx: int) -> int:
        """解析TB CREEP蠕变参数"""
        # TB,CREE, mat_id, ...
        line = lines[start_idx].strip()
        parts = [p.strip() for p in line.split(',')]
        if len(parts) >= 3:
            try:
                mat_id = int(parts[2])
                if mat_id not in self.materials:
                    self.materials[mat_id] = {}
                
                # 查找TBDATA行（可能有多个）
                i = start_idx + 1
                C1 = C2 = C3 = C4 = None
                while i < len(lines) and i < start_idx + 20:
                    next_line = lines[i].strip()
                    if next_line.startswith('TBDATA,'):
                        tbdata_parts = [p.strip() for p in next_line.split(',')]
                        if len(tbdata_parts) >= 2:
                            try:
                                index = int(tbdata_parts[1])
                                # TBDATA,1,C1,C2,C3 或 TBDATA,4,C4
                                if index == 1 and len(tbdata_parts) >= 5:
                                    C1 = float(tbdata_parts[2])
                                    C2 = float(tbdata_parts[3])
                                    C3 = float(tbdata_parts[4])
                                elif index == 4 and len(tbdata_parts) >= 3:
                                    C4 = float(tbdata_parts[2])
                            except (ValueError, IndexError):
                                pass
                    elif next_line.startswith('TB,'):
                        # 遇到下一个TB命令，停止
                        break
                    i += 1
                
                # 保存蠕变参数
                if C1 is not None:
                    self.materials[mat_id]['C1'] = C1
                if C2 is not None:
                    self.materials[mat_id]['C2'] = C2
                if C3 is not None:
                    self.materials[mat_id]['C3'] = C3
                if C4 is not None:
                    self.materials[mat_id]['C4'] = C4
                else:
                    # 如果没有找到C4，默认为0
                    self.materials[mat_id]['C4'] = 0.0
                
                return i - 1
            except (ValueError, IndexError):
                pass
        return start_idx

