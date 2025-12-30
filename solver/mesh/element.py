"""
Element - 单元基类
"""

from typing import List


class Element:
    """单元基类"""
    
    def __init__(self, elem_id: int, elem_type: int, node_ids: List[int]):
        """
        初始化单元
        
        参数:
            elem_id: 单元ID
            elem_type: 单元类型（288=PIPE288, 290=ELBOW290）
            node_ids: 节点ID列表
        """
        self.id = elem_id
        self.type = elem_type  # 288 for PIPE288, 290 for ELBOW290
        self.node_ids = node_ids
    
    def __repr__(self):
        type_name = "PIPE288" if self.type == 288 else "ELBOW290" if self.type == 290 else str(self.type)
        return f"Element({self.id}, {type_name}, nodes={self.node_ids})"
    
    def get_type_name(self):
        """获取单元类型名称"""
        if self.type == 288:
            return "PIPE288"
        elif self.type == 290:
            return "ELBOW290"
        else:
            return f"TYPE{self.type}"

