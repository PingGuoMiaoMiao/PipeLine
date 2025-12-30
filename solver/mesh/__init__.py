"""
Mesh模块：网格相关类
"""

from .node import Node
from .element import Element
from .section import Section

try:
    from .cdb_parser import CDBParser
    __all__ = ['Node', 'Element', 'Section', 'CDBParser']
except ImportError:
    __all__ = ['Node', 'Element', 'Section']

