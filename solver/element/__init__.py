"""
Element模块：单元相关类
"""

from .pipe288_element import PIPE288Element

try:
    from .elbow290_element import ELBOW290Element
    __all__ = ['PIPE288Element', 'ELBOW290Element']
except ImportError:
    __all__ = ['PIPE288Element']

