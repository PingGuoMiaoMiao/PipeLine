"""
Material模块：材料相关类
"""

from .elastic import ElasticMaterial

try:
    from .plastic_biso import BISOPlasticMaterial
    __all__ = ['ElasticMaterial', 'BISOPlasticMaterial']
except ImportError:
    __all__ = ['ElasticMaterial']

try:
    from .creep_strain_hardening import CreepStrainHardeningMaterial
    __all__.append('CreepStrainHardeningMaterial')
except ImportError:
    pass

