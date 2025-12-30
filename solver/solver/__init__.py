"""
Solver模块：求解器相关类
"""

from .nonlinear_static import NonlinearStaticSolver
from .newton_raphson import NewtonRaphsonSolver

try:
    from .optimized_creep_solver import OptimizedCreepSolver
    from .sparse_solver import SparseSolver
    from .adaptive_time_step import AdaptiveTimeStepController
    from .parallel_assembly import ParallelAssembly
    __all__ = [
        'NonlinearStaticSolver', 
        'NewtonRaphsonSolver',
        'OptimizedCreepSolver',
        'SparseSolver',
        'AdaptiveTimeStepController',
        'ParallelAssembly'
    ]
except ImportError:
    __all__ = ['NonlinearStaticSolver', 'NewtonRaphsonSolver']

