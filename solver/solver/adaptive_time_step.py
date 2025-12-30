"""
Adaptive Time Step - 自适应时间步控制
根据收敛性和误差自动调整时间步长
"""

import numpy as np
from typing import Dict, Tuple, Optional


class AdaptiveTimeStepController:
    """自适应时间步控制器"""
    
    def __init__(self,
                 initial_time_step: float,
                 min_time_step: float = 1e-6,
                 max_time_step: float = 1e6,
                 safety_factor: float = 0.8,
                 max_increase_factor: float = 1.5,
                 min_decrease_factor: float = 0.5,
                 target_iterations: int = 5):
        """
        初始化自适应时间步控制器
        
        参数:
            initial_time_step: 初始时间步长
            min_time_step: 最小时间步长
            max_time_step: 最大时间步长
            safety_factor: 安全因子（<1.0）
            max_increase_factor: 最大增长因子（>1.0）
            min_decrease_factor: 最小缩减因子（<1.0）
            target_iterations: 目标迭代次数（用于调整时间步）
        """
        self.initial_time_step = initial_time_step
        self.min_time_step = min_time_step
        self.max_time_step = max_time_step
        self.safety_factor = safety_factor
        self.max_increase_factor = max_increase_factor
        self.min_decrease_factor = min_decrease_factor
        self.target_iterations = target_iterations
        
        self.current_time_step = initial_time_step
        self.previous_time_step = initial_time_step
    
    def compute_next_time_step(self,
                              iterations: int,
                              converged: bool,
                              residual_norm: float,
                              tolerance: float) -> float:
        """
        计算下一个时间步长
        
        参数:
            iterations: 当前迭代次数
            converged: 是否收敛
            residual_norm: 残差范数
            tolerance: 收敛容差
        
        返回:
            下一个时间步长
        """
        if not converged:
            # 未收敛，减小时间步
            new_time_step = self.current_time_step * self.min_decrease_factor
            new_time_step = max(new_time_step, self.min_time_step)
            reason = "未收敛"
        elif iterations < self.target_iterations:
            # 收敛且迭代次数少，可以增大时间步
            # 基于迭代次数调整
            factor = self.target_iterations / max(iterations, 1)
            factor = min(factor, self.max_increase_factor)
            new_time_step = self.current_time_step * factor * self.safety_factor
            new_time_step = min(new_time_step, self.max_time_step)
            reason = f"收敛快（{iterations}次迭代）"
        elif iterations > self.target_iterations * 2:
            # 迭代次数多，减小时间步
            factor = self.target_iterations / iterations
            new_time_step = self.current_time_step * factor * self.safety_factor
            new_time_step = max(new_time_step, self.min_time_step)
            reason = f"收敛慢（{iterations}次迭代）"
        else:
            # 正常情况，保持时间步
            new_time_step = self.current_time_step
            reason = "正常"
        
        # 基于残差范数调整（如果残差接近容差，减小时间步）
        if converged and residual_norm > tolerance * 0.1:
            residual_factor = tolerance / residual_norm
            new_time_step = new_time_step * residual_factor * self.safety_factor
            new_time_step = max(new_time_step, self.min_time_step)
        
        self.previous_time_step = self.current_time_step
        self.current_time_step = new_time_step
        
        return new_time_step
    
    def should_retry_step(self, converged: bool) -> bool:
        """判断是否需要重新尝试当前时间步"""
        return not converged
    
    def get_current_time_step(self) -> float:
        """获取当前时间步长"""
        return self.current_time_step
    
    def reset(self, time_step: Optional[float] = None):
        """重置时间步控制器"""
        if time_step is None:
            self.current_time_step = self.initial_time_step
        else:
            self.current_time_step = time_step
        self.previous_time_step = self.current_time_step


class ErrorEstimator:
    """误差估计器（用于自适应时间步）"""
    
    @staticmethod
    def estimate_time_integration_error(displacement_old: np.ndarray,
                                       displacement_new: np.ndarray,
                                       time_step: float) -> float:
        """
        估计时间积分误差
        
        参数:
            displacement_old: 上一时间步的位移
            displacement_new: 当前时间步的位移
            time_step: 时间步长
        
        返回:
            误差估计值
        """
        # 计算位移增量
        delta_displacement = displacement_new - displacement_old
        
        # 误差估计（基于位移增量的范数）
        error = np.linalg.norm(delta_displacement) / max(time_step, 1e-10)
        
        return error
    
    @staticmethod
    def compute_relative_error(value_old: float, value_new: float) -> float:
        """计算相对误差"""
        if abs(value_old) < 1e-10:
            return abs(value_new)
        return abs((value_new - value_old) / value_old)

