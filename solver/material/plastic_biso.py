"""
Plastic BISO - 双线性等向强化弹塑性材料模型
"""

import numpy as np
from typing import Tuple, Optional


class BISOPlasticMaterial:
    """
    双线性等向强化（Bilinear Isotropic Hardening）弹塑性材料模型
    
    材料参数：
    - E: 弹性模量 (Pa)
    - nu: 泊松比
    - yield_stress: 屈服应力 (Pa)
    - tangent_modulus: 切线模量 (Pa)，用于定义硬化斜率
    """
    
    def __init__(self, E: float, nu: float, yield_stress: float, 
                 tangent_modulus: float):
        """
        初始化BISO材料
        
        参数:
            E: 弹性模量 (Pa)
            nu: 泊松比
            yield_stress: 屈服应力 (Pa)
            tangent_modulus: 切线模量 (Pa)
        """
        self.E = E
        self.nu = nu
        self.yield_stress = yield_stress
        
        # 计算硬化参数
        # 弹性模量和切线模量定义硬化斜率
        # Et = dσ/dεp，其中Et为切线模量
        # H = Et * E / (E - Et) 为塑性模量（等向硬化模量）
        if E > tangent_modulus > 0:
            self.H = tangent_modulus * E / (E - tangent_modulus)
        else:
            # 如果切线模量接近弹性模量，使用很小的硬化
            self.H = E * 1e-6
        
        # 计算剪切模量和拉梅常数
        self.G = E / (2.0 * (1.0 + nu))
        self.lambda_lame = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
        
        # 弹性刚度矩阵（3D，用于平面应变/应力）
        self._compute_elastic_stiffness()
    
    def _compute_elastic_stiffness(self):
        """计算弹性刚度矩阵（3D）"""
        # 3D弹性刚度矩阵（使用拉梅常数）
        # 对于各向同性材料，使用拉梅常数和剪切模量
        # C_ijkl = lambda * delta_ij * delta_kl + G * (delta_ik * delta_jl + delta_il * delta_jk)
        
        # 这里我们使用6×6 Voigt记号的刚度矩阵
        # [σ11, σ22, σ33, σ12, σ13, σ23]
        self.C_elastic = np.zeros((6, 6))
        
        # 对角项
        for i in range(3):
            for j in range(3):
                if i == j:
                    self.C_elastic[i, j] = self.lambda_lame + 2 * self.G
                else:
                    self.C_elastic[i, j] = self.lambda_lame
        
        # 剪切项
        for i in range(3, 6):
            self.C_elastic[i, i] = self.G
    
    def compute_stress_from_strain(self, strain: np.ndarray,
                                   plastic_strain: Optional[np.ndarray] = None,
                                   equivalent_plastic_strain: float = 0.0) -> Tuple[np.ndarray, float]:
        """
        从应变计算应力（弹塑性）
        
        参数:
            strain: 总应变 (6×1向量，Voigt记号)
            plastic_strain: 塑性应变 (6×1向量，可选)
            equivalent_plastic_strain: 等效塑性应变（历史变量）
        
        返回:
            (stress, new_equivalent_plastic_strain): 应力和新的等效塑性应变
        """
        if plastic_strain is None:
            plastic_strain = np.zeros(6)
        
        # 弹性应变
        elastic_strain = strain - plastic_strain
        
        # 弹性应力预测
        stress_trial = self.C_elastic @ elastic_strain
        
        # 计算等效Mises应力
        sigma_mises = self._compute_mises_stress(stress_trial)
        
        # 当前屈服应力（考虑硬化）
        current_yield = self.yield_stress + self.H * equivalent_plastic_strain
        
        # 屈服判断
        if sigma_mises <= current_yield:
            # 弹性状态
            return stress_trial, equivalent_plastic_strain
        else:
            # 塑性状态，使用返回映射算法
            return self._return_mapping(stress_trial, strain, plastic_strain,
                                       equivalent_plastic_strain, current_yield)
    
    def _compute_mises_stress(self, stress: np.ndarray) -> float:
        """
        计算Mises等效应力
        
        参数:
            stress: 应力向量 (6×1，Voigt记号: [σ11, σ22, σ33, σ12, σ13, σ23])
        
        返回:
            Mises等效应力
        """
        # 应力偏量
        s = np.zeros(6)
        s[0:3] = stress[0:3] - np.mean(stress[0:3])  # 偏应力
        s[3:6] = stress[3:6]  # 剪切应力
        
        # Mises应力: σ_mises = sqrt(3/2 * s_ij * s_ij)
        mises_sq = (3.0 / 2.0) * (
            s[0]**2 + s[1]**2 + s[2]**2 +
            2.0 * (s[3]**2 + s[4]**2 + s[5]**2)
        )
        
        return np.sqrt(mises_sq)
    
    def _return_mapping(self, stress_trial: np.ndarray, strain: np.ndarray,
                       plastic_strain_old: np.ndarray,
                       equivalent_plastic_strain_old: float,
                       current_yield: float) -> Tuple[np.ndarray, float]:
        """
        返回映射算法（Return Mapping Algorithm）
        
        参数:
            stress_trial: 弹性应力预测
            strain: 总应变
            plastic_strain_old: 上一增量步的塑性应变
            equivalent_plastic_strain_old: 上一增量步的等效塑性应变
            current_yield: 当前屈服应力
        
        返回:
            (stress, new_equivalent_plastic_strain): 更新后的应力和等效塑性应变
        """
        # 计算应力偏量
        mean_stress = np.mean(stress_trial[0:3])
        s_trial = np.zeros(6)
        s_trial[0:3] = stress_trial[0:3] - mean_stress
        s_trial[3:6] = stress_trial[3:6]
        
        # Mises应力
        mises_trial = self._compute_mises_stress(stress_trial)
        
        # 塑性乘子增量
        delta_lambda = (mises_trial - current_yield) / (3 * self.G + self.H)
        
        # 更新等效塑性应变
        new_equivalent_plastic_strain = equivalent_plastic_strain_old + delta_lambda
        
        # 更新应力（径向返回）
        if mises_trial > 1e-10:
            scale_factor = current_yield / mises_trial
        else:
            scale_factor = 0.0
        
        # 更新应力偏量
        s_new = scale_factor * s_trial
        
        # 应力 = 平均应力 + 偏应力
        stress_new = np.zeros(6)
        stress_new[0:3] = mean_stress + s_new[0:3]
        stress_new[3:6] = s_new[3:6]
        
        # 更新塑性应变（增量形式）
        # dεp = delta_lambda * (3/2) * (s / σ_mises)
        if mises_trial > 1e-10:
            plastic_strain_increment = (3.0 / 2.0) * delta_lambda * s_trial / mises_trial
        else:
            plastic_strain_increment = np.zeros(6)
        
        # 注意：这里返回的是应力，塑性应变由调用者管理
        
        return stress_new, new_equivalent_plastic_strain
    
    def compute_tangent_stiffness(self, stress: np.ndarray,
                                 equivalent_plastic_strain: float) -> np.ndarray:
        """
        计算切线刚度矩阵（用于Newton-Raphson迭代）
        
        参数:
            stress: 当前应力
            equivalent_plastic_strain: 当前等效塑性应变
        
        返回:
            切线刚度矩阵 (6×6)
        """
        # 计算Mises应力
        mises = self._compute_mises_stress(stress)
        
        # 当前屈服应力
        current_yield = self.yield_stress + self.H * equivalent_plastic_strain
        
        if mises <= current_yield:
            # 弹性状态
            return self.C_elastic
        else:
            # 塑性状态，计算弹塑性切线刚度矩阵
            return self._compute_elastoplastic_tangent(stress, mises, current_yield)
    
    def _compute_elastoplastic_tangent(self, stress: np.ndarray,
                                      mises: float, current_yield: float) -> np.ndarray:
        """
        计算弹塑性切线刚度矩阵
        
        使用一致切线模量（Consistent Tangent Modulus）
        """
        # 应力偏量
        mean_stress = np.mean(stress[0:3])
        s = np.zeros(6)
        s[0:3] = stress[0:3] - mean_stress
        s[3:6] = stress[3:6]
        
        # 归一化的应力偏量方向
        if mises > 1e-10:
            n = s / mises
        else:
            n = np.zeros(6)
        
        # 一致性切线模量
        # C_ep = C_el - (9*G^2 / (3*G + H)) * (n ⊗ n)
        
        # 标量因子
        factor = 9.0 * self.G**2 / (3.0 * self.G + self.H)
        
        # 外积 n ⊗ n
        n_outer_n = np.outer(n, n)
        
        # 弹塑性切线刚度矩阵
        C_ep = self.C_elastic - factor * n_outer_n
        
        return C_ep
    
    def is_plastic(self, stress: np.ndarray, 
                   equivalent_plastic_strain: float = 0.0) -> bool:
        """
        判断是否进入塑性状态
        
        参数:
            stress: 应力
            equivalent_plastic_strain: 等效塑性应变
        
        返回:
            True如果进入塑性，False如果弹性
        """
        mises = self._compute_mises_stress(stress)
        current_yield = self.yield_stress + self.H * equivalent_plastic_strain
        return mises > current_yield
