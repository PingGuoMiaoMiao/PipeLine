"""
Ovalization Estimator - 椭圆化幅值估计器
用于在VTK后处理阶段根据弯矩等力学量估计椭圆化幅值
椭圆化模态不参与刚度矩阵，仅用于显示
"""

import numpy as np
from typing import Optional


class OvalizationEstimator:
    """椭圆化幅值估计器"""
    
    @staticmethod
    def estimate_ovalization_from_moment(
        M: float,  # 弯矩 (N·m)
        R_curvature: float,  # 弯管曲率半径 (m)
        R_mean: float,  # 管道平均半径 (m)
        t: float,  # 壁厚 (m)
        E: float,  # 弹性模量 (Pa)
        display_amplification: float = 1.0  # 显示放大系数
    ) -> float:
        """
        根据弯矩经验公式估计椭圆化幅值
        
        经验公式（基于薄壁弯管理论）：
        A ≈ C * (M * R) / (E * t^2 * R_mean)
        
        其中：
        - M: 弯矩
        - R: 曲率半径
        - E: 弹性模量
        - t: 壁厚
        - R_mean: 平均半径
        - C: 经验系数（通常为0.1-0.5，取决于弯管几何）
        
        参数:
            M: 弯矩 (N·m)
            R_curvature: 弯管曲率半径 (m)
            R_mean: 管道平均半径 (m)
            t: 壁厚 (m)
            E: 弹性模量 (Pa)
            display_amplification: 显示放大系数（用于可视化增强）
        
        返回:
            椭圆化幅值 (m)
        """
        if abs(M) < 1e-10 or R_curvature < 1e-6 or t < 1e-6:
            return 0.0
        
        # 经验系数（可根据实际测试调整）
        C_empirical = 0.2
        
        # 椭圆化幅值估计
        A = C_empirical * (abs(M) * R_curvature) / (E * t**2 * R_mean)
        
        # 应用显示放大系数
        A_display = A * display_amplification
        
        return A_display
    
    @staticmethod
    def compute_moment_from_displacements(
        node1_disp: np.ndarray,  # 节点1位移（6 DOF）[ux, uy, uz, rotx, roty, rotz]
        node2_disp: np.ndarray,  # 节点2位移（6 DOF）
        L: float,  # 单元长度 (m)
        E: float,  # 弹性模量 (Pa)
        I: float,  # 惯性矩 (m⁴)
        is_bent: bool = False,  # 是否为弯管
        node3_disp: Optional[np.ndarray] = None  # 节点3位移（6 DOF，仅用于弯管，可选）
    ) -> float:
        """
        从节点位移计算等效弯矩
        
        简化处理：使用节点转角估计弯矩
        M ≈ E * I * (θ2 - θ1) / L
        
        参数:
            node1_disp: 节点1位移（前6个DOF）
            node2_disp: 节点2位移（前6个DOF）
            node3_disp: 节点3位移（前6个DOF，弯管使用）
            L: 单元长度 (m)
            E: 弹性模量 (Pa)
            I: 惯性矩 (m⁴)
            is_bent: 是否为弯管
        
        返回:
            等效弯矩幅值 (N·m)
        """
        # 确保L是标量
        if hasattr(L, '__len__') and not isinstance(L, (str, bytes)):
            L = float(L[0]) if len(L) > 0 else 0.0
        else:
            L = float(L)
        
        if L < 1e-10:
            return 0.0
        
        # 提取转角（DOF 3-5为转动：rotx, roty, rotz）
        if len(node1_disp) >= 6 and len(node2_disp) >= 6:
            # 计算节点间转角差
            rot1 = node1_disp[3:6]  # [rotx, roty, rotz]
            rot2 = node2_disp[3:6]
            
            # 对于弯管，也可以考虑节点3
            if is_bent and len(node3_disp) >= 6:
                rot3 = node3_disp[3:6]
                # 使用平均转角差
                delta_rot = (np.linalg.norm(rot2 - rot1) + np.linalg.norm(rot3 - rot2)) / 2.0
            else:
                delta_rot = np.linalg.norm(rot2 - rot1)
            
            # 估计弯矩：M ≈ E * I * Δθ / L
            M = E * I * delta_rot / L
        else:
            M = 0.0
        
        return abs(M)
    
    @staticmethod
    def estimate_ovalization_for_element(
        node_displacements: list,  # 节点位移列表（每个节点6或10 DOF）
        R_curvature: float,  # 弯管曲率半径 (m)
        R_mean: float,  # 管道平均半径 (m)
        t: float,  # 壁厚 (m)
        E: float,  # 弹性模量 (Pa)
        I: float,  # 惯性矩 (m⁴)
        L: float,  # 单元长度 (m)
        display_amplification: float = 1.0  # 显示放大系数
    ) -> float:
        """
        为单元估计椭圆化幅值
        
        参数:
            node_displacements: 节点位移列表
            R_curvature: 弯管曲率半径 (m)
            R_mean: 管道平均半径 (m)
            t: 壁厚 (m)
            E: 弹性模量 (Pa)
            I: 惯性矩 (m⁴)
            L: 单元长度 (m)
            display_amplification: 显示放大系数
        
        返回:
            椭圆化幅值 (m)
        """
        if len(node_displacements) < 2:
            return 0.0
        
        # 提取标准DOF（前6个）
        node1_disp = node_displacements[0][:6] if len(node_displacements[0]) >= 6 else np.zeros(6)
        node2_disp = node_displacements[1][:6] if len(node_displacements[1]) >= 6 else np.zeros(6)
        node3_disp = node_displacements[2][:6] if len(node_displacements) > 2 and len(node_displacements[2]) >= 6 else np.zeros(6)
        
        # 计算弯矩
        is_bent = len(node_displacements) > 2
        node3_disp_arg = node3_disp if is_bent else None
        M = OvalizationEstimator.compute_moment_from_displacements(
            node1_disp, node2_disp, L, E, I, is_bent, node3_disp_arg
        )
        
        # 估计椭圆化幅值
        A = OvalizationEstimator.estimate_ovalization_from_moment(
            M, R_curvature, R_mean, t, E, display_amplification
        )
        
        return A

