"""
Parallel Transport Frame - 平行传输标架计算模块
用于计算曲线上的平行传输标架（Parallel Transport Frame）
通过最小旋转从前一截面传播法向，避免截面翻转和卷曲
"""

import numpy as np
from typing import Tuple, List, Optional


class ParallelTransportFrame:
    """
    平行传输标架计算类
    
    原理：
    1. 第一个截面法向手动初始化
    2. 后续截面法向通过最小旋转从前一截面传播
    3. 旋转轴 = T_prev × T_cur
    4. 旋转角度 = arccos(T_prev · T_cur)
    """
    
    @staticmethod
    def normalize(v: np.ndarray) -> np.ndarray:
        """归一化向量"""
        norm = np.linalg.norm(v)
        if norm < 1e-10:
            return np.zeros_like(v)
        return v / norm
    
    @staticmethod
    def rotate_vector_around_axis(v: np.ndarray, axis: np.ndarray, angle: float) -> np.ndarray:
        """
        绕轴旋转向量（Rodrigues旋转公式）
        
        参数:
            v: 要旋转的向量
            axis: 旋转轴（需归一化）
            angle: 旋转角度（弧度）
        
        返回:
            旋转后的向量
        """
        if np.linalg.norm(axis) < 1e-10:
            return v
        
        axis = ParallelTransportFrame.normalize(axis)
        
        # Rodrigues旋转公式：v' = v*cos(θ) + (axis×v)*sin(θ) + axis*(axis·v)*(1-cos(θ))
        cos_angle = np.cos(angle)
        sin_angle = np.sin(angle)
        
        v_rot = (v * cos_angle + 
                np.cross(axis, v) * sin_angle + 
                axis * np.dot(axis, v) * (1 - cos_angle))
        
        return v_rot
    
    @staticmethod
    def transport_normal(
        T_prev: np.ndarray,  # 上一段切向
        T_cur: np.ndarray,   # 当前切向
        N_prev: np.ndarray   # 上一段法向
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        平行传输法向（Parallel Transport）
        
        步骤：
        1. 计算旋转轴 = cross(T_prev, T_cur)
        2. 若|axis|很小：N_cur = N_prev
        3. 否则：将N_prev绕axis旋转angle
        4. 归一化N_cur
        5. B_cur = cross(T_cur, N_cur)
        
        参数:
            T_prev: 上一段切向（归一化）
            T_cur: 当前切向（归一化）
            N_prev: 上一段法向（归一化）
        
        返回:
            (N_cur, B_cur): 当前法向和副法向（归一化）
        """
        # 归一化切向
        T_prev = ParallelTransportFrame.normalize(T_prev)
        T_cur = ParallelTransportFrame.normalize(T_cur)
        N_prev = ParallelTransportFrame.normalize(N_prev)
        
        # 计算旋转轴 = T_prev × T_cur
        axis = np.cross(T_prev, T_cur)
        axis_norm = np.linalg.norm(axis)
        
        if axis_norm < 1e-6:
            # 如果切向几乎平行，法向保持不变
            N_cur = N_prev.copy()
        else:
            # 归一化旋转轴
            axis = axis / axis_norm
            
            # 计算旋转角度
            # cos(angle) = T_prev · T_cur
            cos_angle = np.clip(np.dot(T_prev, T_cur), -1.0, 1.0)
            angle = np.arccos(cos_angle)
            
            # 将N_prev绕axis旋转angle
            N_cur = ParallelTransportFrame.rotate_vector_around_axis(N_prev, axis, angle)
        
        # 归一化法向
        N_cur = ParallelTransportFrame.normalize(N_cur)
        
        # 确保N_cur与T_cur垂直
        # 投影到垂直于T_cur的平面
        N_cur = N_cur - np.dot(N_cur, T_cur) * T_cur
        N_cur = ParallelTransportFrame.normalize(N_cur)
        
        # 如果N_cur为零，选择一个与T_cur垂直的方向
        if np.linalg.norm(N_cur) < 1e-10:
            # 选择一个与T_cur垂直的参考方向
            if abs(T_cur[2]) < 0.9:
                ref = np.array([0, 0, 1])
            else:
                ref = np.array([1, 0, 0])
            N_cur = ref - np.dot(ref, T_cur) * T_cur
            N_cur = ParallelTransportFrame.normalize(N_cur)
        
        # 计算副法向 B_cur = T_cur × N_cur
        B_cur = np.cross(T_cur, N_cur)
        B_cur = ParallelTransportFrame.normalize(B_cur)
        
        # 确保N和B正交
        N_cur = np.cross(B_cur, T_cur)
        N_cur = ParallelTransportFrame.normalize(N_cur)
        
        return N_cur, B_cur
    
    @staticmethod
    def compute_initial_frame(
        T: np.ndarray,  # 切向
        reference_direction: Optional[np.ndarray] = None  # 参考方向（用于初始化法向）
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        计算初始标架（第一个截面的标架）
        
        参数:
            T: 切向（归一化）
            reference_direction: 参考方向（可选），用于初始化法向
        
        返回:
            (T, N, B): 切向、法向、副法向（归一化）
        """
        T = ParallelTransportFrame.normalize(T)
        
        # 初始化法向
        if reference_direction is not None:
            ref = ParallelTransportFrame.normalize(reference_direction)
            # 投影到垂直于T的平面
            N = ref - np.dot(ref, T) * T
            N = ParallelTransportFrame.normalize(N)
            if np.linalg.norm(N) < 1e-10:
                # 如果ref与T平行，选择默认方向
                if abs(T[2]) < 0.9:
                    ref = np.array([0, 0, 1])
                else:
                    ref = np.array([1, 0, 0])
                N = ref - np.dot(ref, T) * T
                N = ParallelTransportFrame.normalize(N)
        else:
            # 默认：选择与T垂直的方向
            if abs(T[2]) < 0.9:
                ref = np.array([0, 0, 1])
            else:
                ref = np.array([1, 0, 0])
            N = ref - np.dot(ref, T) * T
            N = ParallelTransportFrame.normalize(N)
        
        # 计算副法向 B = T × N
        B = np.cross(T, N)
        B = ParallelTransportFrame.normalize(B)
        
        # 确保N和B正交
        N = np.cross(B, T)
        N = ParallelTransportFrame.normalize(N)
        
        return T, N, B
    
    @staticmethod
    def compute_frames_for_curve(
        points: List[np.ndarray],  # 点列表
        initial_reference: Optional[np.ndarray] = None  # 初始参考方向
    ) -> List[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """
        为一系列点计算平行传输标架
        
        参数:
            points: 点列表（按顺序）
            initial_reference: 初始参考方向（用于第一个截面的法向初始化）
        
        返回:
            frames: 标架列表，每个元素为(T, N, B)
        """
        n = len(points)
        if n < 2:
            return []
        
        frames = []
        
        # 计算切向列表
        tangents = []
        for i in range(n - 1):
            v = points[i + 1] - points[i]
            T = ParallelTransportFrame.normalize(v)
            tangents.append(T)
        # 最后一个点的切向与前一个相同
        if len(tangents) > 0:
            tangents.append(tangents[-1])
        else:
            tangents.append(np.array([1, 0, 0]))  # 默认切向
        
        # 第一个点的标架（初始化）
        T0, N0, B0 = ParallelTransportFrame.compute_initial_frame(
            tangents[0], initial_reference
        )
        frames.append((T0, N0, B0))
        
        # 后续点的标架（平行传输）
        for i in range(1, n):
            T_prev = frames[i - 1][0]  # 上一段切向（使用上一段的切向）
            T_cur = tangents[i]  # 当前切向
            N_prev = frames[i - 1][1]  # 上一段法向
            
            # 平行传输法向
            N_cur, B_cur = ParallelTransportFrame.transport_normal(T_prev, T_cur, N_prev)
            
            frames.append((T_cur, N_cur, B_cur))
        
        return frames

