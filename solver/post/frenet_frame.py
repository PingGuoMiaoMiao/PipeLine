"""
Frenet Frame - Frenet标架计算模块
用于计算曲线上的Frenet标架（切向T、法向N、副法向B）
"""

import numpy as np
from typing import Tuple, List, Optional


class FrenetFrame:
    """Frenet标架计算类"""
    
    @staticmethod
    def compute_frenet_frame(
        p1: np.ndarray,  # 点1坐标 (m 或 mm)
        p2: np.ndarray,  # 点2坐标 (m 或 mm)
        p3: Optional[np.ndarray] = None  # 点3坐标（用于弯管，可选）
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        计算Frenet标架（T, N, B）
        
        对于两点（直管）：
        - T: 从p1到p2的归一化方向
        - N: 构造垂直方向
        - B: T × N
        
        对于三点（弯管）：
        - T: 在p2处的切向（p1->p2和p2->p3的平均方向）
        - N: 指向曲率中心的方向（垂直于T）
        - B: T × N
        
        参数:
            p1: 点1坐标
            p2: 点2坐标
            p3: 点3坐标（可选，用于弯管）
        
        返回:
            (T, N, B): 切向、法向、副法向单位向量
        """
        if p3 is None:
            # 直管：两点情况
            v = p2 - p1
            T = v / (np.linalg.norm(v) + 1e-10)
            
            # 构造法向N（选择与T垂直的方向）
            # 如果T接近[0,0,1]，使用[1,0,0]作为参考
            if abs(T[2]) < 0.9:
                ref = np.array([0, 0, 1])
            else:
                ref = np.array([1, 0, 0])
            
            N = np.cross(T, ref)
            N_norm = np.linalg.norm(N)
            if N_norm < 1e-10:
                # 如果N为零，选择另一个参考方向
                ref = np.array([0, 1, 0])
                N = np.cross(T, ref)
                N_norm = np.linalg.norm(N)
            
            if N_norm > 1e-10:
                N = N / N_norm
            else:
                N = np.array([1, 0, 0])
            
            # 副法向
            B = np.cross(T, N)
            B_norm = np.linalg.norm(B)
            if B_norm > 1e-10:
                B = B / B_norm
            else:
                B = np.array([0, 1, 0])
            
            # 确保N和B正交
            N = np.cross(B, T)
            N_norm = np.linalg.norm(N)
            if N_norm > 1e-10:
                N = N / N_norm
            else:
                N = np.array([1, 0, 0])
        
        else:
            # 弯管：三点情况
            v1 = p2 - p1  # 第一段方向
            v2 = p3 - p2  # 第二段方向
            
            L1 = np.linalg.norm(v1)
            L2 = np.linalg.norm(v2)
            
            if L1 < 1e-10 or L2 < 1e-10:
                # 如果某段长度为0，退化为两点情况
                return FrenetFrame.compute_frenet_frame(p1, p3)
            
            # 归一化
            v1_norm = v1 / L1
            v2_norm = v2 / L2
            
            # 切向T：两段的平均方向（加权平均，考虑长度）
            T = (v1_norm + v2_norm) / 2.0
            T_norm = np.linalg.norm(T)
            if T_norm > 1e-10:
                T = T / T_norm
            else:
                T = v1_norm  # 退化情况
            
            # 法向N：指向曲率中心的方向
            # 计算垂直于两段的向量（指向曲率中心）
            # N应该垂直于T，并且指向曲率中心
            # 对于弯管，曲率中心在v1和v2的角平分线的垂直方向上
            
            # 计算角平分线方向
            bisector = (v1_norm + v2_norm)
            bisector_norm = np.linalg.norm(bisector)
            if bisector_norm > 1e-10:
                bisector = bisector / bisector_norm
            else:
                # 如果v1和v2反向，使用垂直方向
                bisector = np.cross(v1_norm, np.array([0, 0, 1]))
                bisector_norm = np.linalg.norm(bisector)
                if bisector_norm > 1e-10:
                    bisector = bisector / bisector_norm
                else:
                    bisector = np.array([1, 0, 0])
            
            # N应该垂直于T，并且在v1和v2所在的平面内
            # N = normalize(cross(T, cross(v1_norm, v2_norm)))
            plane_normal = np.cross(v1_norm, v2_norm)
            plane_norm = np.linalg.norm(plane_normal)
            
            if plane_norm > 1e-10:
                plane_normal = plane_normal / plane_norm
                # N垂直于T且在平面内
                N = np.cross(T, plane_normal)
                N_norm = np.linalg.norm(N)
                if N_norm > 1e-10:
                    N = N / N_norm
                else:
                    # 退化情况：T和plane_normal平行
                    N = np.cross(T, np.array([0, 0, 1]))
                    N_norm = np.linalg.norm(N)
                    if N_norm > 1e-10:
                        N = N / N_norm
                    else:
                        N = np.array([1, 0, 0])
            else:
                # v1和v2平行，退化为直管
                return FrenetFrame.compute_frenet_frame(p1, p3)
            
            # 副法向
            B = np.cross(T, N)
            B_norm = np.linalg.norm(B)
            if B_norm > 1e-10:
                B = B / B_norm
            else:
                B = np.array([0, 1, 0])
            
            # 确保正交性
            N = np.cross(B, T)
            N_norm = np.linalg.norm(N)
            if N_norm > 1e-10:
                N = N / N_norm
            else:
                N = np.array([1, 0, 0])
        
        return T, N, B
    
    @staticmethod
    def compute_frenet_frames_for_segments(
        points: List[np.ndarray],  # 点列表
        is_closed: bool = False  # 是否为封闭曲线
    ) -> List[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """
        为一系列点计算Frenet标架
        
        参数:
            points: 点列表
            is_closed: 是否为封闭曲线
        
        返回:
            frenet_frames: Frenet标架列表，每个元素为(T, N, B)
        """
        n = len(points)
        if n < 2:
            return []
        
        frenet_frames = []
        
        if n == 2:
            # 只有两点，使用两点方法
            T, N, B = FrenetFrame.compute_frenet_frame(points[0], points[1])
            frenet_frames.append((T, N, B))
            frenet_frames.append((T, N, B))  # 两点使用相同的标架
        else:
            # 多点情况
            for i in range(n):
                if i == 0:
                    # 第一个点：使用前两点
                    T, N, B = FrenetFrame.compute_frenet_frame(points[0], points[1])
                elif i == n - 1:
                    # 最后一个点：使用后两点
                    T, N, B = FrenetFrame.compute_frenet_frame(points[n-2], points[n-1])
                else:
                    # 中间点：使用三点方法
                    T, N, B = FrenetFrame.compute_frenet_frame(
                        points[i-1], points[i], points[i+1]
                    )
                
                frenet_frames.append((T, N, B))
        
        return frenet_frames
    
    @staticmethod
    def smooth_frenet_frames(
        frenet_frames: List[Tuple[np.ndarray, np.ndarray, np.ndarray]],
        smoothing_factor: float = 0.5
    ) -> List[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """
        平滑Frenet标架，防止翻转和扭曲
        
        使用相邻标架的平均来平滑，确保连续过渡
        
        参数:
            frenet_frames: 原始Frenet标架列表
            smoothing_factor: 平滑因子（0-1），0表示不平滑，1表示完全平均
        
        返回:
            smoothed_frames: 平滑后的Frenet标架列表
        """
        n = len(frenet_frames)
        if n <= 1:
            return frenet_frames
        
        smoothed_frames = []
        
        for i in range(n):
            T, N, B = frenet_frames[i]
            
            if smoothing_factor > 0 and 0 < i < n - 1:
                # 中间点：与相邻标架平均
                T_prev, N_prev, B_prev = frenet_frames[i-1]
                T_next, N_next, B_next = frenet_frames[i+1]
                
                # 平滑切向
                T_smooth = T * (1 - smoothing_factor) + (T_prev + T_next) / 2.0 * smoothing_factor
                T_smooth = T_smooth / (np.linalg.norm(T_smooth) + 1e-10)
                
                # 平滑法向和副法向
                N_smooth = N * (1 - smoothing_factor) + (N_prev + N_next) / 2.0 * smoothing_factor
                N_smooth = N_smooth / (np.linalg.norm(N_smooth) + 1e-10)
                
                # 重新计算副法向以确保正交
                B_smooth = np.cross(T_smooth, N_smooth)
                B_norm = np.linalg.norm(B_smooth)
                if B_norm > 1e-10:
                    B_smooth = B_smooth / B_norm
                else:
                    B_smooth = B
                
                # 重新计算法向以确保正交
                N_smooth = np.cross(B_smooth, T_smooth)
                N_norm = np.linalg.norm(N_smooth)
                if N_norm > 1e-10:
                    N_smooth = N_smooth / N_norm
                else:
                    N_smooth = N
                
                smoothed_frames.append((T_smooth, N_smooth, B_smooth))
            else:
                smoothed_frames.append((T, N, B))
        
        return smoothed_frames

