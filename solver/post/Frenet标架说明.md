# Frenet标架在弯管几何展开中的应用

## 概述

在弯管几何展开中引入Frenet标架，确保管道截面的连续性和自然过渡，防止截面翻转和扭曲。

## Frenet标架

Frenet标架由三个正交单位向量组成：

1. **T（切向 - Tangent）**：沿曲线的切线方向
2. **N（法向 - Normal）**：指向曲率中心的方向
3. **B（副法向 - Binormal）**：B = T × N，垂直于T和N

## 实现方法

### 1. 直管（两点）

对于直管，使用两点计算Frenet标架：

- **T**：从节点1到节点2的归一化方向向量
- **N**：构造垂直于T的方向（选择与T垂直的参考方向）
- **B**：T × N

### 2. 弯管（三点）

对于弯管（ELBOW290），使用三点计算Frenet标架：

- **T**：在中间节点处的切向，取两段的加权平均方向
- **N**：指向曲率中心的方向，通过计算两段所在平面的法向量得到
- **B**：T × N

### 3. 截面生成

管道截面点在N-B平面内生成：

```python
# 环向角度
phi = 0 到 2π

# 截面点坐标
r_local = R * (cos(phi) * N + sin(phi) * B)
section_point = node_pos + r_local
```

其中：
- `R`：当前半径（考虑椭圆化：`R = R_base + ovalization_amp * cos(2*phi)`）
- `N`：法向单位向量
- `B`：副法向单位向量

### 4. 防止翻转和扭曲

为确保沿中心线的连续过渡：

1. **方向一致性检查**：检查相邻节点的N向量是否反向
   ```python
   if np.dot(N1, N2) < 0:
       N2 = -N2
       B2 = -B2
   ```

2. **平滑过渡**：为每个节点单独计算Frenet标架，确保连续性

## 优势

1. **几何连续**：使用Frenet标架确保截面沿中心线连续过渡
2. **防止扭曲**：通过方向一致性检查防止截面翻转
3. **自然显示**：弯管几何更加自然和真实
4. **数学严谨**：基于微分几何理论，保证标架的正交性和连续性

## 代码结构

- `solver/post/frenet_frame.py`：Frenet标架计算模块
  - `compute_frenet_frame()`：计算单个Frenet标架
  - `compute_frenet_frames_for_segments()`：为一系列点计算标架
  - `smooth_frenet_frames()`：平滑标架（可选）

- `solver/post/pipe_geometry_expander.py`：几何展开模块
  - `expand_pipe_segment()`：使用Frenet标架展开管段
  - 为每个节点计算独立的Frenet标架
  - 在N-B平面生成截面点

## 使用示例

```python
from solver.post.frenet_frame import FrenetFrame
import numpy as np

# 直管：两点
p1 = np.array([0, 0, 0])
p2 = np.array([1, 0, 0])
T, N, B = FrenetFrame.compute_frenet_frame(p1, p2)

# 弯管：三点
p1 = np.array([0, 0, 0])
p2 = np.array([1, 1, 0])
p3 = np.array([2, 1, 1])
T, N, B = FrenetFrame.compute_frenet_frame(p1, p2, p3)
```

## 注意事项

1. 对于极端几何（如180度弯曲），可能需要额外的平滑处理
2. 椭圆化变形在N-B平面内应用，保持与标架的一致性
3. 确保N和B的正交性，通过重新计算B = T × N来保证

