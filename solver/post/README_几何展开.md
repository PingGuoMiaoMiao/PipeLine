# 管单元几何展开模块说明

## 概述

本模块实现了**力学求解与几何显示的完全解耦**：

- **求解阶段**：仅计算管中心线的节点位移、转角、椭圆化模态幅值
- **后处理阶段**：独立模块根据中心线结果生成三维管壳曲面

## 设计原则

1. **完全解耦**：求解器不生成任何三维网格，只输出中心线信息
2. **独立后处理**：几何展开模块仅基于节点坐标、位移、截面参数工作
3. **保持精度**：不改变原有求解结果，只改变可视化方式

## 模块结构

### `pipe_geometry_expander.py`

核心几何展开模块，包含：

- `PipeGeometryExpander.expand_pipe_segment()`: 展开单个管段
- `PipeGeometryExpander.expand_pipe_centerline_to_surface()`: 展开整个系统

### `vtk_writer.py`

VTK输出模块，已重构为使用几何展开模块：

- `VTKWriter.write_polylines()`: 统一接口，支持表面或线条输出
- `VTKWriter._write_surface_from_centerline()`: 使用几何展开模块生成表面

## 使用方法

### 基本用法

```python
from solver.post.vtk_writer import VTKWriter

# 求解完成后
VTKWriter.write_polylines(
    output_file,
    elements,  # 单元字典
    nodes,     # 节点字典
    displacements,  # 位移数组（m）
    sections=parser.sections,  # 截面参数
    output_surface=True  # 输出为表面
)
```

### 支持椭圆化（ELBOW290）

```python
# 对于ELBOW290，可以从位移数组提取椭圆化DOF
ovalization_dofs = {}  # {node_id: [oval_cos, oval_sin, ...]}
# 或者让模块自动从位移数组提取（如果dof_per_node=10）

VTKWriter.write_polylines(
    output_file,
    elements, nodes, displacements,
    sections=parser.sections,
    ovalization_dofs=ovalization_dofs,  # 可选
    output_surface=True
)
```

## 椭圆化变形

椭圆化变形公式：`r(θ) = r0 + a * cos(2θ)`

其中：
- `r0`: 原始半径
- `a`: 椭圆化幅值 = sqrt(oval_cos² + oval_sin²)
- `θ`: 环向角度

对于ELBOW290单元：
- DOF 6-7: 椭圆化的cos和sin项
- 模块会自动计算椭圆化幅值并应用

## 输出格式

输出为VTK PolyData格式：
- `POINTS`: 三维点坐标（mm）
- `POLYGONS`: 四边形面片（用于显示管壳曲面）

可在ParaView中：
- 查看三维管壳几何（而非中心线）
- 显示位移、应力等云图
- 查看椭圆化变形效果

## 优势

1. **计算效率**：求解器只处理1D问题，不生成3D网格
2. **内存效率**：存储中心线数据而非完整3D网格
3. **灵活性**：可以调整环向划分点数而不影响求解
4. **可维护性**：求解和可视化完全分离

