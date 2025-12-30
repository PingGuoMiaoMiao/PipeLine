# ELBOW290 弯管单元实现说明

## 1. 数学模型

### 1.1 弯管几何

ELBOW290是3节点弯管单元，用于描述管道系统的弯管段。

**几何参数：**
- 曲率半径 R：弯管的曲率半径（从中心线到曲率中心的距离）
- 弯曲角度 θ：弯管的弯曲角度
- 外径 D：管道外径
- 壁厚 t：管道壁厚

**节点定义：**
- 节点1：起始节点
- 节点2：中间节点（用于定义弯曲）
- 节点3：结束节点

### 1.2 自由度定义

每个节点有**10个自由度**：

1. **标准自由度（6个）**：
   - UX, UY, UZ：平移自由度
   - ROTX, ROTY, ROTZ：转动自由度

2. **截面变形自由度（4个）**：
   - δ_2c：椭圆化cos(2φ)项系数
   - δ_2s：椭圆化sin(2φ)项系数
   - δ_3c：翘曲cos(3φ)项系数
   - δ_3s：翘曲sin(3φ)项系数

### 1.3 环向傅里叶级数形函数

环向位移场用傅里叶级数表示：

\[
\Delta r(\phi) = \delta_2^c \cos(2\phi) + \delta_2^s \sin(2\phi) + 
                \delta_3^c \cos(3\phi) + \delta_3^s \sin(3\phi)
\]

其中：
- φ：环向角度（0 ≤ φ < 2π）
- δ_2^c, δ_2^s：椭圆化系数（n=2项）
- δ_3^c, δ_3^s：翘曲系数（n=3项）

变形后的截面坐标：

\[
x(\phi) = (R_{mean} + \Delta r(\phi)) \cos(\phi)
\]
\[
y(\phi) = (R_{mean} + \Delta r(\phi)) \sin(\phi)
\]

其中R_mean是平均半径。

### 1.4 刚度矩阵

弯管单元的刚度矩阵（30×30）包括：

1. **标准DOF刚度**（类似PIPE288，但考虑曲率修正）
   - 轴向刚度
   - 弯曲刚度
   - 扭转刚度

2. **椭圆化刚度**：
   \[
   K_{oval} \approx \frac{E t^3}{R^2}
   \]
   其中E为弹性模量，t为壁厚，R为曲率半径

3. **翘曲刚度**：
   \[
   K_{warp} \approx 0.1 \frac{E t^3}{R^3}
   \]
   翘曲刚度通常比椭圆化刚度小

## 2. 关键类设计

### 2.1 ELBOW290Element类

**位置：** `solver/element/elbow290_element.py`

**主要方法：**
- `__init__()`: 初始化弯管几何
- `_compute_bend_geometry()`: 计算曲率半径和角度
- `_compute_local_coordinate_system()`: 计算局部坐标系
- `compute_stiffness_matrix()`: 计算单元刚度矩阵（30×30）
- `compute_load_vector()`: 计算单元载荷向量（30×1）
- `get_section_points()`: 生成截面点坐标（用于VTK输出）

**特点：**
- 3节点单元，每个节点10个自由度
- 支持椭圆化和翘曲变形
- 刚度矩阵包括曲率修正

### 2.2 FourierShapeFunction类

**位置：** `solver/element/shape_function_fourier.py`

**主要方法：**
- `compute_shape_function()`: 计算n阶傅里叶形函数值
- `compute_radial_displacement()`: 计算环向径向位移
- `compute_section_coordinates()`: 计算变形后的截面坐标

### 2.3 PipeSurfaceExpander类

**位置：** `solver/post/pipe_surface_expand.py`

**主要方法：**
- `expand_elbow290_to_surface()`: 将ELBOW290单元扩展为三维曲面

**功能：**
- 根据椭圆化和翘曲自由度生成变形后的截面
- 连接相邻截面形成四边形面片
- 输出VTK POLYGON格式

## 3. 实现细节

### 3.1 几何计算

曲率半径通过节点坐标估算：
- 计算两个向量（节点1→节点2，节点2→节点3）
- 计算向量夹角
- 估算曲率半径：R ≈ L / θ

### 3.2 刚度矩阵组装

采用简化策略：
1. 将弯管视为两段直管的组合
2. 添加椭圆化和翘曲的简化刚度项
3. 考虑曲率修正

### 3.3 VTK输出

**曲面生成步骤：**
1. 为每个节点生成截面点（考虑椭圆化和翘曲）
2. 连接相邻节点的截面形成四边形面片
3. 输出为VTK POLYGON格式

**环向划分：**
- 默认20个点（可配置）
- 角度范围：0 到 2π

## 4. 使用示例

### 4.1 基本用法

```python
from solver.mesh.cdb_parser import CDBParser
from solver.element.elbow290_element import ELBOW290Element
from solver.solver.nonlinear_static import NonlinearStaticSolver
from solver.post.vtk_writer import VTKWriter

# 解析CDB文件
parser = CDBParser("ELBOW290.cdb")
parser.parse()

# 创建ELBOW290单元
elbow290_elem = ELBOW290Element(
    elem_id=1,
    node1=parser.nodes[1],
    node2=parser.nodes[2],
    node3=parser.nodes[3],
    material=material,
    section=section
)

# 求解
solver = NonlinearStaticSolver(...)
displacements = solver.solve(load_vector)

# 输出VTK
VTKWriter.write_polylines(
    "output.vtk",
    parser.elements,
    parser.nodes,
    displacements,
    elbow290_elements={1: elbow290_elem}
)
```

### 4.2 运行程序

```bash
python -m solver.main examples/ELBOW290_PLAST.cdb output_elbow290.vtk
```

## 5. 验证说明

### 5.1 验证要点

1. **几何正确性**：
   - 曲率半径计算正确
   - 节点顺序正确

2. **变形形态**：
   - 椭圆化变形趋势与ANSYS一致
   - 翘曲变形合理

3. **数值精度**：
   - 位移量级合理
   - 边界条件正确应用

### 5.2 可视化检查

在ParaView中打开输出的VTK文件，检查：

1. **几何形状**：
   - 弯管形状正确
   - 截面圆形（未变形时）

2. **椭圆化趋势**：
   - 施加弯矩后，截面应呈现椭圆化
   - 椭圆化方向与载荷方向一致

3. **位移云图**：
   - 位移分布合理
   - 最大值位置正确

## 6. 简化说明

当前实现采用了以下简化：

1. **线性化**：采用线性刚度矩阵，忽略几何非线性
2. **简化刚度**：椭圆化和翘曲刚度使用简化公式
3. **均匀截面**：假设截面沿弧长不变
4. **曲率估算**：使用节点坐标估算曲率半径（非精确值）

这些简化对于初步实现和验证是合理的，后续可以逐步完善。

## 7. 后续扩展方向

1. **精确曲率计算**：从CDB文件中读取或精确计算曲率半径
2. **非线性分析**：考虑几何非线性和材料非线性
3. **更精确的刚度**：使用更精确的弯管理论计算刚度矩阵
4. **高阶项**：支持n≥4的高阶傅里叶项
5. **数值积分**：实现环向和轴向的数值积分

