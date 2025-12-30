# ELBOW290 弯管单元实现总结

## 已完成的工作

### 1. 核心模块实现

#### 1.1 ELBOW290Element类
**文件：** `solver/element/elbow290_element.py`

**功能：**
- ✅ 3节点弯管单元定义
- ✅ 弯管几何计算（曲率半径、角度）
- ✅ 局部坐标系计算
- ✅ 30×30刚度矩阵计算（包含椭圆化和翘曲刚度）
- ✅ 载荷向量计算（内压、重力、热膨胀）
- ✅ 截面点坐标生成（用于VTK输出）

**特点：**
- 每个节点10个自由度（6个标准DOF + 4个截面变形DOF）
- 简化的椭圆化和翘曲刚度计算
- 基于PIPE288的刚度矩阵，考虑曲率修正

#### 1.2 FourierShapeFunction类
**文件：** `solver/element/shape_function_fourier.py`

**功能：**
- ✅ 傅里叶级数形函数计算
- ✅ 环向径向位移计算（椭圆化+翘曲）
- ✅ 变形后截面坐标计算
- ✅ 环向角度点生成

#### 1.3 PipeSurfaceExpander类
**文件：** `solver/post/pipe_surface_expand.py`

**功能：**
- ✅ 将ELBOW290单元扩展为三维曲面
- ✅ 四边形面片生成
- ✅ 支持椭圆化和翘曲变形的可视化

#### 1.4 VTKWriter扩展
**文件：** `solver/post/vtk_writer.py`

**功能：**
- ✅ 支持ELBOW290曲面输出（POLYGON格式）
- ✅ 自动检测单元类型并选择输出格式
- ✅ 保留PIPE288的PolyLine输出

#### 1.5 NonlinearStaticSolver更新
**文件：** `solver/solver/nonlinear_static.py`

**功能：**
- ✅ 支持不同DOF数的单元（6或10 DOF/节点）
- ✅ 自动识别单元类型并正确组装
- ✅ 边界条件处理（只约束标准DOF）

### 2. 主程序

**文件：** `solver/main_elbow290.py`

**功能：**
- ✅ 支持PIPE288和ELBOW290混合分析
- ✅ 自动选择DOF数（优先使用ELBOW290的10 DOF/节点）
- ✅ 输出椭圆化和翘曲统计信息
- ✅ 自动选择VTK输出格式（PolyLine或曲面）

## 使用方法

### 基本运行

```bash
# 运行ELBOW290分析
python -m solver.main_elbow290 examples/ELBOW290_PLAST.cdb output_elbow290.vtk

# 运行PIPE288分析（也可以使用）
python -m solver.main_elbow290 examples/PIPE288_PLAST.cdb output_pipe288.vtk
```

### 代码中使用

```python
from solver.element.elbow290_element import ELBOW290Element
from solver.material.elastic import ElasticMaterial

# 创建材料
material = ElasticMaterial(1, E=200e9, nu=0.3, density=7800.0,
                          alpha=1.2e-5, T_ref=25.0)

# 创建ELBOW290单元
section = {'diameter': 90.0, 'wall_thickness': 2.0}
elem = ELBOW290Element(1, node1, node2, node3, material, section)

# 计算刚度矩阵
K = elem.compute_stiffness_matrix()  # 30×30

# 计算载荷向量
F = elem.compute_load_vector(pressure=1e6, gravity=[0,0,9800],
                             temperature=350, ref_temperature=25)  # 30×1
```

## 数学模型

### 自由度定义

每个节点10个自由度：
1. **标准DOF（0-5）**：UX, UY, UZ, ROTX, ROTY, ROTZ
2. **椭圆化DOF（6-7）**：δ_2c (cos(2φ)), δ_2s (sin(2φ))
3. **翘曲DOF（8-9）**：δ_3c (cos(3φ)), δ_3s (sin(3φ))

### 截面变形公式

环向径向位移：
\[
\Delta r(\phi) = \delta_2^c \cos(2\phi) + \delta_2^s \sin(2\phi) + 
                \delta_3^c \cos(3\phi) + \delta_3^s \sin(3\phi)
\]

### 刚度矩阵

- **标准DOF刚度**：基于PIPE288，考虑曲率修正
- **椭圆化刚度**：\( K_{oval} \approx \frac{E t^3}{R^2} \)
- **翘曲刚度**：\( K_{warp} \approx 0.1 \frac{E t^3}{R^3} \)

## 输出格式

### VTK文件

ELBOW290单元输出为**POLYGON格式**，包含：
- 三维曲面点
- 四边形面片（连接相邻截面）
- 单元类型标记

### 在ParaView中查看

1. 打开ParaView
2. 加载VTK文件
3. 应用"Surface"或"Surface with Edges"过滤器
4. 查看椭圆化变形（如果存在）

## 简化说明

当前实现采用了以下简化：

1. **线性化**：采用线性刚度矩阵，忽略几何非线性
2. **简化刚度**：椭圆化和翘曲刚度使用简化公式
3. **均匀截面**：假设截面沿弧长不变
4. **曲率估算**：使用节点坐标估算曲率半径（非精确值）
5. **忽略耦合**：椭圆化和翘曲DOF之间的耦合项被忽略

这些简化对于初步实现和验证是合理的，后续可以逐步完善。

## 验证建议

1. **几何验证**：
   - 检查曲率半径计算是否合理
   - 检查局部坐标系是否正确

2. **变形验证**：
   - 对比ANSYS的位移结果
   - 检查椭圆化变形趋势是否一致

3. **边界条件验证**：
   - 确保边界条件正确应用
   - 检查固定端的位移是否为零

4. **载荷验证**：
   - 检查内压、重力、热膨胀载荷是否正确计算
   - 验证载荷分配是否合理

## 后续扩展方向

1. **精确曲率计算**：从CDB文件读取或精确计算曲率半径
2. **非线性分析**：实现几何非线性和材料非线性
3. **精确刚度矩阵**：使用更精确的弯管理论
4. **耦合项**：考虑椭圆化和翘曲DOF之间的耦合
5. **数值积分**：实现环向和轴向的数值积分
6. **混合单元**：支持PIPE288和ELBOW290的混合模型

