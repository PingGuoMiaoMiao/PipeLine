# 蠕变分析实现说明

## 1. 实现概述

已实现高温蠕变分析功能，支持各向同性应变硬化（Strain Hardening）蠕变模型，与弹性、塑性耦合分析。

## 2. 蠕变模型

### 2.1 蠕变率公式

各向同性应变硬化蠕变模型：

\[
\dot{\varepsilon}_{cr} = C_1 \cdot \sigma^{C_2} \cdot t^{C_3} \cdot \exp(-C_4/T)
\]

其中：
- \(C_1\): 蠕变系数
- \(C_2\): 应力指数（通常≥1）
- \(C_3\): 时间指数（通常<1）
- \(C_4\): 温度相关参数（通常为Q/R，Q为激活能，R为气体常数）
- \(\sigma\): Mises等效应力
- \(t\): 时间（或等效蠕变应变，用于Strain Hardening模型）
- \(T\): 温度（绝对温度，K）

### 2.2 参数说明

从算例说明中，典型参数值为：
- C1 = 1.000000e-014
- C2 = 2
- C3 = -1
- C4 = 0

## 3. 核心模块

### 3.1 CreepStrainHardeningMaterial类

**文件：** `solver/material/creep_strain_hardening.py`

**功能：**
- 蠕变应变率计算
- 蠕变应变增量计算（显式/隐式）
- 蠕变应变方向计算
- Mises等效应力计算

**关键方法：**
- `compute_creep_strain_rate()`: 计算蠕变应变率
- `compute_creep_strain_increment()`: 计算蠕变应变增量
- `compute_creep_strain_direction()`: 计算蠕变应变方向

### 3.2 CreepTimeIntegrationSolver类

**文件：** `solver/solver/creep_time_integration.py`

**功能：**
- 时间步进求解
- 蠕变应变更新
- 与Newton-Raphson迭代耦合
- 时间历史记录

**关键方法：**
- `solve_time_step()`: 求解一个时间步
- `_update_creep_strain()`: 更新蠕变应变

### 3.3 PIPE288ElementCreep类

**文件：** `solver/element/pipe288_element_creep.py`

**功能：**
- 支持蠕变的PIPE288单元
- 积分点蠕变应变管理
- 蠕变状态输出

**关键方法：**
- `update_creep_strain()`: 更新蠕变应变
- `get_creep_state()`: 获取蠕变状态

### 3.4 IntegrationPoint扩展

**文件：** `solver/integration/integration_point.py`

**新增状态变量：**
- `creep_strain`: 蠕变应变向量（6×1）
- `equivalent_creep_strain`: 等效蠕变应变（标量）

## 4. 使用方法

### 4.1 基本用法

```python
from solver.material.creep_strain_hardening import CreepStrainHardeningMaterial
from solver.solver.creep_time_integration import CreepTimeIntegrationSolver

# 创建蠕变材料
creep_material = CreepStrainHardeningMaterial(
    C1=1.0e-14,
    C2=2.0,
    C3=-1.0,
    C4=0.0
)

# 创建时间积分求解器
time_solver = CreepTimeIntegrationSolver(
    nodes, elements, boundary_conditions,
    dof_per_node=6, max_iterations=20, tolerance=1e-6
)

# 时间步进分析
for step in range(num_steps):
    time = step * time_increment
    displacements, info = time_solver.solve_time_step(
        load_vector, time, time_increment, temperature_kelvin
    )
```

### 4.2 运行程序

```bash
# 运行蠕变分析
python -m solver.main_creep examples/PIPE288_CREEP.cdb output_creep 10000 100
```

参数说明：
- `input.cdb`: 输入CDB文件
- `output_prefix`: 输出文件前缀
- `total_time`: 总时间（默认10000）
- `num_steps`: 时间步数（默认100）

## 5. 时间积分方法

### 5.1 显式方法

当前实现使用显式Euler方法：

\[
\varepsilon_{cr}^{(n+1)} = \varepsilon_{cr}^{(n)} + \dot{\varepsilon}_{cr}^{(n)} \Delta t
\]

**优点：**
- 实现简单
- 计算效率高

**缺点：**
- 可能不稳定（需要小的时间步长）
- 精度较低

### 5.2 隐式方法（未来扩展）

可以考虑使用隐式方法：

\[
\varepsilon_{cr}^{(n+1)} = \varepsilon_{cr}^{(n)} + \dot{\varepsilon}_{cr}^{(n+1)} \Delta t
\]

需要迭代求解，但稳定性更好。

## 6. 耦合分析

### 6.1 弹性-蠕变耦合

总应变：

\[
\varepsilon_{total} = \varepsilon_{elastic} + \varepsilon_{creep}
\]

应力-应变关系：

\[
\sigma = E \cdot (\varepsilon_{total} - \varepsilon_{creep})
\]

### 6.2 弹性-塑性-蠕变耦合

总应变：

\[
\varepsilon_{total} = \varepsilon_{elastic} + \varepsilon_{plastic} + \varepsilon_{creep}
\]

需要同时考虑：
- 塑性屈服
- 蠕变变形
- 两者的相互作用

## 7. 输出格式

### 7.1 控制台输出

- 每个时间步的进度
- 收敛信息
- 位移统计
- 蠕变应变统计
- VTK文件输出提示

### 7.2 VTK输出

每个输出时刻生成一个VTK文件，包含：

1. **节点数据**：
   - Displacement: 位移向量
   - DisplacementMagnitude: 位移幅值

2. **单元数据**：
   - ElementType: 单元类型
   - EquivalentCreepStrain: 等效蠕变应变（云图）

在ParaView中：
1. 打开VTK文件序列
2. 选择"EquivalentCreepStrain"
3. 应用颜色映射查看蠕变应变分布
4. 使用时间滑块查看不同时刻的变形

## 8. 验证建议

1. **蠕变率验证**：
   - 检查蠕变应变率是否随应力增加而增加
   - 检查温度依赖性是否正确

2. **时间积分验证**：
   - 检查蠕变应变是否随时间累积
   - 检查时间步长的影响

3. **耦合验证**：
   - 检查弹性-蠕变耦合是否正确
   - 检查塑性-蠕变耦合（如果启用）

4. **与ANSYS对比**：
   - 对比位移-时间曲线
   - 对比蠕变应变分布
   - 检查变形趋势是否一致

## 9. 简化说明

当前实现采用了以下简化：

1. **显式积分**：使用显式Euler方法，可能需要小的时间步长
2. **简化应力计算**：对于梁单元，使用简化的应力计算方法
3. **线性化**：在某些情况下使用线性化的刚度矩阵
4. **温度恒定**：假设温度在整个分析过程中恒定

这些简化对于初步实现和验证是合理的，后续可以逐步完善。

## 10. 后续扩展方向

1. **隐式积分**：实现隐式时间积分方法以提高稳定性
2. **自适应时间步**：根据收敛性和精度自动调整时间步长
3. **温度场**：支持温度场分布（而非均匀温度）
4. **多种蠕变模型**：支持Time Hardening、Generalized Strain Hardening等
5. **精确应力计算**：在积分点处进行精确的应力计算
6. **蠕变疲劳**：考虑蠕变-疲劳交互作用

