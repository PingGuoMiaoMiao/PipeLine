# 弹塑性分析实现说明

## 1. 实现概述

已实现BISO（双线性等向强化）弹塑性材料模型，支持Newton-Raphson非线性迭代求解和增量加载。

## 2. 核心模块

### 2.1 BISOPlasticMaterial类

**文件：** `solver/material/plastic_biso.py`

**功能：**
- 双线性等向强化弹塑性材料模型
- 返回映射算法（Return Mapping Algorithm）进行应力更新
- 一致性切线模量（Consistent Tangent Modulus）计算
- Mises屈服准则

**材料参数：**
- E: 弹性模量 (Pa)
- nu: 泊松比
- yield_stress: 屈服应力 (Pa)
- tangent_modulus: 切线模量 (Pa)

**关键方法：**
- `compute_stress_from_strain()`: 从应变计算应力（弹塑性）
- `_return_mapping()`: 返回映射算法
- `compute_tangent_stiffness()`: 计算切线刚度矩阵
- `is_plastic()`: 判断是否进入塑性

### 2.2 IntegrationPoint类

**文件：** `solver/integration/integration_point.py`

**功能：**
- 存储积分点的状态变量
- 应力、应变、塑性应变
- 等效塑性应变（历史变量）
- 状态保存和更新

### 2.3 NewtonRaphsonSolver类

**文件：** `solver/solver/newton_raphson.py`

**功能：**
- Newton-Raphson非线性迭代求解
- 残差计算和收敛判断
- 切线刚度矩阵组装
- 内部力向量组装
- 增量加载支持

**关键方法：**
- `solve()`: 求解非线性方程组
- `_compute_residual()`: 计算残差向量
- `_assemble_tangent_stiffness()`: 组装切线刚度矩阵
- `_assemble_internal_force()`: 组装内部力向量

### 2.4 PIPE288ElementPlastic类

**文件：** `solver/element/pipe288_element_plastic.py`

**功能：**
- 支持弹塑性的PIPE288单元
- 积分点管理
- 切线刚度矩阵计算
- 内部力向量计算
- 塑性状态输出

**关键方法：**
- `compute_tangent_stiffness()`: 计算切线刚度矩阵
- `compute_internal_force()`: 计算内部力向量
- `get_plastic_state()`: 获取塑性状态

## 3. 使用方法

### 3.1 基本用法

```python
from solver.material.plastic_biso import BISOPlasticMaterial
from solver.solver.newton_raphson import NewtonRaphsonSolver

# 创建弹塑性材料
plastic_material = BISOPlasticMaterial(
    E=200e9,
    nu=0.3,
    yield_stress=25e6,  # 25 MPa
    tangent_modulus=100e9  # 100 GPa
)

# 创建求解器
solver = NewtonRaphsonSolver(
    nodes, elements, boundary_conditions,
    dof_per_node=6, max_iterations=20, tolerance=1e-6
)

# 增量加载
for step in range(num_steps):
    load_factor = (step + 1) / num_steps
    displacements, info = solver.solve(load_vector, incremental=True, 
                                      load_factor=load_factor)
    print(f"步骤 {step+1}: 收敛={info['converged']}, "
          f"迭代={info['iterations']}")
```

### 3.2 运行程序

```bash
# 运行弹塑性分析
python -m solver.main_plastic examples/PIPE288_PLAST.cdb output_plastic.vtk 10
```

参数说明：
- `input.cdb`: 输入CDB文件
- `output.vtk`: 输出VTK文件
- `num_steps`: 增量步数（可选，默认10）

## 4. 输出格式

### 4.1 控制台输出

- 每个增量步的载荷因子
- 每次迭代的残差范数
- 收敛信息
- 最终位移统计
- 塑性单元统计

### 4.2 VTK输出

VTK文件包含以下数据：

1. **节点位移**：
   - Displacement: 位移向量
   - DisplacementMagnitude: 位移幅值

2. **单元数据**：
   - ElementType: 单元类型（288）
   - PlasticState: 塑性状态（0=弹性，1=塑性）
   - EquivalentPlasticStrain: 等效塑性应变

在ParaView中：
1. 打开VTK文件
2. 选择"PlasticState"或"EquivalentPlasticStrain"
3. 应用颜色映射查看塑性区分布

## 5. 算法说明

### 5.1 返回映射算法

返回映射算法（Return Mapping Algorithm）是弹塑性分析中常用的应力更新算法：

1. **弹性预测**：计算弹性应力
2. **屈服判断**：检查是否超过屈服面
3. **塑性修正**：如果进入塑性，进行径向返回

### 5.2 Newton-Raphson迭代

非线性方程求解：

\[
R(U) = F_{ext} - F_{int}(U) = 0
\]

迭代公式：

\[
K_{tan}^{(i)} \Delta U^{(i)} = R^{(i)}
\]
\[
U^{(i+1)} = U^{(i)} + \Delta U^{(i)}
\]

其中：
- \(K_{tan}\): 切线刚度矩阵
- \(F_{int}\): 内部力向量
- \(F_{ext}\): 外部力向量
- \(R\): 残差向量

### 5.3 增量加载

为了处理路径相关的弹塑性问题，采用增量加载：

1. 将总载荷分为多个增量步
2. 每个增量步使用Newton-Raphson迭代求解
3. 更新状态变量（应力、塑性应变等）

## 6. 验证建议

1. **弹性阶段验证**：
   - 载荷较小时，应该为纯弹性
   - 结果应与线性分析一致

2. **屈服验证**：
   - 当应力达到屈服应力时，应该开始进入塑性
   - 检查塑性单元的位置是否合理

3. **硬化验证**：
   - 进入塑性后，应力应该继续增加
   - 检查切线模量的影响

4. **收敛性验证**：
   - 残差应该单调下降
   - 通常应在5-10次迭代内收敛

5. **与ANSYS对比**：
   - 位移趋势应该一致
   - 塑性区分布应该相似
   - 等效塑性应变量级应该合理

## 7. 简化说明

当前实现采用了以下简化：

1. **梁单元简化**：对于PIPE288单元，主要考虑轴向和弯曲的弹塑性，使用简化的刚度降低方法
2. **积分点简化**：使用2个高斯积分点，实际应用中可能需要更多积分点
3. **塑性判断简化**：基于单元平均应变判断，而非精确的积分点应力

这些简化对于初步实现和验证是合理的，后续可以逐步完善。

## 8. 后续扩展方向

1. **精确的积分点应力更新**：在每个积分点处进行精确的应力更新
2. **多种屈服准则**：支持Tresca、Drucker-Prager等
3. **多种硬化模型**：支持随动强化、混合强化等
4. **塑性铰模型**：针对梁单元的塑性铰模型
5. **自适应增量步**：根据收敛性自动调整增量步大小

