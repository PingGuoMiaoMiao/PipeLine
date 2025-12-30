# 算法说明文档

## 1. 概述

本求解器实现了基于管单元的有限元分析，支持线弹性、弹塑性和蠕变分析。核心算法包括单元刚度矩阵计算、非线性迭代求解、时间积分等。

## 2. 单元类型

### 2.1 PIPE288 直管单元

**特性：**
- 2节点单元
- 每个节点6个自由度（UX, UY, UZ, ROTX, ROTY, ROTZ）
- 基于空间梁单元理论

**刚度矩阵：**
采用Euler-Bernoulli梁理论，包括：
- 轴向刚度：\(k_{axial} = \frac{EA}{L}\)
- 弯曲刚度：\(k_{bending} = \frac{EI}{L^3}\)
- 扭转刚度：\(k_{torsion} = \frac{GJ}{L}\)

详细推导见 `pipe288_stiffness_derivation.md`

### 2.2 ELBOW290 弯管单元

**特性：**
- 3节点单元
- 每个节点10个自由度（6个标准DOF + 4个截面变形DOF）
- 支持椭圆化（Ovalization）和翘曲（Warping）

**截面变形：**
使用傅里叶级数展开：
\[
\Delta r(\phi) = \delta_2^c \cos(2\phi) + \delta_2^s \sin(2\phi) + 
                \delta_3^c \cos(3\phi) + \delta_3^s \sin(3\phi)
\]

详细说明见 `ELBOW290_实现说明.md`

## 3. 材料模型

### 3.1 线弹性材料

**本构关系：**
\[
\sigma = E \varepsilon
\]

**参数：**
- E: 弹性模量 (Pa)
- ν: 泊松比
- G: 剪切模量 = \(E / (2(1+\nu))\)

### 3.2 双线性等向强化（BISO）弹塑性材料

**屈服准则：** Mises屈服准则
\[
\sigma_{mises} = \sqrt{\frac{3}{2} s_{ij} s_{ij}} \leq \sigma_y
\]

**硬化规律：**
\[
\sigma_y(\varepsilon_p) = \sigma_{y0} + H \varepsilon_p
\]

**应力更新：** 返回映射算法（Return Mapping Algorithm）

详细说明见 `PLASTIC_实现说明.md`

### 3.3 各向同性应变硬化蠕变材料

**蠕变率公式：**
\[
\dot{\varepsilon}_{cr} = C_1 \cdot \sigma^{C_2} \cdot \varepsilon_{cr}^{C_3} \cdot \exp(-C_4/T)
\]

**参数：**
- C1: 蠕变系数
- C2: 应力指数
- C3: 时间指数（通常 < 0）
- C4: 温度相关参数
- T: 绝对温度 (K)

**时间积分：** 显式Euler方法
\[
\varepsilon_{cr}^{(n+1)} = \varepsilon_{cr}^{(n)} + \dot{\varepsilon}_{cr}^{(n)} \Delta t
\]

详细说明见 `CREEP_实现说明.md`

## 4. 求解算法

### 4.1 线性静力分析

**控制方程：**
\[
K U = F
\]

**求解：**
- 直接方法：`numpy.linalg.solve`
- 稀疏矩阵：`scipy.sparse.linalg.spsolve`（优化版本）

### 4.2 非线性静力分析（Newton-Raphson）

**控制方程：**
\[
R(U) = F_{ext} - F_{int}(U) = 0
\]

**迭代公式：**
\[
K_{tan}^{(i)} \Delta U^{(i)} = R^{(i)}
\]
\[
U^{(i+1)} = U^{(i)} + \Delta U^{(i)}
\]

**收敛准则：**
\[
\|R\| < \varepsilon_{tol}
\]

**切线刚度矩阵：**
- 弹性：\(K_{tan} = K_{el}\)
- 弹塑性：\(K_{tan} = K_{ep}\)（一致性切线模量）

### 4.3 蠕变时间积分

**算法流程：**

```
for each time step:
    1. 保存上一时间步的状态
    2. Newton-Raphson迭代求解：
        - 更新蠕变应变
        - 组装切线刚度矩阵
        - 计算残差
        - 求解线性方程组
        - 更新位移
    3. 检查收敛
    4. 自适应调整时间步（可选）
```

**自适应时间步控制：**
- 收敛快 → 增大时间步
- 收敛慢 → 减小时间步
- 未收敛 → 大幅减小时间步并重试

## 5. 载荷处理

### 5.1 内压载荷

**等效轴向力：**
\[
F_{axial} = p \cdot \pi R_i^2
\]

**分配：** 均匀分配到两个节点

### 5.2 重力载荷

**单位长度重量：**
\[
w = \rho A g
\]

**等效节点力：**
\[
F_{node} = \frac{wL}{2}
\]

### 5.3 热膨胀载荷

**热应变：**
\[
\varepsilon_{th} = \alpha (T - T_{ref})
\]

**等效节点力：**
\[
F_{th} = EA \alpha (T - T_{ref})
\]

## 6. 边界条件

**应用方法：** 置大数法（Penalty Method）

对于固定自由度 \(u_i = 0\)：
\[
K_{ii} = 10^{12}, \quad F_i = 0
\]

**优点：** 不改变矩阵大小，实现简单

## 7. 性能优化

### 7.1 稀疏矩阵

**存储格式：** CSR（Compressed Sparse Row）

**优势：**
- 内存：O(nnz) vs O(n²)
- 求解：针对稀疏矩阵优化的算法

### 7.2 并行计算

**策略：** 多进程并行装配

**适用场景：**
- 单元数 > 100
- 多核CPU

### 7.3 自适应时间步

**策略：**
- 根据迭代次数调整
- 根据收敛性调整
- 根据误差估计调整

详细说明见 `PERFORMANCE_OPTIMIZATION.md`

## 8. 数值精度

### 8.1 精度保证

- 所有优化都是数值等价的
- 稀疏矩阵：结果差异 < 1e-12
- 并行计算：结果差异 < 1e-14
- 自适应时间步：相对误差 < 1e-6

### 8.2 收敛性

- Newton-Raphson：通常5-10次迭代收敛
- 时间积分：显式方法，需要小时间步长保证稳定性

## 9. 参考文档

- `pipe288_stiffness_derivation.md` - PIPE288刚度矩阵推导
- `ELBOW290_实现说明.md` - ELBOW290单元实现
- `PLASTIC_实现说明.md` - 弹塑性分析实现
- `CREEP_实现说明.md` - 蠕变分析实现
- `PERFORMANCE_OPTIMIZATION.md` - 性能优化说明

