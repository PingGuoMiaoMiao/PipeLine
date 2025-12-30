# PIPE288 单元刚度矩阵推导说明

## 1. 单元基本假设

PIPE288是一个2节点空间直管单元，每个节点有6个自由度（DOF）：
- 3个平移自由度：UX, UY, UZ
- 3个转动自由度：ROTX, ROTY, ROTZ

在线弹性分析中，可以将PIPE288简化为**空间梁单元**（Space Beam Element）来处理。

## 2. 局部坐标系与全局坐标系

### 2.1 局部坐标系定义
- X轴：沿管道轴线方向，从节点1指向节点2
- Y轴、Z轴：垂直于X轴，形成右手坐标系

### 2.2 坐标变换
需要将单元刚度矩阵从局部坐标系转换到全局坐标系：
\[
K_{global} = T^T K_{local} T
\]
其中T为坐标变换矩阵。

## 3. 局部坐标系下的单元刚度矩阵

### 3.1 几何参数
- L：单元长度
- A：管道横截面积 \( A = \pi (R_o^2 - R_i^2) \)
- Iy, Iz：截面惯性矩 \( I_y = I_z = \frac{\pi}{4}(R_o^4 - R_i^4) \)
- J：极惯性矩（扭转） \( J = \frac{\pi}{2}(R_o^4 - R_i^4) \)

其中：
- \( R_o = D/2 \)：外半径
- \( R_i = R_o - t \)：内半径
- D：管道直径
- t：壁厚

### 3.2 材料参数
- E：弹性模量
- G：剪切模量 \( G = \frac{E}{2(1+\nu)} \)
- ν：泊松比

### 3.3 局部刚度矩阵（12×12）

局部坐标下的单元刚度矩阵可以分解为几个独立的部分：

#### (1) 轴向刚度
\[
k_{axial} = \frac{EA}{L}
\begin{bmatrix}
1 & -1 \\
-1 & 1
\end{bmatrix}
\]
对应自由度：UX1, UX2

#### (2) 扭转刚度
\[
k_{torsion} = \frac{GJ}{L}
\begin{bmatrix}
1 & -1 \\
-1 & 1
\end{bmatrix}
\]
对应自由度：ROTX1, ROTX2

#### (3) 弯曲刚度（Y方向）
基于Euler-Bernoulli梁理论（忽略剪切变形）：
\[
k_{bending\_y} = \frac{EI_y}{L^3}
\begin{bmatrix}
12 & 6L & -12 & 6L \\
6L & 4L^2 & -6L & 2L^2 \\
-12 & -6L & 12 & -6L \\
6L & 2L^2 & -6L & 4L^2
\end{bmatrix}
\]
对应自由度：UY1, ROTZ1, UY2, ROTZ2

#### (4) 弯曲刚度（Z方向）
\[
k_{bending\_z} = \frac{EI_z}{L^3}
\begin{bmatrix}
12 & -6L & -12 & -6L \\
-6L & 4L^2 & 6L & 2L^2 \\
-12 & 6L & 12 & 6L \\
-6L & 2L^2 & 6L & 4L^2
\end{bmatrix}
\]
对应自由度：UZ1, ROTY1, UZ2, ROTY2

注意：Z方向的符号与Y方向不同，这是因为右手坐标系的关系。

### 3.4 组装完整的局部刚度矩阵

将上述各部分组合，得到12×12的局部刚度矩阵：

\[
K_{local} = \begin{bmatrix}
k_{axial} & 0 & 0 & 0 & 0 & 0 \\
0 & k_{bending\_y} & 0 & 0 & 0 & 0 \\
0 & 0 & k_{bending\_z} & 0 & 0 & 0 \\
0 & 0 & 0 & k_{torsion} & 0 & 0
\end{bmatrix}
\]

自由度顺序：UX1, UY1, UZ1, ROTX1, ROTY1, ROTZ1, UX2, UY2, UZ2, ROTX2, ROTY2, ROTZ2

## 4. 坐标变换矩阵

### 4.1 方向向量
设节点1坐标为\( (x_1, y_1, z_1) \)，节点2坐标为\( (x_2, y_2, z_2) \)

局部X轴方向：
\[
\vec{e_x} = \frac{(x_2-x_1, y_2-y_1, z_2-z_1)}{L}
\]

### 4.2 构建正交基
选择全局Y轴或Z轴作为参考，构造局部Y轴和Z轴：

如果\( |e_{x,z}| < 0.9 \)（X轴不接近Z轴）：
\[
\vec{e_y} = \frac{\vec{e_x} \times \vec{e_z^{global}}}{|\vec{e_x} \times \vec{e_z^{global}}|}
\]
\[
\vec{e_z} = \vec{e_x} \times \vec{e_y}
\]

否则（X轴接近Z轴）：
\[
\vec{e_y} = \frac{\vec{e_x} \times \vec{e_y^{global}}}{|\vec{e_x} \times \vec{e_y^{global}}|}
\]
\[
\vec{e_z} = \vec{e_x} \times \vec{e_y}
\]

### 4.3 变换矩阵
12×12的坐标变换矩阵T：

\[
T = \begin{bmatrix}
R & 0 & 0 & 0 \\
0 & R & 0 & 0 \\
0 & 0 & R & 0 \\
0 & 0 & 0 & R
\end{bmatrix}
\]

其中R是3×3的旋转矩阵：
\[
R = \begin{bmatrix}
e_{x,x} & e_{x,y} & e_{x,z} \\
e_{y,x} & e_{y,y} & e_{y,z} \\
e_{z,x} & e_{z,y} & e_{z,z}
\end{bmatrix}
\]

## 5. 载荷向量

### 5.1 内压等效轴向力
对于内压p，在管道端部产生的等效轴向力：
\[
F_{axial} = p \cdot A_{internal} = p \cdot \pi R_i^2
\]

分配给两个节点，各承担一半（考虑端盖效应）：
- 节点1：\( -F_{axial}/2 \)（沿X轴负方向）
- 节点2：\( +F_{axial}/2 \)（沿X轴正方向）

### 5.2 自重载荷
单位长度的重力：
\[
w = \rho \cdot A \cdot g
\]

其中ρ为密度，g为重力加速度（-Z方向，9.8 m/s²）。

将均布载荷转换为等效节点力（使用形函数积分）：
- 节点1：\( wL/2 \)（Z方向）
- 节点2：\( wL/2 \)（Z方向）

### 5.3 热膨胀载荷
对于温度变化ΔT = T - T_ref，产生的热应变：
\[
\varepsilon_{thermal} = \alpha \cdot \Delta T
\]

等效热载荷（轴向）：
\[
F_{thermal} = EA \cdot \alpha \cdot \Delta T
\]

分配给两个节点：
- 节点1：\( -F_{thermal}/2 \)
- 节点2：\( +F_{thermal}/2 \)

## 6. 边界条件处理

对于固定边界条件（D命令），采用**置大数法**或**对角元素法**：
- 在对应自由度位置，将刚度矩阵对角元素置为大数（如1e12）
- 或直接删除对应行和列（需要重新编号）

## 7. 求解方程

整体平衡方程：
\[
K_{global} \cdot U = F
\]

其中：
- \( K_{global} \)：全局刚度矩阵（组装后）
- U：节点位移向量
- F：节点载荷向量

求解：
\[
U = K_{global}^{-1} \cdot F
\]

## 8. 简化说明

本实现采用了以下简化：
1. **忽略剪切变形**：使用Euler-Bernoulli梁理论，未考虑Timoshenko梁的剪切效应
2. **小变形假设**：假设变形足够小，可以忽略几何非线性
3. **等截面假设**：假设单元内截面不变
4. **均匀材料**：假设单元内材料参数恒定

这些简化对于细长直管结构通常是合理的。
