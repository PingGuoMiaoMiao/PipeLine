# Solver 目录结构说明

本目录包含基于管单元的高温管道蠕变高效求解程序的核心代码，按照模块化设计组织。

## 目录结构

```
solver/
├── __init__.py              # 包初始化
├── main.py                  # 程序入口
│
├── mesh/                    # 网格模块
│   ├── __init__.py
│   ├── node.py             # 节点类
│   ├── element.py          # 单元基类
│   ├── section.py          # 截面类
│   └── cdb_parser.py       # ANSYS CDB文件解析器
│
├── element/                 # 单元模块
│   ├── __init__.py
│   ├── pipe288_element.py  # PIPE288直管单元
│   ├── elbow290_element.py # ELBOW290弯管单元（占位）
│   └── shape_function_fourier.py  # 环向傅里叶形函数（占位）
│
├── material/                # 材料模块
│   ├── __init__.py
│   ├── elastic.py          # 弹性材料
│   ├── plastic_biso.py     # 弹塑性材料（占位）
│   ├── creep_strain_hardening.py  # 蠕变材料（占位）
│   └── temperature_dependent.py   # 温度相关材料（占位）
│
├── integration/             # 积分模块（占位）
│   └── __init__.py
│
├── solver/                  # 求解器模块
│   ├── __init__.py
│   └── nonlinear_static.py # 非线性静力求解器（当前为线性）
│
├── load/                    # 载荷模块
│   ├── __init__.py
│   ├── pressure.py         # 压力载荷
│   ├── gravity.py          # 重力载荷
│   └── thermal_strain.py   # 热应变载荷
│
├── boundary/                # 边界条件模块
│   ├── __init__.py
│   └── displacement.py     # 位移边界条件
│
└── post/                    # 后处理模块
    ├── __init__.py
    └── vtk_writer.py       # VTK文件写入器
```

## 使用方法

### 作为模块运行

```bash
python -m solver.main examples/PIPE288_PLAST.cdb output.vtk
```

### 在代码中导入使用

```python
from solver.mesh.cdb_parser import CDBParser
from solver.material.elastic import ElasticMaterial
from solver.element.pipe288_element import PIPE288Element
from solver.solver.nonlinear_static import NonlinearStaticSolver
from solver.post.vtk_writer import VTKWriter
```

## 模块说明

### mesh/
- **node.py**: 节点类，存储节点坐标和位移
- **element.py**: 单元基类，定义单元的基本接口
- **section.py**: 截面类，存储截面几何参数
- **cdb_parser.py**: ANSYS CDB文件解析器，解析节点、单元、材料、截面、载荷、边界条件

### element/
- **pipe288_element.py**: PIPE288直管单元，实现线弹性分析
- **elbow290_element.py**: ELBOW290弯管单元（待实现）
- **shape_function_fourier.py**: 环向傅里叶级数形函数（待实现）

### material/
- **elastic.py**: 弹性材料，支持弹性模量、泊松比、密度、热膨胀系数
- **plastic_biso.py**: 双线性等向强化弹塑性材料（待实现）
- **creep_strain_hardening.py**: 应变硬化蠕变材料（待实现）
- **temperature_dependent.py**: 温度相关材料（待实现）

### solver/
- **nonlinear_static.py**: 非线性静力求解器（当前实现为线性求解）

### load/
- **pressure.py**: 压力载荷，计算内压等效轴向力
- **gravity.py**: 重力载荷，计算单位长度重力载荷
- **thermal_strain.py**: 热应变载荷，计算热膨胀等效载荷

### boundary/
- **displacement.py**: 位移边界条件，支持置大数法处理边界条件

### post/
- **vtk_writer.py**: VTK文件写入器，输出节点位移和几何信息

## 扩展方向

1. **element/elbow290_element.py**: 实现ELBOW290弯管单元
2. **material/plastic_biso.py**: 实现弹塑性材料模型
3. **material/creep_strain_hardening.py**: 实现蠕变材料模型
4. **integration/**: 实现轴向、环向、厚度方向的数值积分
5. **solver/nonlinear_static.py**: 实现真正的非线性求解（Newton-Raphson迭代）
6. **element/shape_function_fourier.py**: 实现环向傅里叶级数形函数

## 依赖

- Python 3.6+
- NumPy 1.19+

