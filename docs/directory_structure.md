# 目录结构说明

## 顶层目录

```
Piprline/
├── solver/              # 核心求解器代码
├── examples/            # 示例文件目录
├── docs/                # 文档目录
├── README.md            # 项目说明文档
└── 其他文档和脚本
```

## solver/ 目录结构

```
solver/
├── __init__.py          # 包初始化文件
├── main.py              # 线弹性分析主程序
├── main_plastic.py      # 弹塑性分析主程序
├── main_creep.py        # 蠕变分析主程序
│
├── mesh/                # 网格模块
│   ├── __init__.py
│   ├── node.py          # 节点类（坐标、位移）
│   ├── element.py       # 单元基类
│   ├── section.py       # 截面类
│   └── cdb_parser.py    # ANSYS CDB文件解析器
│
├── element/             # 单元模块
│   ├── __init__.py
│   ├── pipe288_element.py          # PIPE288单元（线弹性）
│   ├── pipe288_element_plastic.py  # PIPE288单元（弹塑性）
│   ├── pipe288_element_creep.py    # PIPE288单元（蠕变）
│   ├── elbow290_element.py         # ELBOW290弯管单元
│   └── shape_function_fourier.py   # 傅里叶级数形函数
│
├── material/            # 材料模块
│   ├── __init__.py
│   ├── elastic.py                   # 弹性材料
│   ├── plastic_biso.py              # BISO弹塑性材料
│   ├── creep_strain_hardening.py    # 蠕变材料
│   ├── plastic_biso.py              # （占位：其他材料模型）
│   └── temperature_dependent.py     # （占位：温度相关材料）
│
├── solver/              # 求解器模块
│   ├── __init__.py
│   ├── nonlinear_static.py          # 线性/非线性静力求解器
│   ├── newton_raphson.py            # Newton-Raphson迭代求解器
│   ├── creep_time_integration.py    # 蠕变时间积分求解器
│   ├── optimized_creep_solver.py    # 优化的蠕变求解器（整合所有优化）
│   ├── sparse_solver.py             # 稀疏矩阵求解器
│   ├── adaptive_time_step.py        # 自适应时间步控制
│   └── parallel_assembly.py         # 并行装配
│
├── load/                # 载荷模块
│   ├── __init__.py
│   ├── pressure.py      # 压力载荷计算
│   ├── gravity.py       # 重力载荷计算
│   └── thermal_strain.py # 热应变载荷计算
│
├── boundary/            # 边界条件模块
│   ├── __init__.py
│   └── displacement.py  # 位移边界条件处理
│
├── integration/         # 积分模块
│   ├── __init__.py
│   └── integration_point.py  # 积分点状态变量（应力、应变、塑性、蠕变）
│
└── post/                # 后处理模块
    ├── __init__.py
    ├── vtk_writer.py            # VTK文件写入器
    └── pipe_surface_expand.py   # 管单元曲面展开（ELBOW290）
```

## 模块说明

### mesh/ - 网格模块

负责网格数据的存储和管理：
- `node.py`: 节点类，存储节点坐标和位移
- `element.py`: 单元基类，定义单元的基本接口
- `section.py`: 截面类，存储截面几何参数
- `cdb_parser.py`: ANSYS CDB文件解析器，解析节点、单元、材料、载荷等

### element/ - 单元模块

实现不同类型的有限元单元：
- `pipe288_element.py`: PIPE288直管单元（线弹性）
- `pipe288_element_plastic.py`: PIPE288单元（弹塑性）
- `pipe288_element_creep.py`: PIPE288单元（蠕变）
- `elbow290_element.py`: ELBOW290弯管单元（支持椭圆化和翘曲）
- `shape_function_fourier.py`: 环向傅里叶级数形函数

### material/ - 材料模块

实现不同的材料本构模型：
- `elastic.py`: 线弹性材料
- `plastic_biso.py`: 双线性等向强化弹塑性材料
- `creep_strain_hardening.py`: 各向同性应变硬化蠕变材料

### solver/ - 求解器模块

实现不同的求解算法：
- `nonlinear_static.py`: 线性/非线性静力求解器（基础版）
- `newton_raphson.py`: Newton-Raphson非线性迭代求解器
- `creep_time_integration.py`: 蠕变时间积分求解器
- `optimized_creep_solver.py`: 优化的蠕变求解器（整合稀疏矩阵、并行、自适应时间步）
- `sparse_solver.py`: 稀疏矩阵求解器
- `adaptive_time_step.py`: 自适应时间步控制器
- `parallel_assembly.py`: 并行装配器

### load/ - 载荷模块

处理不同类型的载荷：
- `pressure.py`: 内压载荷（等效轴向力）
- `gravity.py`: 重力载荷（单位长度重量）
- `thermal_strain.py`: 热应变载荷（热膨胀等效力）

### boundary/ - 边界条件模块

处理边界条件：
- `displacement.py`: 位移边界条件（固定位移/转动）

### integration/ - 积分模块

积分点状态管理：
- `integration_point.py`: 积分点类，存储应力、应变、塑性应变、蠕变应变等历史变量

### post/ - 后处理模块

结果输出和可视化：
- `vtk_writer.py`: VTK文件写入器，输出位移、应力、塑性、蠕变等结果
- `pipe_surface_expand.py`: 管单元曲面展开，用于ELBOW290的三维曲面输出

## 主程序文件

### main.py
线弹性分析的主程序入口：
- 解析CDB文件
- 创建单元和材料
- 线性求解
- 输出VTK文件

### main_plastic.py
弹塑性分析的主程序入口：
- 支持BISO弹塑性材料
- Newton-Raphson迭代
- 增量加载
- 输出塑性状态

### main_creep.py
蠕变分析的主程序入口：
- 支持蠕变材料
- 时间步进分析
- 输出多个时刻的VTK文件
- 支持自适应时间步

## 设计原则

### 模块化设计

- 每个模块职责单一
- 模块间低耦合
- 易于扩展和维护

### 可扩展性

- 单元类型：易于添加新单元（如ELBOW290）
- 材料模型：易于添加新材料（如其他蠕变模型）
- 求解算法：易于添加新算法（如隐式时间积分）

### 性能优化

- 稀疏矩阵：大规模问题
- 并行计算：多核CPU
- 自适应控制：提高效率

### 工程化

- 完整的文档
- 清晰的接口
- 易于使用
- 结果可视化

## 扩展方向

### 已实现
- ✅ PIPE288单元（线弹性、弹塑性、蠕变）
- ✅ ELBOW290单元（线弹性）
- ✅ BISO弹塑性材料
- ✅ 应变硬化蠕变材料

### 待扩展（占位文件已创建）
- ⏳ ELBOW290弹塑性和蠕变
- ⏳ 其他材料模型（温度相关、其他蠕变模型等）
- ⏳ 隐式时间积分
- ⏳ GPU加速
- ⏳ 更多单元类型

## 文件命名规范

- **模块文件**：小写字母，下划线分隔（如 `cdb_parser.py`）
- **类名**：驼峰命名（如 `CDBParser`）
- **函数名**：小写字母，下划线分隔（如 `compute_stiffness_matrix`）
- **主程序**：`main*.py` 格式

## 代码组织原则

1. **单一职责**：每个类/模块只负责一个功能
2. **接口清晰**：公共方法有明确的文档字符串
3. **依赖最小化**：模块间依赖关系清晰
4. **可测试性**：关键函数易于单独测试

