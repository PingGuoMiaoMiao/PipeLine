# 5分钟快速上手

## 目标

让您能在5分钟内理解并运行管单元蠕变求解器。

## 第1分钟：环境准备

### 安装Python（如果未安装）

```bash
# 检查Python版本（需要3.7+）
python --version

# 如果没有Python，从 https://www.python.org/ 下载安装
```

### 安装依赖

```bash
# 必需
pip install numpy

# 推荐（性能优化）
pip install scipy
```

### 验证安装

```bash
python -c "import numpy; print('NumPy:', numpy.__version__)"
python -c "import scipy; print('SciPy:', scipy.__version__)"  # 可选
```

## 第2分钟：准备输入文件

### 示例文件位置

确保 `examples/` 目录下有CDB文件：
- `examples/PIPE288_PLAST.cdb` - 弹塑性算例
- `examples/ELBOW290_PLAST.cdb` - 弯管算例
- `examples/PIPE288_CREEP.cdb` - 蠕变算例

### 检查文件

```bash
# Windows
dir examples\*.cdb

# Linux/Mac
ls examples/*.cdb
```

## 第3分钟：运行第一个分析

### 线弹性分析（最简单）

```bash
python -m solver.main examples/PIPE288_PLAST.cdb output.vtk
```

**您将看到：**
```
解析CDB文件: examples/PIPE288_PLAST.cdb
解析完成:
  - 节点数: 31
  - 单元数: 30
开始有限元分析...
最大位移: UX=2.231 mm, UY=2.930 mm, UZ=5.027 mm
写入VTK文件: output.vtk
完成！
```

**关键信息：**
- ✅ 自动解析CDB文件
- ✅ 计算位移
- ✅ 生成VTK文件

## 第4分钟：查看结果

### 打开ParaView

```bash
# 启动ParaView
paraview
```

### 加载结果

1. **File → Open** → 选择 `output.vtk`
2. 点击 **Apply**
3. 在左侧选择 **DisplacementMagnitude**
4. 查看位移云图

### 显示变形

1. **Filters → Warp By Vector**
2. 选择 **Displacement**
3. 设置 **Scale Factor = 100**（放大变形）
4. 点击 **Apply**

**现在您应该能看到：**
- ✅ 管道的变形形状
- ✅ 颜色表示位移大小
- ✅ 变形被放大了100倍以便观察

## 第5分钟：尝试更多功能

### 弹塑性分析

```bash
python -m solver.main_plastic examples/PIPE288_PLAST.cdb output_plastic.vtk 5
```

**在ParaView中：**
- 选择 **PlasticState** 查看哪些单元进入塑性
- 选择 **EquivalentPlasticStrain** 查看塑性应变大小

### 蠕变分析

```bash
python -m solver.main_creep examples/PIPE288_CREEP.cdb output_creep 1000 10
```

**在ParaView中：**
- 加载多个时间步文件（output_creep_t*.vtk）
- 播放动画查看蠕变发展
- 选择 **EquivalentCreepStrain** 查看蠕变应变

## 核心概念（理解要点）

### 1. 单元类型

- **PIPE288**：直管单元（2节点，6 DOF/节点）
- **ELBOW290**：弯管单元（3节点，10 DOF/节点，支持椭圆化）

### 2. 材料模型

- **线弹性**：小变形，线性关系
- **弹塑性**：考虑塑性变形，非线性
- **蠕变**：高温下时间相关的变形

### 3. 分析类型

- **线性静力**：一次求解，快速
- **非线性静力**：迭代求解，考虑非线性
- **蠕变时间积分**：多个时间步，考虑时间效应

### 4. 输出格式

- **VTK格式**：ParaView标准格式
- **包含数据**：位移、应力、塑性应变、蠕变应变等

## 常见命令速查

```bash
# 线弹性分析
python -m solver.main <input.cdb> <output.vtk>

# 弹塑性分析
python -m solver.main_plastic <input.cdb> <output.vtk> <steps>

# 蠕变分析
python -m solver.main_creep <input.cdb> <output_prefix> <total_time> <steps>
```

## 下一步

- 📖 阅读 `README.md` 了解详细功能
- 📖 查看 `docs/algorithm.md` 了解算法原理
- 📖 参考 `docs/visualization.md` 学习ParaView技巧
- 📖 运行 `performance_benchmark.py` 查看性能优化效果

## 遇到问题？

1. **Python未找到**：检查PATH环境变量
2. **模块导入错误**：确认在项目根目录运行
3. **ParaView无法打开**：检查VTK文件是否生成成功
4. **结果异常**：检查CDB文件格式是否正确

**提示**：所有问题都可以在 `README.md` 的"常见问题"部分找到答案。

---

**恭喜！** 您已经掌握了基本使用方法。现在可以：
- ✅ 运行不同类型的分析
- ✅ 在ParaView中查看结果
- ✅ 理解核心概念

开始探索更多功能吧！🚀

