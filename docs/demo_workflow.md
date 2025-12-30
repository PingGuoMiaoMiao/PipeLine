# 示例演示流程

## 目标

在5分钟内展示管单元蠕变求解器的核心功能和使用方法。

## 演示脚本（5分钟版本）

### 第1分钟：项目介绍和快速运行（1分钟）

```bash
# 1. 展示项目结构
tree solver -L 2

# 2. 快速运行线弹性分析
python -m solver.main examples/PIPE288_PLAST.cdb output_demo.vtk

# 输出关键信息：
# - 节点数、单元数
# - 材料参数
# - 最大位移
# - 生成VTK文件
```

**要点：**
- ✅ 简单易用：一条命令完成分析
- ✅ 输出直观：VTK格式，ParaView可视化
- ✅ 功能完整：支持PIPE288和ELBOW290

### 第2分钟：ParaView可视化（1分钟）

```bash
# 在ParaView中打开结果
paraview output_demo.vtk
```

**演示步骤：**
1. 打开VTK文件
2. 选择 "DisplacementMagnitude" 显示位移云图
3. 应用 "Warp By Vector" 显示变形（Scale Factor = 100）
4. 添加 "Contour" 显示等值线

**要点：**
- ✅ 可视化效果直观
- ✅ 操作简单
- ✅ 多种显示方式

### 第3分钟：弹塑性分析（1分钟）

```bash
# 运行弹塑性分析
python -m solver.main_plastic examples/PIPE288_PLAST.cdb output_plastic.vtk 5

# 输出关键信息：
# - 每个增量步的迭代过程
# - 收敛信息
# - 塑性单元统计
```

**演示要点：**
- ✅ Newton-Raphson迭代过程
- ✅ 收敛监控
- ✅ 塑性区分布

**ParaView展示：**
- 选择 "PlasticState" 显示塑性区
- 选择 "EquivalentPlasticStrain" 显示塑性应变云图

### 第4分钟：蠕变分析（1分钟）

```bash
# 运行蠕变分析（简化：10个时间步）
python -m solver.main_creep examples/PIPE288_CREEP.cdb output_creep 1000 10

# 输出关键信息：
# - 时间步进过程
# - 蠕变应变累积
# - 多个时刻的VTK文件
```

**演示要点：**
- ✅ 时间步进分析
- ✅ 蠕变应变累积
- ✅ 时间序列输出

**ParaView展示：**
- 加载时间序列文件
- 播放动画查看蠕变发展
- 选择 "EquivalentCreepStrain" 显示蠕变应变

### 第5分钟：性能优化和总结（1分钟）

```bash
# 展示性能优化（如果时间允许）
python performance_benchmark.py examples/PIPE288_PLAST.cdb 10
```

**总结要点：**

1. **核心功能**
   - ✅ PIPE288和ELBOW290单元
   - ✅ 线弹性、弹塑性、蠕变分析
   - ✅ 完整的材料模型

2. **技术亮点**
   - ✅ 稀疏矩阵优化（内存减少90%+）
   - ✅ 并行计算（加速4-8倍）
   - ✅ 自适应时间步
   - ✅ 与ANSYS结果对比验证

3. **工程化**
   - ✅ 模块化设计
   - ✅ 完整文档
   - ✅ 易于使用
   - ✅ 可视化友好

## 完整演示脚本（详细版）

### 准备阶段（演示前）

```bash
# 1. 检查环境
python --version
python -c "import numpy; print('NumPy OK')"
python -c "import scipy; print('SciPy OK')"  # 可选

# 2. 准备示例文件
ls examples/*.cdb

# 3. 启动ParaView（后台）
paraview &
```

### 演示阶段

#### 步骤1：线弹性分析（2分钟）

```bash
# 运行分析
echo "=== 线弹性分析 ===" 
python -m solver.main examples/PIPE288_PLAST.cdb demo_elastic.vtk

# 在ParaView中展示：
# - 位移云图
# - 变形显示
# - 等值线
```

#### 步骤2：弹塑性分析（2分钟）

```bash
# 运行分析
echo "=== 弹塑性分析 ===" 
python -m solver.main_plastic examples/PIPE288_PLAST.cdb demo_plastic.vtk 5

# 在ParaView中展示：
# - 塑性状态分布
# - 等效塑性应变
# - 变形对比
```

#### 步骤3：蠕变分析（2分钟）

```bash
# 运行分析
echo "=== 蠕变分析 ===" 
python -m solver.main_creep examples/PIPE288_CREEP.cdb demo_creep 1000 10

# 在ParaView中展示：
# - 时间序列动画
# - 蠕变应变云图
# - 位移-时间曲线（如果可能）
```

#### 步骤4：性能展示（1分钟）

```bash
# 性能对比（如果时间允许）
echo "=== 性能优化 ===" 
python performance_benchmark.py examples/PIPE288_PLAST.cdb 10
```

### 总结阶段

**核心成果：**
1. ✅ 完整的管单元分析能力（PIPE288 + ELBOW290）
2. ✅ 多物理场耦合（弹性+塑性+蠕变）
3. ✅ 高性能优化（稀疏矩阵+并行+自适应）
4. ✅ 工程化实现（模块化+文档+可视化）

## 演示技巧

### 1. 控制时间

- **5分钟版本**：重点展示线弹性和弹塑性，快速演示蠕变
- **10分钟版本**：完整展示所有功能
- **15分钟版本**：包含性能优化和算法细节

### 2. 突出亮点

- **功能完整性**：支持多种单元和材料模型
- **性能优化**：稀疏矩阵、并行计算、自适应时间步
- **工程化程度**：模块化设计、完整文档、易于使用
- **可视化**：ParaView友好，结果直观

### 3. 常见问题准备

- **Q: 与ANSYS对比如何？**
  - A: 已通过多个算例验证，位移趋势一致，详见验证文档

- **Q: 性能如何？**
  - A: 优化后加速6-7倍，内存减少90%+，详见性能报告

- **Q: 精度如何保证？**
  - A: 使用标准有限元算法，所有优化都是数值等价的

### 4. 展示建议

**推荐展示顺序：**
1. 快速运行（展示易用性）
2. 可视化结果（展示直观性）
3. 多物理场（展示完整性）
4. 性能优化（展示技术深度）

**不推荐：**
- 深入算法细节（除非评委提问）
- 展示大量代码（重点在功能）
- 过度强调性能数字（重点是功能完整性）

## 演示检查清单

演示前确认：
- [ ] 示例CDB文件已准备好
- [ ] Python环境和依赖已安装
- [ ] ParaView已安装并可运行
- [ ] 输出目录可写入
- [ ] 演示脚本已测试

演示中注意：
- [ ] 控制时间，重点突出
- [ ] 解释关键输出
- [ ] 展示可视化效果
- [ ] 回答问题时引用文档

演示后总结：
- [ ] 强调核心功能
- [ ] 突出技术亮点
- [ ] 说明工程化程度
- [ ] 提供文档位置

