# ParaView 可视化指南

## 1. 打开VTK文件

### 方法1：命令行

```bash
paraview output.vtk
```

### 方法2：图形界面

1. 启动ParaView
2. File → Open
3. 选择 `.vtk` 文件
4. 点击 "Apply"

## 2. 基本可视化操作

### 2.1 显示位移

1. **位移向量**
   - 在左侧面板选择 "Displacement"
   - 类型选择 "Vectors"
   - 调整箭头大小（Glyph → Scale Factor）

2. **位移幅值**
   - 选择 "DisplacementMagnitude"
   - 类型选择 "Surface"
   - 调整颜色映射（Coloring）

### 2.2 显示塑性状态

1. 选择 "PlasticState"
   - 0 = 弹性状态（蓝色）
   - 1 = 塑性状态（红色）

2. 选择 "EquivalentPlasticStrain"
   - 显示等效塑性应变云图
   - 调整颜色范围（Edit Color Map）

### 2.3 显示蠕变应变（蠕变分析）

1. 选择 "EquivalentCreepStrain"
2. 调整颜色映射范围
3. 使用等值线（Contour过滤器）显示等高线

## 3. 高级可视化技巧

### 3.1 变形显示

1. **放大变形**
   - Filters → Warp By Vector
   - 选择 "Displacement"
   - 调整 Scale Factor（例如：100）

2. **原始+变形对比**
   - 复制当前视图
   - 一个显示原始形状，一个显示变形后形状

### 3.2 等值线和切片

1. **等值线（Contour）**
   - Filters → Contour
   - 选择标量场（如DisplacementMagnitude）
   - 设置等值线数量或数值

2. **切片（Slice）**
   - Filters → Slice
   - 选择切平面方向
   - 显示内部状态分布

### 3.3 颜色映射

1. **配色方案**
   - 常用：Rainbow, Cool to Warm, Blue to Red
   - 位移：Cool to Warm（蓝色=小，红色=大）
   - 应变：Rainbow（多色渐变）

2. **调整范围**
   - 右键颜色条 → Edit Color Map
   - 手动设置 Min/Max 值
   - 使用 "Rescale to Data Range" 自动调整

## 4. 蠕变分析动画

### 4.1 加载时间序列

1. 使用通配符加载多个文件：
   ```
   File → Open → output_creep_t*.vtk
   ```

2. 或使用 ParaView 的 "Group Files" 功能

### 4.2 播放动画

1. 在时间控制栏点击 "Play"
2. 调整播放速度（FPS）
3. 使用滑块查看特定时刻

### 4.3 导出动画

1. File → Save Animation
2. 选择输出格式（AVI, PNG序列等）
3. 设置帧率和分辨率

## 5. 对比分析

### 5.1 多结果对比

1. 加载多个VTK文件（不同分析或不同时刻）
2. 使用不同的颜色映射区分
3. 并排显示（Layout → Split View）

### 5.2 与ANSYS结果对比

1. 导入ANSYS结果文件（如 `.rst` 文件）
2. 加载本求解器的VTK文件
3. 使用相同的数据场和颜色映射对比

## 6. 典型可视化场景

### 场景1：线弹性分析结果

**推荐显示：**
- DisplacementMagnitude（位移幅值云图）
- Displacement（位移向量）
- 等值线显示位移分布

**操作步骤：**
1. 加载VTK文件
2. 选择 "DisplacementMagnitude"
3. 应用 Contour 过滤器
4. 调整颜色映射

### 场景2：弹塑性分析结果

**推荐显示：**
- PlasticState（弹性/塑性状态分布）
- EquivalentPlasticStrain（等效塑性应变云图）
- DisplacementMagnitude（变形后的位移）

**操作步骤：**
1. 选择 "PlasticState" 查看塑性区
2. 切换到 "EquivalentPlasticStrain" 查看塑性应变大小
3. 使用 Warp By Vector 显示变形

### 场景3：蠕变分析结果

**推荐显示：**
- EquivalentCreepStrain（蠕变应变云图）
- DisplacementMagnitude（不同时刻的位移）
- 时间序列动画

**操作步骤：**
1. 加载时间序列文件
2. 选择 "EquivalentCreepStrain"
3. 播放动画查看蠕变发展
4. 提取关键时刻进行详细分析

### 场景4：ELBOW290弯管分析

**推荐显示：**
- 三维曲面（使用Surface with Edges）
- 椭圆化变形（通过截面查看）
- 位移云图

**操作步骤：**
1. 确保使用曲面格式输出
2. 选择 "Surface with Edges" 显示网格
3. 使用 Slice 过滤器查看截面
4. 应用变形显示椭圆化

## 7. 导出图像

### 7.1 截图

- File → Save Screenshot
- 选择格式（PNG, JPG等）
- 设置分辨率

### 7.2 保存状态

- File → Save State
- 保存所有可视化设置
- 下次可以直接加载状态文件

## 8. 常见问题

### Q: 颜色条显示不正常？

A: 右键颜色条 → Edit Color Map → 点击 "Rescale to Data Range"

### Q: 变形太小看不清？

A: 使用 Warp By Vector 过滤器，增大 Scale Factor

### Q: 如何同时显示多个标量场？

A: 复制当前视图，在不同视图中显示不同字段

### Q: 动画播放太快/太慢？

A: 调整时间控制栏的播放速度（FPS设置）

## 9. 推荐配置

### 快速预览配置

- 显示方式：Surface
- 颜色：DisplacementMagnitude
- 配色：Cool to Warm
- 自动调整颜色范围

### 详细分析配置

- 显示方式：Surface with Edges
- 添加：Contour（等值线）
- 添加：Warp By Vector（变形显示）
- 多视图对比

### 演示配置

- 高质量渲染（Render View → Render Settings）
- 添加图例（Color Legend）
- 添加标量条（Scalar Bar）
- 导出高分辨率图像

