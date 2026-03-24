# B1 校正整合流程说明（Rat dead muscle）

本项目用于批量执行 CEST 数据的 B1 校正，并自动完成验证实验、B1采样数量优化及结果导出。

核心主脚本：`run_B1_all_tasks.m`
核心校正函数：`Z_B1_correction.m`（5D Z 谱）与 `contrast_B1_correction.m`（2D/3D 对比图像）

---

##1. 两张图片表示的含义

###1) `allinone_diag_PCr.png`
该图是**单像素诊断图**（脚本中默认取图像中心像素和中间 offset）。

- 横轴：`B1_actual (μT)`（实际 B1，计算方式为 `B1_input * rel_B1map(pixel)`）
-纵轴：Z 值
- 红色圆点（orig）：该像素在不同 B1 下的原始 Z 值
- 方形连线（B1-corrected）：经过 `Z_B1_correction` 校正后的结果

用途：
- 快速看该像素在 B1维度上的响应是否平滑；
-观察校正前后曲线变化，判断校正是否合理。

---

###2) `validation_experiments_PCr.png`
该图是**三个验证实验 + RMSE 对比**的汇总图（2x2 子图）：

- 左上（实验1）：删除一个偏频点后再校正，与完整数据校正结果比较；
-右上（实验2）：删除一个 B1采样点，用其余 B1 插值重建，再与原始比较；
- 左下（实验3）：在实验2基础上，把大误差点替换回原始值（混合策略）；
-右下：实验1/2/3 的 RMSE 柱状对比。

是否基于单体素：
- 左上/右上/左下：是，均为同一单像素（中心像素）可视化；
-右下 RMSE：不是单像素，而是对整体数据误差 `diff(:)`统计得到。

用途：
-评估流程对“缺失偏频”与“缺失 B1 点”的鲁棒性；
- 比较纯插值与混合替代策略的误差表现。

---

##2.代码流程（对应 `run_B1_all_tasks.m`）

1.选择数据系列：`series='PCr'` 或 `series='Cr'`；
2.设定拟合方式：`fit_type`（默认 `smoothingspline`）；
3. 自动选择数据子目录：`poly-fit-40-water` 或 `poly-fit-water`；
4.读取 `B1-map.mat`，并生成 `rel_B1map`；
5.读取各 B1 条件下 `zspect_dataavecorrect.mat`，构建 `Z_stack(y,x,z,offset,B1)`；
6. 执行基准 B1 校正，输出 `Z_stack_corr_full`；
7.生成单像素诊断图（`allinone_diag_*.png`）；
8.运行三个验证实验并导出 RMSE、重建数据及验证图；
9. 执行 B1 数量优化探索并保存策略结果。

---

##3. 如何使用

###3.1 数据准备
将以下内容放在项目根目录与对应采集子目录下：

- 根目录：`B1-map.mat`
- 每个 B1采集文件夹（例如 `053_CESTEPI_dead_PCr_0.35_...`）下的：
 - `poly-fit-40-water/zspect_dataavecorrect.mat` 或
 - `poly-fit-water/zspect_dataavecorrect.mat`
 - `offset.mat`

> 脚本会自动探测使用哪个子目录。

###3.2 参数设置
在 `run_B1_all_tasks.m` 顶部修改：

- `series`：`'PCr'` 或 `'Cr'`
- `fit_type`：`'smoothingspline' | 'linear' | 'spline' | 'poly2'~'poly5'`
- `run_experiments`：是否执行三个验证实验
- `run_optimization`：是否执行 B1 数量优化

###3.3运行
在 MATLAB 当前目录切到项目根目录后运行：

```matlab
run_B1_all_tasks
```

---

##4.主要输出文件

按 `series` 自动命名（如 `PCr`）：

- `allinone_baseline_PCr.mat`：基准校正结果（含 `Z_stack_corr_full` 等）
- `allinone_diag_PCr.png`：单像素诊断图
- `exp1_results_PCr.mat`：实验1结果
- `exp2_results_PCr.mat`：实验2结果（含插值与重建）
- `exp3_results_PCr.mat`：实验3结果（混合替代）
- `validation_experiments_PCr.png`：三个实验可视化与 RMSE 对比
- `optimization_results_PCr.mat`：B1 数量优化策略结果

---

##5. B1 校正函数说明

### `Z_B1_correction.m`
输入5D Z 谱数据，对每个像素、每个 offset 在 B1维度进行拟合/插值，输出目标 `B1_output` 下的校正值。

- 支持 `smoothingspline`、多项式拟合与线性/样条插值；
- 自动处理 NaN/Inf 与分割掩膜；
- 内部会对传入 B1 map 再做 `/2` 缩放（主脚本已按该逻辑调用）。

### `contrast_B1_correction.m`
用于2D/3D 图像堆栈的同类 B1 校正，思路与上面一致。

---

##6.结果解读建议

-先看 `allinone_diag_*.png`：确认单像素曲线校正后趋势合理；
- 再看 `validation_experiments_*.png`：重点比较实验2和实验3的 RMSE；
- 最后结合 `optimization_results_*.mat`选择可接受误差下的最少 B1采样策略。
