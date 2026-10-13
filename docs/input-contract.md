# 输入契约（v0.1）

每个工况由一对同名文件组成：`case_<整数编号>_<名称>.mat` 和 `case_<整数编号>_<名称>_meta.json`。例如 `case_1001_example.mat` 对应 `case_1001_example_meta.json`。MAT 使用 MATLAB v5/v7 格式，包含：

- `data`：`N × C` 数值矩阵，至少两行，每列是一条信号。
- `header`：`C × 2` 元胞数组，第一列是与 `data` 列顺序对应的信号名；第二列可留空。

所有通道先转换成下表单位。`t_s` 必须是有限且严格递增的秒数；必需信号不允许 NaN 或 Inf。缺少必需字段时，`pcs-postprocess validate` 会指出文件和字段名。

| 信号名 | 单位与含义 | 要求 |
|---|---|---|
| `t_s` | s，输入录波时间轴 | 七类模式必需 |
| `Va_V`、`Vb_V`、`Vc_V` | V，三相瞬时相电压 | 七类模式必需 |
| `Ia_A`、`Ib_A`、`Ic_A` | A，三相瞬时相电流 | 七类模式必需 |
| `P_W` | W，有功功率 | 七类模式必需 |
| `Q_var` | var，无功功率 | PWR 必需，其他模式可选 |
| `f_hz` | Hz，频率 | 调频、惯量、频率保护必需 |
| `Vpos_pu` | pu，正序电压幅值 | 可选；未提供时由三相电压计算 |
| `breaker` | 0/1，断路器状态 | 可选；保护跳闸时间会使用 |

## JSON 时间基准

JSON 必须包含 `mode`、非负整数 `case_index`、非空 `case_name` 和 `time_basis`。文件名中的编号必须与 `case_index` 一致。

- `"time_basis": "recorded"`：`t_s` 是录波器的相对时间。处理器根据实测事件沿对齐元数据中的事件时间；若无法可靠检测，可使用项目配置中的 `record_start_nominal_s` 回退。
- `"time_basis": "absolute"`：`t_s` 与事件元数据使用同一仿真时钟，首个采样点可以不是零。原始时间轴不平移；各模式的指标起点按模式定义检测，显示时间归一化到 1 s。

`record_pre_s`、`record_post_s` 和 `sim_time_s` 是可选的观察窗口信息。每类测试还需要下列数值事件字段：

| `mode` | 事件开始 / 持续时间 | 其他必需字段 |
|---|---|---|
| `frt` | `fault_start_s` / `fault_duration_s` | `V_fault_pu`、`P_ref_pu` |
| `pwr` | `step_start_s` / `step_dur_s` | `P_base_pu`、`P_step_pu`、`Q_base_pu`、`Q_step_pu` |
| `freq_reg` | `event_start_s` / `event_hold_s` | `P_ref_pu`、`freq_deviation_hz` |
| `inertia` | `event_start_s` / `event_hold_s` | `P_ref_pu`、`freq_deviation_hz`、`freq_ramp_rate_hz_per_s` |
| `volt_prot` | `event_start_s` / `event_hold_s` | `P_ref_pu`、`V_target_pu` |
| `freq_prot` | `event_start_s` / `event_hold_s` | `P_ref_pu`、`f_target_hz` |
| `scr` | `disturbance_start_s` / `disturbance_duration_s` | `P_ref_pu`、`SCR_value` |

PWR 有第二次阶跃时，再添加 `step2_start_s`、`step2_dur_s`、`P_step2_pu` 和 `Q_step2_pu`。`test_type_code=4` 表示 Q 阶跃；其他 PWR 代码按 P 阶跃计算。

## 项目额定值与检查

项目配置 JSON 必须提供正数 `power_base_w`、`voltage_base_peak_v`、`current_base_peak_a` 和 `nominal_frequency_hz`；`record_start_nominal_s` 可选。它们分别用于功率、电压、电流的 pu 换算、频率基准和录波对时回退，不应把示例额定值直接用于真实项目。

```sh
pcs-postprocess validate --mode frt --input examples/generated/frt
pcs-postprocess run --mode frt --input examples/generated/frt --output outputs/frt --config examples/project.json
```

`validate` 检查输入文件；`run` 才生成处理后 MAT、PNG 和 `summary_<mode>.csv`。完整的合成输入可由 `python examples/make_synthetic.py` 生成。

## inertia 与 volt_prot 的指标时间轴

`inertia_event_start_input_s`、`inertia_metric_end_input_s` 和
`protection_event_start_input_s`、`protection_trip_input_s` 使用输入时间轴：
recorded 是从首样本计时的录波相对秒数，absolute 是原始仿真秒数。
绘图时统一按 `t_normalized = t_input - t_event_input + 1` 转换。
volt_prot 的 `t1_s=1` 只在实测电压沿已确认时成立；`t2_s` 使用同一变换，
`delta_t_s=t2_s-t1_s`。未确认沿时，名义时间仅用于展示，不产生 Δt。
`breaker` 保持可选，但缺失会明确标记 `missing_breaker`，不能视为正常未跳闸。

## SCR 扰动、记录窗口和稳定边界

新对比需要 `V_target_pu`、`disturbance_start_s`、`disturbance_duration_s`、
`sim_time_s`、`record_pre_s` 和 `record_post_s`。本项目新表为 10/0.2/12.2/10/2 s。
手动录波的 `scr_event_search_center_input_s` 是输入轴上的搜索先验，不能冒充实测沿。
`scr_configured_start_input_s` / `scr_configured_end_input_s` 是设定扰动窗，
`scr_event_start_input_s` / `scr_recovery_input_s` 是实测电压沿，
`scr_display_end_input_s` 控制实际时间轴的图片范围。两类时刻分别记录。
`scr_metrics_version=scr_common_horizon_v2`；compare 要求双方和汇总的时序配置一致。
不要把旧录波的 2 s 起点改成 10 s；必须重新测试、导出和分析。
单源/对比图保留启动和提前脱网过程，稳定判断采用双方共同观察窗。

路径配置可选 `"compare_cases": {"scr": [7001, 7033, 7034]}` 指定代表例，
未配置则使用 HIL 已有编号。共享包不硬编码项目编号。

## FREQ_PROT 实测线路扰动计时

`f_target_hz` 表示试验频率目标。主时钟由 `Va_V/Vb_V/Vc_V` 的相位变化估算，
不由 `f_hz` 或保护门槛选取；`f_hz` 仍按原输入契约提供，用于辅助曲线和基线检查。
已有 `frequency_protection_reference` 元数据可保留，但不参与新计时。

`protection_event_start_input_s` 为已确认的实测线路扰动起点，
`protection_trip_input_s` 为其后的持续 breaker 下降沿。
`delta_t_s` 为两者之差；归一化后 `t1_s=1`、`t2_s=1+delta_t_s`。
`protection_display_start_input_s` 在检测失败时可使用名义时刻作展示，不能产生 Δt。
`freq_fit_*` 字段记录拟合依据与质量，`freq_input_sampling_s` 记录输入采样周期。
新版本 `voltage_disturbance_breaker_v2` 与旧 PLL 门槛口径不兼容，必须重新分析。
未跳闸、缺通道、起点未确认、扰动前已断开、观察不足分别记录状态和原因。
