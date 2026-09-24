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
- `"time_basis": "absolute"`：`t_s` 与事件元数据使用同一仿真时钟。处理器直接使用原始时间轴，不做事件沿修正；首个采样点可以不是零。

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
