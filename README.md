# PCS Postprocess

面向不同录波来源的 PCS 离线指标计算与绘图工具。支持 `frt`、`pwr`、`freq_reg`、`inertia`、`volt_prot`、`freq_prot`、`scr` 七类测试。伙伴先把设备或仿真数据转换为[统一 MAT/JSON 格式](docs/input-contract.md)，再运行同一条 `pcs-postprocess` 命令。本包不连接硬件，也不生成验收报告。

## 安装并运行合成示例

需要 Python 3.10 或更新版本。在本仓库根目录执行：

```sh
python -m pip install -e .
python examples/make_synthetic.py
pcs-postprocess validate --mode frt --input examples/generated/frt
pcs-postprocess run --mode frt --input examples/generated/frt --output outputs/frt --config examples/project.json
```

第二条命令在 `examples/generated/<mode>/` 为七类模式各生成一对合成 MAT/JSON。`validate` 只检查输入；`run` 读取项目额定值、计算指标并写出结果。上例的额定值只供合成数据使用。

### Windows / VS Code 终端找不到 `pcs-postprocess`

如果 PowerShell 提示“无法将 `pcs-postprocess` 项识别为 cmdlet、函数、脚本文件或可运行程序的名称”，说明终端没有找到命令入口。常见原因是包尚未安装，或包已安装在某个 Python 环境中，但当前终端没有激活该环境。仅进入项目目录不会激活环境。

先检查当前终端使用的 Python 和命令路径：

```powershell
python -c "import sys; print(sys.executable)"
Get-Command pcs-postprocess
```

如果 Python 指向 `Microsoft\WindowsApps\python.exe`，或指向其他环境，而包安装在 Conda 环境中，应先激活安装包的环境。以下以 `claudecode-env` 为环境名；请替换为自己的环境名和 Miniconda 安装目录：

```powershell
# 若当前终端已能使用 conda，可直接执行下一行 activate。
& "C:\path\to\miniconda3\shell\condabin\conda-hook.ps1"
conda activate claudecode-env
pcs-postprocess --help
```

激活后再次检查路径：`python` 应指向该环境的 `python.exe`，`pcs-postprocess` 应指向该环境的 `Scripts\pcs-postprocess.exe`。如果该环境尚未安装包，在本包仓库根目录执行 `python -m pip install -e .`；若本包作为父项目的 `external/pcs-postprocess` 子模块使用，则在父项目根目录执行 `python -m pip install -e external/pcs-postprocess`。

也可以直接使用已安装命令的完整路径，无需激活环境。PowerShell 中需保留路径前的 `&`；以下路径同样需按本机环境替换：

```powershell
& "C:\path\to\miniconda3\envs\claudecode-env\Scripts\pcs-postprocess.exe" --help
& "C:\path\to\miniconda3\envs\claudecode-env\Scripts\pcs-postprocess.exe" run --mode freq_reg --source sim
```

`--source sim` 需要先完成下面的个人配置，并准备该模式的标准 MAT/JSON 输入。

在 VS Code 中，按 `Ctrl+Shift+P`，执行 **Python: Select Interpreter** 并选择安装包的环境，然后关闭旧终端并新建终端。选择解释器通常不会改变已经打开的终端；若新终端仍未自动激活环境，执行上面的激活命令。仅更改 VS Code 的解释器选择，也不会把包安装到所选环境。

## 每个人配置自己的默认目录

在本仓库根目录复制 `postprocess_paths.example.json` 为 `postprocess_paths.local.json`，再复制 `postprocess_project.example.json` 为 `postprocess_project.json`。两份实际配置已由 Git 忽略；填写自己的目录和额定值。例如：

```json
{
  "project_config": "postprocess_project.json",
  "sources": {
    "hil": {"results_root": "data/results_hil", "raw_subdir": "hil_raw", "input_subdir": "canonical", "output_subdir": "postprocess"},
    "sim": {"results_root": "data/results_sim", "input_subdir": "canonical", "output_subdir": "postprocess"}
  }
}
```

相对路径以此配置文件所在目录为准，也可以填写绝对路径。来源名称可自行增加。`results_root` 下每类数据使用 `results_<mode>/`：输入在 `input_subdir`，分析结果在 `output_subdir`。`project_config` 指向额定值配置，`raw_subdir` 供来源适配器定位原始文件；公开包只读取标准输入。三个子目录名均可省略，默认依次为 `canonical`、`postprocess`、`hil_raw`。

从当前目录或父项目目录运行均可：命令先查找当前目录的个人配置，再查找 editable 安装所在项目的配置。安装为普通 wheel 时，可用 `--paths-config` 显式指定配置。

```powershell
pcs-postprocess validate --mode frt --source hil --case 1001
pcs-postprocess run --mode frt --source hil --case 1001
pcs-postprocess run --mode frt --source sim --case 1097
```

也可使用 `--paths-config <文件>` 指定其他位置的个人配置。`--input`、`--output`、`--config` 可以分别覆盖配置中的路径。没有个人配置时，原有显式路径命令继续可用；`--source` 指定了未配置的来源时会报错。

`run --source sim` 写入 `sim_summary_<mode>.csv`，`run --source hil` 写入
`hil_summary_<mode>.csv`，均位于对应输出目录。未指定来源或指定其他来源名称时，
仍写入 `summary_<mode>.csv`。旧汇总文件不会自动删除；读取方应优先使用来源明确的新名称。

## 数据经过哪些后处理

`validate` 只核对每例 MAT/JSON 的配对、模式与事件字段、必需信号及时间轴，不计算指标，也不生成图片。`run` 对每例依次执行：

1. **读取标准数据和额定值。**从 MAT 的 `data/header` 读取时间、三相电压电流、P/Q 等信号，从同名 JSON 读取事件参数，从项目配置读取功率、电压、电流基准值和额定频率。按基准值换算标幺量；正序电压优先使用输入的 `Vpos_pu`，否则由三相电压计算。
2. **对时。**`time_basis=absolute` 的仿真数据保留输入 `t_s`，直接与 JSON 中的事件时刻对应；`time_basis=recorded` 的录波先以首样本为相对时间，再依据实测事件沿求与事件时刻之间的偏移。事件沿无法可靠检出时，使用批次中位偏移或录波起点先验；保护模式按其元数据时序处理。汇总 CSV 的 `rec_offset_s` 和 `offset_source` 记录所用偏移及来源。
3. **检查与裁剪。**按各模式的事件定义检查扰动前稳态等状态，计算相应统计量；将首个事件统一放在输出时间轴的 **1 s**，保留事件前约 1 s 至末次事件后配置的记录窗口。三相瞬时波形在此阶段不做 20 ms 平滑或 5 ms 降采样。
4. **计算指标。**七类模式各自使用对时后的信号及事件窗口计算汇总。一般指标使用原采样数据，绘图专用的平滑与降采样不会反过来改变它们。**PWR 例外：**P/Q 阶跃指标沿用 20 ms 平滑、5 ms 采样后的序列；响应沿检测在降采样前的原采样时间轴上单独进行，检测时可能使用平滑。PWR 使用的降采样 P/Q 同时导出为 CSV。
5. **绘图。**P/Q、正序电压、频率、Ip/Iq 等标量曲线仅在绘图时做 **20 ms 平滑，再映射到 5 ms 时间网格**；三相电压电流瞬时曲线使用裁剪后的原采样。无效采样和原时间缺口保留为空值，不跨缺口补造曲线。`--no-fig` 跳过图片生成。

`run` 默认只保存汇总 CSV、PNG，以及 PWR 的降采样 CSV。需要复查对时和裁剪后的数组时，可加 `--save-processed`，额外保存 `processed/*.mat`；其中是绘图平滑**之前**的数据，不是另一份原始录波或标准输入。

| 输出位置 | 内容 |
|---|---|
| `outputs/frt/processed/` | 对时、裁剪后的 MAT；仅 `--save-processed` 时写入 |
| `outputs/frt/figures/` | PNG 曲线；使用 `--no-fig` 时不生成 |
| `outputs/frt/summary_frt.csv` | 每例一行的汇总指标 |
| `outputs/pwr/downsampled/` | PWR 指标计算所用的降采样 P/Q CSV |

## 七类模式和常用参数

### PWR SIM/HIL 对比

完成两种来源的 canonical 转换和单源分析后，在项目根目录运行：

```powershell
pcs-postprocess compare --mode pwr --case 2001
pcs-postprocess compare --mode pwr --all
```

实现位于包内 `pcs_postprocess.compare_pwr`；也可运行
`python -m pcs_postprocess.compare_pwr --all`。支持 `--paths-config <文件>` 和
`--output <新目录>`。默认输出为当前工作目录的 `result_compare/pwr/<运行时间>/`，
已有输出目录不覆盖。读取配置中 `sim`、`hil` 两种来源的 canonical 和后处理目录。

所有叠加图集中在 `figures/`，每例包含三相电压、三相电流、P、Q。显示范围从
扰动前 1 s 开始，横轴从 0 s 开始、第一事件位于 1 s，保留恢复段和第二阶跃。
HIL 三相统一红色实线，SIM 三相统一蓝色虚线；P/Q 使用同样的来源配色。
三相量仅抽点，P/Q 复用 20 ms 平滑和 5 ms 采样。

`compare_summary_pwr.csv` 每例输出 event1 的四行指标（超调量、上升时间、峰值时间、
稳定时间），直接读取原位置的 `sim_summary_pwr.csv` 和 `hil_summary_pwr.csv`，
不复制汇总、不重新计算指标或检测响应沿。新名称不存在时读取旧 `summary_pwr.csv`，
并兼容已有带来源后缀的指标列。原始录波或分析设置变化后，应先重新运行单源分析。
差值为 SIM − HIL，超调量差值单位为百分点、时间差值单位为秒；缺失值留空。
双方方法、原因、人工确认栏和有效性标记一并导出，估算值不作为正式验收结论。

`manifest.json` 引用实际使用的 summary、canonical 和图片路径，路径相对于输出目录；
跨盘引用使用绝对路径。报告脚本以该目录解析 CSV 的 `figure_png`。重复编号、缺失数据
或参数不一致记录为异常，批处理继续其他工况，有异常时返回非零状态。PWR 对比不生成
Word 报告；FRT 的图和分段偏差计算见下一节。

### FREQ_REG SIM/HIL 对比

先重新生成双方单源汇总，再运行对比；`--no-fig` 可用于只更新指标：

```powershell
pcs-postprocess run --mode freq_reg --source sim --no-fig
pcs-postprocess run --mode freq_reg --source hil --no-fig
pcs-postprocess compare --mode freq_reg --case 3019
pcs-postprocess compare --mode freq_reg --all
```

支持 `--paths-config` 和 `--output <新目录>`。默认输出为
`result_compare/freq_reg/<运行时间>/`，只含 `figures/`、
`compare_summary_freq_reg.csv` 和 `manifest.json`，不覆盖已有目录，不生成 Word。
两份来源汇总优先取 `sim_summary_freq_reg.csv` / `hil_summary_freq_reg.csv`，
只有前缀文件不存在时才回退 `summary_freq_reg.csv`。普通工况的旧汇总缺少
第一阶跃指标时会要求重新运行单源分析。对比只引用汇总，不重新计算指标或检测响应沿。

普通工况每例输出第一段有功阶跃的 `M_p`、`T_up`、`T_p`、`T_s` 四行。
频率事件仍用于录波对齐；计时起点独立检测有功响应，复用 PWR 的响应沿及五点辅助方法，
标记为估算依据，不能当作频率指令时刻。阶跃检测幅值由频率扰动前 0.5 s 与第一保持段
结束前 0.5 s 的实测均值差估算，预期有功方向与频率偏差相反；没有足够数据或方向错误
时给出原因。测量窗口止于恢复前或录波结束，不纳入恢复阶跃。

P 使用 20 ms 平滑和 5 ms 采样。超调量、90% 上升时间与 ±5% 稳定带，分别按各来源
实测起值和末段均值定义，采用 PWR 同一套阶跃算法；`T_s` 优先取第 4 次进入稳定带，
否则回退持续稳定判据，并保留原因及人工确认栏。末段漂移时保留的暂定值也带说明。
`deviation = SIM - HIL`，超调差值单位为百分点，时间差值单位为秒。无效值留空，
状态为 `valid`、`estimated` 或 `unavailable`，不输出正式验收通过结论。

每例两张图沿用单源 HIL 布局：频率／三相电压／三相电流，以及频率／P／Ip／Q／Iq。
HIL 为红色实线、SIM 为蓝色虚线，图例仅列来源；事件前 1 s 对应横轴 0 s，频率事件
对应 1 s，保留恢复段及有功响应延迟。三相仅作绘图抽点，标量不外推、不跨采样缺口补值；
P/Q 使用汇总记录的有功起点分段平滑，Ip/Iq 沿用现有瞬时方向电流推导及显示处理。
CSV 的 `figure_png` 与 `power_current_figure_png` 分别引用两张图。

死区工况不绘图，每例一行 `P_change_abs`，比较已有 `P_change_abs_pu`，同时保留
`sim_5pct_pass` / `hil_5pct_pass`；判据为绝对功率变化 `<= 0.05 pu`。现有事件末段
1 s 均值口径保持不变。缺一方录波时保留另一方死区结果，缺失值为空，清单记录原因。
重复编号、参数不一致及缺失来源均记录异常，批次继续其他工况并返回非零状态。
清单的 `completed` 统计完整配对，`partial` 统计保留单方死区结果的工况。
所有输入及图片引用相对于清单目录，跨盘引用保留绝对路径；CSV 使用 UTF-8 BOM。

#### 普通调频的频率对时与人工复核

普通 HIL 调频在原始录波上确定频率起止时刻，再应用
`t_plot = t_recorder - freq_event_start_recorder_s + 1`。
三相量、频率、功率和方向电流共用这一平移量；有功响应起点独立计算，保留相对频率扰动
的延迟。SIM 绝对时间保持原模型时钟，死区、惯量、PWR 和 FRT 沿用原流程。

PLL 频率做 20 ms 平滑、5 ms 采样，以扰动前中位数和 MAD 估计噪声。
阈值仅确认变化区间，不把阈值越过时刻当成扰动开始；起止沿分别拟合
“基线—线性过渡—平台”，在单例名义事件前 0.3 s、后 0.35 s 内以 5 ms 间隔搜索，
过渡时长为 10–350 ms。检查方向、幅值、起止间隔、残差、窗口敏感性和边界。
普通调频不再使用跨工况中位数回退；死区的历史回退与计算行为保留。

频率对时优先使用同步指令标记，其次使用人工确认值，最后使用通过检查的自动估计。
检查失败时，图暂以该案例自己的名义起点归一化，明确标记 `review_required`。
实际有功起点始终在完整录波和名义起点前 0.3 s 至后 0.35 s 的范围内独立检测；
功率基线先取搜索区间之前的样本，末段测量止于恢复前。功率明显早于估计频率起点时
触发复核，不能通过改用功率起点掩盖对时问题。

单源汇总增加 `freq_event_start_recorder_s`、`freq_event_end_recorder_s`、
`freq_timing_method`、`freq_timing_confidence`、`freq_timing_review_required`、
`freq_timing_reason`、候选起止时刻、噪声和残差、有功原始起点、功率测量结束时刻及复核图路径。
置信状态为 `absolute`、`synchronized`、`manual_confirmed`、`accepted_estimate` 或
`review_required`。自动或人工频率估计均不是同步指令时刻。

正常绘图运行会为待确认案例输出 `postprocess/timing_review/*_timing_review.png`，
包含原始 PLL 频率、平滑曲线、拟合过渡、当前锚点、拟合候选及独立有功响应起点。
`--no-fig` 禁止包括复核图在内的全部绘图。对比 CSV 保留双方对时字段，清单中的
`timing_review_cases` 列出待确认编号；存在待确认值时 `run_status=needs_timing_review`，
该批结果只能作为待复核结果，不能当作已确认的对时结果。

伙伴适配器可将人工确认的 `freq_event1_start_recorder_s`、`freq_event1_end_recorder_s`
及 `freq_timing_override_source=manual_confirmed` 写入 canonical 元数据。
起止值必须在当前录波时间轴内且开始早于结束；同步指令使用现有
`hil_command_event1_recorder_s` / `hil_command_time_source=synchronized`，可额外提供
`hil_command_event1_end_recorder_s`。不要把主机调度时间或 PLL 拟合标成同步指令。

### FRT SIM/HIL 对比与 F1/F2/F3/FG 偏差

```powershell
pcs-postprocess compare --mode frt --case 1001
pcs-postprocess compare --mode frt --all
```

普通 `LVRT_SINGLE/HVRT_SINGLE` 工况在 `result_compare/frt/<运行时间>/` 输出两张叠加图：
`figures/case_<编号>_sim_hil_overlay.png` 包含正序电压、三相电压和三相电流；
`*_overlay_power_current.png` 包含正序电压、P、Ip、Q、Iq。扰动位于 1 s、横轴从
0 s 开始；HIL 红色实线、SIM 蓝色虚线。电流方向计算和平滑复用单源流程，低于
0.05 pu 电压的无效方向电流保留缺失值。两边实际恢复时刻分别标记，图不拉伸时间。

`compare_summary_frt.csv` 每例五行（正序电压、P、Q、Ip、Iq），单位均为 pu；沿用参考偏差表
的九个字段。单源 summary 只读取现有对时信息，不复制或重写，其他偏差由成对波形计算。
优先引用 `sim_summary_frt.csv` / `hil_summary_frt.csv`，兼容旧 `summary_frt.csv`。

| 段 | 指标窗口 |
| --- | --- |
| A 扰动前稳态 | 0–0.9 s，即扰动前 1 s 至前 0.1 s |
| B1 故障开始暂态 | 从 1 s 开始，长度为双方较短实测故障时长的 1/4，最多 0.2 s |
| B2 故障中稳态 | B1 之后至较早恢复时刻前；保护间隔为较短故障时长的 10%，最多 0.1 s |
| C1 恢复暂态 | 各自恢复起点后的相同时间段，长度与 B1 一致 |
| C2 恢复后稳态 | C1 之后至各自恢复后 2 s，只用共同数据覆盖的部分 |

恢复时刻复用 FRT 恢复沿检测，检测失败回退元数据并标记估算。比较使用共同 5 ms
网格和半开窗口，标量沿用 20 ms 平滑，不外推或跨缺失采样补值。
F1 为 A/B2/C2 的均值差绝对值，F2 为 B1/C1 的平均绝对误差，F3 为 A/B2/C2 的最大
绝对误差，FG 为五段按有效采样数加权的平均绝对误差。缺失采样或空窗口对应指标留空；
任一段不可计算时 FG 留空。`status` 是数据可用性说明，不表示正式验收通过。
Ip/Iq 指标使用功率电流叠加图中的同一条曲线：三相电压定向瞬时电流除以配置的峰值
电流基准，再经 20 ms 平滑和 5 ms 降采样。低于 0.05 pu 电压的无效电流保留 NaN，
受影响窗口的指标及该电流的 FG 留空，其他信号仍独立计算。

故障持续时间配置不同时仍绘图、计算。清单的 `duration_assessment` 比较恢复分别对齐与
只对齐起点的九项指标，并用较长实测故障时长构造暂态/保护间隔作窗口敏感性检查。
记录指标变化和采样数变化；`affected` 表示观察到影响，`no_observed_effect` 表示未见
影响（数值比较容差 1e-12 pu，不是验收限值）。CSV 的 `duration_effect` 保留该结论。

`LVRT_BND/HVRT_BND` 扫描工况和缺少 SIM 的工况写入 `manifest.json` 的 `skips`，
不算运行失败；没有永久编号上限。重复编号、参数冲突、数据损坏等进入 `errors`，
整批继续其他工况并返回非零状态。清单还记录每例窗口、有效采样数、对时依据和输入引用。

合成数据可逐类运行。以下 PowerShell 示例将七类结果分别放入自己的目录：

```powershell
$modes = @('frt', 'pwr', 'freq_reg', 'inertia', 'volt_prot', 'freq_prot', 'scr')
foreach ($mode in $modes) {
    pcs-postprocess validate --mode $mode --input "examples/generated/$mode"
    pcs-postprocess run --mode $mode --input "examples/generated/$mode" --output "outputs/$mode" --config examples/project.json
}
```

只处理一例时用 `--case`，编号来自 JSON 的 `case_index`；可用英文或中文逗号填写多个编号。只需要数据和汇总、不需要图片时加 `--no-fig`：

```sh
pcs-postprocess run --mode frt --input examples/generated/frt --output outputs/one-frt --config examples/project.json --case 1000
pcs-postprocess run --mode pwr --input examples/generated/pwr --output outputs/pwr-metrics --config examples/project.json --case 1001 --no-fig
```

SIM 和 HIL 的 `freq_reg` 死区案例（`FREQ_REG_DEADBAND`）以及 FRT 边界案例
（`LVRT_BND` / `HVRT_BND`）默认不生成任何单案例图片，但仍计算全部指标并写入
summary；`--save-processed` 仍可保存处理后的 MAT。识别依据是元数据的类型或
案例名称前缀，不依据编号或频差。单独指定 `--case` 不会恢复绘图；需要这些案例的
图片时显式加 `--include-summary-only-figures`：

```sh
pcs-postprocess run --mode frt --source sim --case 1097 --include-summary-only-figures
```

`--no-fig` 优先于上述开关，始终禁止全部绘图。默认跳过不会删除已有图片。
这些规则适用于本包的 `run` 入口，旧 Python/MATLAB 后处理入口保持原行为。

绝对仿真时间输入允许不同文件名共用 `case_index`；两份文件按上述规则分别决定是否绘图，并在汇总中保留两行。录波时间输入仍要求编号唯一，以免事件沿对时混淆。

## 接入自己的录波

1. 按[输入契约](docs/input-contract.md)确认时间轴、三相电压电流、P/Q、频率及事件字段的名称和单位。每例输出 `case_<编号>_<名称>.mat` 与同名 `_meta.json`。
2. 在 JSON 中设置 `mode` 和 `time_basis`。录波器相对时间使用 `recorded`，由实测事件沿对时；与事件元数据共用仿真时钟时使用 `absolute`，原始 `t_s` 保持不变。
3. 把项目额定值写进自己的配置 JSON，先执行 `validate`，再用 `--case` 跑一例并检查汇总与图片，最后处理整批。

下例是 FRT 的最小事件元数据；MAT 仍须包含输入契约列出的必需信号：

```json
{
  "mode": "frt",
  "case_index": 1000,
  "case_name": "example",
  "time_basis": "absolute",
  "fault_start_s": 2.0,
  "fault_duration_s": 0.2,
  "V_fault_pu": 0.5,
  "P_ref_pu": 0.3
}
```

若已有列名和单位均符合契约的 CSV，可从示例适配器开始：

```sh
python adapters/csv_to_canonical.py --csv path/to/wave.csv --meta path/to/case_meta.json --output outputs/csv-input
pcs-postprocess validate --mode frt --input outputs/csv-input
pcs-postprocess run --mode frt --input outputs/csv-input --output outputs/csv-result --config path/to/project.json
```

CSV 第一行应是 `t_s,Va_V,Vb_V,Vc_V,Ia_A,Ib_A,Ic_A,P_W,...` 等标准列名。适配器读取 `--meta` 指定的 JSON，并按其中编号和名称生成 MAT/JSON 文件。原始列名或单位不同的伙伴，应在自己的项目中改写适配器，不要修改共享指标算法。

FRTtestScript 项目的 HIL 字段和 Simulink `SimulationOutput` 由其本地 `postprocess_adapters/` 转为该格式；新包只读取转换结果。伙伴使用 Codex 接入其他来源时，可按 [AGENTS.md](AGENTS.md) 的步骤开发和验证适配器。

## 测试与边界

```sh
python -m unittest discover -s tests -v
```

处理算法迁移自 FRTtestScript 的现有离线流程。汇总是数据和统计结果，不构成产品验收判定。SIM/HIL 对比读取离线录波和单源汇总；连续 FRT 手动流程、Word 报告、RT-LAB 控制脚本及真实录波留在伙伴仓库。生成的 `examples/generated/` 和 `outputs/` 已被 Git 忽略。

## INERTIA / VOLT_PROT SIM/HIL 对比

```powershell
pcs-postprocess run --mode inertia --source sim
pcs-postprocess run --mode inertia --source hil
pcs-postprocess compare --mode inertia --case 4002
pcs-postprocess compare --mode inertia --all

pcs-postprocess run --mode volt_prot --source sim
pcs-postprocess run --mode volt_prot --source hil
pcs-postprocess compare --mode volt_prot --case 5005
pcs-postprocess compare --mode volt_prot --all
```

个人路径仍由 `--paths-config` 指定；默认输出为
`result_compare/<mode>/<运行时间>/`，含 `figures/`、
`compare_summary_<mode>.csv` 和 `manifest.json`。每次使用新目录。
compare 只读取单源指标，不重新计算；优先读取 `<source>_summary_<mode>.csv`，
缺失时读取 `summary_<mode>.csv`。旧字段或旧指标版本会要求重算，不能静默使用旧结果。
单例 run 会重写来源汇总，整批 compare 前必须完整 run。

INERTIA 每例两张叠加图：频率/Vabc/Iabc，以及频率/P/Ip/Q/Iq。
图和指标共用实测频率斜坡起点，归一化为 1 s。HIL 使用频率斜坡拟合；
SIM 同样检测斜坡起点，但保留输入的绝对时间轴。不能确认拟合时，仅用名义起点显示，
四指标留空。P 先分段进行 20 ms 平滑，再映射至 5 ms 网格，避免恢复响应污染首次响应。
观察段包含斜坡与保持；结束取计划边界、实测频率恢复前边界、录波结束的最早值。

INERTIA 汇总的 `M_p_pct` 为相对实测变化量的超调百分比，
`T_up_s` 为从频率斜坡起点到实测变化量 90% 的时间，
`T_p_s` 为响应方向峰值时间，`T_s_s` 为持续进入末值 ±5% 实测变化幅值带的时间。
基线与末值分别取起点前及观察段末尾 0.5 s；末值需通过漂移和半窗均值检查。
这里描述频率斜坡下的有功响应，不是功率指令阶跃验收。
稳定时间不采用第 4 次入带估算。缺样、变化过小、方向异常、未确认末值、
未持续稳定或录波不足均保留空值和原因；已有可靠峰值时间独立保留。
每例四行，双方数值来自单源汇总，`deviation=SIM-HIL`。
超调差值单位为 `percentage_points`，其余为 `s`。

VOLT_PROT 每例一张 V+/P/Q 三行叠加图、一行 breaker 响应时间对比。
不绘制 breaker 子图；跳闸时刻以竖线标出，第一行分别标注双方 t2 和 Δt。
电压沿沿用 FRT 的 20 ms 平滑、目标方向、自适应阈值和持续越限判据，
不是同步指令时刻；双方均标记 `measured_voltage_edge` 估算。
`t1_s=1`，在原始 breaker 时间轴查找此后的 1→0 下降沿，
要求至少 10 个连续原始采样低于 0.5；短毛刺和已经断开的状态不能作为本次响应。
`delta_t_s=t2_s-t1_s`，缺少实测起点时不计算 Δt。

`protection_status` 区分 `tripped`、`no_trip`、`missing_breaker`、
`invalid_breaker`、`disturbance_unconfirmed`、`already_open`、`insufficient_observation`。
完整观察窗口内未跳闸是正常可记录结果；一方跳闸而另一方未跳闸时，
`response_match=trip_mismatch`，差值为空，不伪造零值，也不判为处理失败。
缺通道等输入/检测问题写入 manifest 并返回非零状态。
原始输入时刻、归一化转换、指标版本及来源文件保存在汇总和清单中。
CSV 使用 UTF-8 BOM；重复编号、缺失来源/汇总行及参数不一致记录为异常，批次继续。

## SCR：10 s 扰动与边界代表例

本项目用例表已改为扰动起点 10 s、持续 0.2 s、总时长 12.2 s，
record_pre=10 s、record_post=2 s。新 SIM 数据需重新测试、重新导出；
不能把旧 2 s 录波改写元数据后当作新数据。
路径配置可设置 `"compare_cases": {"scr": [7001, 7033, 7034]}`。

```powershell
pcs-postprocess run --mode scr --source sim --case 7001,7033,7034
pcs-postprocess run --mode scr --source hil --case 7001,7033,7034
pcs-postprocess compare --mode scr --all
```

每例一张 Vabc/Iabc/P 对比图，并分别保存双方状态及其原因。
单源和 compare 都保留实际录波时间轴，按配置显示 0～12.2 s、标注设定扰动 10 s；
实测电压下降、恢复沿另行保存在汇总中，完整启动和提前脱网不被裁掉。
HIL 边界代表例为 P100 的 7001 不正常、Pm100 的 7033 正常和 7034 不正常；
SIM 根据新录波重新判断，不预设状态。
新 `scr_common_horizon_v2` 汇总保存设定时序和实测沿；compare 检查双方设定起点、
持续时间、总时长、记录窗口以及单源汇总的一致性，不允许混用新旧时序。

稳定判据使用相同的恢复后观察时长，受配置终点和较短录波限制，上限 8 s；
跳过恢复后 0.3 s，剩余至少 1 s。P 经 20 ms 平滑和 5 ms 采样，
后/前半段去均值 RMS ≥1.3 且后段 RMS ≥max(0.01 pu, 3×基线标准差)，
或 |P|≥2 pu 持续 20 ms，判为振荡失稳。持续 10 个原始低电平采样的
breaker 下降沿记为脱网；扰动前脱网是适应性失败，不证明振荡失稳。
观察或检测不足记为无法判定。`compare_boundary_scr.csv` 只报告已测区间，
不跨未知点推断，不给稀疏或非单调结果虚构临界 SCR。

## FREQ_PROT：实测线路扰动至脱网总响应时间

```powershell
pcs-postprocess run --mode freq_prot --source sim
pcs-postprocess run --mode freq_prot --source hil
pcs-postprocess compare --mode freq_prot --case 6017
pcs-postprocess compare --mode freq_prot --all
```

**Δt = 持续 breaker 下降沿时刻 − 实测线路频率扰动起点。**
线路起点来自三相电压旋转相量、额定周期的对称平滑及连续分段相位拟合，
校验方向、基线、残差、搜索边界和观察长度，并排除跳闸对拟合的污染。
不再用 PLL 达到保护门槛的时间，也不依赖 `frequency_protection_reference`。
新版本为 `voltage_disturbance_breaker_v2`，旧门槛计时汇总必须重算。

双方图上实测起点归一化为 1 s，breaker 保持原始采样检测；
每例一张 f/P/Q 三行图、一行主 CSV，差值为 SIM−HIL。
双方跳闸也保留完整扰动及扰动后录波，不自动缩短横轴。当前 60 s 扰动从 1 s 到 61 s，
含 2 s 扰动后录波的完整数据图显示 0～63 s；录波不足时只展示已有数据。
不绘制 breaker 子图；跳闸时刻以竖线标出，第一行分别标注双方 t2 和 Δt。
线路频率为主曲线，录波 PLL 可作辅助。总响应包含测频动态、保护计时及输出链路，
不解释为控制器内部定时器的配置延时。
缺少实测起点、缺 breaker、扰动前已断开、观察不足或确实未跳闸时留空并说明原因。
PLL 晚越限不再造成空值；正常未跳闸和一方跳闸不算处理错误。
诊断 CSV 保存输入起止时刻和拟合质量，不重复输出固件参数 JSON。

## 对比归档与单源覆盖

compare 默认输出 `result_compare/<mode>/<运行时间>/`，每次新建时间戳目录。
单源 run 继续覆盖各自 `postprocess/` 中的同名图片与汇总，不增加时间戳。
compare 不修改单源结果。显式 `--output` 可指定一个新的对比目录。
日常查看 `figures/`、`compare_summary_<mode>.csv` 和边界/计时诊断 CSV。
每次对比最多一份 `manifest.json`，仅用于输入来源、版本、对时依据和异常追溯；
不生成额外调试/验收 JSON，不复制原始 metadata JSON，不在说明中转储参数 JSON。
旧历史结果保留，按新的运行时间戳识别本次输出。

## Python 对比接口：筛选与展示配置

各模式的 `compare_<mode>.run(...)` 可选传入 `case_ids=[...]` 与
`presentation={case_index: {...}}`。不传时保持既有命令的行为；`case_ids` 与单个
`case_id` 互斥，必须是非空、唯一的整数列表。

展示项可指定 `display_case_name`、`title`、`reported_P_ref_pu`、
`power_axis_limits`、`power_tick_pu` 和 `broad_q_axis`。新名称用于安全的 PNG 文件名；
标题保留原图中的状态说明。功率坐标扩大到容纳全部显示曲线，竖直事件标记不参与取值。
绘图使用原始参数和既有信号处理结果，展示功率不参与配对校验、参考值或指标计算。

指定展示配置时，对比 CSV 与逐案例清单追加 `display_case_name`、
`reported_P_ref_pu`、`source_P_ref_pu` 和 `approximation`，原有字段保留原值。
项目专有的案例编号、近似许可和连续录波配对由伙伴仓库管理，本包不硬编码这些选择。
