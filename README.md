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

绝对仿真时间输入允许不同文件名共用 `case_index`；两份文件分别绘图，并在汇总中保留两行。录波时间输入仍要求编号唯一，以免事件沿对时混淆。

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

处理算法迁移自 FRTtestScript 的现有离线流程。汇总是数据和统计结果，不构成产品验收判定。SIM/HIL 对比、连续 FRT 手动流程、Word 报告、RT-LAB 控制脚本及真实录波不在本仓库中。生成的 `examples/generated/` 和 `outputs/` 已被 Git 忽略。
