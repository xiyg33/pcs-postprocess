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

图中的 P/Q、正序电压、频率、Ip/Iq 等标量曲线使用 20 ms 平滑和 5 ms 采样；三相瞬时波形保留原始采样。此处理只改变绘图，汇总指标仍依据原始数据。PWR 的 P/Q 指标继续使用原有的同参数降采样数据。默认不保存对时、裁剪后的额外 MAT；需要检查中间数组时为 `run` 加 `--save-processed`。

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
