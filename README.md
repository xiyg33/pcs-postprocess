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

| 输出位置 | 内容 |
|---|---|
| `outputs/frt/processed/` | 对时、裁剪后的 MAT |
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
