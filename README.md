# PCS Postprocess

电力变流器测试录波的数据处理与绘图。支持 `frt`、`pwr`、`freq_reg`、`inertia`、`volt_prot`、`freq_prot`、`scr` 七种模式。输入来自任何设备或仿真器，只需先转换为[统一 MAT/JSON 格式](docs/input-contract.md)。本仓库不连接硬件，也不生成验收报告。

## 快速开始

```sh
python -m pip install -e .
python examples/make_synthetic.py
pcs-postprocess validate --mode frt --input examples/generated/frt
pcs-postprocess run --mode frt --input examples/generated/frt --output outputs/frt --config examples/project.json
```

`outputs/frt/` 中会有 `processed/`、`figures/` 和 `summary_frt.csv`。PWR 另有 `downsampled/`。可用 `--case 1000,1002` 筛选工况，`--no-fig` 只生成处理数据与汇总。

## 接入自己的数据

1. 将原始数据的信号与单位映射到[输入契约](docs/input-contract.md)；CSV 数据可从 `adapters/csv_to_canonical.py` 开始修改。
2. 为每例生成 `case_<编号>_<名称>.mat` 和同名 `_meta.json`。记录相对录波时间设 `time_basis=recorded`；准确的仿真绝对时间设 `absolute`。
3. 用 `pcs-postprocess validate` 检查数据，再用 `run` 处理。将自己的额定值写入项目配置 JSON，不要使用示例额定值处理真实项目。

伙伴使用 Codex 时，可以按仓库根目录的 [AGENTS.md](AGENTS.md) 开发自己的适配器。运行测试：`python -m unittest discover -s tests -v`。

已有 FRTtestScript 项目的 HIL 字段和 Simulink `SimulationOutput` 均由该项目自己的 `postprocess_adapters/` 转成上述公开格式；本包的运行命令只读取转换后的文件。

## 项目边界

分析计算与图片生成来自 FRTtestScript 的现有离线后处理流程。输出是数据和统计结果，不是产品验收判定。`compare/`、Word 报告、RT-LAB 控制脚本及真实录波不在本仓库中。适配器可放在自己的项目仓库，通过输入契约与本包连接。
