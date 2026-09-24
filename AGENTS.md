# 伙伴数据源接入步骤

本仓库只负责离线指标计算和绘图。硬件控制、真实录波和项目专有路径留在伙伴自己的仓库。

使用 Codex 接入新来源时：

1. 选一例有代表性的原始录波，列出源字段名、单位、采样时间、时间基准、事件参数和缺失通道。先阅读 `docs/input-contract.md`，明确每个源字段对应的标准列。
2. 在伙伴仓库编写独立适配器；通用 CSV 示例位于 `adapters/csv_to_canonical.py`。为每例输出 `data/header` MAT 和同名 `_meta.json`。只读原始文件，转换结果写入另一个目录。
3. 运行 `pcs-postprocess validate --mode <模式> --input <转换目录>`。缺字段或单位错误应在适配器中修正，不要在共享分析代码里猜测数值。
4. 用 `pcs-postprocess run ... --case <编号>` 先处理一例，检查 `summary_<mode>.csv` 和图片；确认无误后再跑整批。为新适配器增加一例合成或脱敏回归数据，并运行 `python -m unittest discover -s tests -v`。

注释请解释单位换算、时间对齐、字段取舍的原因，避免只重复代码表面行为。不要提交生成结果、私有数据，也不要为了适配某一种来源而修改共享指标算法。
