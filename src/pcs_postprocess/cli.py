"""Public command line workflow for canonical PCS recordings."""

import argparse
import importlib
import json
from pathlib import Path
import re

from . import processing
from .config import ProjectConfig
from .contract import validate_case
from .paths import load_paths_config

MODES = tuple(processing.TEST_MODE_CONFIG)


def _case_filter(value: str | None) -> set[int] | None:
    if value is None:
        return None
    try:
        result = {int(part.strip()) for part in value.replace("，", ",").split(",")}
    except ValueError as exc:
        raise ValueError("--case must contain comma-separated integers") from exc
    return result


def _cases(input_dir: Path, mode: str, filter_ids: set[int] | None):
    if not input_dir.is_dir():
        raise ValueError(f"Input directory does not exist: {input_dir}")
    cases = []
    ids = {}
    for path in sorted(input_dir.glob("case_*.mat")):
        # 单例模式先按文件名编号筛选，目录内其他旧文件不影响指定用例。
        match = re.match(r"^case_(\d+)(?:_|$)", path.stem)
        if filter_ids is not None and (match is None or int(match.group(1)) not in filter_ids):
            continue
        meta = validate_case(path, mode)
        case_id = meta["case_index"]
        if filter_ids is not None and case_id not in filter_ids:
            continue
        if case_id in ids and (meta["time_basis"] != "absolute"
                               or ids[case_id] != "absolute"):
            raise ValueError(f"Duplicate recorded case_index {case_id} in {input_dir}")
        # 绝对时间 SIM 不需要按编号求录波偏移；允许同编号、不同名称的原始工况。
        ids[case_id] = meta["time_basis"]
        cases.append((case_id, str(path), meta))
    if not cases:
        raise ValueError(f"No matching case_*.mat files in {input_dir}")
    return cases


def run(mode: str, input_dir: Path, output_dir: Path, config_path: Path,
        filter_ids: set[int] | None = None, no_fig: bool = False,
        save_processed: bool = False) -> int:
    config = ProjectConfig.load(config_path)
    cases = _cases(input_dir, mode, filter_ids)
    processing.configure_project(config)
    recorded = [case for case in cases if case[2]["time_basis"] == "recorded"]
    offsets = ({ci: (0.0, "absolute_time") for ci, _, meta in cases
                if meta["time_basis"] == "absolute"})
    if recorded:
        detected, _, _ = processing.resolve_case_offsets(recorded)
        offsets.update(detected)
    renderer = None
    if not no_fig:
        from . import plot_common
        plot_common.configure_nominal_frequency(config.nominal_frequency_hz)
        renderer = importlib.import_module(f".plot_{mode}", __package__)
        if hasattr(renderer, "F_NOMINAL"):
            renderer.F_NOMINAL = config.nominal_frequency_hz
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for ci, path, meta in cases:
        offset, source = offsets[ci]
        case = processing.prepare_case(path, meta, offset, source)
        # 处理数据供当前运行的指标和绘图使用；仅明确要求时写入额外 MAT。
        if save_processed:
            processing.save_bridged(case, str(output_dir / "processed"))
        if mode == "pwr":
            processing.save_pwr_downsampled(case, str(output_dir / "downsampled"))
        if renderer:
            figures = output_dir / "figures"
            figures.mkdir(exist_ok=True)
            renderer.render_case(case, str(figures))
        rows.append(case.summary)
        print(f"case {ci}: processed")
    rows.sort(key=lambda row: int(row["case_index"]))
    summary = output_dir / f"summary_{mode}.csv"
    processing.write_summary_csv(rows, str(summary), mode)
    print(f"summary: {summary}")
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="pcs-postprocess")
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("validate", "run"):
        command = commands.add_parser(name)
        command.add_argument("--mode", required=True, choices=MODES)
        command.add_argument("--source", help="配置中的数据来源名称，例如 hil 或 sim")
        command.add_argument("--paths-config", type=Path,
                             help="个人目录配置；默认查找当前目录和包项目目录")
        command.add_argument("--input", type=Path)
        command.add_argument("--case")
        if name == "run":
            command.add_argument("--output", type=Path)
            command.add_argument("--config", type=Path)
            command.add_argument("--no-fig", action="store_true")
            command.add_argument("--save-processed", action="store_true",
                                 help="额外保存对时、裁剪后的 MAT；默认不写入")
    args = parser.parse_args(argv)
    try:
        filter_ids = _case_filter(args.case)
        # 显式参数优先；只有缺少目录或额定值路径时才读取来源配置。
        paths = None
        if args.source is not None:
            paths = load_paths_config(args.paths_config)
            source_paths = paths.source(args.source)
        elif args.paths_config is not None:
            raise ValueError("--paths-config requires --source")
        input_dir = args.input or (source_paths.input_dir(args.mode) if paths else None)
        if input_dir is None:
            raise ValueError("Specify --input or --source")
        if args.command == "validate":
            cases = _cases(input_dir, args.mode, filter_ids)
            print(json.dumps({"mode": args.mode, "valid_cases": len(cases)}))
            return 0
        output_dir = args.output or (source_paths.output_dir(args.mode) if paths else None)
        config_path = args.config or (paths.project_config if paths else None)
        if output_dir is None or config_path is None:
            raise ValueError("Specify --output and --config, or use --source")
        return run(args.mode, input_dir, output_dir, config_path,
                   filter_ids, args.no_fig, args.save_processed)
    except (ValueError, KeyError, OSError) as exc:
        parser.exit(2, f"error: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
