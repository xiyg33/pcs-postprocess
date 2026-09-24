"""Public command line workflow for canonical PCS recordings."""

import argparse
import importlib
import json
from pathlib import Path

from . import processing
from .config import ProjectConfig
from .contract import validate_case

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
    ids = set()
    for path in sorted(input_dir.glob("case_*.mat")):
        meta = validate_case(path, mode)
        case_id = meta["case_index"]
        if filter_ids is not None and case_id not in filter_ids:
            continue
        if case_id in ids:
            raise ValueError(f"Duplicate case_index {case_id} in {input_dir}")
        ids.add(case_id)
        cases.append((case_id, str(path), meta))
    if not cases:
        raise ValueError(f"No matching case_*.mat files in {input_dir}")
    return cases


def run(mode: str, input_dir: Path, output_dir: Path, config_path: Path,
        filter_ids: set[int] | None = None, no_fig: bool = False) -> int:
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
        command.add_argument("--input", required=True, type=Path)
        command.add_argument("--case")
        if name == "run":
            command.add_argument("--output", required=True, type=Path)
            command.add_argument("--config", required=True, type=Path)
            command.add_argument("--no-fig", action="store_true")
    args = parser.parse_args(argv)
    try:
        filter_ids = _case_filter(args.case)
        if args.command == "validate":
            cases = _cases(args.input, args.mode, filter_ids)
            print(json.dumps({"mode": args.mode, "valid_cases": len(cases)}))
            return 0
        return run(args.mode, args.input, args.output, args.config,
                   filter_ids, args.no_fig)
    except (ValueError, KeyError, OSError) as exc:
        parser.exit(2, f"error: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
