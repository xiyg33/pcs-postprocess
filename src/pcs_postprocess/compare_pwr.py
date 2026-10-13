"""PWR overlays from one second before the event and referenced report metrics."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
import json
import math
import os
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pcs_postprocess import processing
from pcs_postprocess.config import ProjectConfig
from pcs_postprocess.contract import validate_case
from pcs_postprocess.paths import load_paths_config
from pcs_postprocess.plot_common import pwr_power_axis_limits, PWR_POWER_TICK_STEP_PU
from matplotlib.ticker import MultipleLocator
from pcs_postprocess.compare_presentation import (PRESENTATION_FIELDS, apply_presentation,
                                                 figure_path, presentation_record, render_plot, select_case_ids)


METRICS = (
    ("M_p", "超调量", "M_p_pct", "%", "percentage_points"),
    ("T_up", "上升时间", "T_up_s", "s", "s"),
    ("T_p", "峰值时间", "T_p_s", "s", "s"),
    ("T_s", "稳定时间", "T_s_s", "s", "s"),
)
FIELDS = ["case_index", "case_name", "test_type", "event", "signal", "metric",
          "parameter", "unit", "sim_value", "hil_value", "deviation",
          "absolute_deviation", "deviation_unit", "status"]
for _side in ("sim", "hil"):
    FIELDS += [f"{_side}_{key}" for key in (
        "status", "reason", "step_timing_basis", "step_response_start_s",
        "T_s_method", "metric_sampling_ms", "manual_T_s")]
FIELDS += ["figure_png"]


def relative(path: Path, output: Path) -> str:
    try:
        return Path(os.path.relpath(path.resolve(), output.resolve())).as_posix()
    except ValueError:
        # Windows 的不同盘符之间无法表示相对路径。
        return path.resolve().as_posix()


def finite(value, label: str, optional: bool = False) -> float | None:
    if optional and (value is None or str(value).strip().lower() in ("", "nan", "none")):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid {label}: {value!r}") from exc
    if not math.isfinite(number):
        raise ValueError(f"Non-finite {label}: {value!r}")
    return number


def summary_path(directory: Path, side: str) -> Path:
    preferred = directory / f"{side}_summary_pwr.csv"
    path = preferred if preferred.is_file() else directory / "summary_pwr.csv"
    if not path.is_file():
        raise FileNotFoundError(f"Missing {preferred} (legacy summary_pwr.csv also absent)")
    return path


def read_summary(path: Path, side: str) -> dict[int, list[dict]]:
    rows = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"case_index", "case_name", "test_type", "rec_offset_s",
                    "offset_source", "step_timing_basis", "step_response_start_s",
                    "T_s_method", "metric_sampling_ms", "人工确认Ts"}
        columns = set(reader.fieldnames or [])
        if required - columns:
            raise ValueError(f"{path}: missing columns {sorted(required - columns)}")
        aliases = {}
        for metric, _, column, _, _ in METRICS:
            # 已有的来源汇总曾把来源写进指标列名；读取时兼容，不重写文件。
            legacy = column.replace(metric + "_", metric + f"_{side}_", 1)
            aliases[column] = column if column in columns else legacy
            if aliases[column] not in columns or metric + "_reason" not in columns:
                raise ValueError(f"{path}: missing {column} / {metric}_reason")
        for row in reader:
            case_id = int(row["case_index"])
            for column, alias in aliases.items():
                row[column] = row[alias]
            rows.setdefault(case_id, []).append(row)
    return rows


def recordings(directory: Path) -> dict[int, list[Path]]:
    if not directory.is_dir():
        raise FileNotFoundError(directory)
    result = {}
    for path in sorted(directory.glob("case_*.mat")):
        match = re.match(r"case_(\d+)(?:_|$)", path.stem)
        if match:
            result.setdefault(int(match.group(1)), []).append(path)
    return result


def one(mapping, case_id: int, label: str):
    found = mapping.get(case_id, [])
    if len(found) != 1:
        raise ValueError(f"{label}: case {case_id} requires one entry, found {len(found)}")
    return found[0]


def validate_pair(sim: dict, hil: dict) -> None:
    for key in ("case_index", "test_type", "test_type_code"):
        if sim.get(key) != hil.get(key):
            raise ValueError(f"SIM/HIL parameter mismatch: {key}")
    keys = ["P_base_pu", "P_step_pu", "Q_base_pu", "Q_step_pu", "step_start_s", "step_dur_s"]
    if int(sim["test_type_code"]) == 5:
        keys += ["step2_start_s", "step2_dur_s", "P_step2_pu", "Q_step2_pu"]
    if int(sim["test_type_code"]) not in (3, 4, 5):
        raise ValueError("Unsupported PWR test_type_code")
    for key in keys:
        a, b = finite(sim.get(key), key), finite(hil.get(key), key)
        if not np.isclose(a, b, atol=1e-8, rtol=0):
            raise ValueError(f"SIM/HIL parameter mismatch: {key} ({a} vs {b})")


def full_signals(path: Path, meta: dict, summary: dict) -> dict:
    """Use saved alignment, never calculate or redetect source summary metrics."""
    offset = finite(summary["rec_offset_s"], "rec_offset_s")
    ev = processing.event_times_from_meta(meta)
    timing = processing.event_timing_metadata(meta, ev, offset, summary["offset_source"])
    raw_time, signals = processing.load_case_file(str(path), absolute=meta["time_basis"] == "absolute")
    t = raw_time + timing["recorder_to_normalized_shift_s"]
    events = {key: value - ev["event1_start"] + 1.0 for key, value in ev.items()
              if key in ("event1_start", "event1_end", "event2_start", "event2_end") and value is not None}
    start = finite(summary["step_response_start_s"], "step_response_start_s", optional=True)
    command = timing["command_event1_normalized_s"]
    if summary["step_timing_basis"] == "synchronized_command" and np.isfinite(command):
        end = command + float(meta["step_dur_s"])
    else:
        end = events["event1_end"]
    if "event2_start" in events:
        end = min(end, events["event2_start"])
    end = min(end, t[-1])
    data = {"t_s": t, "P_pu": signals["P"] / processing.HIL_P_BASE_W,
            "Q_pu": signals["Q"] / processing.HIL_P_BASE_W}
    # 先用全部输入平滑，再裁掉显示范围外的样本，避免改变 0 s 附近的滤波结果。
    shift = max(0.0, -math.floor(t[0] / processing.PWR_METRIC_SAMPLE_S) * processing.PWR_METRIC_SAMPLE_S)
    data["t_s"] = t + shift
    power = processing.pwr_downsampled_signals(
        data, None if start is None else start + shift, end + shift)
    power["t_s"] -= shift
    keep = t >= 0.0
    if not np.any(keep):
        raise ValueError("No recording samples from one second before the first event")
    power_keep = power["t_s"] >= 0.0
    power = {key: values[power_keep] for key, values in power.items()}
    return {"t_s": t[keep],
            "Vabc_pu": np.column_stack([signals[f"V{p}"] for p in "abc"])[keep] / processing.HIL_V_BASE_PEAK_V,
            "Iabc_pu": np.column_stack([signals[f"I{p}"] for p in "abc"])[keep] / processing.HIL_I_BASE_PEAK_A,
            "power": power, "events": events, "timing": {
                "offset_s": offset, "offset_source": summary["offset_source"],
                "recorder_to_normalized_shift_s": timing["recorder_to_normalized_shift_s"]}}


def plot_pair(pair: dict, metas: dict, path: Path, case_id: int, presentation=None) -> None:
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    fig, axes = plt.subplots(4, 1, figsize=(12, 10), sharex=True)
    try:
        for side, style in (("hil", "-"), ("sim", "--")):
            data = pair[side]
            color = "tab:red" if side == "hil" else "tab:blue"
            stride = max(1, math.ceil(len(data["t_s"]) / 40000))
            indexes = np.unique(np.r_[np.arange(0, len(data["t_s"]), stride), len(data["t_s"]) - 1])
            for ax, key in zip(axes[:2], ("Vabc_pu", "Iabc_pu")):
                for i in range(3):
                    ax.plot(data["t_s"][indexes], data[key][indexes, i], color=color,
                            ls=style, lw=.45, alpha=.8, label=side.upper() if i == 0 else "_nolegend_")
            for ax, key in zip(axes[2:], ("P_pu", "Q_pu")):
                ax.plot(data["power"]["t_s"], data["power"][key], ls=style,
                        color=color, lw=.9, label=side.upper())
        for ax, label in zip(axes, ("三相电压 / pu", "三相电流 / pu", "有功功率 / pu", "无功功率 / pu")):
            ax.set_ylabel(label)
            ax.legend(loc="upper right", ncol=2, fontsize=8)
            ax.grid(True, alpha=.3)
            for event in pair["hil"]["events"].values():
                ax.axvline(event, color="gray", ls=":", lw=.8)
        for ax, channel in zip(axes[2:], ("P", "Q")):
            values = np.concatenate([pair[s]["power"][channel + "_pu"] for s in ("sim", "hil")])
            ax.set_ylim(pwr_power_axis_limits(values, metas["hil"], channel))
            ax.yaxis.set_major_locator(MultipleLocator(PWR_POWER_TICK_STEP_PU))
        axes[-1].set_xlim(0.0, max(pair[s]["t_s"][-1] for s in pair))
        axes[-1].set_xlabel("时间 / s")
        fig.suptitle(f"Case {case_id}: {metas['hil']['case_name']} — SIM / HIL")
        fig.tight_layout(rect=(0, 0, 1, .96))
        apply_presentation(fig, case_id, presentation)
        fig.savefig(path, dpi=150, facecolor="white")
    finally:
        plt.close(fig)


def comparison_rows(case_id: int, meta: dict, summaries: dict, figure: str) -> list[dict]:
    rows = []
    for metric, parameter, column, unit, deviation_unit in METRICS:
        row = dict(case_index=case_id, case_name=meta["case_name"], test_type=meta["test_type"],
                   event="event1", signal="Q_pu" if int(meta["test_type_code"]) == 4 else "P_pu",
                   metric=metric, parameter=parameter, unit=unit, deviation_unit=deviation_unit,
                   figure_png=figure)
        states = []
        for side in ("sim", "hil"):
            source = summaries[side]
            value = finite(source[column], column, optional=True)
            reason = source[metric + "_reason"]
            basis = source["step_timing_basis"]
            estimated = (bool(reason.strip()) or (metric != "M_p" and basis not in
                         ("sim_command", "synchronized_command")) or
                         (metric == "T_s" and source["T_s_method"] != "continuous_5pct"))
            status = "unavailable" if value is None else "estimated" if estimated else "valid"
            states.append(status)
            row.update({f"{side}_value": value, f"{side}_status": status, f"{side}_reason": reason,
                        f"{side}_manual_T_s": source["人工确认Ts"]})
            for key in ("step_timing_basis", "step_response_start_s", "T_s_method", "metric_sampling_ms"):
                row[f"{side}_{key}"] = source[key]
        difference = None if any(row[s + "_value"] is None for s in ("sim", "hil")) else row["sim_value"] - row["hil_value"]
        row.update(deviation=difference, absolute_deviation=None if difference is None else abs(difference),
                   status="unavailable" if "unavailable" in states else "estimated" if "estimated" in states else "valid")
        rows.append(row)
    return rows


def run(paths_config: Path | None, output: Path, case_id: int | None = None,
        *, case_ids=None, presentation=None) -> dict:
    paths = load_paths_config(paths_config)
    processing.configure_project(ProjectConfig.load(paths.project_config))
    inventories, summaries, references = {}, {}, {}
    for side in ("sim", "hil"):
        profile = paths.source(side)
        inventories[side] = recordings(profile.input_dir("pwr"))
        path = summary_path(profile.output_dir("pwr"), side)
        summaries[side] = read_summary(path, side)
        references[side] = relative(path, output)
    ids = select_case_ids(case_id, case_ids, set(inventories["sim"]) | set(inventories["hil"]))
    if not ids:
        raise ValueError("No PWR canonical recordings")
    output.mkdir(parents=True, exist_ok=False)
    figures = output / "figures"
    figures.mkdir()
    rows, completed, errors = [], [], []
    for ci in ids:
        try:
            files = {s: one(inventories[s], ci, s + " canonical") for s in inventories}
            sources = {s: one(summaries[s], ci, s + " summary") for s in summaries}
            metas = {s: validate_case(files[s], "pwr") for s in files}
            for side in files:
                if metas[side]["case_index"] != ci:
                    raise ValueError(f"{side}: filename and metadata index disagree")
                for key in ("case_name", "test_type"):
                    if sources[side][key] != metas[side][key]:
                        raise ValueError(f"{side}: summary / canonical mismatch: {key}")
            validate_pair(metas["sim"], metas["hil"])
            image_path = figure_path(figures, ci, "overlay", presentation)
            image_ref = relative(image_path, output)
            metrics = comparison_rows(ci, metas["hil"], sources, image_ref)
            pair = {s: full_signals(files[s], metas[s], sources[s]) for s in files}
            render_plot(plot_pair, pair, metas, image_path, ci, presentation=presentation)
            for row in metrics:
                row.update(presentation_record(metas["hil"], ci, presentation))
            rows.extend(metrics)
            completed.append({"case_index": ci, "figure_png": image_ref,
                              **presentation_record(metas["hil"], ci, presentation),
                              "canonical": {s: relative(files[s], output) for s in files},
                              "metadata": {s: relative(files[s].with_name(files[s].stem + "_meta.json"), output) for s in files},
                              "alignment": {s: pair[s]["timing"] for s in pair},
                              "recording_range_s": {s: [float(pair[s]["t_s"][0]), float(pair[s]["t_s"][-1])] for s in pair}})
            print(f"case {ci}: compared")
        except (ValueError, KeyError, OSError, TypeError) as exc:
            errors.append({"case_index": ci, "reason": str(exc)})
            print(f"case {ci}: error: {exc}")
    with (output / "compare_summary_pwr.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[*FIELDS, *PRESENTATION_FIELDS] if presentation is not None else FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    manifest = {"selected_case_ids": ids, "schema_version": 1, "mode": "pwr", "metric_event": "event1",
                "summaries": references, "project_config": relative(paths.project_config, output),
                "compare_csv": "compare_summary_pwr.csv", "completed": len(completed),
                "exceptions": len(errors), "cases": completed, "errors": errors}
    (output / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    return manifest


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--mode", default="pwr", choices=("pwr",))
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--case", type=int)
    selection.add_argument("--all", action="store_true")
    parser.add_argument("--paths-config", type=Path)
    parser.add_argument("--output", type=Path)


def execute(args, parser: argparse.ArgumentParser) -> int:
    output = args.output or Path.cwd() / "result_compare" / "pwr" / datetime.now().strftime("%Y-%m-%d_%H%M%S_%f")
    try:
        result = run(args.paths_config, output, args.case)
    except (ValueError, KeyError, OSError) as exc:
        parser.exit(2, f"error: {exc}\n")
    print(json.dumps({"output": str(output.resolve()), "completed": result["completed"], "exceptions": result["exceptions"]}))
    return int(bool(result["exceptions"]))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    return execute(parser.parse_args(argv), parser)


if __name__ == "__main__":
    raise SystemExit(main())
