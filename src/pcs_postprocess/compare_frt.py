"""FRT overlays and F1/F2/F3/FG deviations with independently aligned recovery."""

import csv
import json
import math
from pathlib import Path

import numpy as np
from matplotlib.ticker import MultipleLocator

from . import processing
from .compare_pwr import finite, one, plt, recordings, relative
from .config import ProjectConfig
from .contract import validate_case
from .paths import load_paths_config
from .compare_presentation import (PRESENTATION_FIELDS, apply_presentation, figure_path,
                                   presentation_record, render_plot, select_case_ids)
from .plot_common import pwr_power_axis_limits, PWR_POWER_TICK_STEP_PU


SIGNALS = {"Vpos_pu": "正序电压", "P_pu": "有功功率", "Q_pu": "无功功率",
           "Ip_pu": "有功电流", "Iq_pu": "无功电流"}
METRIC_COLUMNS = {
    "f1_a_pre_steady_mean": ("A", "mean_deviation"),
    "f1_b2_disturbed_steady_mean": ("B2", "mean_deviation"),
    "f1_c2_post_steady_mean": ("C2", "mean_deviation"),
    "f2_b1_rise_transient_mae": ("B1", "mae"),
    "f2_c1_recovery_transient_mae": ("C1", "mae"),
    "f3_a_pre_steady_max": ("A", "max_abs_error"),
    "f3_b2_disturbed_steady_max": ("B2", "max_abs_error"),
    "f3_c2_post_steady_max": ("C2", "max_abs_error"),
}
ALL_METRICS = [*METRIC_COLUMNS, "fg_weighted_overall_mae"]
FIELDS = ["case_index", "case_name", "test_type", "signal", "signal_label", "unit",
          *ALL_METRICS, "status", "reason", "sim_recovery_basis", "hil_recovery_basis",
          "configured_duration_difference_s", "measured_recovery_difference_s",
          "duration_effect", "figure_png", "power_current_figure_png"]
DT = processing.PWR_METRIC_SAMPLE_S


def summary_path(directory: Path, side: str) -> Path:
    path = directory / f"{side}_summary_frt.csv"
    if not path.is_file():
        path = directory / "summary_frt.csv"
    if not path.is_file():
        raise FileNotFoundError(f"Missing {side} FRT summary in {directory}")
    return path


def read_summary(path: Path) -> dict:
    result = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"case_index", "case_name", "test_type", "rec_offset_s", "offset_source"}
        if required - set(reader.fieldnames or []):
            raise ValueError(f"{path}: missing columns {sorted(required - set(reader.fieldnames or []))}")
        for row in reader:
            result.setdefault(int(row["case_index"]), []).append(row)
    return result


def load_metadata(path: Path) -> dict:
    return json.loads(path.with_name(path.stem + "_meta.json").read_text(encoding="utf-8-sig"))


def is_scan(meta: dict) -> bool:
    return str(meta.get("test_type", "")).strip().upper() in ("LVRT_BND", "HVRT_BND")


def validate_pair(sim: dict, hil: dict) -> None:
    for key in ("case_index", "case_name", "test_type"):
        if sim.get(key) != hil.get(key):
            raise ValueError(f"SIM/HIL parameter mismatch: {key}")
    if sim["test_type"] not in ("LVRT_SINGLE", "HVRT_SINGLE"):
        raise ValueError(f"Unsupported ordinary FRT type: {sim['test_type']}")
    for key in ("P_ref_pu", "V_fault_pu", "fault_start_s"):
        a, b = finite(sim.get(key), key), finite(hil.get(key), key)
        if not np.isclose(a, b, atol=1e-8, rtol=0):
            raise ValueError(f"SIM/HIL parameter mismatch: {key} ({a} vs {b})")
    # 旧 SIM 元数据未保存 fault_type，但 case_name 保留 3P/AB；两边都有字段时再校验。
    if sim.get("fault_type") is not None and hil.get("fault_type") is not None:
        if sim["fault_type"] != hil["fault_type"]:
            raise ValueError("SIM/HIL parameter mismatch: fault_type")
    for side, meta in (("sim", sim), ("hil", hil)):
        if finite(meta.get("fault_duration_s"), side + " fault_duration_s") <= 0:
            raise ValueError(f"{side}: nonpositive fault duration")


def prepare_signals(path: Path, meta: dict, summary: dict) -> dict:
    """Reuse source offsets and FRT processing without recalculating summaries."""
    offset = finite(summary["rec_offset_s"], "rec_offset_s")
    ev = processing.event_times_from_meta(meta)
    raw_t, sig = processing.load_case_file(str(path), absolute=meta["time_basis"] == "absolute")
    origin = ev["event1_start"] - 1.0
    recovery = processing.detect_frt_recovery_time(raw_t, sig["Vp_pu"], meta, ev, offset)
    basis = "detected" if recovery is not None else "metadata_estimate"
    end = ev["event1_end"] if recovery is None else recovery
    ev = dict(ev, event1_end=end, last_event_end=end)
    time, events, _ = processing.normalize_case_time(raw_t + offset, ev)
    keep = time >= 0.0
    if np.count_nonzero(keep) < 3:
        raise ValueError("Insufficient recording after normalized 0 s")
    vabc = np.column_stack([sig[f"V{p}"] for p in "abc"])
    iabc = np.column_stack([sig[f"I{p}"] for p in "abc"])
    ip, iq, valid = processing.calculate_instantaneous_current(vabc, iabc)
    # scalar_plot_data 的 Ip/Iq 同样基于瞬时方向电流，再统一平滑和降采样。
    bridged = {"t_s": time[keep], "Vp_pu": sig["Vp_pu"][keep],
               "P_pu": sig["P"][keep] / processing.HIL_P_BASE_W,
               "Q_pu": (sig["Q"][keep] / processing.HIL_P_BASE_W if sig["Q"] is not None
                        else np.full(np.count_nonzero(keep), np.nan)),
               "Ip_pu": ip[keep] / processing.HIL_I_BASE_PEAK_A,
               "Iq_pu": iq[keep] / processing.HIL_I_BASE_PEAK_A}
    display = processing.downsample_scalar_signals(
        bridged, ("Vp_pu", "P_pu", "Q_pu", "Ip_pu", "Iq_pu"), 1.0, end - origin)
    display["Vpos_pu"] = display["Vp_pu"]
    return {"t_s": time[keep], "Vabc_pu": vabc[keep] / processing.HIL_V_BASE_PEAK_V,
            "Iabc_pu": iabc[keep] / processing.HIL_I_BASE_PEAK_A,
            "current_valid": valid[keep], "display": display,
            "recovery_s": end - origin, "configured_recovery_s": 1.0 + float(meta["fault_duration_s"]),
            "recovery_basis": basis, "events": events,
            "timing": {"offset_s": offset, "offset_source": summary["offset_source"],
                       "recorder_to_normalized_shift_s": offset - origin,
                       "recovery_basis": basis, "recovery_s": end - origin}}


def windows(pair: dict, recovery_aligned: bool = True, longer_windows: bool = False) -> list[dict]:
    ends = {side: float(data["recovery_s"]) for side, data in pair.items()}
    durations = [value - 1.0 for value in ends.values()]
    if min(durations) <= 0:
        raise ValueError("Fault recovery must be after fault onset")
    duration = max(durations) if longer_windows else min(durations)
    transient = min(duration / 4, .2)
    guard = min(duration / 10, .1)
    recovery_zero = ends if recovery_aligned else {side: min(ends.values()) for side in ends}
    return [
        {"code": "A", "start_s": 0., "end_s": .9, "origins_s": {s: 0. for s in ends}},
        {"code": "B1", "start_s": 1., "end_s": 1. + transient, "origins_s": {s: 0. for s in ends}},
        {"code": "B2", "start_s": 1. + transient, "end_s": min(ends.values()) - guard,
         "origins_s": {s: 0. for s in ends}},
        {"code": "C1", "start_s": 0., "end_s": transient, "origins_s": recovery_zero},
        {"code": "C2", "start_s": transient, "end_s": 2., "origins_s": recovery_zero},
    ]


def segment_metrics(pair: dict, signal: str, window: dict) -> dict:
    start, stop = window["start_s"], window["end_s"]
    for side, data in pair.items():
        origin = window["origins_s"][side]
        start = max(start, float(data["display"]["t_s"][0]) - origin,
                    float(data["t_s"][0]) - origin)
        stop = min(stop, float(data["display"]["t_s"][-1]) - origin)
    result = {"window": window, "common_start_s": start, "common_end_s": stop,
              "sample_count": 0, "valid_sample_count": 0, "reason": ""}
    if stop <= start:
        return dict(result, reason="empty_common_window")
    # 半开区间：同一采样不能同时计入相邻暂态和稳态窗口。
    grid = window["start_s"] + np.arange(max(0, int(math.ceil((stop - window["start_s"]) / DT)))) * DT
    grid = grid[(grid >= start - 1e-12) & (grid < stop - 1e-12)]
    if not grid.size:
        return dict(result, reason="empty_common_window")
    values = {}
    for side, data in pair.items():
        display = data["display"]
        values[side] = np.interp(grid + window["origins_s"][side], display["t_s"],
                                 display[signal], left=np.nan, right=np.nan)
    valid = np.isfinite(values["sim"]) & np.isfinite(values["hil"])
    result.update(sample_count=int(grid.size), valid_sample_count=int(valid.sum()))
    if not valid.all():
        return dict(result, reason="missing_samples")
    error = values["sim"] - values["hil"]
    result.update(mean_deviation=float(abs(np.mean(error))), mae=float(np.mean(abs(error))),
                  max_abs_error=float(np.max(abs(error))))
    return result


def deviations(pair: dict, recovery_aligned: bool = True, longer_windows: bool = False) -> tuple[list, dict]:
    definitions = windows(pair, recovery_aligned, longer_windows)
    estimated = any(data["recovery_basis"] != "detected" or
                    data["timing"]["offset_source"] not in ("edge", "absolute_time") for data in pair.values())
    return deviations_for_windows(pair, definitions, estimated=estimated)


def deviations_for_windows(pair: dict, definitions: list[dict], *, estimated: bool = False) -> tuple[list, dict]:
    """Summarize explicit A/B1/B2/C1/C2 windows on per-side time origins.

    Continuous events supply their own adjacent-stage windows while retaining
    the single-fault formulas, missing-data policy and per-unit signals.
    """
    rows, details = [], {}
    for signal, label in SIGNALS.items():
        segments = {window["code"]: segment_metrics(pair, signal, window) for window in definitions}
        row = {"signal": signal, "signal_label": label, "unit": "pu"}
        for column, (code, metric) in METRIC_COLUMNS.items():
            row[column] = segments[code].get(metric)
        failures = [f"{code}: {data['reason']}" for code, data in segments.items() if data["reason"]]
        count = sum(data["sample_count"] for data in segments.values())
        row["fg_weighted_overall_mae"] = (sum(data["mae"] * data["sample_count"] for data in segments.values()) / count
                                         if count and not failures else None)
        row.update(status="unavailable" if failures else "estimated" if estimated else "valid",
                   reason="; ".join(failures + (["timing_estimated"] if estimated else [])))
        rows.append(row)
        details[signal] = segments
    return rows, details


def metric_changes(baseline: list, alternate: list) -> dict:
    result = {}
    for before, after in zip(baseline, alternate):
        result[before["signal"]] = {}
        for column in ALL_METRICS:
            a, b = before[column], after[column]
            result[before["signal"]][column] = {
                "baseline": a, "alternate": b, "delta": None if a is None or b is None else b - a,
                "changed": (a is None) != (b is None) or (a is not None and b is not None and abs(a - b) > 1e-12)}
    return result


def duration_assessment(pair: dict, meta: dict, rows: list, details: dict) -> dict:
    configured = float(meta["sim"]["fault_duration_s"]) - float(meta["hil"]["fault_duration_s"])
    measured = pair["sim"]["recovery_s"] - pair["hil"]["recovery_s"]
    result = {"configured_duration_difference_s": configured, "measured_recovery_difference_s": measured,
              "status": "same_configuration", "effect": "not_applicable"}
    if abs(configured) <= 1e-8:
        return result
    onset, onset_details = deviations(pair, recovery_aligned=False)
    longer, longer_details = deviations(pair, longer_windows=True)
    onset_changes, window_changes = metric_changes(rows, onset), metric_changes(rows, longer)
    sample_changes = {s: {code: {"baseline": details[s][code]["sample_count"],
                               "alternate": longer_details[s][code]["sample_count"]}
                         for code in details[s] if details[s][code]["sample_count"] != longer_details[s][code]["sample_count"]}
                      for s in SIGNALS}
    affected = any(v["changed"] for changes in (onset_changes, window_changes) for signal in changes.values() for v in signal.values()) or any(sample_changes.values())
    result.update(status="evaluated", effect="affected" if affected else "no_observed_effect",
                  explanation="Comparison checks numerical changes, not acceptance limits; recovery is aligned separately in the reported metrics.",
                  onset_only_metric_changes=onset_changes, longer_duration_window_metric_changes=window_changes,
                  window_sample_count_changes=sample_changes, onset_only_windows=onset_details,
                  longer_duration_windows=longer_details)
    return result


def plot_pair(pair: dict, paths: tuple[Path, Path], case_id: int, case_name: str, presentation=None) -> None:
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    layouts = (("Vpos_pu", "Vabc_pu", "Iabc_pu"), ("Vpos_pu", "P_pu", "Ip_pu", "Q_pu", "Iq_pu"))
    labels = {**SIGNALS, "Vabc_pu": "三相电压", "Iabc_pu": "三相电流", "Ip_pu": "有功电流", "Iq_pu": "无功电流"}
    for panels, path in zip(layouts, paths):
        fig, axes = plt.subplots(len(panels), 1, figsize=(12, 9), sharex=True)
        try:
            for ax, signal in zip(axes, panels):
                for side, color, style in (("hil", "tab:red", "-"), ("sim", "tab:blue", "--")):
                    data = pair[side]
                    if signal in ("Vabc_pu", "Iabc_pu"):
                        stride = max(1, math.ceil(len(data["t_s"]) / 40000))
                        for phase in range(3):
                            ax.plot(data["t_s"][::stride], data[signal][::stride, phase], color=color,
                                    ls=style, lw=.45, alpha=.8, label=side.upper() if phase == 0 else "_nolegend_")
                    else:
                        ax.plot(data["display"]["t_s"], data["display"][signal], color=color,
                                ls=style, lw=.85, label=side.upper())
                    ax.axvline(data["recovery_s"], color=color, ls=":", lw=.7, alpha=.7)
                    if signal in ("Ip_pu", "Iq_pu"):
                        ax.fill_between(data["t_s"], 0, 1, where=~data["current_valid"],
                                        transform=ax.get_xaxis_transform(), color="gray", alpha=.08, step="mid")
                ax.axvline(1., color="gray", ls=":", lw=.8)
                ax.set_ylabel(labels[signal] + " / pu")
                ax.legend(loc="upper right", ncol=2, fontsize=8)
                ax.grid(True, alpha=.3)
                if signal == "Vpos_pu":
                    values = np.concatenate([pair[s]["display"][signal] for s in pair])
                    values = values[np.isfinite(values)]
                    ax.set_ylim(min(-.1, float(values.min()) - .05) if values.size else -.1,
                                max(1.3, float(values.max()) + .05) if values.size else 1.3)
                elif signal in ("P_pu", "Q_pu", "Ip_pu", "Iq_pu"):
                    values = np.concatenate([pair[s]["display"][signal] for s in pair])
                    ax.set_ylim(pwr_power_axis_limits(values, {}, "Q" if signal.startswith("Q") or signal == "Iq_pu" else "P"))
                    ax.yaxis.set_major_locator(MultipleLocator(PWR_POWER_TICK_STEP_PU))
            axes[-1].set_xlim(0., max(data["t_s"][-1] for data in pair.values()))
            axes[-1].set_xlabel("时间 / s")
            fig.suptitle(f"Case {case_id}: {case_name} — SIM / HIL")
            fig.tight_layout(rect=(0, 0, 1, .96))
            apply_presentation(fig, case_id, presentation)
            fig.savefig(path, dpi=150, facecolor="white")
        finally:
            plt.close(fig)


def run(paths_config: Path | None, output: Path, case_id: int | None = None,
        *, case_ids=None, presentation=None) -> dict:
    config = load_paths_config(paths_config)
    processing.configure_project(ProjectConfig.load(config.project_config))
    inventories, summaries, references = {}, {}, {}
    for side in ("sim", "hil"):
        profile = config.source(side)
        inventories[side] = recordings(profile.input_dir("frt"))
        path = summary_path(profile.output_dir("frt"), side)
        summaries[side] = read_summary(path)
        references[side] = relative(path, output)
    ids = select_case_ids(case_id, case_ids, set(inventories["sim"]) | set(inventories["hil"]))
    if not ids:
        raise ValueError("No FRT canonical recordings")
    output.mkdir(parents=True, exist_ok=False)
    figures = output / "figures"
    figures.mkdir()
    rows, completed, skipped, errors = [], [], [], []
    for ci in ids:
        try:
            files = {side: one(data, ci, side + " canonical") for side, data in inventories.items() if ci in data}
            meta = {side: load_metadata(path) for side, path in files.items()}
            if any(is_scan(item) for item in meta.values()):
                skipped.append({"case_index": ci, "reason": "boundary_scan_excluded"})
                continue
            if "sim" not in files and "hil" in files:
                skipped.append({"case_index": ci, "reason": "missing_sim"})
                continue
            if set(files) != {"sim", "hil"}:
                raise ValueError(f"case {ci}: missing canonical source")
            sources = {side: one(summaries[side], ci, side + " summary") for side in files}
            meta = {side: validate_case(path, "frt") for side, path in files.items()}
            for side in files:
                if meta[side]["case_index"] != ci:
                    raise ValueError(f"{side}: filename / metadata index mismatch")
                for key in ("case_name", "test_type"):
                    if sources[side][key] != meta[side][key]:
                        raise ValueError(f"{side}: summary / canonical mismatch: {key}")
            validate_pair(meta["sim"], meta["hil"])
            pair = {side: prepare_signals(files[side], meta[side], sources[side]) for side in files}
            metrics, details = deviations(pair)
            assessment = duration_assessment(pair, meta, metrics, details)
            image_paths = (figure_path(figures, ci, "overlay", presentation),
                           figure_path(figures, ci, "overlay_power_current", presentation))
            images = [relative(path, output) for path in image_paths]
            render_plot(plot_pair, pair, image_paths, ci, meta["hil"]["case_name"], presentation=presentation)
            for row in metrics:
                row.update(presentation_record(meta["hil"], ci, presentation))
                row.update(case_index=ci, case_name=meta["hil"]["case_name"], test_type=meta["hil"]["test_type"],
                           sim_recovery_basis=pair["sim"]["recovery_basis"], hil_recovery_basis=pair["hil"]["recovery_basis"],
                           configured_duration_difference_s=assessment["configured_duration_difference_s"],
                           measured_recovery_difference_s=assessment["measured_recovery_difference_s"],
                           duration_effect=assessment["effect"], figure_png=images[0], power_current_figure_png=images[1])
                if assessment["effect"] == "affected":
                    row["reason"] = "; ".join(filter(None, (row["reason"], "duration_difference_has_observed_effect; recovery_aligned_separately")))
            rows.extend(metrics)
            completed.append({"case_index": ci, "figures": images,
                              **presentation_record(meta["hil"], ci, presentation),
                              "canonical": {s: relative(p, output) for s, p in files.items()},
                              "metadata": {s: relative(p.with_name(p.stem + "_meta.json"), output) for s, p in files.items()},
                              "alignment": {s: data["timing"] for s, data in pair.items()},
                              "windows": details, "duration_assessment": assessment})
            print(f"case {ci}: compared")
        except (ValueError, KeyError, OSError, TypeError) as exc:
            errors.append({"case_index": ci, "reason": str(exc)})
            print(f"case {ci}: error: {exc}")
    with (output / "compare_summary_frt.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[*FIELDS, *PRESENTATION_FIELDS] if presentation is not None else FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    manifest = {"selected_case_ids": ids, "schema_version": 1, "mode": "frt", "metric_sampling_s": DT,
                "summary_definition": "F1/F2/F3/FG; onset aligned, recovery independently aligned",
                "summaries": references, "project_config": relative(config.project_config, output),
                "compare_csv": "compare_summary_frt.csv", "completed": len(completed),
                "skipped": len(skipped), "exceptions": len(errors), "cases": completed, "skips": skipped, "errors": errors}
    (output / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    return manifest
