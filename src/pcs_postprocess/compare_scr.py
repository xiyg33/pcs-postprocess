"""SCR comparisons on HIL-present cases and common observation horizons."""

import csv
import json
import re
from pathlib import Path

import numpy as np

from . import processing as p, scr_stability as stability
from .compare_common import read_summary, signals, validate_pair
from .compare_presentation import (PRESENTATION_FIELDS, apply_presentation, figure_path,
                                   presentation_record, render_plot, select_case_ids)
from .compare_pwr import finite, one, plt, recordings, relative
from .config import ProjectConfig
from .contract import validate_case
from .paths import load_paths_config


EVIDENCE_FIELDS = ["state", "reason", "rms_early_pu", "rms_late_pu", "growth_ratio",
                   "growth_noise_floor_pu", "baseline_std_pu", "trip_input_s",
                   "observation_start_after_recovery_s", "observation_end_after_recovery_s"]
FIELDS = ["case_index", "case_name", "P_ref_pu", "SCR_value", "comparison_version",
          "common_observation_duration_s", "state_match", "figure_png"]
for _side in ("sim", "hil"):
    FIELDS += [f"{_side}_{key}" for key in (*EVIDENCE_FIELDS, *stability.SUMMARY_FIELDS)]
BOUNDARY_FIELDS = ["P_ref_pu", "power_direction", "paired_case_count", "boundary_match"]
BOUNDARY_PRESENTATION_FIELDS = ["source_P_ref_pu", "reported_P_ref_pu", "approximation"]
for _side in ("sim", "hil"):
    BOUNDARY_FIELDS += [f"{_side}_{key}" for key in (
        "status", "lower_scr", "upper_scr", "minimum_tested_stable_scr",
        "maximum_tested_failed_scr", "transitions", "points")]


def write_csv(path, fields, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def boundary_rows(rows):
    groups = sorted(set(float(row["P_ref_pu"]) for row in rows))
    result = []
    for power in groups:
        selected = [row for row in rows if float(row["P_ref_pu"]) == power]
        row = {"P_ref_pu": power, "power_direction": "charge" if power < 0 else "discharge",
               "paired_case_count": len(selected)}
        states = {}
        for side in ("sim", "hil"):
            states[side] = stability.boundary([{"scr": float(x["SCR_value"]), "case_index": x["case_index"],
                                               "state": x[side + "_state"]} for x in selected])
            row.update({f"{side}_{k}": json.dumps(v, ensure_ascii=False) if isinstance(v, list) else v
                        for k, v in states[side].items()})
        sim, hil = states["sim"], states["hil"]
        row["boundary_match"] = ("same_tested_states" if all(x["sim_state"] == x["hil_state"] for x in selected)
                                 and all(x["sim_state"] != "UNDETERMINED" for x in selected) else
                                 "overlapping_brackets" if sim["status"] == hil["status"] == "bracketed"
                                 and max(sim["lower_scr"], hil["lower_scr"]) < min(sim["upper_scr"], hil["upper_scr"]) else
                                 "different_tested_states")
        if all("reported_P_ref_pu" in item for item in selected):
            targets = {float(item["reported_P_ref_pu"]) for item in selected}
            row.update(source_P_ref_pu=power,
                       reported_P_ref_pu=next(iter(targets)) if len(targets) == 1 else None,
                       approximation=any(not np.isclose(target, power, atol=1e-8, rtol=0.) for target in targets))
        result.append(row)
    return result


def plot_pair(pair, results, meta, path, horizon, presentation=None):
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
    try:
        for side, color, style in (("hil", "tab:red", "-"), ("sim", "tab:blue", "--")):
            data = pair[side]["recording"]
            stride = max(1, len(data["t_s"]) // 40000)
            for ax, key in zip(axes[:2], ("Vabc_pu", "Iabc_pu")):
                for phase in range(3):
                    ax.plot(data["t_s"][::stride], data[key][::stride, phase], color=color, ls=style,
                            lw=.45, alpha=.75, label=side.upper() if phase == 0 else "_nolegend_")
            axes[2].plot(data["display"]["t_s"], data["display"]["P_pu"], color=color, ls=style, lw=.85, label=side.upper())
            for ax in axes:
                ax.axvline(data["event_s"], color=color, ls="--", lw=.75, alpha=.7)
            recovery = data["recovery_s"]
            if recovery is not None:
                for ax in axes:
                    ax.axvline(recovery, color=color, ls=":", lw=.8)
                if horizon is not None:
                    axes[2].axvspan(recovery + .3, recovery + horizon, color=color, alpha=.045)
            trip = results[side]["trip_input_s"]
            if trip is not None:
                capture_trip = trip - data["input_origin_s"]
                if capture_trip >= 0:
                    for ax in axes:
                        ax.axvline(capture_trip, color=color, ls="-.", lw=.9)
                    axes[2].annotate(f"{side.upper()} 脱网 {capture_trip:.3f} s", (capture_trip, 0.),
                                     xytext=(6, 14), textcoords="offset points", color=color, fontsize=8)
        for ax, label in zip(axes, ("三相电压 / pu", "三相电流 / pu", "有功功率 / pu")):
            ax.set_ylabel(label)
            ax.grid(True, alpha=.3)
            ax.legend(loc="upper right", ncol=2, fontsize=8)
        end = min(float(data["recording"]["display_end_s"]) for data in pair.values())
        axes[-1].set_xlim(0., end)
        axes[-1].set_xlabel("各来源录波相对时间 / s；扰动和脱网分别标注")
        labels = {"STABLE": "正常", "UNSTABLE": "不正常（振荡）", "TRIP": "不正常（脱网）", "UNDETERMINED": "无法判定"}
        notes = " / ".join(f"{s.upper()}: {labels[r['state']]}" for s, r in results.items())
        power = float(meta["P_ref_pu"])
        power_tag = f"P{'m' if power < 0 else ''}{abs(power) * 100:03.0f}"
        fig.suptitle(f"Case {meta['case_index']}: {power_tag}, SCR={meta['SCR_value']} — HIL {labels[results['hil']['state']]}\n{notes}", fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, .94))
        apply_presentation(fig, meta["case_index"], presentation)
        fig.savefig(path, dpi=150, facecolor="white")
    finally:
        plt.close(fig)


def run(paths_config: Path | None, output: Path, case_id: int | None = None, overwrite=False,
        *, case_ids=None, presentation=None):
    paths = load_paths_config(paths_config)
    p.configure_project(ProjectConfig.load(paths.project_config))
    output.mkdir(parents=True, exist_ok=overwrite)
    (output / "figures").mkdir(exist_ok=overwrite)
    inventories, summaries, references, source_issues = {}, {}, {}, {}
    for side in ("sim", "hil"):
        source = paths.source(side)
        try:
            inventories[side] = recordings(source.input_dir("scr"))
        except OSError as exc:
            inventories[side] = {}
            source_issues[side] = str(exc)
        try:
            summary, summaries[side] = read_summary(source.output_dir("scr"), side, "scr",
                                                    ["case_index", "case_name", "test_type", *stability.SUMMARY_FIELDS])
            references[side] = relative(summary, output)
        except (ValueError, OSError) as exc:
            summaries[side] = {}
            source_issues[side] = str(exc)
    ids = select_case_ids(case_id, case_ids, paths.compare_cases.get("scr", sorted(inventories["hil"])))
    if not ids and case_ids is None:
        ids = sorted(summaries["hil"])
    rows, cases, errors = [], [], []
    if not ids:
        errors.append({"reason": source_issues.get("hil", "No HIL SCR recordings")})
    for ci in ids:
        try:
            if source_issues:
                raise ValueError(str(source_issues))
            files = {s: one(inventories[s], ci, s + " canonical") for s in inventories}
            sources = {s: one(summaries[s], ci, s + " summary") for s in files}
            metas = {s: validate_case(files[s], "scr") for s in files}
            for side in files:
                for key in ("case_index", "case_name", "test_type"):
                    if str(sources[side][key]) != str(metas[side][key]):
                        raise ValueError(f"{side}: summary / canonical mismatch: {key}")
                if sources[side]["scr_metrics_version"] != stability.VERSION:
                    raise ValueError(f"{side}: outdated SCR summary; rerun single-source analysis")
                for key in ("P_ref_pu", "SCR_value", "V_target_pu", "disturbance_start_s", "disturbance_duration_s", "sim_time_s", "record_pre_s", "record_post_s"):
                    if finite(sources[side][key], key) != finite(metas[side][key], key):
                        raise ValueError(f"{side}: summary / canonical parameter mismatch: {key}")
            validate_pair(metas, ("P_ref_pu", "SCR_value", "V_target_pu", "disturbance_start_s", "disturbance_duration_s", "sim_time_s", "record_pre_s", "record_post_s"))
            pair, available = {}, []
            for side in files:
                start = finite(sources[side]["scr_display_start_input_s"], "SCR display start")
                recovery = stability.optional_float(sources[side]["scr_recovery_input_s"])
                end = recovery if recovery is not None else start + float(metas[side]["disturbance_duration_s"])
                # Stability uses up to eight seconds of the available capture,
                # independent of the original display-only record_post request.
                pair[side] = signals(files[side], dict(metas[side], record_post_s=8.), start, end)
                pair[side]["input_t_s"] = pair[side]["t_s"] + start - 1.
                pair[side]["recovery_s"] = None if recovery is None else recovery - start + 1.
                pair[side]["timing"].update({key: sources[side][key] for key in stability.SUMMARY_FIELDS})
                raw_t, raw_sig = p.load_case_file(str(files[side]), absolute=metas[side]["time_basis"] == "absolute")
                origin = float(raw_t[0])
                capture_t = raw_t - origin
                capture_power = p.downsample_scalar_signals(
                    {"t_s": capture_t, "P_pu": raw_sig["P"] / p.HIL_P_BASE_W}, ("P_pu",))
                # SCR boundary evidence may occur well before the voltage
                # pulse; preserve the complete capture instead of clipping it.
                pair[side]["recording"] = {
                    "t_s": capture_t, "input_origin_s": origin,
                    "Vabc_pu": np.column_stack([raw_sig["V" + phase] for phase in "abc"]) / p.HIL_V_BASE_PEAK_V,
                    "Iabc_pu": np.column_stack([raw_sig["I" + phase] for phase in "abc"]) / p.HIL_I_BASE_PEAK_A,
                    "display": capture_power,
                    "event_s": finite(sources[side]["scr_configured_start_input_s"], "configured SCR start") - origin,
                    "display_end_s": finite(sources[side]["scr_display_end_input_s"], "SCR display end") - origin,
                    "recovery_s": None if recovery is None else recovery - origin}
                if recovery is not None:
                    available.append(min(float(sources[side]["scr_input_end_s"]),
                                         float(sources[side]["scr_display_end_input_s"])) - recovery - .01)
            horizon = min(8., *available) if len(available) == 2 else None
            results = {s: stability.classify(pair[s], sources[s], horizon) for s in files}
            image_path = figure_path(output / "figures", ci, "overlay", presentation)
            image_ref = relative(image_path, output)
            render_plot(plot_pair, pair, results, metas["hil"], image_path, horizon, presentation=presentation)
            row = {"case_index": ci, "case_name": metas["hil"]["case_name"],
                   "P_ref_pu": metas["hil"]["P_ref_pu"], "SCR_value": metas["hil"]["SCR_value"],
                   "comparison_version": stability.VERSION,
                   "common_observation_duration_s": None if horizon is None else max(0., horizon - .3),
                   "state_match": "unavailable" if any(r["state"] == "UNDETERMINED" for r in results.values()) else
                                  "same" if results["sim"]["state"] == results["hil"]["state"] else "different",
                   "figure_png": image_ref}
            for side in files:
                row.update({f"{side}_{key}": value for key, value in results[side].items()})
                row.update({f"{side}_{key}": sources[side][key] for key in stability.SUMMARY_FIELDS})
                if results[side]["state"] == "UNDETERMINED":
                    errors.append({"case_index": ci, "source": side, "reason": results[side]["reason"]})
            row.update(presentation_record(metas["hil"], ci, presentation))
            rows.append(row)
            cases.append({"case_index": ci, "figure": image_ref,
                          **presentation_record(metas["hil"], ci, presentation),
                          "canonical": {s: relative(path, output) for s, path in files.items()},
                          "metadata": {s: relative(path.with_name(path.stem + "_meta.json"), output) for s, path in files.items()},
                          "alignment": {s: pair[s]["timing"] for s in files},
                          "common_horizon_after_recovery_s": horizon})
            print(f"case {ci}: compared")
        except (ValueError, KeyError, TypeError, OSError) as exc:
            errors.append({"case_index": ci, "reason": str(exc)})
            print(f"case {ci}: error: {exc}")
    removed = []
    if overwrite and case_id is None and not errors:
        # Keep the active SCR gallery limited to the selected representative
        # cases. Only recognized generated case PNGs in this mode are removed.
        active = {case["figure"] for case in cases}
        figure_dir = (output / "figures").resolve()
        for path in sorted(figure_dir.glob("case_*.png")):
            match = re.match(r"case_(\d+)_", path.name)
            if (match and int(match.group(1)) in inventories["hil"]
                    and relative(path, output) not in active and path.resolve().parent == figure_dir):
                removed.append(relative(path, output))
                path.unlink()
    write_csv(output / "compare_summary_scr.csv",
              [*FIELDS, *PRESENTATION_FIELDS] if presentation is not None else FIELDS, rows)
    boundaries = boundary_rows(rows)
    write_csv(output / "compare_boundary_scr.csv",
              [*BOUNDARY_FIELDS, *BOUNDARY_PRESENTATION_FIELDS] if presentation is not None else BOUNDARY_FIELDS,
              boundaries)
    manifest = {"schema_version": 1, "mode": "scr", "metric_version": stability.VERSION,
                "selection": "configured representative cases" if "scr" in paths.compare_cases else "HIL-present cases only",
                "selected_case_ids": ids, "removed_generated_figures": removed,
                "output_policy": "overwrite_fixed_compare_directory" if overwrite else "explicit_output_directory",
                "summaries": references, "project_config": relative(paths.project_config, output),
                "compare_csv": "compare_summary_scr.csv", "boundary_csv": "compare_boundary_scr.csv",
                "completed": len(cases), "exceptions": len(errors), "cases": cases, "errors": errors,
                "stability_definition": {"smooth_s": .020, "sample_s": .005, "skip_after_recovery_s": .3,
                                         "maximum_horizon_after_recovery_s": 8., "minimum_observation_s": 1.,
                                         "growth_ratio": 1.3, "growth_floor_pu": .01, "baseline_noise_multiple": 3.,
                                         "power_limit_pu": 2., "power_limit_hold_s": .020,
                                         "breaker_hold_samples": 10}}
    manifest_name = "manifest.json"
    manifest["manifest_file"] = manifest_name
    (output / manifest_name).write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    return manifest
