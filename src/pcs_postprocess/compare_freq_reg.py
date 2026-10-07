"""Referenced first-P-step metrics and frequency-aligned regulation overlays."""

import csv
import json
import math
from pathlib import Path

import numpy as np
from matplotlib.ticker import MultipleLocator

from . import processing
from .freq_reg_timing import SUMMARY_FIELDS as TIMING_FIELDS
from .compare_pwr import (FIELDS as PWR_FIELDS, METRICS, comparison_rows as step_rows,
                          finite, one, plt, recordings, relative)
from .config import ProjectConfig
from .contract import validate_case
from .paths import load_paths_config
from .plot_common import (configure_nominal_frequency, freq_reg_display_limits,
                          pwr_power_axis_limits, FREQ_REG_DISPLAY_TICK_HZ,
                          PWR_POWER_TICK_STEP_PU)


FIELDS = [*PWR_FIELDS, "power_current_figure_png", "P_ref_pu", "freq_deviation_hz",
          "limit_pu", "sim_5pct_pass", "hil_5pct_pass",
          "sim_offset_source", "hil_offset_source"]
for _side in ("sim", "hil"):
    FIELDS += [f"{_side}_{key}" for key in TIMING_FIELDS]


def summary_path(directory: Path, side: str) -> Path:
    path = directory / f"{side}_summary_freq_reg.csv"
    if not path.is_file():
        path = directory / "summary_freq_reg.csv"
    if not path.is_file():
        raise FileNotFoundError(f"Missing {side} freq_reg summary in {directory}")
    return path


def read_summary(path: Path, side: str) -> dict:
    rows = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        columns = set(reader.fieldnames or [])
        required = {"case_index", "case_name", "test_type", "rec_offset_s", "offset_source",
                    "P_change_abs_pu", "P_change_5pct_pass"}
        if required - columns:
            raise ValueError(f"{path}: missing columns {sorted(required - columns)}; rerun single-source analysis")
        aliases = {column: column if column in columns else column.replace(metric + "_", metric + f"_{side}_", 1)
                   for metric, _, column, _, _ in METRICS}
        dynamic = (set(processing.PWR_STEP_SUMMARY_COLS) - set(aliases)) | set(TIMING_FIELDS)
        for row in reader:
            if not processing.is_freq_reg_deadband(row):
                missing = (dynamic - columns) | {column for column, alias in aliases.items() if alias not in columns}
                if missing:
                    raise ValueError(f"{path}: missing first-step columns {sorted(missing)}; "
                                     "rerun pcs-postprocess run --mode freq_reg --source " + side)
                for column, alias in aliases.items():
                    row[column] = row[alias]
            rows.setdefault(int(row["case_index"]), []).append(row)
    return rows


def validate_pair(sim: dict, hil: dict) -> None:
    for key in ("case_index", "case_name", "test_type", "test_type_code"):
        if sim.get(key) != hil.get(key):
            raise ValueError(f"SIM/HIL parameter mismatch: {key}")
    for key in ("P_ref_pu", "freq_deviation_hz", "event_start_s", "event_hold_s"):
        if not np.isclose(finite(sim.get(key), key), finite(hil.get(key), key), atol=1e-8, rtol=0):
            raise ValueError(f"SIM/HIL parameter mismatch: {key}")


def prepare_signals(path: Path, meta: dict, summary: dict) -> dict:
    """Read saved frequency alignment and P timing; never redetect or measure."""
    offset = finite(summary["rec_offset_s"], "rec_offset_s")
    ev = processing.event_times_from_meta(meta)
    selected_start = finite(summary["freq_event_start_recorder_s"], "freq_event_start_recorder_s")
    selected_end = finite(summary["freq_event_end_recorder_s"], "freq_event_end_recorder_s")
    if not np.isclose(ev["event1_start"] - offset, selected_start, atol=1e-8, rtol=0):
        raise ValueError("Frequency alignment summary is inconsistent; rerun single-source analysis")
    ev["event1_end"] = selected_end + offset
    ev["last_event_end"] = ev["event1_end"]
    raw_t, sig = processing.load_case_file(str(path), absolute=meta["time_basis"] == "absolute")
    time, events, keep = processing.normalize_case_time(raw_t + offset, ev)
    vabc = np.column_stack([sig[f"V{p}"] for p in "abc"])
    iabc = np.column_stack([sig[f"I{p}"] for p in "abc"])
    ip, iq, valid = processing.calculate_instantaneous_current(vabc, iabc)
    data = {"t_s": time[keep], "P_pu": sig["P"][keep] / processing.HIL_P_BASE_W,
            "Q_pu": (sig["Q"][keep] / processing.HIL_P_BASE_W if sig["Q"] is not None
                     else np.full(int(keep.sum()), np.nan)), "f_hz": sig["f"][keep],
            "Ip_pu": ip[keep] / processing.HIL_I_BASE_PEAK_A,
            "Iq_pu": iq[keep] / processing.HIL_I_BASE_PEAK_A}
    display = processing.downsample_scalar_signals(data, ("f_hz", "Ip_pu", "Iq_pu"),
                                                  events["event1_start"], events["event1_end"])
    start = finite(summary.get("step_response_start_s"), "step_response_start_s", optional=True)
    metric_end = finite(summary["freq_metric_end_recorder_s"], "freq_metric_end_recorder_s")
    end = min(metric_end - selected_start + 1., float(data["t_s"][-1]))
    power = processing.pwr_downsampled_signals(data, start, end)
    display.update(P_pu=power["P_pu"], Q_pu=power["Q_pu"])
    return {"t_s": data["t_s"], "Vabc_pu": vabc[keep] / processing.HIL_V_BASE_PEAK_V,
            "Iabc_pu": iabc[keep] / processing.HIL_I_BASE_PEAK_A, "display": display,
            "current_valid": valid[keep], "events": events,
            "timing": {"offset_s": offset, "offset_source": summary["offset_source"],
                       "recorder_to_normalized_shift_s": offset - ev["event1_start"] + 1.,
                       "step_response_start_s": start,
                       "step_timing_basis": summary["step_timing_basis"],
                       **{key: summary[key] for key in TIMING_FIELDS}}}


def plot_pair(pair: dict, metas: dict, paths: tuple[Path, Path], case_id: int) -> None:
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    layouts = (("f_hz", "Vabc_pu", "Iabc_pu"), ("f_hz", "P_pu", "Ip_pu", "Q_pu", "Iq_pu"))
    labels = {"f_hz": "频率 / Hz", "Vabc_pu": "三相电压 / pu", "Iabc_pu": "三相电流 / pu",
              "P_pu": "有功功率 / pu", "Q_pu": "无功功率 / pu",
              "Ip_pu": "有功电流 / pu", "Iq_pu": "无功电流 / pu"}
    for keys, path in zip(layouts, paths):
        fig, axes = plt.subplots(len(keys), 1, figsize=(12, 9), sharex=True)
        try:
            for side, color, style in (("hil", "tab:red", "-"), ("sim", "tab:blue", "--")):
                data = pair[side]
                for ax, key in zip(axes, keys):
                    if key in ("Vabc_pu", "Iabc_pu"):
                        stride = max(1, math.ceil(len(data["t_s"]) / 40000))
                        indexes = np.unique(np.r_[np.arange(0, len(data["t_s"]), stride), len(data["t_s"]) - 1])
                        for phase in range(3):
                            ax.plot(data["t_s"][indexes], data[key][indexes, phase], color=color,
                                    ls=style, lw=.45, alpha=.8, label=side.upper() if phase == 0 else "_nolegend_")
                    else:
                        ax.plot(data["display"]["t_s"], data["display"][key],
                                color=color, ls=style, lw=.9, label=side.upper())
            for ax, key in zip(axes, keys):
                ax.set_ylabel(labels[key])
                ax.legend(loc="upper right", ncol=2, fontsize=8)
                ax.grid(True, alpha=.3)
                for event in (pair["hil"]["events"]["event1_start"], pair["hil"]["events"]["event1_end"]):
                    ax.axvline(event, color="gray", ls=":", lw=.8)
                if key == "f_hz":
                    values = np.concatenate([pair[s]["display"][key] for s in pair])
                    ax.set_ylim(freq_reg_display_limits(values))
                    ax.yaxis.set_major_locator(MultipleLocator(FREQ_REG_DISPLAY_TICK_HZ))
                    ax.axhline(processing.NOMINAL_FREQUENCY_HZ, color="gray", ls="--", lw=.6)
                elif key in ("P_pu", "Ip_pu", "Q_pu", "Iq_pu"):
                    values = np.concatenate([pair[s]["display"][key] for s in pair])
                    channel = "P" if key in ("P_pu", "Ip_pu") else "Q"
                    limits = pwr_power_axis_limits(values, {}, channel)
                    if channel == "Q":
                        limits = (min(limits[0], -.4), max(limits[1], .4))
                    ax.set_ylim(limits)
                    ax.yaxis.set_major_locator(MultipleLocator(PWR_POWER_TICK_STEP_PU))
                    ax.axhline(0., color="gray", ls="--", lw=.6)
            axes[-1].set_xlim(0., max(pair[s]["t_s"][-1] for s in pair))
            axes[-1].set_xlabel("时间 / s")
            pending = any(str(pair[s]["timing"].get("freq_timing_review_required", "0")) == "1" for s in pair)
            title_note = "\n频率对时待人工确认：当前使用单例名义起点" if pending else ""
            fig.suptitle(f"Case {case_id}: {metas['hil']['case_name']} — SIM / HIL{title_note}", fontsize=10)
            fig.tight_layout(rect=(0, 0, 1, .96))
            fig.savefig(path, dpi=150, facecolor="white")
        finally:
            plt.close(fig)


def comparison_rows(case_id: int, meta: dict, summaries: dict, figures: tuple[str, str]) -> list[dict]:
    if not processing.is_freq_reg_deadband(meta):
        rows = step_rows(case_id, dict(meta, test_type_code=3), summaries, figures[0])
        for row in rows:
            row["power_current_figure_png"] = figures[1]
            for side in ("sim", "hil"):
                for key in TIMING_FIELDS:
                    row[f"{side}_{key}"] = summaries[side].get(key, "")
                if str(summaries[side].get("freq_timing_review_required")) == "1":
                    row[f"{side}_reason"] = "; ".join(filter(None, (
                        row[f"{side}_reason"], "频率扰动对时待人工确认", summaries[side].get("freq_timing_reason", ""))))
                    if row[f"{side}_status"] == "valid":
                        row[f"{side}_status"] = "estimated"
                    if row["status"] == "valid":
                        row["status"] = "estimated"
    else:
        row = dict(case_index=case_id, case_name=meta["case_name"], test_type="FREQ_REG_DEADBAND",
                   event="event1", signal="P_pu", metric="P_change_abs", parameter="死区绝对功率变化",
                   unit="pu", deviation_unit="pu", limit_pu=processing.FREQ_REG_DEADBAND_P_LIMIT_PU,
                   figure_png="", power_current_figure_png="")
        states = []
        for side in ("sim", "hil"):
            source = summaries[side]
            value = finite(source.get("P_change_abs_pu"), "P_change_abs_pu", optional=True)
            flag = finite(source.get("P_change_5pct_pass"), "P_change_5pct_pass", optional=True)
            reason = source.get("reason", "")
            if value is not None and flag not in (0., 1.):
                raise ValueError(f"{side}: missing or invalid P_change_5pct_pass; rerun single-source analysis")
            if value is not None and flag != int(value <= processing.FREQ_REG_DEADBAND_P_LIMIT_PU + 1e-12):
                # Summary values are rounded to four decimal places; allow a
                # boundary-straddling value without changing the saved verdict.
                if abs(value - processing.FREQ_REG_DEADBAND_P_LIMIT_PU) > .00005:
                    raise ValueError(f"{side}: inconsistent P_change_5pct_pass")
            status = ("unavailable" if value is None else "estimated"
                      if source.get("offset_source") not in ("edge", "absolute_time") else "valid")
            states.append(status)
            row.update({f"{side}_value": value, f"{side}_5pct_pass": flag,
                        f"{side}_status": status, f"{side}_reason": reason,
                        f"{side}_offset_source": source.get("offset_source", "")})
        difference = None if any(row[s + "_value"] is None for s in ("sim", "hil")) else row["sim_value"] - row["hil_value"]
        row.update(deviation=difference, absolute_deviation=None if difference is None else abs(difference),
                   status="unavailable" if "unavailable" in states else "estimated" if "estimated" in states else "valid")
        rows = [row]
    for row in rows:
        row.update(P_ref_pu=meta["P_ref_pu"], freq_deviation_hz=meta["freq_deviation_hz"])
        for side in summaries:
            row[f"{side}_offset_source"] = summaries[side].get("offset_source", "")
    return rows


def missing_summary(reason: str) -> dict:
    row = {key: "" for key in processing.PWR_STEP_SUMMARY_COLS}
    row.update(offset_source="missing_recording", reason=reason)
    for metric, *_ in METRICS:
        row[metric + "_reason"] = reason
    return row


def run(paths_config: Path | None, output: Path, case_id: int | None = None) -> dict:
    paths = load_paths_config(paths_config)
    config = ProjectConfig.load(paths.project_config)
    processing.configure_project(config)
    configure_nominal_frequency(config.nominal_frequency_hz)
    inventories, summaries, references = {}, {}, {}
    for side in ("sim", "hil"):
        source = paths.source(side)
        inventories[side] = recordings(source.input_dir("freq_reg"))
        summary = summary_path(source.output_dir("freq_reg"), side)
        summaries[side] = read_summary(summary, side)
        references[side] = relative(summary, output)
    ids = [case_id] if case_id is not None else sorted(set(inventories["sim"]) | set(inventories["hil"]) |
                                                      set(summaries["sim"]) | set(summaries["hil"]))
    if not ids:
        raise ValueError("No freq_reg cases")
    output.mkdir(parents=True, exist_ok=False)
    (output / "figures").mkdir()
    rows, cases, errors = [], [], []
    for ci in ids:
        try:
            files, metas, source_rows, missing = {}, {}, {}, []
            for side in ("sim", "hil"):
                if not inventories[side].get(ci):
                    reason = f"{side}: missing canonical recording for case {ci}"
                    missing.append(reason)
                    source_rows[side] = missing_summary(reason)
                    continue
                files[side] = one(inventories[side], ci, side + " canonical")
                metas[side] = validate_case(files[side], "freq_reg")
                source_rows[side] = one(summaries[side], ci, side + " summary")
                for key in ("case_name", "test_type"):
                    if source_rows[side][key] != metas[side][key]:
                        raise ValueError(f"{side}: summary / canonical mismatch: {key}")
            if not metas:
                raise ValueError("Missing both canonical recordings")
            meta = metas.get("hil", metas.get("sim"))
            if len(metas) == 2:
                validate_pair(metas["sim"], metas["hil"])
            if missing and not processing.is_freq_reg_deadband(meta):
                raise ValueError("; ".join(missing))
            images = ("", "")
            pair = {}
            if not processing.is_freq_reg_deadband(meta):
                image_paths = (output / "figures" / f"case_{ci:04d}_sim_hil_overlay.png",
                               output / "figures" / f"case_{ci:04d}_sim_hil_power_current.png")
                images = tuple(relative(p, output) for p in image_paths)
                pair = {side: prepare_signals(files[side], metas[side], source_rows[side]) for side in files}
                plot_pair(pair, metas, image_paths, ci)
            rows.extend(comparison_rows(ci, meta, source_rows, images))
            cases.append({"case_index": ci, "test_type": meta["test_type"], "figure_png": images[0],
                          "power_current_figure_png": images[1], "missing_sources": missing,
                          "canonical": {s: relative(p, output) for s, p in files.items()},
                          "metadata": {s: relative(p.with_name(p.stem + "_meta.json"), output) for s, p in files.items()},
                          "alignment": {s: pair[s]["timing"] for s in pair}})
            if missing:
                errors.append({"case_index": ci, "reason": "; ".join(missing)})
            print(f"case {ci}: {'partial, missing source' if missing else 'compared'}")
        except (ValueError, KeyError, OSError, TypeError) as exc:
            errors.append({"case_index": ci, "reason": str(exc)})
            print(f"case {ci}: error: {exc}")
    with (output / "compare_summary_freq_reg.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    manifest = {"schema_version": 1, "mode": "freq_reg", "metric_event": "event1",
                "summaries": references, "project_config": relative(paths.project_config, output),
                "compare_csv": "compare_summary_freq_reg.csv",
                "completed": sum(not c["missing_sources"] for c in cases),
                "partial": sum(bool(c["missing_sources"]) for c in cases),
                "exceptions": len(errors), "cases": cases, "errors": errors}
    manifest["timing_review_cases"] = sorted({int(row["case_index"]) for row in rows
                                              if any(str(row.get(s + "_freq_timing_review_required")) == "1"
                                                     for s in ("sim", "hil"))})
    manifest["run_status"] = "needs_timing_review" if manifest["timing_review_cases"] else "complete"
    (output / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    return manifest
