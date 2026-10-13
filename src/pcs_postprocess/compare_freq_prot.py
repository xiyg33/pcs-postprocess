"""Measured line-frequency disturbance to breaker total response comparison."""

import csv
import json

import numpy as np

from . import processing as p, frequency_protection as timing
from .compare_common import run_comparison, signals
from .compare_presentation import PRESENTATION_FIELDS, apply_presentation
from .compare_pwr import finite, plt
from .compare_volt_prot import comparison_rows as protection_rows, response_markers, full_recording_stop
from .compare_scr import write_csv


FIELDS = ["case_index", "case_name", "test_type", "P_ref_pu", "f_target_hz",
          "event", "signal", "metric", "parameter", "unit", "sim_value", "hil_value",
          "deviation", "absolute_deviation", "deviation_unit", "status", "response_match", "figure_png"]
for _side in ("sim", "hil"):
    FIELDS += [f"{_side}_{key}" for key in ("status", "reason", "t1_s", "t2_s", *timing.SUMMARY_FIELDS)]
DIAGNOSTIC_FIELDS = ["case_index", "case_name", "source", "f_target_hz", "delta_t_s",
                     "protection_event_start_input_s", "protection_trip_input_s",
                     "protection_status", "protection_reason", *timing.DIAGNOSTIC_FIELDS]


def prepare_signals(path, meta, summary):
    start = finite(summary["protection_display_start_input_s"], "protection_display_start_input_s")
    data = signals(path, meta, start, start + float(meta["event_hold_s"]))
    raw_t, sig = p.load_case_file(str(path), absolute=meta["time_basis"] == "absolute")
    vabc = np.column_stack([sig["V" + phase] for phase in "abc"])
    _, frequency, _, half_cycle = timing.voltage_phase(raw_t, vabc, p.NOMINAL_FREQUENCY_HZ)
    if sig["breaker"] is not None:
        for i in np.flatnonzero(np.diff(sig["breaker"]) != 0) + 1:
            # The voltage estimator is auxiliary; switching impulses should
            # not look like measured grid-frequency excursions on this plot.
            frequency[np.abs(raw_t - raw_t[i]) <= half_cycle + .01] = np.nan
    # Reconstruct the displayed frequency only; saved event/trip metrics are
    # never refitted by compare, even when the waveform estimator changes.
    data["display"]["voltage_frequency_hz"] = np.interp(data["display"]["t_s"] + start - 1., raw_t, frequency)
    data["timing"].update({key: summary[key] for key in timing.SUMMARY_FIELDS})
    data["trip_s"] = finite(summary["t2_s"], "t2_s", optional=True)
    data["delta_t_s"] = finite(summary["delta_t_s"], "delta_t_s", optional=True)
    raw_trip = finite(summary["freq_breaker_trip_input_s"], "freq_breaker_trip_input_s", optional=True)
    if data["trip_s"] is None and raw_trip is not None:
        data["trip_s"] = raw_trip - start + 1.
    return data


def comparison_rows(case_id, meta, summaries, figures):
    rows = protection_rows(case_id, meta, summaries, figures, "f_target_hz", timing.SUMMARY_FIELDS,
                           timing.UNAVAILABLE_STATUSES)
    for row in rows:
        row["parameter"] = "频率扰动至脱网总响应时间"
    return rows


def source_error(source):
    return source["protection_reason"] if source["protection_status"] in timing.ERROR_STATUSES else None


def plot_pair(pair, metas, paths, case_id, presentation=None):
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
    try:
        for side, color, style in (("hil", "tab:red", "-"), ("sim", "tab:blue", "--")):
            data = pair[side]
            display = data["display"]
            axes[0].plot(display["t_s"], display["voltage_frequency_hz"], color=color, ls=style,
                         lw=1., label=side.upper() + " 实测线路频率")
            if "f_hz" in display:
                axes[0].plot(display["t_s"], display["f_hz"], color=color, ls=":", lw=.7,
                             alpha=.7, label=side.upper() + " 录波 PLL（辅助）")
            for ax, key in zip(axes[1:3], ("P_pu", "Q_pu")):
                ax.plot(display["t_s"], display[key], color=color, ls=style, lw=.9, label=side.upper())
            response_markers(axes, data, side, color)
        for ax, label in zip(axes, ("频率 / Hz", "有功功率 / pu", "无功功率 / pu")):
            ax.set_ylabel(label)
            ax.axvline(1., color="gray", ls="--", lw=.8)
            ax.grid(True, alpha=.3)
            if ax.get_legend_handles_labels()[0]:
                ax.legend(loc="upper right", ncol=2, fontsize=8)
        stop = full_recording_stop(pair, metas)
        axes[-1].set_xlim(0., stop)
        values = []
        for data in pair.values():
            d = data["display"]
            keep = (d["t_s"] >= 0.) & (d["t_s"] <= stop)
            measured = d["voltage_frequency_hz"][keep]
            measured = measured[np.isfinite(measured)]
            if measured.size:
                values.extend(np.percentile(measured, [1., 99.]).tolist())
            if "f_hz" in d:
                values.extend(d["f_hz"][keep].tolist())
        values.extend(float(meta["f_target_hz"]) for meta in metas.values())
        values.append(p.NOMINAL_FREQUENCY_HZ)
        valid = np.asarray(values)[np.isfinite(values)]
        if valid.size:
            margin = max(.2, .1 * float(np.ptp(valid)))
            axes[0].set_ylim(float(valid.min()) - margin, float(valid.max()) + margin)
        axes[-1].set_xlabel("归一化时间 / s（实测线路扰动开始 = 1 s；未确认时按名义起点展示）")
        notes = " / ".join(f"{s.upper()}: {pair[s]['timing']['protection_status']}" for s in pair)
        fig.suptitle(f"Case {case_id}: {metas['hil']['case_name']}\n频率扰动至脱网总响应时间 — {notes}", fontsize=10)
        fig.tight_layout(rect=(0, 0, 1, .95))
        apply_presentation(fig, case_id, presentation)
        fig.savefig(paths[0], dpi=150, facecolor="white")
    finally:
        plt.close(fig)


def run(paths_config, output, case_id=None, overwrite=False, *, case_ids=None, presentation=None):
    manifest = run_comparison("freq_prot", paths_config, output, case_id,
                              ["t1_s", "t2_s", "delta_t_s", *timing.SUMMARY_FIELDS],
                              ["P_ref_pu", "f_target_hz", "event_start_s", "event_hold_s"],
                              prepare_signals, plot_pair, comparison_rows, FIELDS,
                              "protection_metrics_version", timing.VERSION, 1, source_error,
                              tolerate_source_errors=True, overwrite=overwrite,
                              case_ids=case_ids, presentation=presentation)
    with (output / "compare_summary_freq_prot.csv").open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    diagnostics = []
    for row in rows:
        for side in ("sim", "hil"):
            record = {key: row[key] for key in ("case_index", "case_name", "f_target_hz")}
            record.update(source=side, delta_t_s=row[side + "_value"])
            record.update({key: row[side + "_" + key] for key in DIAGNOSTIC_FIELDS if side + "_" + key in row})
            if presentation is not None:
                record.update({key: row[key] for key in PRESENTATION_FIELDS})
            diagnostics.append(record)
    write_csv(output / "timing_diagnostics_freq_prot.csv",
              [*DIAGNOSTIC_FIELDS, *PRESENTATION_FIELDS] if presentation is not None else DIAGNOSTIC_FIELDS,
              diagnostics)
    manifest["diagnostics_csv"] = "timing_diagnostics_freq_prot.csv"
    manifest["timing_interpretation"] = "Total response from measured line-frequency disturbance onset to sustained breaker falling edge; PLL thresholds do not select timing."
    (output / manifest["manifest_file"]).write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    notes = ("# 频率保护计时解释\n\n"
             "Δt = breaker 持续下降沿时刻 − 实测线路频率扰动开始时刻。图上扰动起点为 1 s。\n\n"
             "线路起点由三相电压旋转相量和连续分段相位拟合估算，包含测频动态、保护计时及输出链路的总响应。"
             "它不等于控制器内部保护定时器延时。录波 PLL 及其门槛不参与主指标。\n\n"
             "实测起点未确认、缺 breaker 或未跳闸时，指标留空并说明原因。"
             "日常查看 figures 和 compare_summary_freq_prot.csv；诊断 CSV 记录拟合质量，"
             "唯一 manifest.json 用于输入来源、版本、对时依据和异常追溯。\n")
    (output / "timing_notes.md").write_text(notes, encoding="utf-8")
    return manifest
