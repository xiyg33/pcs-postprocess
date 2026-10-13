"""Measured voltage-to-breaker response-time differences; no power metrics."""

import numpy as np

from .compare_common import run_comparison, signals
from .compare_presentation import apply_presentation
from .compare_pwr import finite, plt
from .protection_timing import SUMMARY_FIELDS, VERSION, ERROR_STATUSES


FIELDS = ["case_index", "case_name", "test_type", "P_ref_pu", "V_target_pu",
          "event", "signal", "metric", "parameter", "unit", "sim_value", "hil_value",
          "deviation", "absolute_deviation", "deviation_unit", "status", "response_match",
          "figure_png"]
for _side in ("sim", "hil"):
    FIELDS += [f"{_side}_{key}" for key in ("status", "reason", "t1_s", "t2_s", *SUMMARY_FIELDS)]


def prepare_signals(path, meta, summary):
    start = finite(summary["protection_display_start_input_s"], "protection_display_start_input_s")
    data = signals(path, meta, start, start + float(meta["event_hold_s"]))
    data["timing"].update({key: summary[key] for key in SUMMARY_FIELDS})
    data["trip_s"] = finite(summary["t2_s"], "t2_s", optional=True)
    data["delta_t_s"] = finite(summary["delta_t_s"], "delta_t_s", optional=True)
    return data


def comparison_rows(case_id, meta, summaries, figures, target_key="V_target_pu", summary_fields=SUMMARY_FIELDS,
                    unavailable_statuses=ERROR_STATUSES):
    row = dict(case_index=case_id, case_name=meta["case_name"], test_type=meta["test_type"],
               P_ref_pu=meta["P_ref_pu"], event="event1",
               signal="breaker", metric="delta_t", parameter="breaker响应时间",
               unit="s", deviation_unit="s", figure_png=figures[0])
    row[target_key] = meta[target_key]
    states = []
    for side, source in summaries.items():
        state = source["protection_status"]
        value = finite(source["delta_t_s"], "delta_t_s", optional=True)
        start = finite(source["t1_s"], "t1_s", optional=True)
        trip = finite(source["t2_s"], "t2_s", optional=True)
        if state == "tripped":
            if value is None or value < 0 or start != 1. or trip is None or not np.isclose(trip - start, value, atol=1e-8, rtol=0):
                raise ValueError(f"{side}: inconsistent protection timing; rerun single-source analysis")
        elif value is not None or trip is not None:
            raise ValueError(f"{side}: unavailable protection result contains trip metrics")
        if state not in {"tripped", "no_trip", *unavailable_statuses}:
            raise ValueError(f"{side}: unknown protection_status {state!r}")
        states.append(state)
        row.update({f"{side}_value": value, f"{side}_status": "estimated" if value is not None else "unavailable",
                    f"{side}_reason": source["protection_reason"], f"{side}_t1_s": start, f"{side}_t2_s": trip})
        for key in summary_fields:
            row[f"{side}_{key}"] = source[key]
    difference = None if any(row[s + "_value"] is None for s in summaries) else row["sim_value"] - row["hil_value"]
    row.update(deviation=difference, absolute_deviation=None if difference is None else abs(difference),
               status="estimated" if difference is not None else "unavailable",
               response_match="same" if states[0] == states[1] and states[0] in ("tripped", "no_trip") else
                              "trip_mismatch" if set(states) == {"tripped", "no_trip"} else "unavailable")
    return [row]


def response_markers(axes, data, side, color):
    """Show saved timing independently of the signal's scale and trip zoom."""
    trip = data["trip_s"]
    if trip is not None:
        for ax in axes:
            ax.axvline(trip, color=color, ls="-.", lw=.9)
    label = side.upper() + ": "
    if trip is not None:
        label += f"t2={trip:.4f} s"
        if data["delta_t_s"] is not None:
            label += f"; Δt={data['delta_t_s']:.4f} s"
    else:
        label += data["timing"]["protection_status"]
    axes[0].text(.02, .96 if side == "hil" else .79, label,
                 transform=axes[0].transAxes, va="top", color=color, fontsize=8,
                 bbox=dict(facecolor="white", alpha=.85, edgecolor="none"))


def full_recording_stop(pair, metas):
    """Keep the configured event and post-window without extrapolating data."""
    stops = []
    for side, data in pair.items():
        expected = 1. + float(metas[side]["event_hold_s"]) + float(metas[side].get("record_post_s", 2.))
        available = float(data["t_s"][-1])
        dt = float(np.median(np.diff(data["t_s"]))) if len(data["t_s"]) > 1 else 0.
        # A sampled window commonly ends one sample before its requested end.
        stops.append(expected if available + dt + 1e-8 >= expected else available)
    return max(stops)


def plot_pair(pair, metas, paths, case_id, presentation=None):
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
    try:
        for side, color, style in (("hil", "tab:red", "-"), ("sim", "tab:blue", "--")):
            data = pair[side]
            for ax, key in zip(axes[:3], ("Vp_pu", "P_pu", "Q_pu")):
                ax.plot(data["display"]["t_s"], data["display"][key], color=color, ls=style, lw=.9, label=side.upper())
            response_markers(axes, data, side, color)
        for ax, label in zip(axes, ("正序电压 / pu", "有功功率 / pu", "无功功率 / pu")):
            ax.set_ylabel(label)
            ax.axvline(1., color="gray", ls="--", lw=.8)
            ax.grid(True, alpha=.3)
            if ax.get_legend_handles_labels()[0]:
                ax.legend(loc="upper right", ncol=2, fontsize=8)
        axes[-1].set_xlim(0., full_recording_stop(pair, metas))
        axes[-1].set_xlabel("归一化时间 / s；实测电压扰动起点 = 1 s")
        notes = " / ".join(f"{s.upper()}: {pair[s]['timing']['protection_status']}" for s in pair)
        fig.suptitle(f"Case {case_id}: {metas['hil']['case_name']} — SIM / HIL\n{notes}", fontsize=10)
        fig.tight_layout(rect=(0, 0, 1, .95))
        apply_presentation(fig, case_id, presentation)
        fig.savefig(paths[0], dpi=150, facecolor="white")
    finally:
        plt.close(fig)


def source_error(source):
    return source["protection_reason"] if source["protection_status"] in ERROR_STATUSES else None


def run(paths_config, output, case_id=None, *, case_ids=None, presentation=None):
    return run_comparison("volt_prot", paths_config, output, case_id,
                          ["t1_s", "t2_s", "delta_t_s", *SUMMARY_FIELDS],
                          ["P_ref_pu", "V_target_pu", "event_start_s", "event_hold_s"],
                          prepare_signals, plot_pair, comparison_rows, FIELDS,
                          "protection_metrics_version", VERSION, 1, source_error,
                          case_ids=case_ids, presentation=presentation)
