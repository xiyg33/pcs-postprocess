"""Frequency-ramp-aligned SIM/HIL overlays and referenced P response metrics."""

from . import processing as p
from .compare_common import run_comparison, signals
from .compare_freq_reg import plot_pair
from .compare_pwr import FIELDS as PWR_FIELDS, comparison_rows as step_rows, finite
from .inertia_response import SUMMARY_FIELDS, VERSION


FIELDS = [*PWR_FIELDS, "power_current_figure_png", "P_ref_pu", "freq_deviation_hz",
          "freq_ramp_rate_hz_per_s", "interpretation"]
for _side in ("sim", "hil"):
    FIELDS += [f"{_side}_{key}" for key in SUMMARY_FIELDS]


def prepare_signals(path, meta, summary):
    start = finite(summary["inertia_event_start_input_s"], "inertia_event_start_input_s")
    end = finite(summary["inertia_metric_end_input_s"], "inertia_metric_end_input_s")
    data = signals(path, meta, start, end, 1., end - start + 1., currents=True)
    data["timing"].update({key: summary[key] for key in SUMMARY_FIELDS})
    return data


def comparison_rows(case_id, meta, summaries, figures):
    rows = step_rows(case_id, dict(meta, test_type_code=3), summaries, figures[0])
    for row in rows:
        if row["metric"] == "T_up":
            row["parameter"] = "90%响应时间"
        row.update(power_current_figure_png=figures[1], P_ref_pu=meta["P_ref_pu"],
                   freq_deviation_hz=meta["freq_deviation_hz"],
                   freq_ramp_rate_hz_per_s=meta["freq_ramp_rate_hz_per_s"],
                   interpretation="频率斜坡下的有功响应指标，非功率指令阶跃验收指标")
        for side in summaries:
            for key in SUMMARY_FIELDS:
                row[f"{side}_{key}"] = summaries[side][key]
    return rows


def run(paths_config, output, case_id=None, *, case_ids=None, presentation=None):
    return run_comparison("inertia", paths_config, output, case_id,
                          [*p.PWR_STEP_SUMMARY_COLS, *SUMMARY_FIELDS],
                          ["P_ref_pu", "freq_deviation_hz", "freq_ramp_rate_hz_per_s",
                           "event_start_s", "event_hold_s"],
                          prepare_signals, plot_pair, comparison_rows, FIELDS,
                          "inertia_metrics_version", VERSION, 2,
                          case_ids=case_ids, presentation=presentation)
