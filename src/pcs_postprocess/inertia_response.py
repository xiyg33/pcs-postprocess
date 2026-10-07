"""Descriptive P response metrics timed from the measured frequency ramp."""

import numpy as np

from . import processing as p
from .step_response import measure_step, smooth_native


SUMMARY_FIELDS = [
    "inertia_metrics_version", "inertia_event_start_input_s", "inertia_metric_end_input_s",
    "inertia_timing_method", "inertia_end_method", "inertia_initial_P_pu",
    "inertia_final_P_pu", "inertia_P_change_pu",
]
VERSION = "frequency_ramp_response_v1"


def recovery_start(t, frequency, expected, deviation, ramp_end):
    """Conservatively exclude the return-frequency edge from the first response."""
    t = np.asarray(t, float)
    frequency = np.asarray(frequency, float)
    local = (t >= expected - .8) & (t <= expected + .6)
    tt, ff = t[local], frequency[local]
    if len(tt) < 3 or not np.isfinite(ff).all():
        return None
    dt = float(np.median(np.diff(tt)))
    ff = smooth_native(tt, ff, 20.)
    before = (tt >= max(ramp_end, expected - .5)) & (tt < expected - .1)
    if before.sum() < 3:
        return None
    baseline = float(np.median(ff[before]))
    threshold = max(.01, .05 * abs(deviation), 3 * float(np.std(ff[before])))
    if threshold >= .5 * abs(deviation):
        return None
    crossed = ((tt >= max(ramp_end, expected - .5)) & (tt <= expected + .5)
               & (-np.sign(deviation) * (ff - baseline) > threshold))
    hold = max(2, int(np.ceil(.02 / dt)))
    if hold > len(tt):
        return None
    hits = np.flatnonzero(np.convolve(crossed.astype(int), np.ones(hold, int), "valid") == hold)
    # The centered filter contains later samples: exclude its entire half-window
    # as well, rather than allowing the restoration to contaminate the P tail.
    return float(tt[hits[0]] - .01) if hits.size else None


def timing(t, sig, meta, ev, offset, offset_source):
    start = float(ev["event1_start"] - offset)
    method = "frequency_ramp_fit" if offset_source == "edge" else "timing_unavailable"
    if meta["time_basis"] == "absolute":
        fitted = p._legacy_frequency_edge_offset(t, sig, meta, ev, 0.)
        if fitted is not None:
            start = float(ev["event1_start"] - fitted)
            method = "frequency_ramp_fit"
    deviation = float(meta["freq_deviation_hz"])
    rate = float(meta["freq_ramp_rate_hz_per_s"])
    if abs(rate) < 1e-12 or deviation * rate <= 0:
        method = "timing_unavailable"
    ramp = abs(deviation / rate) if abs(rate) > 1e-12 else 0.
    expected_end = start + ramp + float(meta["event_hold_s"])
    restoration = recovery_start(t, sig["f"], expected_end, deviation, start + ramp)
    end = min(expected_end, float(t[-1]))
    end_method = "scheduled_boundary" if restoration is None else "measured_frequency_recovery"
    if restoration is not None:
        end = min(end, restoration)
    if t[-1] < expected_end - .02:
        end_method = "recording_truncated"
    return {"inertia_metrics_version": VERSION, "inertia_event_start_input_s": start,
            "inertia_metric_end_input_s": end, "inertia_timing_method": method,
            "inertia_end_method": end_method}


def summary(t, power_pu, reactive_pu, meta, clock):
    start_input = clock["inertia_event_start_input_s"]
    end = clock["inertia_metric_end_input_s"] - start_input + 1.
    data = {"t_s": np.asarray(t) - start_input + 1.,
            "P_pu": np.asarray(power_pu), "Q_pu": np.asarray(reactive_pu)}
    # Use the complete recording before resampling. Baseline/response are
    # smoothed separately so the ramp and restoration cannot leak across them.
    reduced = p.pwr_downsampled_signals(data, 1., end)
    tt, yy = reduced["t_s"], reduced["P_pu"]
    pre, tail = (tt >= .5) & (tt < 1.), (tt >= end - .5) & (tt < end)
    result = {key: None for key in p.PWR_STEP_SUMMARY_COLS}
    result.update(clock, step_timing_basis=clock["inertia_timing_method"],
                  step_response_start_s=1., T_s_method="unavailable", metric_sampling_ms=5,
                  人工确认Ts="待确认", inertia_initial_P_pu=None,
                  inertia_final_P_pu=None, inertia_P_change_pu=None)
    reason = None
    if clock["inertia_timing_method"] == "timing_unavailable":
        reason = "未确认实测频率斜坡起点"
    elif clock["inertia_end_method"] == "recording_truncated" or end < 1.5:
        reason = "观察时长不足"
    elif min(pre.sum(), tail.sum()) < 3:
        reason = "可用采样点不足"
    elif not np.isfinite(yy[pre | tail]).all():
        reason = "录波存在缺失采样"
    if reason is None:
        initial, final = float(np.mean(yy[pre])), float(np.mean(yy[tail]))
        delta = final - initial
        result.update(inertia_initial_P_pu=initial, inertia_final_P_pu=final,
                      inertia_P_change_pu=delta)
        if abs(delta) < 1e-6:
            reason = "实测有功变化幅值过小"
        else:
            measured = measure_step(tt, yy, 1., end,
                                    -np.sign(meta["freq_deviation_hz"]) * abs(delta),
                                    smooth_ms=0., tail_window_s=.5)
            for metric, column in (("M_p", "M_p_pct"), ("T_up", "T_up_s"),
                                   ("T_p", "T_p_s"), ("T_s", "T_s_s")):
                result[column] = measured["values"][metric]
                code = measured["reasons"].get(metric, "")
                result[metric + "_reason"] = p.PWR_STEP_REASON_TEXT.get(code, code).replace("阶跃", "响应")
            if measured["values"]["T_s"] is not None:
                result.update(T_s_method="continuous_5pct", 人工确认Ts="")
    if reason is not None:
        for metric in ("M_p", "T_up", "T_p", "T_s"):
            result[metric + "_reason"] = reason
    return result
