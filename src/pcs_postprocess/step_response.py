"""Native-time first-step measurements for table D.1 (P or Q)."""
from __future__ import annotations

import numpy as np

METRICS = (("M_p", "超调量", "%", 10.0),
           ("T_up", "上升时间", "s", 0.1),
           ("T_p", "峰值时间", "s", 0.2),
           ("T_s", "稳定时间", "s", 1.0))


def smooth_native(t, y, smooth_ms=20.0):
    """Centered boxcar with edge padding; never silently fill missing samples."""
    t, y = np.asarray(t, float), np.asarray(y, float)
    if t.size < 2:
        return y.copy()
    n = max(1, int(round(smooth_ms / 1000.0 / np.median(np.diff(t)))))
    n = min(n, len(y))
    if n % 2 == 0:
        n = max(1, n - 1)  # symmetric filter, no half-sample timing shift
    return np.convolve(np.pad(y, (n // 2, n // 2), mode="edge"),
                       np.ones(n) / n, mode="valid")


def measure_step(t, y, start_s, end_s, command_delta, smooth_ms=20.0):
    """Measure one interval [start, end), excluding the next command.

    Baseline needs 0.5 s; the tail needs at least 0.1 s. Confirm a stationary
    mean using half-window mean difference and fitted drift, each <=5% of the
    measured step. Ripple is retained: Ts independently requires every sample
    in its remaining observation to lie inside the 5% band for at least a tail
    window. The end is clipped by the caller to recorded observations.
    """
    out = {"values": {key: None for key, *_ in METRICS}, "reasons": {},
           "start_s": float(start_s), "end_s": float(end_s),
           "baseline_window_s": [float(start_s - .5), float(start_s)],
           "tail_window_s": None, "initial_pu": None, "final_pu": None,
           "baseline_stable": False, "tail_stable": False,
           "smooth_ms": smooth_ms, "stability_band_fraction": .05}

    def invalid(reason):
        out["reasons"] = {key: reason for key, value in out["values"].items() if value is None}
        return out

    t, y = np.asarray(t, float).reshape(-1), np.asarray(y, float).reshape(-1)
    if t.size != y.size or t.size < 3 or not np.all(np.isfinite(t)) or np.any(np.diff(t) <= 0):
        return invalid("invalid_time_axis_or_signal_length")
    if not np.isfinite(start_s):
        return invalid("timing_unavailable")
    if not np.isfinite(command_delta) or abs(command_delta) < 1e-6:
        return invalid("zero_command_step")
    dt = float(np.median(np.diff(t)))
    tail_s = min(.5, (end_s - start_s) * .2)
    out["tail_window_s"] = [float(end_s - tail_s), float(end_s)]
    if tail_s < .1 or t[0] > start_s - .5 + dt * 1.1 or t[-1] < end_s - dt * 1.1:
        return invalid("insufficient_observation")
    # Filter each side separately so the next step cannot leak backwards into
    # the steady tail, nor can the first step leak into the baseline.
    pre = (t >= start_s - .5) & (t < start_s)
    event = (t >= start_s) & (t < end_s)
    if np.count_nonzero(pre) < 3 or np.count_nonzero(event) < 3:
        return invalid("insufficient_samples")
    if (np.any(~np.isfinite(y[pre | event]))
            or np.any(np.diff(t[pre | event]) > 3 * dt)):
        return invalid("missing_samples")
    pre_y = smooth_native(t[pre], y[pre], smooth_ms)
    tt = t[event]
    yy = smooth_native(tt, y[event], smooth_ms)
    tail = tt >= end_s - tail_s
    if np.count_nonzero(tail) < 3:
        return invalid("insufficient_samples")
    initial, final = float(np.mean(pre_y)), float(np.mean(yy[tail]))
    delta = final - initial
    out.update(initial_pu=initial, final_pu=final, measured_delta_pu=delta,
               baseline_sample_count=int(pre.sum()), tail_sample_count=int(tail.sum()))
    if abs(delta) < max(1e-6, abs(command_delta) * .01):
        return invalid("near_zero_measured_step")
    if np.sign(delta) != np.sign(command_delta):
        return invalid("response_direction_mismatch")
    band = .05 * abs(delta)
    pre_deviation = float(np.max(np.abs(pre_y - initial)))
    tail_deviation = float(np.max(np.abs(yy[tail] - final)))
    def stationarity(time, values):
        halfway = len(values) // 2
        half_difference = abs(float(np.mean(values[:halfway]) - np.mean(values[halfway:])))
        drift = abs(float(np.polyfit(time - time[0], values, 1)[0] * (time[-1] - time[0])))
        return dict(half_mean_difference_pu=half_difference, fitted_drift_pu=drift,
                    stable=max(half_difference, drift) <= band)

    pre_stationary = stationarity(t[pre], pre_y)
    tail_stationary = stationarity(tt[tail], yy[tail])
    out.update(baseline_max_deviation_pu=pre_deviation,
               tail_max_deviation_pu=tail_deviation,
               baseline_stationarity=pre_stationary, tail_stationarity=tail_stationary,
               baseline_stable=pre_stationary["stable"], tail_stable=tail_stationary["stable"])
    directed = np.sign(delta) * (yy - initial)
    peak_i = int(np.argmax(directed))
    out["values"]["T_p"] = float(tt[peak_i] - start_s)
    out["peak_time_s"] = float(tt[peak_i])
    if not out["baseline_stable"]:
        return invalid("baseline_unstable")
    if not out["tail_stable"]:
        return invalid("final_steady_state_unconfirmed")
    out["values"]["M_p"] = max(0.0, float((directed[peak_i] / abs(delta) - 1) * 100))
    hits = np.flatnonzero(directed >= .9 * abs(delta))
    if hits.size:
        i = int(hits[0])
        crossing = tt[i]
        if i > 0 and directed[i] != directed[i - 1]:
            crossing = tt[i - 1] + ((.9 * abs(delta) - directed[i - 1])
                       / (directed[i] - directed[i - 1]) * (tt[i] - tt[i - 1]))
        out["values"]["T_up"] = float(crossing - start_s)
    else:
        out["reasons"]["T_up"] = "ninety_percent_not_reached"
    outside = np.flatnonzero(np.abs(yy - final) > band + 1e-12 * max(1., abs(delta)))
    settle_i = int(outside[-1] + 1) if outside.size else 0
    if settle_i < len(tt) and end_s - tt[settle_i] >= tail_s - dt:
        out["values"]["T_s"] = float(tt[settle_i] - start_s)
    else:
        out["reasons"]["T_s"] = "not_settled_within_observation"
    return out


def comparison_rows(sim, hil, timing_reliable, signal):
    rows = []
    for key, label, unit, limit in METRICS:
        sv, hv = sim["values"][key], hil["values"][key]
        deviation = None if sv is None or hv is None else sv - hv
        reasons = [f"{side}: {item['reasons'][key]}" for side, item in
                   (("SIM", sim), ("HIL", hil)) if key in item["reasons"]]
        if deviation is None:
            status = "unavailable"
        elif key != "M_p" and not timing_reliable:
            status = "estimated"
            reasons.append("response_edge_estimate_not_formal_acceptance")
        else:
            status = "pass" if abs(deviation) <= limit else "fail"
        rows.append(dict(signal=signal, metric=key, parameter=label, unit=unit,
                         hil_value=hv, sim_value=sv, deviation=deviation,
                         allowed_deviation=limit,
                         deviation_unit="percentage_points" if key == "M_p" else "s",
                         status=status, reason="; ".join(reasons)))
    return rows
