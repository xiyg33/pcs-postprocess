"""SCR event provenance and common-horizon operating stability evidence."""

import json

import numpy as np
from scipy.ndimage import uniform_filter1d

from . import processing as p
from .breaker_events import evidence


VERSION = "scr_common_horizon_v2"
SUMMARY_FIELDS = ["scr_metrics_version", "scr_event_start_input_s", "scr_recovery_input_s",
                  "scr_display_start_input_s", "scr_timing_status", "scr_timing_reason",
                  "scr_input_end_s", "scr_input_sampling_s", "scr_breaker_status",
                  "scr_breaker_reason", "scr_trip_times_input_s", "scr_incomplete_trip_times_input_s",
                  "scr_pre_disconnected", "V_target_pu", "scr_configured_start_input_s",
                  "scr_configured_end_input_s", "scr_display_end_input_s",
                  "disturbance_start_s", "disturbance_duration_s", "sim_time_s", "record_pre_s", "record_post_s"]


def local_voltage_edge(t, voltage, hint, target, absolute):
    """Weak-grid baselines can drift; estimate their local trend, not 1 pu.

    This fallback still requires a measured directed and sustained departure.
    It never assigns a command timestamp to an indistinguishable voltage edge.
    """
    dt = float(np.median(np.diff(t)))
    smooth = uniform_filter1d(np.asarray(voltage, float), max(1, int(round(.020 / dt)) | 1), mode="nearest")
    baseline = (t >= hint - .30) & (t <= hint - .03)
    if baseline.sum() < 10:
        return None
    design = np.column_stack((np.ones(baseline.sum()), t[baseline] - hint))
    coef = np.linalg.lstsq(design, smooth[baseline], rcond=None)[0]
    residual = smooth[baseline] - design @ coef
    noise = 1.4826 * float(np.median(np.abs(residual - np.median(residual))))
    threshold = max(.008, .10 * abs(target - 1.), 3 * noise)
    if threshold >= .8 * abs(target - 1.):
        return None
    trend = coef[0] + coef[1] * (t - hint)
    search = (t >= hint - (0. if absolute else .08)) & (t <= hint + .10)
    crossed = search & (np.sign(target - 1.) * (smooth - trend) > threshold)
    hold = max(2, int(np.ceil(.020 / dt)))
    hits = np.flatnonzero(np.convolve(crossed.astype(int), np.ones(hold, int), "valid") >= hold)
    return float(t[hits[0]]) if hits.size else None


def timing(t, sig, meta, ev, prior):
    t = np.asarray(t, float)
    hint = float(meta.get("scr_event_search_center_input_s", ev["event1_start"] - prior))
    reference = dict(meta, V_fault_pu=meta["V_target_pu"])
    offset = p.detect_frt_edge_offset(t, sig["Vp_pu"], reference, ev, ev["event1_start"] - hint)
    start = None if offset is None else float(ev["event1_start"] - offset)
    if start is None:
        start = local_voltage_edge(t, sig["Vp_pu"], hint, float(meta["V_target_pu"]), meta["time_basis"] == "absolute")
    # Search clearing near the independent scheduled/hinted input time. A
    # delayed threshold crossing must not shift the clearing search window.
    nominal_offset = ev["event1_start"] - hint
    recovery = None if start is None else p.detect_frt_recovery_time(t, sig["Vp_pu"], reference, ev, nominal_offset)
    recovery = None if recovery is None else float(recovery - nominal_offset)
    if (recovery is not None and meta["time_basis"] == "absolute"
            and recovery < hint + float(meta["disturbance_duration_s"]) - .010):
        # Voltage regulation can restore PCC magnitude while the source pulse
        # is still active. That response is not evidence of fault clearing.
        recovery = None
    br = evidence(t, sig.get("breaker"))
    display = hint if start is None else start
    reason = "disturbance_unconfirmed" if start is None else "recovery_unconfirmed" if recovery is None else ""
    index = int(np.searchsorted(t, display))
    disconnected = None
    if br["status"] == "valid" and 0 < index <= len(t):
        disconnected = int(sig["breaker"][index - 1] < .5)
    return {"scr_metrics_version": VERSION, "scr_event_start_input_s": start,
            "scr_recovery_input_s": recovery, "scr_display_start_input_s": display,
            "scr_timing_status": "confirmed" if not reason else "unconfirmed", "scr_timing_reason": reason,
            "scr_input_end_s": float(t[-1]), "scr_input_sampling_s": float(np.median(np.diff(t))),
            "scr_breaker_status": br["status"], "scr_breaker_reason": br["reason"],
            "scr_trip_times_input_s": json.dumps(br["falls"]),
            "scr_incomplete_trip_times_input_s": json.dumps(br["incomplete_falls"]),
            "scr_pre_disconnected": disconnected, "V_target_pu": meta["V_target_pu"],
            "scr_configured_start_input_s": hint,
            "scr_configured_end_input_s": hint + float(meta["disturbance_duration_s"]),
            "scr_display_end_input_s": hint + float(meta["disturbance_duration_s"]) + float(ev["record_post"]),
            **{key: meta.get(key) for key in ("disturbance_start_s", "disturbance_duration_s", "sim_time_s", "record_pre_s", "record_post_s")}}


def classify(data, summary, horizon):
    """Trip remains evidence even when the oscillation window is unavailable."""
    result = {"state": "UNDETERMINED", "reason": "", "rms_early_pu": None,
              "rms_late_pu": None, "growth_ratio": None, "growth_noise_floor_pu": None,
              "baseline_std_pu": None, "trip_input_s": None,
              "observation_start_after_recovery_s": .3, "observation_end_after_recovery_s": horizon}
    falls = json.loads(summary["scr_trip_times_input_s"])
    start = optional_float(summary["scr_event_start_input_s"])
    recovery = optional_float(summary["scr_recovery_input_s"])
    # Only events within the observation horizon count; a later failure is
    # not allowed to penalize a source whose partner ended recording earlier.
    observation_end = None if recovery is None or horizon is None else recovery + horizon
    relevant = [x for x in falls if observation_end is None or x <= observation_end]
    if relevant:
        result.update(state="TRIP", trip_input_s=relevant[0],
                      reason="pre_disturbance_trip" if start is not None and relevant[0] < start else "breaker_trip")
        return result
    if summary["scr_breaker_status"] != "valid":
        result["reason"] = summary["scr_breaker_reason"]
        return result
    if str(summary["scr_pre_disconnected"]) in ("1", "1.0"):
        result["reason"] = "already_open_without_confirmed_falling_edge"
        return result
    if summary["scr_timing_status"] != "confirmed":
        result["reason"] = summary["scr_timing_reason"]
        return result
    if horizon is None or horizon - .3 < 1. - 1e-9:
        result["reason"] = "insufficient_common_observation"
        return result
    if any(start <= x <= observation_end for x in json.loads(summary["scr_incomplete_trip_times_input_s"])):
        result["reason"] = "incomplete_breaker_transition"
        return result
    display = data["display"]
    onset_norm = float(data["timing"]["event_start_input_s"])
    recovery_norm = recovery - onset_norm + 1.
    tt, power = display["t_s"], display["P_pu"]
    grid = np.arange(.3, horizon + 1e-10, p.PWR_METRIC_SAMPLE_S)
    native_t = data["input_t_s"]
    inside = (native_t >= recovery + .3) & (native_t <= recovery + horizon)
    dt = float(summary["scr_input_sampling_s"])
    if np.any(np.diff(native_t[inside]) > 3 * dt):
        result["reason"] = "sampling_gap_in_observation"
        return result
    values = np.interp(recovery_norm + grid, tt, power, left=np.nan, right=np.nan)
    baseline = power[(tt >= .5) & (tt < .95)]
    if not np.isfinite(values).all() or baseline.size < 3 or not np.isfinite(baseline).all():
        result["reason"] = "missing_power_or_baseline"
        return result
    half = len(values) // 2
    early, late = values[:half], values[half:]
    rms_early = float(np.std(early))
    rms_late = float(np.std(late))
    ratio = rms_late / max(rms_early, 1e-6)
    noise = float(np.std(baseline))
    floor = max(.01, 3 * noise)
    # A five-sample span on the 5 ms grid covers twenty milliseconds.
    high = np.abs(late) >= 2.
    overload = len(high) >= 5 and np.any(np.convolve(high.astype(int), np.ones(5, int), "valid") >= 5)
    growth = ratio >= 1.3 and rms_late >= floor
    result.update(state="UNSTABLE" if overload or growth else "STABLE",
                  reason="power_overlimit" if overload else "oscillation_growth" if growth else "no_failure_in_common_window",
                  rms_early_pu=rms_early, rms_late_pu=rms_late, growth_ratio=ratio,
                  growth_noise_floor_pu=floor, baseline_std_pu=noise)
    return result


def optional_float(value):
    return None if value is None or str(value).strip().lower() in ("", "nan", "none") else float(value)


def boundary(points):
    """Report tested brackets only; an unknown point interrupts adjacency."""
    ordered = sorted(points, key=lambda x: x["scr"], reverse=True)
    transitions = []
    for high, low in zip(ordered, ordered[1:]):
        if high["state"] == "UNDETERMINED" or low["state"] == "UNDETERMINED":
            continue
        high_ok, low_ok = high["state"] == "STABLE", low["state"] == "STABLE"
        if high_ok != low_ok:
            transitions.append({"lower_scr": low["scr"], "upper_scr": high["scr"],
                                "higher_state": high["state"], "lower_state": low["state"]})
    known = [x for x in ordered if x["state"] != "UNDETERMINED"]
    # A failed higher SCR followed by a stable lower SCR contradicts the
    # usual single-boundary ordering, even across an unknown intermediate point.
    nonmonotonic = any(a["state"] != "STABLE" and b["state"] == "STABLE"
                      for i, a in enumerate(known) for b in known[i + 1:])
    stable = [x["scr"] for x in known if x["state"] == "STABLE"]
    failed = [x["scr"] for x in known if x["state"] != "STABLE"]
    status = ("nonmonotonic" if nonmonotonic else "bracketed" if len(transitions) == 1
              else "incomplete_boundary" if stable and failed else "stable_side_only" if stable
              else "failed_side_only" if failed else "undetermined")
    bracket = transitions[0] if status == "bracketed" else {}
    return {"status": status, "lower_scr": bracket.get("lower_scr"), "upper_scr": bracket.get("upper_scr"),
            "minimum_tested_stable_scr": min(stable) if stable else None,
            "maximum_tested_failed_scr": max(failed) if failed else None,
            "transitions": transitions, "points": ordered}
