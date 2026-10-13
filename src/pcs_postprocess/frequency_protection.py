"""Measured line-frequency disturbance to native breaker total response."""

import numpy as np
from scipy.ndimage import uniform_filter1d

from . import processing as p
from .breaker_events import evidence
from .protection_timing import SUMMARY_FIELDS as PROTECTION_FIELDS


VERSION = "voltage_disturbance_breaker_v2"
ERROR_STATUSES = {"missing_breaker", "invalid_breaker", "already_open", "insufficient_observation", "disturbance_unconfirmed"}
UNAVAILABLE_STATUSES = ERROR_STATUSES
DIAGNOSTIC_FIELDS = [
    "freq_nominal_start_input_s", "freq_input_sampling_s", "freq_fit_phase_rms_rad",
    "freq_fit_baseline_hz", "freq_fit_after_hz", "freq_fit_window_start_input_s",
    "freq_fit_window_end_input_s", "freq_fit_reason", "freq_breaker_trip_input_s",
    "freq_voltage_start_input_s", "freq_voltage_to_breaker_s", "freq_nominal_delta_t_s",
]
SUMMARY_FIELDS = [*PROTECTION_FIELDS, *DIAGNOSTIC_FIELDS]


def voltage_phase(t, vabc, nominal_hz):
    """Demodulate abc, remove one nominal cycle, and retain phase timestamps.

    The centered filter is zero-phase. Its half-cycle support is excluded at
    recording boundaries and around breaker transitions before fitting.
    """
    t, v = np.asarray(t, float), np.asarray(vabc, float)
    dt = float(np.median(np.diff(t)))
    n = max(3, int(round(1. / nominal_hz / dt)) | 1)
    z = (2 * v[:, 0] - v[:, 1] - v[:, 2]) / 3 + 1j * (v[:, 1] - v[:, 2]) / np.sqrt(3)
    rotating = z * np.exp(-2j * np.pi * nominal_hz * (t - t[0]))
    filtered = (uniform_filter1d(rotating.real, n, mode="nearest")
                + 1j * uniform_filter1d(rotating.imag, n, mode="nearest"))
    phase = np.unwrap(np.angle(filtered))
    amplitude = np.abs(filtered)
    f = nominal_hz + np.gradient(phase, t) / (2 * np.pi)
    valid = amplitude >= .1 * p.HIL_V_BASE_PEAK_V
    valid[:n // 2] = False
    valid[-n // 2:] = False
    f[~valid] = np.nan
    return phase, f, valid, n // 2 * dt


def fit_start(t, vabc, nominal, target_hz, breaker_info):
    t = np.asarray(t, float)
    dt = float(np.median(np.diff(t)))
    phase, frequency, valid, half = voltage_phase(t, vabc, p.NOMINAL_FREQUENCY_HZ)
    start_window = max(float(t[0]) + half, nominal - .5)
    end_window = min(float(t[-1]) - half, nominal + .5)
    # Fit only the connected waveform preceding the first relevant trip.
    falls = [x for x in breaker_info["falls"] if x >= nominal - .35]
    if falls:
        end_window = min(end_window, falls[0] - half - dt)
    result = {"start": None, "frequency": frequency, "freq_fit_phase_rms_rad": None,
              "freq_fit_baseline_hz": None, "freq_fit_after_hz": None,
              "freq_fit_window_start_input_s": start_window,
              "freq_fit_window_end_input_s": end_window, "freq_fit_reason": ""}

    def reject(reason):
        result["freq_fit_reason"] = reason
        return result

    if end_window - start_window < .2:
        return reject("phase_fit_insufficient_observation")
    keep = (t >= start_window) & (t <= end_window)
    if not valid[keep].all() or np.any(np.diff(t[keep]) > 3 * dt):
        return reject("phase_fit_invalid_voltage_or_sampling_gap")
    stride = max(1, int(round(.001 / dt)))
    tt, yy = t[keep][::stride], phase[keep][::stride]
    lo, hi = max(nominal - .35, start_window + .12), min(nominal + .35, end_window - .04)
    if hi <= lo:
        return reject("phase_fit_insufficient_sides")
    step = max(dt, .001)
    candidates = np.arange(lo, hi, step)
    scores, coefficients = [], []
    # Subtract a local origin so translation of absolute clocks cannot change
    # the conditioning of the least-squares design matrix.
    for candidate in candidates:
        design = np.column_stack((np.ones(tt.shape), tt - nominal, np.maximum(tt - candidate, 0.)))
        coef = np.linalg.lstsq(design, yy, rcond=None)[0]
        scores.append(float(np.mean((yy - design @ coef) ** 2)))
        coefficients.append(coef)
    k = int(np.argmin(scores))
    candidate, coef = float(candidates[k]), coefficients[k]
    baseline = p.NOMINAL_FREQUENCY_HZ + float(coef[1]) / (2 * np.pi)
    after = baseline + float(coef[2]) / (2 * np.pi)
    residual = float(np.sqrt(scores[k]))
    result.update(freq_fit_phase_rms_rad=residual, freq_fit_baseline_hz=baseline, freq_fit_after_hz=after)
    change = target_hz - p.NOMINAL_FREQUENCY_HZ
    if k <= 1 or k >= len(candidates) - 2:
        return reject("phase_fit_search_boundary")
    if abs(baseline - p.NOMINAL_FREQUENCY_HZ) > .2:
        return reject("phase_fit_unstable_baseline")
    if abs(change) < .03 or np.sign(change) * (after - baseline) < .5 * abs(change):
        return reject("phase_fit_wrong_or_small_change")
    if abs(after - target_hz) > max(.2, .3 * abs(change)):
        return reject("phase_fit_target_mismatch")
    if residual > .10:
        return reject("phase_fit_large_residual")
    linear = np.column_stack((np.ones(tt.shape), tt - nominal))
    straight = np.linalg.lstsq(linear, yy, rcond=None)[0]
    null_loss = float(np.mean((yy - linear @ straight) ** 2))
    if scores[k] >= .5 * null_loss:
        return reject("phase_fit_no_confirmed_breakpoint")
    result["start"] = candidate
    return result


def measure(t, sig, meta, ev, prior):
    t = np.asarray(t, float)
    dt = float(np.median(np.diff(t)))
    nominal = float(ev["event1_start"] - prior)
    br = evidence(t, sig.get("breaker"))
    vabc = np.column_stack([sig["V" + phase] for phase in "abc"])
    fitted = fit_start(t, vabc, nominal, float(meta["f_target_hz"]), br)
    start = fitted["start"]
    display = nominal if start is None else start
    end = display + float(meta["event_hold_s"])
    falls = [x for x in br["falls"] if display <= x <= end]
    raw_trip = falls[0] if falls else None
    result = {key: None for key in SUMMARY_FIELDS}
    result.update(protection_metrics_version=VERSION, protection_event_start_input_s=start,
                  protection_display_start_input_s=display,
                  protection_timing_method="voltage_phase_breakpoint_estimate" if start is not None else "nominal_display_fallback",
                  protection_status="disturbance_unconfirmed", protection_reason=fitted["freq_fit_reason"],
                  protection_observed_duration_s=max(0., float(t[-1]) - display),
                  t1_s=1. if start is not None else None, t2_s=None, delta_t_s=None,
                  freq_nominal_start_input_s=nominal, freq_input_sampling_s=dt,
                  freq_breaker_trip_input_s=raw_trip, freq_voltage_start_input_s=start,
                  freq_nominal_delta_t_s=None if raw_trip is None else raw_trip - nominal,
                  **{key: fitted[key] for key in fitted if key.startswith("freq_")})

    def status(name, reason):
        result.update(protection_status=name, protection_reason=reason)

    index = int(np.searchsorted(t, display))
    if br["status"] != "valid":
        status(br["status"], br["reason"])
    elif index == 0 or index >= len(t):
        status("insufficient_observation", "缺少扰动前状态")
    elif sig["breaker"][index - 1] < .5:
        status("already_open", "扰动前 breaker 已断开，不能作为本次响应")
    elif start is not None:
        if raw_trip is not None:
            delay = raw_trip - start
            result.update(protection_trip_input_s=raw_trip, t2_s=1. + delay,
                          delta_t_s=delay, freq_voltage_to_breaker_s=delay)
            status("tripped", "")
        elif (t[-1] < end - max(.02, 2 * dt)
              or np.any(np.diff(t[(t >= start) & (t <= end)]) > 3 * dt)
              or any(start <= x <= end for x in br["incomplete_falls"])):
            status("insufficient_observation", "保护观察窗口不足或存在采样缺口")
        else:
            status("no_trip", "完整观察窗口内未检测到持续 breaker 下降沿")
    display_frequency = fitted["frequency"].copy()
    if br["status"] == "valid":
        half_support = .5 / p.NOMINAL_FREQUENCY_HZ + .01
        for i in np.flatnonzero(np.diff(sig["breaker"]) != 0) + 1:
            display_frequency[np.abs(t - t[i]) <= half_support] = np.nan
    return result, display_frequency
