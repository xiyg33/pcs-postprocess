"""Voltage protection timing on the measured voltage and native breaker clock."""

import numpy as np

from . import processing as p


VERSION = "measured_voltage_breaker_v1"
SUMMARY_FIELDS = [
    "protection_metrics_version", "protection_event_start_input_s",
    "protection_trip_input_s", "protection_display_start_input_s",
    "protection_timing_method", "protection_status", "protection_reason",
    "protection_observed_duration_s",
]
ERROR_STATUSES = {"missing_breaker", "invalid_breaker", "disturbance_unconfirmed",
                  "already_open", "insufficient_observation"}


def measure(t, voltage_pu, breaker, meta, ev, prior):
    t = np.asarray(t, float)
    fitted = p.detect_frt_edge_offset(t, voltage_pu,
                                    dict(meta, V_fault_pu=meta["V_target_pu"]), ev, prior)
    start = None if fitted is None else float(ev["event1_start"] - fitted)
    display_start = float(ev["event1_start"] - prior) if start is None else start
    result = {"protection_metrics_version": VERSION,
              "protection_event_start_input_s": start, "protection_trip_input_s": None,
              "protection_display_start_input_s": display_start,
              "protection_timing_method": "nominal_display_fallback" if start is None else "measured_voltage_edge",
              "protection_status": "disturbance_unconfirmed",
              "protection_reason": "未确认实测电压扰动起点",
              "protection_observed_duration_s": max(0., float(t[-1]) - display_start),
              "t1_s": None if start is None else 1., "t2_s": None, "delta_t_s": None}

    def unavailable(status, reason):
        result.update(protection_status=status, protection_reason=reason)
        return result

    if breaker is None:
        return unavailable("missing_breaker", "缺少 breaker 通道；重新导出 canonical")
    breaker = np.asarray(breaker, float)
    if (breaker.shape != t.shape or not np.isfinite(breaker).all()
            or not np.all(np.isclose(breaker, 0.) | np.isclose(breaker, 1.))):
        return unavailable("invalid_breaker", "breaker 必须是有限的 0/1 状态")
    if start is None:
        return result
    index = int(np.searchsorted(t, start))
    if index == 0 or index >= len(t):
        return unavailable("insufficient_observation", "缺少扰动前 breaker 状态")
    if breaker[index - 1] < .5:
        return unavailable("already_open", "扰动前 breaker 已断开，不能作为本次响应")
    # Require an actual high-to-low transition and ten uninterrupted native
    # samples; do not reinterpret an already-open state or a short tail as trip.
    falls = np.flatnonzero((breaker[:-1] >= .5) & (breaker[1:] < .5)) + 1
    dt = float(np.median(np.diff(t)))
    observation_end = start + float(meta["event_hold_s"])
    for i in falls[falls >= index]:
        if t[i] > observation_end:
            break
        if (i + 10 <= len(t) and np.all(breaker[i:i + 10] < .5)
                and np.all(np.diff(t[i - 1:i + 10]) <= 3 * dt)):
            trip = float(t[i])
            result.update(protection_trip_input_s=trip, protection_status="tripped",
                          protection_reason="", t2_s=1. + trip - start, delta_t_s=trip - start)
            return result
    if t[-1] < observation_end - max(2 * dt, .02):
        return unavailable("insufficient_observation", "录波未覆盖完整保护观察窗口，未确认跳闸")
    if np.any(np.diff(t[index:np.searchsorted(t, observation_end)]) > 3 * dt):
        return unavailable("insufficient_observation", "保护观察窗口存在采样缺口，未确认跳闸")
    return unavailable("no_trip", "完整观察窗口内未检测到持续 breaker 下降沿")
