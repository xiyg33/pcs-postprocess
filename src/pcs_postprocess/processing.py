from dataclasses import dataclass
import csv
import glob
import json
import math
import os

import numpy as np
import scipy.io as sio

# Configured by the CLI before processing a batch.
HIL_P_BASE_W = None
HIL_V_BASE_PEAK_V = None
HIL_I_BASE_PEAK_A = None
REC_START_NOMINAL_S = None
NOMINAL_FREQUENCY_HZ = None
TEST_MODE_CONFIG = {name: {} for name in
                    ("frt", "pwr", "freq_reg", "inertia", "volt_prot", "freq_prot", "scr")}


def configure_project(config):
    global HIL_P_BASE_W, HIL_V_BASE_PEAK_V, HIL_I_BASE_PEAK_A
    global REC_START_NOMINAL_S, NOMINAL_FREQUENCY_HZ
    HIL_P_BASE_W = config.power_base_w
    HIL_V_BASE_PEAK_V = config.voltage_base_peak_v
    HIL_I_BASE_PEAK_A = config.current_base_peak_a
    REC_START_NOMINAL_S = config.record_start_nominal_s
    NOMINAL_FREQUENCY_HZ = config.nominal_frequency_hz
    for mode in ("freq_reg", "inertia", "freq_prot"):
        _PROFILES[mode]["edge_baseline"] = NOMINAL_FREQUENCY_HZ


PWR_STEP_SUMMARY_COLS = [
    "M_p_pct", "T_up_s", "T_p_s", "T_s_s",
    "M_p_reason", "T_up_reason", "T_p_reason", "T_s_reason",
    "step_timing_basis", "step_response_start_s", "T_s_method",
    "metric_sampling_ms", "人工确认Ts",
]

PWR_METRIC_SAMPLE_S = .005
PWR_STEP_REASON_TEXT = {
    "timing_unavailable": "未找到可用于计算的阶跃起点",
    "baseline_unstable": "扰动前功率仍在变化",
    "final_steady_state_unconfirmed": "阶跃末段仍在变化",
    "not_settled_within_observation": "观察时间内未持续稳定，稳定时间无法确定",
    "ninety_percent_not_reached": "响应未达到阶跃幅值的90%",
    "near_zero_measured_step": "实测阶跃幅值过小",
    "response_direction_mismatch": "响应方向与指令不符",
    "missing_samples": "录波存在缺失采样",
    "insufficient_samples": "可用采样点不足",
    "insufficient_observation": "观察时长不足",
    "invalid_time_axis_or_signal_length": "时间轴或信号长度异常",
    "zero_command_step": "指令阶跃幅值为零",
    "five_point_invalid_input": "5点辅助检测输入无效",
    "five_point_insufficient_samples": "5点辅助检测采样不足",
    "five_point_missing_samples": "5点辅助检测遇到缺失采样",
    "five_point_baseline_noise_exceeds_step": "扰动前波动大于阶跃，5点起点无法确认",
    "five_point_edge_unconfirmed": "5点辅助检测未确认阶跃起点",
}


def pwr_five_point_edge(t, y, expected_s, command_delta, smooth_ms=20.0):
    """Estimate a PWR response edge from five points after the expected event.

    This fallback is for the HIL summary only. A five-sample crossing must be
    followed by a persistent directional response; it is never a command time.
    """
    from .step_response import smooth_native

    t = np.asarray(t, dtype=float).reshape(-1)
    y = np.asarray(y, dtype=float).reshape(-1)
    if (t.size != y.size or t.size < 6 or not np.isfinite(expected_s)
            or not np.isfinite(command_delta) or abs(command_delta) < 1e-6
            or not np.all(np.isfinite(t)) or np.any(np.diff(t) <= 0)):
        return None, "five_point_invalid_input"
    dt = float(np.median(np.diff(t)))
    baseline = (t >= expected_s - .5) & (t < expected_s - .05)
    search = np.flatnonzero((t >= expected_s) & (t <= expected_s + EDGE_MAX_DEV_S))
    if np.count_nonzero(baseline) < 3 or search.size < 5:
        return None, "five_point_insufficient_samples"
    if (np.any(~np.isfinite(y[baseline])) or np.any(~np.isfinite(y[search]))
            or np.any(np.diff(t[baseline]) > 3 * dt)
            or np.any(np.diff(t[search]) > 3 * dt)):
        return None, "five_point_missing_samples"

    # Smooth the baseline and response separately, as measure_step does.
    base_y = smooth_native(t[baseline], y[baseline], smooth_ms)
    event = (t >= expected_s) & (t <= expected_s + EDGE_MAX_DEV_S + .05)
    event_t = t[event]
    if np.any(~np.isfinite(y[event])) or np.any(np.diff(event_t) > 3 * dt):
        return None, "five_point_missing_samples"
    event_y = smooth_native(event_t, y[event], smooth_ms)
    base = float(np.median(base_y))
    threshold = max(.002, .05 * abs(command_delta), 3 * float(np.std(base_y)))
    if threshold >= .8 * abs(command_delta):
        return None, "five_point_baseline_noise_exceeds_step"
    directed = np.sign(command_delta) * (event_y - base)
    search_count = int(np.count_nonzero(event_t <= expected_s + EDGE_MAX_DEV_S))
    confirm_count = max(5, int(np.ceil(.05 / dt)))
    for i in range(search_count - 4):
        if not np.all(directed[i:i + 5] >= threshold):
            continue
        # Reject a spike or ripple burst that lasts only five samples.
        later = directed[i:i + confirm_count]
        if (later.size < confirm_count or np.count_nonzero(later >= threshold) < .8 * confirm_count
                or float(np.median(later)) < threshold):
            continue
        return float(event_t[i]), ""
    return None, "five_point_edge_unconfirmed"


def pwr_estimate_unstable_tail(t, y, start, end, measured, smooth_ms=20.0):
    """Retain marked, provisional amplitude metrics when the tail still drifts."""
    from .step_response import smooth_native

    if (measured["reasons"].get("M_p") != "final_steady_state_unconfirmed"
            or not measured["baseline_stable"]):
        return
    initial, final = measured["initial_pu"], measured["final_pu"]
    delta = final - initial
    t = np.asarray(t, float).reshape(-1)
    y = np.asarray(y, float).reshape(-1)
    event = (t >= start) & (t < end)
    tt, yy = t[event], smooth_native(t[event], y[event], smooth_ms)
    directed = np.sign(delta) * (yy - initial)
    measured["values"]["M_p"] = max(0.0, float((np.max(directed) / abs(delta) - 1) * 100))
    hits = np.flatnonzero(directed >= .9 * abs(delta))
    if hits.size:
        i = int(hits[0])
        crossing = tt[i]
        if i > 0 and directed[i] != directed[i - 1]:
            crossing = tt[i - 1] + ((.9 * abs(delta) - directed[i - 1])
                       / (directed[i] - directed[i - 1]) * (tt[i] - tt[i - 1]))
        measured["values"]["T_up"] = float(crossing - start)


def pwr_fourth_band_entry(t, y, start, end, measured, smooth_ms=20.0):
    """Return the fourth outside-to-inside entry into the final-value 5% band.

    This is an estimated HIL display time, not proof that later samples stay
    inside the band. Count entries on the same 20 ms smoothed first-step signal
    used by the other metrics; do not inspect the return step after `end`.
    """
    from .step_response import smooth_native

    initial = measured.get("initial_pu")
    final = measured.get("final_pu")
    if (initial is None or final is None or not np.isfinite(initial + final)
            or abs(final - initial) < 1e-6 or not np.isfinite(start + end)):
        return None
    t = np.asarray(t, float).reshape(-1)
    y = np.asarray(y, float).reshape(-1)
    event = (t >= start) & (t < end)
    tt = t[event]
    if (tt.size < 5 or np.any(~np.isfinite(y[event]))
            or np.any(np.diff(tt) > 3 * float(np.median(np.diff(t))))):
        return None
    yy = smooth_native(tt, y[event], smooth_ms)
    inside = np.abs(yy - final) <= .05 * abs(final - initial)
    # A point already inside at the start is not an outside-to-inside entry.
    entries = np.flatnonzero(inside[1:] & ~inside[:-1]) + 1
    return float(tt[entries[3]] - start) if entries.size >= 4 else None


def pwr_downsampled_signals(bridged, start=None, end=None):
    """20 ms smooth P/Q, then interpolate to a 5 ms normalized time grid.

    Keep raw 5-point edge detection separate. A missing sample or raw time gap
    remains NaN on the reduced grid, never an interpolated measurement.
    """
    from .step_response import smooth_native

    t = np.asarray(bridged["t_s"], dtype=float).reshape(-1)
    if t.size < 2 or not np.all(np.isfinite(t)) or np.any(np.diff(t) <= 0):
        raise ValueError("invalid PWR time axis for 5 ms downsampling")
    count = int(np.floor((t[-1] + 1e-12) / PWR_METRIC_SAMPLE_S)) + 1
    grid = np.arange(count, dtype=float) * PWR_METRIC_SAMPLE_S
    result = {"t_s": grid}
    raw_dt = float(np.median(np.diff(t)))
    gaps = np.flatnonzero(np.diff(t) > 3 * raw_dt)
    segments = ()
    if start is not None and end is not None and np.isfinite(start + end):
        segments = ((t >= start - .5) & (t < start),
                    (t >= start) & (t < end))

    def smooth_runs(raw, mask, output):
        indices = np.flatnonzero(mask & np.isfinite(raw))
        if not indices.size:
            return
        breaks = np.flatnonzero((np.diff(indices) > 1)
                                | (np.diff(t[indices]) > 3 * raw_dt)) + 1
        for run in np.split(indices, breaks):
            output[run] = smooth_native(t[run], raw[run], 20.0)

    for name in ("P_pu", "Q_pu"):
        raw = np.asarray(bridged.get(name, np.full(t.shape, np.nan)), dtype=float).reshape(-1)
        if raw.size != t.size:
            raise ValueError("PWR signal length differs from time axis: " + name)
        smoothed = np.full(t.shape, np.nan)
        smooth_runs(raw, np.ones(t.shape, dtype=bool), smoothed)
        for segment in segments:
            smooth_runs(raw, segment, smoothed)
        reduced = np.interp(grid, t, smoothed, left=np.nan, right=np.nan)
        for i in gaps:
            reduced[(grid > t[i]) & (grid < t[i + 1])] = np.nan
        result[name] = reduced
    return result


def pwr_step_summary(bridged, meta, return_downsampled=False):
    """HIL first-step metrics on 5 ms P/Q; detect the edge on raw samples."""
    from .step_response import measure_step
    code = int(safe_num(meta.get("test_type_code")))
    signal = "Q_pu" if code == 4 else "P_pu"
    prefix = signal[0]
    delta = safe_num(meta.get(prefix + "_step_pu")) - safe_num(meta.get(prefix + "_base_pu"))
    t = np.asarray(bridged["t_s"], dtype=float).reshape(-1)
    y = np.asarray(bridged.get(signal, np.full(t.shape, np.nan)), dtype=float).reshape(-1)
    reliable = (int(bridged.get("command_time_reliable", 0)) == 1
                and str(bridged.get("command_time_source")) == "synchronized"
                and np.isfinite(bridged.get("command_event1_normalized_s", np.nan)))
    edge = int(bridged.get("response_edge_detected", 0)) == 1
    start = (float(bridged["command_event1_normalized_s"]) if reliable else
             float(bridged["response_edge_normalized_s"]) if edge else float("nan"))
    basis = ("synchronized_command" if reliable else
             "response_edge_estimate" if edge else "timing_unavailable")
    fallback_reason = ""
    if not reliable:
        five_point_start, fallback_reason = pwr_five_point_edge(
            t, y, float(bridged["event1_start_s"]), delta)
        if five_point_start is not None:
            start = five_point_start
            basis = "five_point_response_estimate"
    command_start = (float(bridged["command_event1_normalized_s"]) if reliable
                     else float(bridged["event1_start_s"]))
    end = command_start + safe_num(meta.get("step_dur_s"))
    second = safe_num(meta.get("step2_start_s"), float("nan"))
    first = safe_num(meta.get("step_start_s"))
    if np.isfinite(second) and second > first:
        end = min(end, command_start + second - first)
    end = min(end, float(t[-1]))
    downsampled = pwr_downsampled_signals(bridged, start, end)
    metric_t = downsampled["t_s"]
    metric_y = downsampled[signal]
    measured = measure_step(metric_t, metric_y, start, end, delta, smooth_ms=0.0)
    pwr_estimate_unstable_tail(metric_t, metric_y, start, end, measured, smooth_ms=0.0)
    fourth_entry = pwr_fourth_band_entry(metric_t, metric_y, start, end, measured,
                                         smooth_ms=0.0)
    if fourth_entry is not None:
        measured["values"]["T_s"] = fourth_entry
        measured["reasons"]["T_s"] = "fourth_band_entry_estimate"
        settling_method = "fourth_band_entry_estimate"
    elif measured["values"]["T_s"] is not None:
        settling_method = "continuous_5pct"
        measured["reasons"]["T_s"] = "continuous_5pct_without_four_entries"
    else:
        settling_method = "unavailable"
        if (measured.get("initial_pu") is not None
                and measured.get("final_pu") is not None):
            measured["reasons"]["T_s"] = "fewer_than_four_band_entries_and_not_settled"
    columns = {"M_p": "M_p_pct", "T_up": "T_up_s",
               "T_p": "T_p_s", "T_s": "T_s_s"}
    result = {"step_timing_basis": basis,
              "step_response_start_s": start if np.isfinite(start) else None,
              "T_s_method": settling_method,
              "metric_sampling_ms": 5,
              "人工确认Ts": "待确认" if settling_method == "unavailable" else ""}
    for metric, column in columns.items():
        result[column] = measured["values"][metric]
        reasons = []
        reason = measured["reasons"].get(metric, "")
        if reason == "timing_unavailable" and fallback_reason:
            reasons.extend((PWR_STEP_REASON_TEXT[reason],
                            PWR_STEP_REASON_TEXT.get(fallback_reason, fallback_reason)))
        elif reason == "fourth_band_entry_estimate":
            reasons.append("第4次进入目标值±5%范围的估算值，之后可能再次离开")
        elif reason == "continuous_5pct_without_four_entries":
            reasons.append("未出现第4次进入5%范围，沿用持续稳定判据")
        elif reason == "fewer_than_four_band_entries_and_not_settled":
            reasons.append("不足4次进入5%范围，且观察期内未持续稳定")
        elif reason:
            reasons.append(PWR_STEP_REASON_TEXT.get(reason, reason))
        if basis == "five_point_response_estimate":
            reasons.append("使用扰动后连续5点辅助起点估算")
            if not measured["baseline_stable"] and reason != "baseline_unstable":
                reasons.append("扰动前功率仍在变化")
        if (metric == "T_s" and settling_method == "fourth_band_entry_estimate"
                and not measured["tail_stable"]):
            reasons.append("阶跃末段仍在变化，数值仅供参考")
        if (measured["values"][metric] is not None
                and reason == "final_steady_state_unconfirmed"):
            reasons.append("数值仅供参考")
        result[metric + "_reason"] = "；".join(reasons)
    return (result, downsampled) if return_downsampled else result

# =============================================================================
#  汇总判据阈值（基于 2026-07-17 批次 23 case 实测整定）
# =============================================================================

PRE_P_TOL_PU   = 0.05   # 事件前 P 均值允许偏差（控制器稳态误差实测 ~3%）
PRE_PP_TOL_PU  = 0.08   # 事件前 P 峰峰值上限（正常纹波实测 ~0.04 pu）
PRE_Q_TOL_PU   = 0.40   # 事件前 Q 偏置上限（正常 VSG 无功交换实测 ~0.18 pu，
                        #   2014 错误平衡点 ~0.73 pu）
UNSTABLE_PP_PU = 0.50   # 判失稳的 P 峰峰值阈值（2005-2008 实测 ~2.0 pu）

# FRT 阈值（注意：监测点在 PCC，大功率时线路阻抗压降会使 V+ 明显偏离 1.0 pu）
FRT_PRE_VP_TOL_PU  = 0.25  # 事件前 V+ 均值允许偏差 (|mean-1.0|)，PCC 处 ±1.0pu 放电实测 ~0.93
FRT_PRE_VPP_TOL_PU = 0.27  # 事件前 V+ pp 上限（±0.3 pu 正常 ~0.25）
FRT_UNSTABLE_VPP   = 0.35  # 判失稳的 V+ pp 阈值（远超正常纹波，大功率功角摆动）

EDGE_SEARCH_AHEAD_S = 1.0   # 事件沿搜索窗：先验位置之后的长度
EDGE_SEARCH_BACK_S  = 0.3   # 事件沿搜索窗：先验位置之前的长度
EDGE_SUSTAIN_S      = 0.03  # 越过阈值后需保持的时长（防纹波/振荡误检）
EDGE_MAX_DEV_S      = 0.35  # 检测偏移与先验的最大允许偏差（超出视为误检，走兜底）

# 边沿检测内部阈值（从 detect_edge_offset 提取为命名常量）
EDGE_BASELINE_BACK_S   = 1.2   # 基线测量窗：先验位置之前的长度
EDGE_MIN_AMP_FRAC      = 0.03  # 最小阶跃幅值（相对 Pbase）
EDGE_NOISE_TOL_FRT     = 0.25  # FRT/VOLT_PROT 噪声容差（pu 域）
EDGE_NOISE_TOL_PWR     = 0.15  # PWR 噪声容差（W 域比例）
EDGE_BASE_TOL          = 0.20  # 基线偏移容差（W 域比例 / pu 域）
EDGE_DELTA_MIN_FRAC    = 0.04  # 离开基线阈值下限（相对 Pbase）
EDGE_DELTA_NOISE_MULT  = 3.0   # 离开基线阈值噪声倍数
EDGE_THRESH_RATIO_CAP  = 0.90  # 阈值相对阶跃最大比例（超过则放弃）
EDGE_PRE_GAP_S         = 1.2   # 基线窗结束到先验的间隔
EDGE_POST_WINDOW_S     = 0.5   # 事件后统计起始延迟
NORMALIZED_PRE_EVENT_S = 1.0   # 输出统一保留 1 s 扰动前稳态
FREQ_REG_EVENT_END_WINDOW_S = 1.0
# Kept as an alias for compatibility with existing deadband terminology.
FREQ_REG_DEADBAND_END_WINDOW_S = FREQ_REG_EVENT_END_WINDOW_S
FREQ_REG_DEADBAND_P_LIMIT_PU = 0.05
FREQ_REG_LIMIT_DEVIATIONS_HZ = (0.25, 0.35, 0.5)
FREQ_REG_FREQ_MATCH_TOL_HZ = 1e-9
CURRENT_MIN_VOLTAGE_PU = 0.05
CURRENT_FILTER_WINDOW_S = 0.020


def supports_power_current(meta, profile_name):
    """Add response data for ordinary FRT, regulation and inertia cases."""
    if profile_name not in ("frt", "freq_reg", "inertia"):
        return False
    if profile_name == "inertia":
        return True
    test_type = str(meta.get("test_type") or "").strip().upper()
    excluded = ("LVRT_BND", "HVRT_BND", "FREQ_REG_DEADBAND")
    if test_type:
        return test_type not in excluded
    name = str(meta.get("case_name") or "").upper()
    return not any(name.startswith(prefix) for prefix in excluded)


def calculate_instantaneous_current(vabc, iabc):
    """Return voltage-oriented Ip/Iq in peak A and a per-sample valid mask.

    Inputs are matching (N, 3) arrays. No sequence extraction or smoothing;
    positive Iq follows q = 1.5 * (v_beta*i_alpha - v_alpha*i_beta).
    """
    voltage = np.asarray(vabc, dtype=float)
    current = np.asarray(iabc, dtype=float)
    if voltage.ndim != 2 or voltage.shape[1] != 3 or current.shape != voltage.shape:
        raise ValueError("vabc and iabc must have matching (N, 3) shapes")
    finite = np.isfinite(voltage).all(axis=1) & np.isfinite(current).all(axis=1)
    ip = np.full(voltage.shape[0], np.nan)
    iq = np.full(voltage.shape[0], np.nan)
    valid = np.zeros(voltage.shape[0], dtype=bool)
    indices = np.flatnonzero(finite)
    v = voltage[indices]
    i = current[indices]
    va = (2 * v[:, 0] - v[:, 1] - v[:, 2]) / 3
    vb = (v[:, 1] - v[:, 2]) / np.sqrt(3)
    ia = (2 * i[:, 0] - i[:, 1] - i[:, 2]) / 3
    ib = (i[:, 1] - i[:, 2]) / np.sqrt(3)
    amplitude = np.hypot(va, vb)
    keep = amplitude / HIL_V_BASE_PEAK_V >= CURRENT_MIN_VOLTAGE_PU
    kept_indices = indices[keep]
    # Normalize the voltage first, avoiding division by small/zero voltages.
    ua, ub = va[keep] / amplitude[keep], vb[keep] / amplitude[keep]
    ip[kept_indices] = ua * ia[keep] + ub * ib[keep]
    iq[kept_indices] = ub * ia[keep] - ua * ib[keep]
    valid[kept_indices] = True
    return ip, iq, valid


def filter_current_for_display(time, values, valid):
    """20 ms centered mean, normalized at edges and never bridging gaps.

    Recorder segment boundaries may repeat a timestamp; the sample spacing is
    estimated from positive intervals while retaining both recorded samples.
    """
    time = np.asarray(time, dtype=float)
    values = np.asarray(values, dtype=float)
    valid = np.asarray(valid, dtype=bool)
    if time.ndim != 1 or values.shape != time.shape or valid.shape != time.shape:
        raise ValueError("time, values and valid must be matching 1-D arrays")
    result = np.full(values.shape, np.nan)
    if not time.size:
        return result
    intervals = np.diff(time)
    if (not np.isfinite(time).all() or np.any(intervals < 0)
            or (time.size > 1 and not np.any(intervals > 0))):
        raise ValueError("time must be finite and nondecreasing")
    dt = (float(np.median(intervals[intervals > 0])) if time.size > 1
          else CURRENT_FILTER_WINDOW_S)
    target = CURRENT_FILTER_WINDOW_S / dt
    lower = max(1, 2 * int(math.floor((target - 1) / 2)) + 1)
    window = min((lower, lower + 2), key=lambda n: abs(n - target))
    half = window // 2
    usable = valid & np.isfinite(values)
    edges = np.diff(np.r_[False, usable, False].astype(int))
    for start, stop in zip(np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)):
        segment = values[start:stop]
        indices = np.arange(segment.size)
        lo = np.maximum(0, indices - half)
        hi = np.minimum(segment.size, indices + half + 1)
        sums = np.r_[0., np.cumsum(segment)]
        result[start:stop] = (sums[hi] - sums[lo]) / (hi - lo)
    return result


# =============================================================================
#  测试类型 Profile 注册表
#
#  每种测试类型的行为差异集中在此处，新增工况只需添加一个条目。
#  各函数通过 _get_profile(code) 查找 profile，消除 is_frt 硬编码分支。
# =============================================================================

# FRT 的 test_type 字符串 → code 映射（meta 中可能没有 test_type_code）
_STR_TO_CODE = {
    "LVRT_SINGLE": 0, "HVRT_SINGLE": 1, "LVRT_CONTINUOUS": 2,
    "FREQ_REG": 6, "FREQ_REG_DEADBAND": 6,
}

_PROFILES = {
    # =========================================================================
    # FRT: Fault Ride-Through (codes 0-2)
    # =========================================================================
    "frt": {
        "codes": [0, 1, 2],
        # ---- 事件时间提取 ----
        "event1_field":     "fault_start_s",
        "event1_dur_field": "fault_duration_s",
        "event2_field":     None,
        "event2_dur_field": None,
        # ---- 边沿检测 ----
        "edge_signal":      "Vp_pu",          # 用 V+ 检测电压跌落/升高
        "edge_baseline":    1.0,              # V+ 额定 = 1.0 pu
        "edge_step_field":  "V_fault_pu",
        "edge_is_Vp_scale": True,             # Vp 信号 ×Pbase 走统一 W 域阈值
        "edge_max_dev_s":   0.50,             # FRT 录波起点漂移略大
        # Drawing layouts live in plot_common.py.
        # ---- 事件前稳态检查 ----
        "pre_check": "frt",
        # ---- 汇总 CSV 额外列 ----
        "summary_cols": [
            "V_fault_pu", "fault_type", "P_ref_pu",
            "pre_Vp_mean_pu", "pre_Vp_pp_pu",
            "event_Vp_mean_pu", "post_Vp_mean_pu", "pre_P_mean_pu",
            "pre_Q_mean_pu",
        ],
    },

    # =========================================================================
    # PWR: Power Step (codes 3-5)
    # =========================================================================
    "pwr": {
        "codes": [3, 4, 5],
        "event1_field":     "step_start_s",
        "event1_dur_field": "step_dur_s",
        "event2_field":     "step2_start_s",
        "event2_dur_field": "step2_dur_s",
        # code 4 = Q 阶跃，其他 = P 阶跃
        "edge_signal":      "auto",
        "edge_baseline_field": "auto",
        "edge_step_field":  "auto",
        "edge_is_Vp_scale": False,
        "edge_max_dev_s":   0.35,
        "pre_check": "pwr",
        "summary_cols": [
            "pre_P_mean_pu", "pre_P_pp_pu", "pre_Q_mean_pu", "step_P_mean_pu",
        ],
    },

    # =========================================================================
    # FREQ_REG: Frequency Regulation (codes 6-7)
    # =========================================================================
    "freq_reg": {
        "codes": [6, 7],
        "event1_field":     "event_start_s",
        "event1_dur_field": "event_hold_s",
        "event2_field":     None,
        "event2_dur_field": None,
        "edge_signal":      "f",
        "edge_baseline":    50.0,
        "edge_step_field":  "freq_deviation_hz",
        "edge_is_freq":     True,
        "edge_max_dev_s":   0.35,
        "pre_check": "freq",
        "summary_cols": [
            "P_ref_pu", "freq_deviation_hz",
            "pre_P_mean_pu", "pre_P_pp_pu", "pre_f_mean_hz",
            "event_end_P_mean_pu", "P_change_signed_pu",
            "freq_reg_upward_limit_pu", "freq_reg_downward_limit_pu",
            "P_change_abs_pu", "P_change_5pct_pass",
        ],
    },

    # =========================================================================
    # INERTIA: Inertia Response (codes 8-9) — 同 freq_reg 布局
    # =========================================================================
    "inertia": {
        "codes": [8, 9],
        "event1_field":     "event_start_s",
        "event1_dur_field": "event_hold_s",
        "event2_field":     None,
        "event2_dur_field": None,
        # INERTIA: event1 时长 = ramp_dur + event_hold_s（含频率斜坡）
        "event1_has_ramp":          True,
        "event1_ramp_dev_field":    "freq_deviation_hz",
        "event1_ramp_rate_field":   "freq_ramp_rate_hz_per_s",
        "edge_signal":      "f",
        "edge_baseline":    50.0,
        "edge_step_field":  "freq_deviation_hz",
        "edge_is_freq":     True,
        "edge_max_dev_s":   0.50,             # 慢 ramp 检测延迟较大
        "pre_check": "freq",
        "summary_cols": [
            "P_ref_pu", "freq_deviation_hz", "freq_ramp_rate_hz_per_s",
            "pre_P_mean_pu", "pre_P_pp_pu", "pre_f_mean_hz",
        ],
    },

    # =========================================================================
    # VOLT_PROT: Voltage Protection (codes 10-11)
    # =========================================================================
    "volt_prot": {
        "codes": [10, 11],
        "event1_field":     "event_start_s",
        "event1_dur_field": "event_hold_s",
        "event2_field":     None,
        "event2_dur_field": None,
        "edge_signal":      "Vp_pu",
        "edge_baseline":    1.0,
        "edge_step_field":  "V_target_pu",
        "edge_is_Vp_scale": True,
        "edge_max_dev_s":   0.35,
        "pre_check": "frt",                # V+ 基线检查（同 FRT）
        "summary_cols": [
            "P_ref_pu", "V_target_pu",
            "pre_Vp_mean_pu", "pre_P_mean_pu", "pre_P_pp_pu",
            "t1_s", "t2_s", "delta_t_s",
        ],
    },

    # =========================================================================
    # FREQ_PROT: Frequency Protection (codes 12-13)
    # =========================================================================
    "freq_prot": {
        "codes": [12, 13],
        "event1_field":     "event_start_s",
        "event1_dur_field": "event_hold_s",
        "event2_field":     None,
        "event2_dur_field": None,
        "edge_signal":      "f",
        "edge_baseline":    50.0,
        "edge_step_field":  "f_target_hz",
        "edge_is_freq":     True,
        "edge_max_dev_s":   0.35,
        "pre_check": "freq",
        "summary_cols": [
            "P_ref_pu", "f_target_hz",
            "pre_f_mean_hz", "pre_P_mean_pu", "pre_P_pp_pu",
            "t1_s", "t2_s", "delta_t_s",
        ],
    },

    # =========================================================================
    # SCR: SCR Adaptability (codes 14-15) — placeholder
    # =========================================================================
    "scr": {
        "codes": [14, 15],
        "event1_field":     "disturbance_start_s",
        "event1_dur_field": "disturbance_duration_s",
        "event2_field":     None,
        "event2_dur_field": None,
        "edge_signal":      "P",
        "edge_baseline_field": "P_ref_pu",
        "edge_step_field":  None,           # SCR 无明确阶跃目标
        "edge_is_Vp_scale": False,
        "edge_max_dev_s":   0.50,
        "pre_check": "pwr",
        "summary_cols": [
            "P_ref_pu", "SCR_value",
            "pre_P_mean_pu", "pre_P_pp_pu",
        ],
    },
}


# =============================================================================
#  基础工具
# =============================================================================

def safe_num(val, default=0.0):
    """meta 中的 NaN/None/缺失 → default"""
    if val is None:
        return default
    try:
        f = float(val)
        if math.isnan(f):
            return default
        return f
    except (TypeError, ValueError):
        return default


def _fmt(v, digits=4):
    """安全格式化统计值：round 或空字符串（NaN 时）"""
    if isinstance(v, float) and math.isnan(v):
        return ""
    try:
        return round(float(v), digits)
    except (TypeError, ValueError):
        return ""


def is_freq_reg_deadband(meta):
    """Return True for new deadband metadata and already-recorded legacy data."""
    test_type = str(meta.get("test_type", "")).strip().upper()
    case_name = str(meta.get("case_name", "")).strip().upper()
    return (
        test_type == "FREQ_REG_DEADBAND"
        or case_name.startswith("FREQ_REG_DEADBAND_")
    )


def _power_change_5pct_pass(pre_P_mean_pu, event_end_P_mean_pu):
    """Return 1 when absolute power change is within 0.05 pu, else 0/NaN."""
    try:
        pre_value = float(pre_P_mean_pu)
        event_value = float(event_end_P_mean_pu)
    except (TypeError, ValueError):
        return float("nan")
    if not (math.isfinite(pre_value) and math.isfinite(event_value)):
        return float("nan")
    delta = abs(event_value - pre_value)
    return int(delta <= FREQ_REG_DEADBAND_P_LIMIT_PU + 1e-12)


def _is_freq_reg_limit_case(meta):
    """Return True for normal FREQ_REG cases at the requested deviations."""
    if is_freq_reg_deadband(meta):
        return False

    test_type = str(meta.get("test_type", "")).strip().upper()
    case_name = str(meta.get("case_name", "")).strip().upper()
    if (
        test_type not in ("FREQ_REG", "FREQ_REG_DEADBAND")
        and not case_name.startswith("FREQ_REG_")
    ):
        return False

    try:
        deviation = abs(float(meta.get("freq_deviation_hz")))
    except (TypeError, ValueError):
        return False
    if not math.isfinite(deviation):
        return False

    return any(
        math.isclose(deviation, target, rel_tol=0.0,
                     abs_tol=FREQ_REG_FREQ_MATCH_TOL_HZ)
        for target in FREQ_REG_LIMIT_DEVIATIONS_HZ
    )


def _event_end_power_change(stats, t_sim, P_pu, ev):
    """Return (event-end mean, signed event-end minus pre-event delta)."""
    ev1 = ev["event1_start"]
    ev1_end = ev["event1_end"]
    end_window_start = max(ev1, ev1_end - FREQ_REG_EVENT_END_WINDOW_S)
    end_m = (t_sim >= end_window_start) & (t_sim <= ev1_end)
    event_end_P_mean, _ = _win_stats(P_pu, end_m)
    pre_P_mean = stats.get("pre_P_mean_pu", float("nan"))
    try:
        pre_value = float(pre_P_mean)
    except (TypeError, ValueError):
        pre_value = float("nan")
    if not (math.isfinite(pre_value) and math.isfinite(event_end_P_mean)):
        P_change_signed = float("nan")
    else:
        P_change_signed = event_end_P_mean - pre_value
    return event_end_P_mean, P_change_signed


def add_freq_reg_power_stats(stats, meta, t_sim, P_pu, ev):
    """Add event-end power metrics for deadband and limit-test cases."""
    deadband_case = is_freq_reg_deadband(meta)
    limit_case = _is_freq_reg_limit_case(meta)
    if not (deadband_case or limit_case):
        return stats

    event_end_P_mean, P_change_signed = _event_end_power_change(
        stats, t_sim, P_pu, ev,
    )
    if math.isfinite(P_change_signed):
        P_change_abs = abs(P_change_signed)
    else:
        P_change_abs = float("nan")

    stats.update({
        "event_end_P_mean_pu": event_end_P_mean,
        "P_change_abs_pu": P_change_abs,
    })

    if deadband_case:
        stats["P_change_5pct_pass"] = _power_change_5pct_pass(
            stats.get("pre_P_mean_pu", float("nan")), event_end_P_mean,
        )
    else:
        stats["P_change_signed_pu"] = P_change_signed
        stats["freq_reg_upward_limit_pu"] = (
            P_change_signed if math.isfinite(P_change_signed) and P_change_signed > 0.0
            else float("nan")
        )
        stats["freq_reg_downward_limit_pu"] = (
            P_change_signed if math.isfinite(P_change_signed) and P_change_signed < 0.0
            else float("nan")
        )
    return stats


def add_freq_reg_deadband_stats(stats, meta, t_sim, P_pu, ev):
    """Backward-compatible wrapper for deadband-only callers."""
    if not is_freq_reg_deadband(meta):
        return stats
    return add_freq_reg_power_stats(stats, meta, t_sim, P_pu, ev)


def _get_profile(code):
    """根据 test_type_code 查找对应 profile name 和 dict。

    可用于 meta 中的 test_type_code（数值）或 test_type（字符串）。
    """
    if isinstance(code, str):
        code = _STR_TO_CODE.get(code, 0)
    else:
        code = int(safe_num(code, 3))

    for name, prof in _PROFILES.items():
        if code in prof["codes"]:
            return name, prof

    # fallback: 按 code 区间
    if code <= 2:
        return "frt", _PROFILES["frt"]
    elif code <= 5:
        return "pwr", _PROFILES["pwr"]
    elif code <= 7:
        return "freq_reg", _PROFILES["freq_reg"]
    elif code <= 9:
        return "inertia", _PROFILES["inertia"]
    elif code <= 11:
        return "volt_prot", _PROFILES["volt_prot"]
    elif code <= 13:
        return "freq_prot", _PROFILES["freq_prot"]
    else:
        return "scr", _PROFILES["scr"]


def _win_stats(x, mask):
    """计算信号 x 在布尔 mask 窗口内的 (均值, 峰峰值)。"""
    if x is None or not np.any(mask):
        return float("nan"), float("nan")
    seg = x[mask]
    return float(seg.mean()), float(seg.max() - seg.min())


def load_case_file(mat_path):
    """读取一个 hil_raw .mat，返回 (t_rel, sig dict)。

    sig 键：Va Vb Vc Ia Ib Ic P Q f Vp_pu（缺失键报 KeyError 并列出 header）。
    """
    d = sio.loadmat(mat_path)
    data = d["data"]
    header = d["header"]

    names = []
    for row in header:
        cell = row[0]
        names.append(str(cell[0]) if getattr(cell, "size", 0) else "")

    def _find_col(pred, label, required=True):
        for i, n in enumerate(names):
            if pred(n):
                return i
        if required:
            raise KeyError("signal '%s' not found in header: %s" % (label, names))
        return None

    col = {
        "t":  _find_col(lambda n: n in ("t_s", "timestamps"), "t_s"),
        "Va": _find_col(lambda n: n == "Va_V" or "Vabc" in n and "[0]" in n, "Va_V"),
        "Vb": _find_col(lambda n: n == "Vb_V" or "Vabc" in n and "[1]" in n, "Vb_V"),
        "Vc": _find_col(lambda n: n == "Vc_V" or "Vabc" in n and "[2]" in n, "Vc_V"),
        "Ia": _find_col(lambda n: n == "Ia_A" or "Iabc" in n and "[0]" in n, "Ia_A"),
        "Ib": _find_col(lambda n: n == "Ib_A" or "Iabc" in n and "[1]" in n, "Ib_A"),
        "Ic": _find_col(lambda n: n == "Ic_A" or "Iabc" in n and "[2]" in n, "Ic_A"),
        "P":  _find_col(lambda n: n == "P_W" or "P_cal" in n, "P_W"),
        "Q":  _find_col(lambda n: n == "Q_var" or "Q_cal" in n, "Q_var", required=False),
        "f":  _find_col(lambda n: n == "f_hz" or "PLL" in n, "f_hz", required=False),
        "breaker": _find_col(lambda n: "breaker" in n, "breaker", required=False),
        "Vp_logged": _find_col(lambda n: n in ("Vpos_pu", "SIM/Up_pu"), "Vpos_pu", required=False),
    }

    t_raw = data[:, col["t"]]
    t_rel = t_raw - t_raw[0]
    sig = {k: (data[:, c] if c is not None else None) for k, c in col.items() if k != "t"}

    # 计算正序电压幅值 V+ (pu)
    alpha = np.exp(2j * np.pi / 3.0)
    Va_c = sig["Va"].astype(np.complex128)
    Vb_c = sig["Vb"].astype(np.complex128)
    Vc_c = sig["Vc"].astype(np.complex128)
    Vp = 2.0 * (Va_c + alpha * Vb_c + alpha * alpha * Vc_c) / 3.0
    sig["Vp_pu"] = (sig["Vp_logged"] if sig["Vp_logged"] is not None
                    else np.abs(Vp) / HIL_V_BASE_PEAK_V)

    return t_rel, sig


# =============================================================================
#  保护测试专用：扰动 / 脱网时间检测
# =============================================================================

def detect_disturbance_time(t, signal):
    """从信号突变中检测扰动开始时刻 t1。

    建立稳态基线，自适应阈值，找首次持续越限点。

    Parameters
    ----------
    t : ndarray
        时间向量。
    signal : ndarray
        信号值（如 f_Hz 或 Vp_pu）。

    Returns
    -------
    float
        扰动开始时刻，未找到则返回 NaN。
    """
    n = len(t)
    if n < 50:
        return np.nan

    t_start = t[0] + 0.1
    t_end = t_start + 0.5
    idx_bl_start = np.searchsorted(t, t_start)
    idx_bl_end = np.searchsorted(t, t_end)
    if idx_bl_end <= idx_bl_start:
        idx_bl_start = 0
        idx_bl_end = min(n, 5000)

    sig_bl = signal[idx_bl_start:idx_bl_end]
    ref = np.median(sig_bl)
    noise = np.std(sig_bl)
    threshold = max(5.0 * noise, 0.02)

    sig_tail = signal[idx_bl_end:]
    t_tail = t[idx_bl_end:]
    dev = np.abs(sig_tail - ref)

    min_hold = 20
    run = 0
    for i in range(len(dev)):
        if dev[i] > threshold:
            run += 1
            if run >= min_hold:
                return t_tail[i - run + 1]
        else:
            run = 0

    return np.nan


def detect_trip_time(t_brk, breaker, t1):
    """从 breaker_flag 下降沿检测脱网时刻 t2。

    从 t1 开始向后搜索连续 10 个样本低于 0.5 的第一个位置。

    Parameters
    ----------
    t_brk : ndarray
        断路器信号的时间向量。
    breaker : ndarray or None
        断路器信号（1=并网，0=脱网）。
    t1 : float
        扰动开始时刻（可为 NaN）。

    Returns
    -------
    float
        脱网时刻，未找到则返回 NaN。
    """
    if breaker is None:
        return np.nan

    if np.isnan(t1):
        idx_start = 0
    else:
        idx_start = np.searchsorted(t_brk, t1)
    if idx_start >= len(breaker):
        return np.nan

    brk_seg = breaker[idx_start:]
    t_seg = t_brk[idx_start:]

    if len(brk_seg) < 2:
        return np.nan

    if brk_seg[0] < 0.5:
        return t_seg[0]

    min_hold = 10
    for i in range(len(brk_seg) - min_hold + 1):
        if np.all(brk_seg[i:i + min_hold] < 0.5):
            return t_seg[i]

    for tail in range(min_hold, 1, -1):
        if len(brk_seg) >= tail and np.all(brk_seg[-tail:] < 0.5):
            return t_seg[-tail]

    return np.nan


# =============================================================================
#  事件时间提取（Profile 驱动）
# =============================================================================

def event_times_from_meta(meta):
    """由 meta 计算实效事件时刻（仿真时间轴）。

    使用 profile 的 event1_field / event1_dur_field / event2_field，
    不再硬编码 is_frt 分支。

    返回 dict: profile_name / code / event1_start / event1_end /
              event2_start / event2_end / last_event_end /
              record_pre / record_post / sim_time_eff
    """
    shift = safe_num(meta.get("hil_settle_shift_s"), 0.0)

    pname = meta["mode"]
    prof = _PROFILES[pname]
    tc = meta.get("test_type_code", meta.get("test_type"))
    code = (_STR_TO_CODE.get(tc, prof["codes"][0]) if isinstance(tc, str)
            else int(safe_num(tc, prof["codes"][0])))

    # event1：优先用 hil_event1_start_eff_s（settle 后的实效值），否则从 CSV 字段 + shift
    ev1 = safe_num(meta.get("hil_event1_start_eff_s"),
                   safe_num(meta.get(prof["event1_field"])) + shift)
    dur = safe_num(meta.get(prof["event1_dur_field"]))
    # INERTIA 含频率斜坡，需额外加上 ramp_dur
    if prof.get("event1_has_ramp"):
        dev = safe_num(meta.get(prof.get("event1_ramp_dev_field", "")))
        rate = safe_num(meta.get(prof.get("event1_ramp_rate_field", "")))
        if abs(rate) > 1e-12:
            ramp_dur = abs(dev / rate)
        else:
            ramp_dur = 0.0
        dur = ramp_dur + dur
    ev1_end = ev1 + dur
    last_end = ev1_end

    # event2（可选，仅 PWR DOUBLE 等有第二个事件的类型）
    ev2 = ev2_end = None
    ev2_field = prof.get("event2_field")
    if ev2_field:
        ev2_raw = safe_num(meta.get(ev2_field))
        if ev2_raw + shift > ev1:  # compare on the same effective time axis
            ev2 = ev2_raw + shift
            ev2_end = ev2 + safe_num(meta.get(prof.get("event2_dur_field", "")))
            last_end = ev2_end

    return {
        "profile_name": pname,
        "code": code,
        "event1_start": ev1,
        "event1_end": ev1_end,
        "event2_start": ev2,
        "event2_end": ev2_end,
        "last_event_end": last_end,
        "record_pre": safe_num(meta.get("record_pre_s"), 2.0),
        "record_post": safe_num(meta.get("record_post_s"), 2.0),
        "sim_time_eff": safe_num(meta.get("hil_sim_time_eff_s"),
                                 safe_num(meta.get("sim_time_s")) + shift),
    }


def normalize_case_time(t_sim, ev, disturbance_time=None,
                        pre_event_s=NORMALIZED_PRE_EVENT_S):
    """将时间轴统一为从 0 s 开始、首次扰动位于 pre_event_s。"""
    if disturbance_time is None or not np.isfinite(disturbance_time):
        disturbance_time = ev["event1_start"]

    origin = float(disturbance_time) - float(pre_event_s)
    t_norm = np.asarray(t_sim, dtype=float) - origin
    ev_norm = dict(ev)
    event1_start = float(ev["event1_start"])

    # 事件时序以首次扰动为锚点，持续时间及第二事件间隔保持不变。
    for key in ("event1_start", "event1_end", "event2_start",
                "event2_end", "last_event_end"):
        value = ev.get(key)
        if value is not None:
            ev_norm[key] = float(pre_event_s) + float(value) - event1_start

    w1 = ev_norm["last_event_end"] + ev_norm["record_post"]
    wm = (t_norm >= 0.0) & (t_norm <= w1)
    if not np.any(wm):
        wm = np.ones_like(t_norm, dtype=bool)
    else:
        # 裁剪边界通常落在两个采样点之间；将首个保留点规范为精确 0 s。
        first = int(np.flatnonzero(wm)[0])
        if len(t_norm) > 1:
            sample_step = float(np.median(np.diff(t_norm)))
            if 0.0 <= t_norm[first] <= 2.0 * sample_step:
                t_norm[first] = 0.0
    return t_norm, ev_norm, wm


# =============================================================================
#  对时：事件沿检测（Profile 驱动）
# =============================================================================

def detect_frequency_edge_offset(t_rel, sig, meta, ev, offset_prior):
    """Detect the actual FREQ_REG/INERTIA frequency event edge.

    FREQ_REG uses a sustained departure from the measured 50 Hz baseline.
    INERTIA uses a piecewise-linear ramp fit because its initial deviation is
    smaller than the recorder noise and cannot be found with a step threshold.
    """
    f = sig.get("f")
    if f is None:
        return None

    t_rel = np.asarray(t_rel, dtype=float)
    f = np.asarray(f, dtype=float)
    if t_rel.size < 3 or f.size != t_rel.size:
        return None

    te0 = ev["event1_start"] - offset_prior
    base_m = (t_rel >= te0 - EDGE_BASELINE_BACK_S) & \
             (t_rel <= te0 - EDGE_SEARCH_BACK_S)
    if not np.any(base_m):
        return None

    base_meas = float(np.median(f[base_m]))
    noise_half = float(np.max(f[base_m]) - np.min(f[base_m])) / 2.0
    base = NOMINAL_FREQUENCY_HZ
    if abs(base_meas - base) > 0.20 * base or noise_half > 0.10 * base:
        return None

    dev = safe_num(meta.get("freq_deviation_hz"))
    rate = safe_num(meta.get("freq_ramp_rate_hz_per_s"))
    max_dev = ev.get("edge_max_dev_s", EDGE_MAX_DEV_S)

    # INERTIA: fit the known frequency ramp over the complete event window.
    # This estimates the ramp onset instead of waiting for a large deviation.
    if abs(rate) > 1e-12:
        ramp_dur = abs(dev / rate)
        if ramp_dur <= 0.0:
            return None

        fit_m = (t_rel >= te0 - max_dev) & \
                (t_rel <= min(te0 + ramp_dur + max_dev, t_rel[-1]))
        if np.count_nonzero(fit_m) < 20:
            return None

        dt = float(np.median(np.diff(t_rel)))
        if not np.isfinite(dt) or dt <= 0.0:
            return None

        # Downsample the fit only for speed; preserve the signal shape.
        sample_step = max(1, int(round(0.01 / dt)))
        tt = t_rel[fit_m][::sample_step]
        xx = f[fit_m][::sample_step]
        candidates = np.arange(te0 - max_dev, te0 + max_dev + 0.002, 0.002)
        scores = []
        for candidate in candidates:
            ramp_progress = np.clip(tt - candidate, 0.0, ramp_dur)
            model = base_meas + rate * ramp_progress
            residual = xx - model
            scores.append(float(np.mean(residual * residual)))

        if not scores:
            return None
        edge_time = float(candidates[int(np.argmin(scores))])
        offset = ev["event1_start"] - edge_time
        return offset if abs(offset - offset_prior) <= max_dev else None

    # FREQ_REG: detect the first sustained departure from the measured
    # baseline.  Thresholds scale with the requested step and noise instead
    # of using the 50 Hz value as an absolute deviation threshold.
    if abs(dev) < 0.03:
        return None
    delta = max(4.0 * noise_half, 0.05 * abs(dev), 0.01)
    if delta >= 0.80 * abs(dev):
        return None

    search_m = (t_rel >= te0 - EDGE_SEARCH_BACK_S) & \
               (t_rel <= te0 + EDGE_SEARCH_AHEAD_S)
    if not np.any(search_m):
        return None

    tt = t_rel[search_m]
    xx = f[search_m]
    crossed = (xx - base_meas) * np.sign(dev) > delta
    dt = float(np.median(np.diff(tt)))
    if not np.isfinite(dt) or dt <= 0.0:
        return None
    n_sustain = max(1, int(round(EDGE_SUSTAIN_S / dt)))

    for idx in np.flatnonzero(crossed):
        sustained = crossed[idx:idx + n_sustain]
        if sustained.size == n_sustain and np.all(sustained):
            offset = ev["event1_start"] - float(tt[idx])
            return offset if abs(offset - offset_prior) <= max_dev else None

    return None

def detect_pwr_edge_offset(t_rel, sig, meta, ev, offset_prior):
    """P/Q response estimate on 20 ms filtered data, shared with SIM.

    This is not a command timestamp. Noise sets the threshold; indistinguishable
    steps return None rather than a nominal timestamp presented as detection.
    """
    name = "Q" if ev["code"] == 4 else "P"
    raw = sig.get(name)
    if raw is None:
        return None
    t = np.asarray(t_rel, dtype=float).reshape(-1)
    x = np.asarray(raw, dtype=float).reshape(-1) / HIL_P_BASE_W
    if len(t) < 3 or x.size != t.size or np.any(np.diff(t) <= 0):
        return None
    dt = float(np.median(np.diff(t)))
    n = min(len(x), max(1, int(round(.020 / dt))))
    if n % 2 == 0:
        n = max(1, n - 1)
    x = np.convolve(np.pad(x, (n // 2, n // 2), mode="edge"), np.ones(n) / n, mode="valid")
    expected = ev["event1_start"] - offset_prior
    baseline = (t >= expected - .5) & (t < expected - .05)
    if np.count_nonzero(baseline) < 3 or not np.all(np.isfinite(x[baseline])):
        return None
    base = float(np.median(x[baseline]))
    amp = safe_num(meta.get(name + "_step_pu")) - safe_num(meta.get(name + "_base_pu"))
    threshold = max(.002, .05 * abs(amp), 3 * float(np.std(x[baseline])))
    if abs(amp) < 1e-6 or threshold >= .8 * abs(amp):
        return None
    search = ((t >= expected - EDGE_SEARCH_BACK_S)
              & (t <= expected + EDGE_MAX_DEV_S))
    crossed = search & (np.sign(amp) * (x - base) >= threshold)
    sustain = max(2, int(np.ceil(.020 / dt)))
    if sustain > len(t):
        return None
    counts = np.convolve(crossed.astype(float), np.ones(sustain), mode="valid")
    hits = np.flatnonzero(counts >= sustain)
    return ev["event1_start"] - float(t[hits[0]]) if hits.size else None


def detect_frt_edge_offset(t_rel, voltage_pu, meta, ev, offset_prior):
    """Locate an ordinary FRT voltage step after removing 100 Hz ripple."""
    t = np.asarray(t_rel, dtype=float)
    voltage = np.asarray(voltage_pu, dtype=float)
    if t.size < 3 or voltage.shape != t.shape or not np.isfinite(voltage).all():
        return None
    dt = float(np.median(np.diff(t)))
    if not np.isfinite(dt) or dt <= 0:
        return None
    window = max(1, int(round(0.020 / dt)))
    if window % 2 == 0:
        window += 1
    window = min(window, t.size if t.size % 2 else t.size - 1)
    half = window // 2
    smooth = np.convolve(np.pad(voltage, (half, half), mode="edge"),
                         np.ones(window) / window, mode="valid")
    expected = ev["event1_start"] - offset_prior
    baseline_mask = ((t >= expected - EDGE_BASELINE_BACK_S)
                     & (t <= expected - EDGE_SEARCH_BACK_S))
    if np.count_nonzero(baseline_mask) < 3:
        return None
    baseline_values = smooth[baseline_mask]
    baseline = float(np.median(baseline_values))
    noise = float(np.ptp(baseline_values)) / 2
    amplitude = safe_num(meta.get("V_fault_pu")) - 1.0
    if (not np.isfinite(amplitude) or abs(amplitude) < EDGE_MIN_AMP_FRAC
            or abs(baseline - 1.0) > EDGE_BASE_TOL or noise > EDGE_NOISE_TOL_FRT):
        return None
    threshold = max(0.020, 0.15 * abs(amplitude), 3 * noise)
    if threshold >= EDGE_THRESH_RATIO_CAP * abs(amplitude):
        return None
    search = ((t >= expected - EDGE_SEARCH_BACK_S)
              & (t <= expected + EDGE_SEARCH_AHEAD_S))
    crossed = search & (np.sign(amplitude) * (smooth - baseline) > threshold)
    hold = max(1, int(math.ceil(EDGE_SUSTAIN_S / dt)))
    if hold > t.size:
        return None
    sustained = np.convolve(crossed.astype(int), np.ones(hold, dtype=int), mode="valid")
    candidates = np.flatnonzero(sustained == hold)
    if not candidates.size:
        return None
    offset = ev["event1_start"] - float(t[candidates[0]])
    if abs(offset - offset_prior) > _PROFILES["frt"]["edge_max_dev_s"]:
        return None
    return offset


def detect_frt_recovery_time(t_rel, voltage_pu, meta, ev, offset):
    """Find the measured fault clearing edge near the scheduled end."""
    t = np.asarray(t_rel, dtype=float)
    voltage = np.asarray(voltage_pu, dtype=float)
    if t.size < 3 or voltage.shape != t.shape or not np.isfinite(voltage).all():
        return None
    intervals = np.diff(t)
    positive = intervals[intervals > 0]
    if not positive.size or np.any(intervals < 0):
        return None
    dt = float(np.median(positive))
    window = max(1, int(round(0.020 / dt)))
    if window % 2 == 0:
        window += 1
    window = min(window, t.size if t.size % 2 else t.size - 1)
    half = window // 2
    smooth = np.convolve(np.pad(voltage, (half, half), mode="edge"),
                         np.ones(window) / window, mode="valid")
    expected = ev["event1_end"] - offset
    duration = ev["event1_end"] - ev["event1_start"]
    pre_width = min(0.08, 0.40 * duration)
    before = smooth[(t >= expected - pre_width) & (t <= expected - 0.02)]
    after = smooth[(t >= expected + 0.03) & (t <= expected + 0.15)]
    if before.size < 3 or after.size < 3:
        return None
    low = float(np.median(before))
    high = float(np.median(after))
    direction = -np.sign(safe_num(meta.get("V_fault_pu")) - 1.0)
    change = direction * (high - low)
    if not np.isfinite(change) or change < 0.025:
        return None
    threshold = 0.5 * (low + high)
    search_start = max(expected - 0.08,
                       ev["event1_start"] - offset + min(0.04, duration / 3))
    search = (t >= search_start) & (t <= expected + 0.20)
    crossed = search & (direction * (smooth - threshold) >= 0)
    hold = max(1, int(math.ceil(0.010 / dt)))
    if hold > t.size:
        return None
    sustained = np.convolve(crossed.astype(int), np.ones(hold, dtype=int), mode="valid")
    crossings = np.flatnonzero(sustained == hold)
    if not crossings.size:
        return None
    measured = float(t[crossings[0]])
    if measured <= ev["event1_start"] - offset:
        return None
    return measured + offset


def detect_edge_offset(t_rel, sig, meta, ev, offset_prior):
    """在阶跃信号上检测事件"离开基线"沿，返回 offset（录制起点的仿真时刻）或 None。

    offset 定义：t_sim = t_rel + offset；检测成功时 offset = event1_start - t_edge_rel。

    用实测基线（先验窗中位数）而非 CSV 基准值定阈值——控制器稳态误差
    可达 ~3% pu，与小阶跃幅值同量级，直接用 CSV 中点会误检。
    检测"离开基线"而非"越过中点"，避免大阶跃斜坡传播带来的滞后。

    信号选择、基线、阶跃目标均由 profile 驱动，不再硬编码 is_frt 分支。
    """
    pname = ev.get("profile_name", "pwr")
    prof = _PROFILES.get(pname, _PROFILES["pwr"])
    code = ev.get("code", 3)

    if pname == "pwr":
        return detect_pwr_edge_offset(t_rel, sig, meta, ev, offset_prior)

    if pname == "frt":
        return detect_frt_edge_offset(t_rel, sig["Vp_pu"], meta, ev, offset_prior)

    if pname in ("freq_reg", "inertia"):
        return detect_frequency_edge_offset(t_rel, sig, meta, ev, offset_prior)

    # ---- 选择检测信号、基线、阶跃目标 ----
    edge_sig = prof["edge_signal"]

    if edge_sig == "auto":
        # PWR: code 4 → Q, 其他 → P
        edge_sig = "Q" if code == 4 else "P"

    if prof.get("edge_is_Vp_scale"):
        # FRT / VOLT_PROT: V+ 信号，转换到 "等效 W" 域走统一阈值逻辑
        x = sig["Vp_pu"] * HIL_P_BASE_W
        base = 1.0 * HIL_P_BASE_W
        step = safe_num(meta.get(prof["edge_step_field"])) * HIL_P_BASE_W
        _noise_tol = EDGE_NOISE_TOL_FRT * HIL_P_BASE_W
        _base_tol = EDGE_BASE_TOL * HIL_P_BASE_W
        _min_amp = EDGE_MIN_AMP_FRAC * HIL_P_BASE_W
        _delta_min = EDGE_DELTA_MIN_FRAC * HIL_P_BASE_W
    elif prof.get("edge_is_freq"):
        # FREQ_REG / INERTIA / FREQ_PROT: 频率信号，Hz 域
        x = sig["f"]
        if x is None:
            return None
        base = NOMINAL_FREQUENCY_HZ
        step = base + safe_num(meta.get(prof["edge_step_field"]))
        _noise_tol = 0.10 * base     # ~5 Hz — 保守
        _base_tol = 0.20 * base      # ~10 Hz
        _min_amp = 0.02 * base       # ~1 Hz
        _delta_min = 0.04 * base     # ~2 Hz
    else:
        # PWR / SCR: 功率信号，W 域
        x = sig[edge_sig]
        if x is None:
            return None
        base_field = prof["edge_baseline_field"]
        if base_field == "auto":
            base_field = "Q_base_pu" if code == 4 else "P_base_pu"
        step_field = prof["edge_step_field"]
        if step_field == "auto":
            step_field = "Q_step_pu" if code == 4 else "P_step_pu"
        base = safe_num(meta.get(base_field)) * HIL_P_BASE_W
        _min_amp = EDGE_MIN_AMP_FRAC * HIL_P_BASE_W
        if step_field is not None:
            step = safe_num(meta.get(step_field)) * HIL_P_BASE_W
        else:
            # SCR 等无明确阶跃目标：用基线 ± 最小检测幅值做占位
            step = base + _min_amp
        _noise_tol = EDGE_NOISE_TOL_PWR * HIL_P_BASE_W
        _base_tol = EDGE_BASE_TOL * HIL_P_BASE_W
        _delta_min = EDGE_DELTA_MIN_FRAC * HIL_P_BASE_W

    amp = step - base
    if abs(amp) < _min_amp:
        return None

    te0 = ev["event1_start"] - offset_prior  # 事件沿的先验相对时刻

    # ---- 实测基线：先验窗 [te0-1.2, te0-0.3] ----
    base_m = (t_rel >= te0 - EDGE_BASELINE_BACK_S) & (t_rel <= te0 - EDGE_SEARCH_BACK_S)
    if not np.any(base_m):
        return None
    base_seg = x[base_m]
    base_meas = float(np.median(base_seg))
    noise_half = float(base_seg.max() - base_seg.min()) / 2.0

    # ---- 基线有效性校验（失稳/漂移 case 在此被拒，走批中位数兜底） ----
    if abs(base_meas - base) > _base_tol or noise_half > _noise_tol:
        return None

    # ---- 离开基线阈值 ----
    delta = max(_delta_min, EDGE_DELTA_NOISE_MULT * noise_half)
    if delta >= EDGE_THRESH_RATIO_CAP * abs(amp):
        return None  # 阈值相对阶跃过大（小阶跃+纹波/稳态误差），放弃检测，走兜底
    thresh = base_meas + np.sign(amp) * delta

    # ---- 搜索窗 [te0-0.3, te0+1.0]，首次越过阈值且保持 EDGE_SUSTAIN_S ----
    m = (t_rel >= te0 - EDGE_SEARCH_BACK_S) & (t_rel <= te0 + EDGE_SEARCH_AHEAD_S)
    if not np.any(m):
        return None
    tt, xx = t_rel[m], x[m]
    dt = float(np.median(np.diff(tt)))
    crossed = (xx - thresh) * np.sign(amp) > 0
    n_sustain = max(1, int(EDGE_SUSTAIN_S / dt))
    idx = None
    for i in np.flatnonzero(crossed):
        seg = crossed[i:i + n_sustain]
        if seg.size == n_sustain and np.all(seg):
            idx = i
            break
    if idx is None:
        return None

    offset = ev["event1_start"] - float(tt[idx])
    _max_dev = prof.get("edge_max_dev_s", EDGE_MAX_DEV_S)
    if abs(offset - offset_prior) > _max_dev:
        return None
    return offset


def detect_bnd_edge_offset(t_rel, sig, meta, ev, offset_prior):
    """检测 FRT_BND 的短、小幅电压扰动沿。"""
    signal = np.asarray(sig["Vp_pu"], dtype=float)
    if len(t_rel) < 3 or len(signal) != len(t_rel):
        return None

    dt = float(np.median(np.diff(t_rel)))
    if not np.isfinite(dt) or dt <= 0.0:
        return None

    te0 = ev["event1_start"] - offset_prior
    smooth_n = max(3, int(round(0.03 / dt)))
    smooth = np.convolve(signal, np.ones(smooth_n) / smooth_n, mode="same")
    baseline_mask = (t_rel >= te0 - 0.8) & (t_rel <= te0 - 0.1)
    if not np.any(baseline_mask):
        return None

    baseline = float(np.median(smooth[baseline_mask]))
    noise = float(np.std(smooth[baseline_mask]))
    target = safe_num(meta.get("V_fault_pu"), baseline)
    direction = np.sign(target - 1.0)
    if direction == 0.0:
        return None

    threshold = max(0.2 * abs(target - 1.0), 4.0 * noise, 0.008)
    search_mask = (t_rel >= te0 - 0.1) & (t_rel <= te0 + 0.5)
    crossed = search_mask & (direction * (smooth - baseline) > threshold)
    hold_n = max(1, int(round(0.015 / dt)))
    held = np.convolve(crossed.astype(int), np.ones(hold_n, dtype=int), mode="same")
    candidates = np.flatnonzero(held >= hold_n)
    if len(candidates) == 0:
        return None

    # 卷积窗口居中，回退半个窗口得到持续越限的起点。
    idx = max(0, int(candidates[0]) - hold_n // 2)
    offset = ev["event1_start"] - float(t_rel[idx])
    if abs(offset - offset_prior) > ev.get("edge_max_dev_s", EDGE_MAX_DEV_S):
        return None
    return offset


def resolve_case_offsets(cases):
    """为批次 case 解析录波偏移；检测失败时使用批次中位数。"""
    detected = {}
    fallback_samples = []
    resolved = {}

    for ci, mat_path, meta in cases:
        ev = event_times_from_meta(meta)
        prior = safe_num(meta.get("hil_rec_start_target_s"), REC_START_NOMINAL_S)
        # BND 是短、小幅扰动，使用平滑后的专用检测；失败时回退自身
        # 元数据，不能套用普通 FRT 批次中位数。
        if "BND" in str(meta.get("test_type", "")).upper():
            try:
                t_rel, sig = load_case_file(mat_path)
                offset = detect_bnd_edge_offset(t_rel, sig, meta, ev, prior)
            except Exception:
                offset = None
            resolved[ci] = ((offset, "edge") if offset is not None
                            else (prior, "nominal"))
            if offset is not None:
                detected[ci] = offset
            continue
        if ev.get("profile_name") not in ("frt", "pwr", "freq_reg", "inertia"):
            resolved[ci] = (prior, "nominal")
            continue
        try:
            t_rel, sig = load_case_file(mat_path)
            offset = detect_edge_offset(t_rel, sig, meta, ev, prior)
            if offset is not None:
                detected[ci] = offset
                fallback_samples.append(offset)
        except Exception:
            pass

    # BND 使用不同录波窗口，其偏移不能参与普通 FRT/PWR 的兜底中位数。
    median_offset = (float(np.median(fallback_samples))
                     if fallback_samples else None)
    for ci, _mat_path, meta in cases:
        if ci in resolved:
            continue
        if ci in detected:
            resolved[ci] = (detected[ci], "edge")
        elif median_offset is not None:
            resolved[ci] = (median_offset, "batch_median")
        else:
            prior = safe_num(meta.get("hil_rec_start_target_s"), REC_START_NOMINAL_S)
            resolved[ci] = (prior, "nominal")
    return resolved, detected, median_offset


# =============================================================================
#  事件前稳态检查（Profile 驱动）
# =============================================================================

def _check_pre_event(meta, t_sim, sig, ev, P_pu, Q_pu, Vp_pu, f_sig):
    """根据 profile 的 pre_check 类型做事件前稳态检查。

    返回 (pre_event_ok, pre_unstable, post_unstable, stats_dict)。
    """
    pname = ev.get("profile_name", "pwr")
    prof = _PROFILES.get(pname, _PROFILES["pwr"])
    check_type = prof.get("pre_check", "pwr")

    ev1 = ev["event1_start"]
    ev1_end = ev["event1_end"]
    # SIM starts at model initialization; use the displayed final second before
    # the fault so the startup ramp is not mistaken for an unstable baseline.
    pre_window = (min(ev["record_pre"], NORMALIZED_PRE_EVENT_S)
                  if meta.get("data_source") == "sim" else ev["record_pre"])
    pre_m = (t_sim >= ev1 - pre_window) & (t_sim <= ev1 - EDGE_SEARCH_BACK_S)
    event_m = (t_sim >= ev1) & (t_sim <= ev1_end)
    post_m = (t_sim >= ev["last_event_end"] + EDGE_POST_WINDOW_S)

    if check_type == "frt":
        # FRT / VOLT_PROT: V+ 基线检查
        Vf = safe_num(meta.get("V_fault_pu"))
        pre_Vp_mean, pre_Vp_pp = _win_stats(Vp_pu, pre_m)
        event_Vp_mean, _       = _win_stats(Vp_pu, event_m)
        post_Vp_mean, _        = _win_stats(Vp_pu, post_m)
        pre_P_mean, pre_P_pp   = _win_stats(P_pu, pre_m)
        pre_Q_mean, _          = _win_stats(Q_pu, pre_m)

        pre_unstable = (not math.isnan(pre_Vp_pp)) and pre_Vp_pp > FRT_UNSTABLE_VPP
        post_unstable = False
        pre_event_ok = (
            not pre_unstable
            and not math.isnan(pre_Vp_mean)
            and abs(pre_Vp_mean - 1.0) <= FRT_PRE_VP_TOL_PU
            and pre_Vp_pp <= FRT_PRE_VPP_TOL_PU
        )
        stats = {
            "pre_Vp_mean_pu": pre_Vp_mean, "pre_Vp_pp_pu": pre_Vp_pp,
            "event_Vp_mean_pu": event_Vp_mean, "post_Vp_mean_pu": post_Vp_mean,
            "pre_P_mean_pu": pre_P_mean,
            "pre_Q_mean_pu": pre_Q_mean,
        }
        return pre_event_ok, pre_unstable, post_unstable, stats

    elif check_type == "pwr":
        # PWR / SCR: P/Q 基线检查
        P_base_pu = safe_num(meta.get("P_base_pu"))
        Q_base_pu = safe_num(meta.get("Q_base_pu"))
        pre_P_mean, pre_P_pp = _win_stats(P_pu, pre_m)
        pre_Q_mean, pre_Q_pp = _win_stats(Q_pu, pre_m)
        step_P_mean, _       = _win_stats(P_pu, event_m)
        _, post_P_pp         = _win_stats(P_pu, post_m)

        pre_unstable = (not math.isnan(pre_P_pp)) and pre_P_pp > UNSTABLE_PP_PU
        post_unstable = (not math.isnan(post_P_pp)) and post_P_pp > UNSTABLE_PP_PU
        pre_event_ok = (
            not pre_unstable
            and not math.isnan(pre_P_mean)
            and abs(pre_P_mean - P_base_pu) <= PRE_P_TOL_PU
            and pre_P_pp <= PRE_PP_TOL_PU
            and (math.isnan(Q_base_pu) or abs(pre_Q_mean - Q_base_pu) <= PRE_Q_TOL_PU)
        )
        stats = {
            "pre_P_mean_pu": pre_P_mean, "pre_P_pp_pu": pre_P_pp,
            "pre_Q_mean_pu": pre_Q_mean,
            "step_P_mean_pu": step_P_mean,
        }
        return pre_event_ok, pre_unstable, post_unstable, stats

    else:  # "freq"
        # FREQ_REG / INERTIA / FREQ_PROT: P/Q/f 基线检查
        P_ref_pu = safe_num(meta.get("P_ref_pu"))
        Q_ref_pu = safe_num(meta.get("Q_ref_pu"))
        pre_P_mean, pre_P_pp = _win_stats(P_pu, pre_m)
        pre_Q_mean, _        = _win_stats(Q_pu, pre_m)
        pre_f_mean, pre_f_pp = _win_stats(f_sig, pre_m)
        event_P_mean, _      = _win_stats(P_pu, event_m)
        _, post_P_pp         = _win_stats(P_pu, post_m)

        pre_unstable = (not math.isnan(pre_P_pp)) and pre_P_pp > UNSTABLE_PP_PU
        post_unstable = (not math.isnan(post_P_pp)) and post_P_pp > UNSTABLE_PP_PU

        f_ok = (math.isnan(pre_f_mean)
                or (not math.isnan(pre_f_mean) and abs(pre_f_mean - NOMINAL_FREQUENCY_HZ) <= 0.2))
        pre_event_ok = (
            not pre_unstable
            and not math.isnan(pre_P_mean)
            and abs(pre_P_mean - P_ref_pu) <= PRE_P_TOL_PU
            and pre_P_pp <= PRE_PP_TOL_PU
            and (math.isnan(Q_ref_pu) or abs(pre_Q_mean - Q_ref_pu) <= PRE_Q_TOL_PU)
            and f_ok
        )
        stats = {
            "pre_P_mean_pu": pre_P_mean, "pre_P_pp_pu": pre_P_pp,
            "pre_Q_mean_pu": pre_Q_mean,
            "step_P_mean_pu": event_P_mean,
            "pre_f_mean_hz": pre_f_mean,
        }
        return pre_event_ok, pre_unstable, post_unstable, stats


# =============================================================================
#  汇总 CSV 构建（共享数据处理工具）
# =============================================================================

def build_summary_row(meta, stem, offset, offset_source,
                      pre_event_ok, pre_unstable, post_unstable,
                      stats, profile_name):
    """构建 模式汇总 CSV 的一行 dict。

    Parameters
    ----------
    profile_name : str
        _PROFILES 的 key（如 "pwr", "frt", "freq_reg"），用于查找 summary_cols。
    """
    prof = _PROFILES.get(profile_name, _PROFILES["pwr"])
    row = {
        "case_index":    int(safe_num(meta.get("case_index"))),
        "case_name":     str(meta.get("case_name", stem)),
        "test_type":     (
            "FREQ_REG_DEADBAND"
            if profile_name == "freq_reg" and is_freq_reg_deadband(meta)
            else str(meta.get("test_type", ""))
        ),
        "rec_offset_s":  round(offset, 4),
        "offset_source": offset_source,
        "pre_unstable":  int(pre_unstable),
        "post_unstable": int(post_unstable),
        "pre_event_ok":  int(pre_event_ok),
    }
    for col in prof["summary_cols"]:
        if col in stats:
            formatted = _fmt(stats[col])
            row[col] = (
                int(formatted)
                if col == "P_change_5pct_pass" and formatted != ""
                else formatted
            )
        elif col in meta:
            row[col] = safe_num(meta.get(col))
        else:
            row[col] = ""
        # pre_Q_mean_pu 缺失时填 0（如 HIL 未录 Q 信号）
        if col == "pre_Q_mean_pu" and (
            row[col] == "" or (isinstance(row[col], float) and math.isnan(row[col]))
        ):
            row[col] = 0
    return row


def write_summary_csv(rows, output_path, profile_name="pwr"):
    """写入模式汇总 CSV。

    Parameters
    ----------
    profile_name : str
        _PROFILES 的 key，用于获取 summary_cols 列顺序。
    """
    common_cols = ["case_index", "case_name", "test_type", "rec_offset_s",
                   "offset_source", "pre_unstable", "post_unstable", "pre_event_ok"]
    prof = _PROFILES.get(profile_name, _PROFILES["pwr"])
    fieldnames = common_cols + prof["summary_cols"]
    if profile_name == "pwr":
        fieldnames += PWR_STEP_SUMMARY_COLS
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def add_missing_freq_reg_summary_rows(rows, case_csv, case_filter=None):
    """Add explicit placeholder rows for expected FREQ_REG recordings not found."""
    present_ids = {int(row["case_index"]) for row in rows}
    with open(case_csv, "r", encoding="utf-8-sig", newline="") as f:
        expected_rows = list(csv.DictReader(f))

    common_cols = [
        "case_index", "case_name", "test_type", "rec_offset_s",
        "offset_source", "pre_unstable", "post_unstable", "pre_event_ok",
    ]
    summary_cols = _PROFILES["freq_reg"]["summary_cols"]
    for expected in expected_rows:
        case_index = int(safe_num(expected.get("case_index")))
        if case_filter and case_index not in case_filter:
            continue
        if case_index in present_ids:
            continue

        placeholder = {col: "" for col in common_cols + summary_cols}
        placeholder.update({
            "case_index": case_index,
            "case_name": str(expected.get("case_name", "")),
            "test_type": str(expected.get("test_type", "")),
            "offset_source": "missing_recording",
            "pre_event_ok": 0,
        })
        for col in summary_cols:
            if col in expected:
                placeholder[col] = safe_num(expected.get(col))
        rows.append(placeholder)

    rows.sort(key=lambda row: int(row["case_index"]))
    return rows


# =============================================================================
#  单 case 处理（Profile 驱动）
# =============================================================================

@dataclass
class ProcessedCase:
    stem: str
    meta: dict
    profile_name: str
    bridged: dict
    events: dict
    summary: dict
    pwr_downsampled: dict | None = None


def event_timing_metadata(meta, ev, offset, offset_source):
    """Describe timing provenance without calling a response edge a command.

    Optional synchronized marker contract: hil_command_event1_recorder_s is
    seconds on load_case_file's recorder axis, with hil_command_time_source
    explicitly set to 'synchronized'. Host scheduling timestamps do not qualify.
    """
    origin = ev["event1_start"] - NORMALIZED_PRE_EVENT_S
    detected = offset_source == "edge"
    recorder_edge = ev["event1_start"] - offset if detected else float("nan")
    command = safe_num(meta.get("hil_command_event1_recorder_s"), float("nan"))
    is_sim = meta.get("data_source") == "sim"
    if is_sim:
        command = ev["event1_start"]
    reliable = (is_sim or meta.get("hil_command_time_source") == "synchronized") and np.isfinite(command)
    if not reliable:
        command = float("nan")
    prof = _PROFILES[ev["profile_name"]]
    signal = prof["edge_signal"]
    if signal == "auto":
        signal = "Q" if ev["code"] == 4 else "P"
    return {
        "timing_schema_version": 1,
        "offset_source": offset_source,
        "response_edge_detected": int(detected),
        "response_edge_recorder_s": recorder_edge,
        "response_edge_normalized_s": recorder_edge + offset - origin,
        "edge_detection_signal": signal,
        "edge_detection_method": ("pwr_smoothed_noise_sustained_v1"
                                  if ev["profile_name"] == "pwr" else "profile_detector"),
        "edge_smooth_ms": 20.0 if ev["profile_name"] == "pwr" else float("nan"),
        "response_edge_kind": "response_edge_estimate" if detected else "unavailable",
        "command_event1_effective_s": ev["event1_start"],
        "command_effective_time_axis": "effective_case_time_metadata",
        "command_event1_recorder_s": command,
        "command_event1_normalized_s": command + offset - origin,
        "command_time_reliable": int(reliable),
        "command_time_source": ("sim_command" if is_sim else
                                "synchronized" if reliable else "metadata_only"),
        "recorder_to_normalized_shift_s": offset - origin,
        "timing_basis": ("sim_command" if is_sim else
                         "synchronized_command" if reliable else
                         "response_edge_estimate" if detected else "unavailable"),
    }


def prepare_case(mat_path, meta, offset, offset_source):
    """Prepare aligned, normalized data and statistics without rendering or writing."""
    base = os.path.basename(mat_path)
    stem = os.path.splitext(base)[0]          # case_XXXX_<name>
    t_rel, sig = load_case_file(mat_path)
    ev = event_times_from_meta(meta)
    t_sim = t_rel + offset

    pname = ev.get("profile_name", "pwr")
    if pname == "frt" and "BND" not in str(meta.get("test_type", "")).upper():
        recovery = detect_frt_recovery_time(t_rel, sig["Vp_pu"], meta, ev, offset)
        if recovery is not None:
            ev["event1_end"] = recovery
            ev["last_event_end"] = recovery
    add_power_current = supports_power_current(meta, pname)
    if add_power_current:
        ip, iq, current_valid = calculate_instantaneous_current(
            np.column_stack([sig["Va"], sig["Vb"], sig["Vc"]]),
            np.column_stack([sig["Ia"], sig["Ib"], sig["Ic"]]))
        ip_filtered = filter_current_for_display(t_rel, ip, current_valid)
        iq_filtered = filter_current_for_display(t_rel, iq, current_valid)

    P_pu = sig["P"] / HIL_P_BASE_W
    Q_pu = sig["Q"] / HIL_P_BASE_W if sig["Q"] is not None else np.full_like(P_pu, np.nan)
    Vp_pu = sig["Vp_pu"]
    f_sig = sig["f"]

    # ---- 事件前稳态检查 ----
    pre_event_ok, pre_unstable, post_unstable, stats = _check_pre_event(
        meta, t_sim, sig, ev, P_pu, Q_pu, Vp_pu, f_sig)
    add_freq_reg_power_stats(stats, meta, t_sim, P_pu, ev)

    # offset 由实际事件沿校准后，事件沿与元数据 event1_start 重合；
    # 保护测试不做信号沿猜测，直接以元数据事件作为 t1。
    t_norm, ev_norm, wm = normalize_case_time(t_sim, ev, ev["event1_start"])

    if pname in ("volt_prot", "freq_prot"):
        t1_raw = ev["event1_start"]
        t2_raw = detect_trip_time(t_sim, sig["breaker"], t1_raw)
        origin = t1_raw - NORMALIZED_PRE_EVENT_S
        stats["t1_s"] = NORMALIZED_PRE_EVENT_S
        stats["t2_s"] = t2_raw - origin if np.isfinite(t2_raw) else np.nan
        stats["delta_t_s"] = (t2_raw - t1_raw
                              if np.isfinite(t2_raw) else np.nan)

    # ---- 桥接输出（裁剪窗口） ----
    bridged = {
        "t_s":       t_norm[wm],
        "Vabc_V":    np.column_stack([sig["Va"][wm], sig["Vb"][wm], sig["Vc"][wm]]),
        "Iabc_A":    np.column_stack([sig["Ia"][wm], sig["Ib"][wm], sig["Ic"][wm]]),
        "P_W":       sig["P"][wm],
        "Vabc_pu":   np.column_stack([sig["Va"][wm], sig["Vb"][wm], sig["Vc"][wm]]) / HIL_V_BASE_PEAK_V,
        "Iabc_pu":   np.column_stack([sig["Ia"][wm], sig["Ib"][wm], sig["Ic"][wm]]) / HIL_I_BASE_PEAK_A,
        "P_pu":      P_pu[wm],
        "Vp_pu":     Vp_pu[wm],
        "event1_start_s": ev_norm["event1_start"],
        "event1_end_s":   ev_norm["event1_end"],
        "rec_offset_s":   offset,
        "case_index":     int(safe_num(meta.get("case_index"))),
        "meta_json":      json.dumps(meta, ensure_ascii=False),
    }
    if sig.get("breaker") is not None:
        bridged["breaker"] = sig["breaker"][wm]
    bridged.update(event_timing_metadata(meta, ev, offset, offset_source))
    if sig["Q"] is not None:
        bridged["Q_var"] = sig["Q"][wm]
        bridged["Q_pu"]  = Q_pu[wm]
    if add_power_current:
        bridged.update({
            "Ip_A": ip[wm], "Iq_A": iq[wm],
            "Ip_pu": ip[wm] / HIL_I_BASE_PEAK_A,
            "Iq_pu": iq[wm] / HIL_I_BASE_PEAK_A,
            "current_valid": current_valid[wm],
            "Ip_filtered_A": ip_filtered[wm], "Iq_filtered_A": iq_filtered[wm],
            "Ip_filtered_pu": ip_filtered[wm] / HIL_I_BASE_PEAK_A,
            "Iq_filtered_pu": iq_filtered[wm] / HIL_I_BASE_PEAK_A,
        })
    if f_sig is not None:
        bridged["f_hz"] = f_sig[wm]
    if pname in ("volt_prot", "freq_prot"):
        bridged["t1_s"] = stats["t1_s"]
        bridged["t2_s"] = stats["t2_s"]
        bridged["delta_t_s"] = stats["delta_t_s"]
    if ev_norm["event2_start"] is not None:
        bridged["event2_start_s"] = ev_norm["event2_start"]
        bridged["event2_end_s"] = ev_norm["event2_end"]
    # ---- 汇总行 ----
    row = build_summary_row(meta, stem, offset, offset_source,
                            pre_event_ok, pre_unstable, post_unstable,
                            stats, pname)
    downsampled = None
    if pname == "pwr":
        step_summary, downsampled = pwr_step_summary(bridged, meta, return_downsampled=True)
        row.update(step_summary)

    return ProcessedCase(stem, dict(meta), pname, bridged, ev_norm, row, downsampled)


def save_bridged(case, bridged_dir):
    """Write a prepared case in the existing MATLAB-compatible format."""
    os.makedirs(bridged_dir, exist_ok=True)
    path = os.path.join(bridged_dir, case.stem + ".mat")
    sio.savemat(path, case.bridged, do_compression=True)
    return path


def save_pwr_downsampled(case, output_dir):
    """Export the exact 5 ms P/Q series used by HIL PWR metrics and plots."""
    if case.profile_name != "pwr" or case.pwr_downsampled is None:
        return None
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, case.stem + "_5ms.csv")
    data = case.pwr_downsampled
    with open(path, "w", newline="", encoding="utf-8-sig") as file:
        writer = csv.writer(file)
        writer.writerow(("t_s", "P_pu_20ms", "Q_pu_20ms"))
        for t, p, q in zip(data["t_s"], data["P_pu"], data["Q_pu"]):
            writer.writerow((t, p if np.isfinite(p) else "",
                             q if np.isfinite(q) else ""))
    return path


# =============================================================================
#  批处理主流程
# =============================================================================

def select_latest_cases(cases):
    """Select one recording per case index before resolving batch offsets.

    Imported batches can leave older case names in hil_raw. File modification
    time selects the newest recording; path order breaks equal-time ties.
    Original recordings are never removed.
    """
    selected = {}
    for case in cases:
        ci, path, _meta = case
        previous = selected.get(ci)
        if previous is None or (os.stat(path).st_mtime_ns, path) > (
                os.stat(previous[1]).st_mtime_ns, previous[1]):
            selected[ci] = case
    result = [selected[ci] for ci in sorted(selected)]
    for ci, path, _meta in cases:
        if path != selected[ci][1]:
            print("SKIP older recording: %s; selected: %s" %
                  (os.path.basename(path), os.path.basename(selected[ci][1])))
    return result


def summary_path_for_mode(result_dir, mode, for_read=False):
    """New writes use mode-specific names; reads support existing legacy files."""
    if mode not in TEST_MODE_CONFIG:
        raise ValueError("Unknown mode: %s" % mode)
    path = os.path.join(result_dir, "hil_summary_%s.csv" % mode)
    if for_read and not os.path.exists(path):
        legacy = os.path.join(result_dir, "hil_summary.csv")
        if os.path.exists(legacy):
            return legacy
    return path


def load_cases(raw_dir, mode, case_filter=None):
    """Load paired metadata, filter cases, and apply the mode's selection rule."""
    cases = []
    for path in sorted(glob.glob(os.path.join(raw_dir, "case_*.mat"))):
        meta_path = os.path.splitext(path)[0] + "_meta.json"
        if not os.path.exists(meta_path):
            print("WARN: no meta for %s, skipped" % os.path.basename(path))
            continue
        with open(meta_path, encoding="utf-8") as fh:
            meta = json.load(fh)
        ci = int(safe_num(meta.get("case_index")))
        if case_filter is not None and ci not in case_filter:
            continue
        cases.append((ci, path, meta))
    return select_latest_cases(cases) if mode == "frt" else cases
