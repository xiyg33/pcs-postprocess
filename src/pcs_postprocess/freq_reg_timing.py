"""Per-recording frequency alignment, independent of the power metric clock."""

from pathlib import Path

import numpy as np


SUMMARY_FIELDS = [
    "freq_event_start_recorder_s", "freq_event_end_recorder_s",
    "freq_timing_method", "freq_timing_confidence", "freq_timing_review_required",
    "freq_timing_reason", "freq_candidate_start_recorder_s", "freq_candidate_end_recorder_s",
    "freq_fit_noise_hz", "freq_fit_rmse_hz", "freq_power_start_recorder_s",
    "freq_timing_review_figure",
    "freq_metric_end_recorder_s",
]

REASON_LABELS = {
    "insufficient_frequency_baseline": "扰动前频率基线不足",
    "frequency_samples_missing": "频率录波存在缺失采样",
    "frequency_transition_unconfirmed": "未确认持续频率变化",
    "frequency_fit_unavailable": "起始沿或恢复沿无法可靠拟合",
    "fit_at_search_boundary": "拟合贴近搜索边界",
    "frequency_amplitude_mismatch": "实测频率变化幅值与工况不一致",
    "frequency_hold_duration_mismatch": "起止沿间隔与保持时长不一致",
    "fit_sensitive_to_window": "拟合起点随观察窗口变化",
    "frequency_fit_residual_high": "拟合残差偏大",
    "frequency_fit_outliers_excessive": "频率波动与拟合过渡不一致",
    "frequency_step_low_signal_to_noise": "扰动幅值较小，频率纹波影响明显",
    "power_response_precedes_frequency_estimate": "有功响应明显早于估计频率起点",
}


def is_ordinary(meta):
    from .processing import is_freq_reg_deadband
    return meta.get("mode") == "freq_reg" and not is_freq_reg_deadband(meta)


def _fit(grid, values, nominal, deviation, noise, padding=.30):
    """Fit baseline/ramp/plateau; the threshold never becomes the time anchor."""
    keep = (grid >= nominal - padding) & (grid <= nominal + .75)
    t, y = grid[keep], values[keep]
    if len(t) < 30 or np.any(~np.isfinite(y)):
        return None
    starts = nominal + np.arange(-60, 71) * .005
    durations = np.arange(2, 71) * .005
    onset = np.repeat(starts, len(durations))
    duration = np.tile(durations, len(starts))
    ramps = np.clip((t[None, :] - onset[:, None]) / duration[:, None], 0., 1.)
    # Robust endpoint levels keep a spike in the search window from moving the
    # regression baseline or shrinking the fitted event amplitude.
    before = values[(grid >= nominal - 1.2) & (grid <= nominal - .3)]
    before = before[np.isfinite(before)]
    after = y[t >= nominal + .50]
    if min(len(before), len(after)) < 10:
        return None
    baseline = float(np.median(before))
    amplitude = float(np.median(after)) - baseline
    if not .6 <= amplitude / deviation <= 1.4:
        return None
    residual = y[None, :] - baseline - amplitude * ramps
    # Clipped loss reduces the influence of isolated PLL outliers. Direction
    # and amplitude constraints prevent fitting a small ripple as an event.
    scale = max(3. * noise, .003)
    score = np.mean(np.minimum(residual ** 2, scale ** 2), axis=1)
    i = int(np.argmin(score))
    if not np.isfinite(score[i]):
        return None
    return {"start_s": float(onset[i]), "duration_s": float(duration[i]),
            "amplitude_hz": amplitude, "baseline_hz": baseline,
            "rmse_hz": float(np.sqrt(score[i])),
            "outlier_fraction": float(np.mean(np.abs(residual[i]) > scale)),
            "boundary": bool(i // len(durations) in (0, len(starts) - 1)
                             or i % len(durations) in (0, len(durations) - 1))}


def power_onset(t, power_pu, nominal, deviation, end, search_before_s=.3):
    """Use a clean metadata baseline on the complete recorder axis."""
    from .processing import pwr_five_point_edge, detect_pwr_edge_offset, HIL_P_BASE_W
    t, y = np.asarray(t, float), np.asarray(power_pu, float)
    search_begin = nominal - search_before_s
    baseline = (t >= search_begin - .5) & (t < search_begin - .05)
    tail = (t >= end - .5) & (t < end)
    if min(int(baseline.sum()), int(tail.sum())) < 3:
        return None
    initial, final = np.mean(y[baseline]), np.mean(y[tail])
    if not np.isfinite(initial + final) or abs(final - initial) < 1e-6:
        return None
    delta = -np.sign(deviation) * abs(final - initial)
    edge, _ = pwr_five_point_edge(t, y, search_begin, delta, max_search_s=.35 + search_before_s)
    if edge is not None:
        return float(edge)
    meta = {"P_base_pu": float(initial), "P_step_pu": float(initial + delta)}
    offset = detect_pwr_edge_offset(t, {"P": y * HIL_P_BASE_W}, meta,
                                   {"code": 3, "event1_start": nominal}, 0.)
    return None if offset is None else float(nominal - offset)


def resolve(t, sig, meta, ev, prior, power_base_w):
    """Resolve a complete pair of event times; unresolved fits remain provisional."""
    from .processing import downsample_scalar_signals
    t = np.asarray(t, float)
    absolute = meta.get("time_basis") == "absolute"
    nominal = float(ev["event1_start"] if absolute else ev["event1_start"] - prior)
    hold = float(meta["event_hold_s"])
    selected = nominal if not absolute else float(ev["event1_start"])
    result = {"freq_event_start_recorder_s": selected,
              "freq_event_end_recorder_s": selected + hold,
              "freq_timing_method": "absolute_time" if absolute else "frequency_nominal_unconfirmed",
              "freq_timing_confidence": "absolute" if absolute else "review_required",
              "freq_timing_review_required": int(not absolute), "freq_timing_reason": "",
              "freq_candidate_start_recorder_s": None, "freq_candidate_end_recorder_s": None,
              "freq_fit_noise_hz": None, "freq_fit_rmse_hz": None,
              "freq_power_start_recorder_s": None, "freq_timing_review_figure": "",
              "freq_metric_end_recorder_s": None}
    deviation = float(meta["freq_deviation_hz"])
    reasons = []
    if not absolute:
        display = downsample_scalar_signals({"t_s": t, "f_hz": sig["f"]}, ("f_hz",))
        grid, frequency = display["t_s"], display["f_hz"]
        baseline = (grid >= nominal - 1.2) & (grid <= nominal - .3)
        clean = frequency[baseline & np.isfinite(frequency)]
        noise = float(1.4826 * np.median(np.abs(clean - np.median(clean)))) if len(clean) >= 10 else None
        result["freq_fit_noise_hz"] = noise
        if noise is None or abs(deviation) < 1e-12:
            reasons.append("insufficient_frequency_baseline")
        else:
            threshold = max(3. * noise, .05 * abs(deviation), .002)
            region = (grid >= nominal - .3) & (grid <= nominal + .70)
            excursion = (frequency[region] - np.median(clean)) * np.sign(deviation)
            if np.any(~np.isfinite(excursion)):
                reasons.append("frequency_samples_missing")
            elif np.count_nonzero(excursion > threshold) < 5:
                reasons.append("frequency_transition_unconfirmed")
            else:
                start_fit = _fit(grid, frequency, nominal, deviation, noise)
                end_fit = _fit(grid, frequency, nominal + hold, -deviation, noise)
                result["start_fit"], result["end_fit"] = start_fit, end_fit
                if start_fit is None or end_fit is None:
                    reasons.append("frequency_fit_unavailable")
                else:
                    start, end = start_fit["start_s"], end_fit["start_s"]
                    result.update(freq_candidate_start_recorder_s=start,
                                  freq_candidate_end_recorder_s=end,
                                  freq_fit_rmse_hz=max(start_fit["rmse_hz"], end_fit["rmse_hz"]))
                    if start_fit["boundary"] or end_fit["boundary"]:
                        reasons.append("fit_at_search_boundary")
                    for expected, amplitude in ((deviation, start_fit["amplitude_hz"]),
                                                (-deviation, end_fit["amplitude_hz"])):
                        if abs(amplitude - expected) > max(.2 * abs(deviation), 3. * noise):
                            reasons.append("frequency_amplitude_mismatch")
                    if abs(end - start - hold) > .05:
                        reasons.append("frequency_hold_duration_mismatch")
                    alternate_start = _fit(grid, frequency, nominal, deviation, noise, padding=.20)
                    alternate_end = _fit(grid, frequency, nominal + hold, -deviation, noise, padding=.20)
                    if (alternate_start is None or alternate_end is None
                            or abs(alternate_start["start_s"] - start) > .0200001
                            or abs(alternate_end["start_s"] - end) > .0200001):
                        reasons.append("fit_sensitive_to_window")
                    if result["freq_fit_rmse_hz"] > max(2. * noise, .12 * abs(deviation), .004):
                        reasons.append("frequency_fit_residual_high")
                    if max(start_fit["outlier_fraction"], end_fit["outlier_fraction"]) > .15:
                        reasons.append("frequency_fit_outliers_excessive")
                    if abs(deviation) < 8. * noise:
                        reasons.append("frequency_step_low_signal_to_noise")
                    if not reasons:
                        result.update(freq_event_start_recorder_s=start, freq_event_end_recorder_s=end,
                                      freq_timing_method="frequency_piecewise_fit",
                                      freq_timing_confidence="accepted_estimate", freq_timing_review_required=0)
    # Synchronized command markers outrank manually reviewed frequency estimates.
    if (not absolute and meta.get("hil_command_time_source") == "synchronized"
            and "hil_command_event1_recorder_s" in meta):
        start = float(meta["hil_command_event1_recorder_s"])
        end = float(meta.get("hil_command_event1_end_recorder_s", start + hold))
        _validate_times(start, end, t)
        result.update(freq_event_start_recorder_s=start, freq_event_end_recorder_s=end,
                      freq_timing_method="synchronized_command", freq_timing_confidence="synchronized",
                      freq_timing_review_required=0)
        reasons = []
    elif not absolute and meta.get("freq_timing_override_source") == "manual_confirmed":
        start = float(meta["freq_event1_start_recorder_s"])
        end = float(meta["freq_event1_end_recorder_s"])
        _validate_times(start, end, t)
        result.update(freq_event_start_recorder_s=start, freq_event_end_recorder_s=end,
                      freq_timing_method="manual_frequency_estimate", freq_timing_confidence="manual_confirmed",
                      freq_timing_review_required=0)
        reasons = []
    # P has its own axis and baseline, even when the display must use a nominal
    # frequency start pending review. It must never replace the frequency anchor.
    metric_end = result["freq_candidate_end_recorder_s"] if result["freq_timing_review_required"] else result["freq_event_end_recorder_s"]
    if metric_end is None:
        metric_end = nominal + hold
    p_start = power_onset(t, sig["P"] / power_base_w, nominal, deviation, min(metric_end, float(t[-1])),
                          search_before_s=0. if absolute else .3)
    result["freq_power_start_recorder_s"] = p_start
    candidate = result["freq_candidate_start_recorder_s"]
    if (not absolute and result["freq_timing_confidence"] == "accepted_estimate"
            and p_start is not None and candidate is not None and p_start < candidate - .015):
        reasons.append("power_response_precedes_frequency_estimate")
        result.update(freq_event_start_recorder_s=nominal, freq_event_end_recorder_s=nominal + hold,
                      freq_timing_method="frequency_nominal_unconfirmed",
                      freq_timing_confidence="review_required", freq_timing_review_required=1)
    result["freq_timing_reason"] = "; ".join(dict.fromkeys(reasons))
    return result


def _validate_times(start, end, t):
    if not np.isfinite(start + end) or not t[0] <= start < end <= t[-1]:
        raise ValueError("Frequency timing override must satisfy recorder first <= start < end <= last")


def render_review(t, sig, meta, details, output):
    """Show the full-resolution observations and candidates for human review."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .processing import downsample_scalar_signals, HIL_P_BASE_W
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    display = downsample_scalar_signals({"t_s": t, "f_hz": sig["f"],
                                        "P_pu": sig["P"] / HIL_P_BASE_W}, ("f_hz", "P_pu"))
    fig, axes = plt.subplots(2, 2, figsize=(12, 7), sharex="col")
    try:
        for col, edge in enumerate(("start", "end")):
            selected = details[f"freq_event_{edge}_recorder_s"]
            candidate = details[f"freq_candidate_{edge}_recorder_s"]
            center = selected if candidate is None else candidate
            window = (np.asarray(t) >= center - .4) & (np.asarray(t) <= center + .65)
            axes[0, col].plot(np.asarray(t)[window][::10], sig["f"][window][::10],
                              color="gray", alpha=.35, lw=.5, label="原始 PLL 频率")
            fitted = details.get(edge + "_fit")
            if fitted is not None:
                model_t = np.linspace(center - .4, center + .65, 400)
                model_f = fitted["baseline_hz"] + fitted["amplitude_hz"] * np.clip(
                    (model_t - fitted["start_s"]) / fitted["duration_s"], 0., 1.)
                axes[0, col].plot(model_t, model_f, color="tab:purple", ls=":", lw=1., label="拟合过渡")
            for ax, name, label in ((axes[0, col], "f_hz", "频率 / Hz"),
                                     (axes[1, col], "P_pu", "有功 / pu")):
                ax.plot(display["t_s"], display[name], color="tab:blue", lw=1., label="20 ms 平滑")
                ax.axvline(selected, color="gray", ls="--", label="当前显示锚点")
                if candidate is not None:
                    ax.axvline(candidate, color="tab:orange", ls="--", label="待复核拟合候选")
                proposal = details.get("manual_proposed_" + edge + "_recorder_s")
                if proposal is not None:
                    ax.axvline(proposal, color="black", ls="-.", lw=.9, label="人工复核建议值（未确认）")
                p_start = details["freq_power_start_recorder_s"]
                if edge == "start" and p_start is not None:
                    ax.axvline(p_start, color="tab:green", ls=":", label="独立有功响应起点")
                ax.set_xlim(center - .4, center + .65)
                values = display[name][(display["t_s"] >= center - .4) & (display["t_s"] <= center + .65)]
                values = values[np.isfinite(values)]
                if len(values):
                    pad = max(float(np.ptp(values)) * .12, .005 if name == "f_hz" else .01)
                    ax.set_ylim(values.min() - pad, values.max() + pad)
                ax.set_ylabel(label)
                ax.grid(True, alpha=.3)
                ax.legend(fontsize=8, loc="best")
            axes[0, col].set_title("扰动开始" if edge == "start" else "扰动恢复")
            axes[1, col].set_xlabel("原始录波时间 / s")
        readable = "；".join(REASON_LABELS.get(reason, reason) for reason in details['freq_timing_reason'].split("; "))
        fig.suptitle(f"Case {meta['case_index']}：对时待人工确认\n{readable}", fontsize=10)
        fig.tight_layout(rect=(0, 0, 1, .93))
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=150)
    finally:
        plt.close(fig)
