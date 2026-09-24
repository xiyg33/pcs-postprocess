import math
import os
import traceback
from copy import deepcopy

import matplotlib
matplotlib.use("Agg")
import numpy as np

from .processing import safe_num, supports_power_current


def draw_processed_case(case, fig_dir, include_power_current=False):
    """Render only prepared arrays; the caller chooses response figures."""
    data = case.bridged
    profile = deepcopy(PLOT_PROFILES[case.profile_name])
    if case.profile_name == "freq_reg" and not supports_power_current(case.meta, "freq_reg"):
        profile["fig_panels"][2] = {
            "type": "P_pu", "ylabel": "P (pu)", "xlabel": True,
            "ref_lines": [{"field": "P_ref_pu"}],
        }
    _draw_figure(profile, data["t_s"], data, case.meta, data["P_pu"],
                 data.get("Q_pu", np.full_like(data["t_s"], np.nan)),
                 data["Vp_pu"], data.get("f_hz"), case.events, case.stem,
                 data["rec_offset_s"], case.summary["pre_event_ok"], fig_dir)
    if include_power_current:
        _draw_power_current_figure(case.profile_name, data["t_s"], data,
                                   case.meta, case.events, case.stem,
                                   data["rec_offset_s"], case.summary["pre_event_ok"], fig_dir)
    return case.stem


PWR_POWER_TICK_STEP_PU = 0.2
PWR_POWER_MIN_SPAN_PU = 0.8
FREQ_REG_DISPLAY_LIMITS_HZ = (49.4, 50.6)
FREQ_REG_DISPLAY_TICK_HZ = 0.2
NOMINAL_FREQUENCY_HZ = 50.0


def configure_nominal_frequency(nominal_hz):
    """Adjust frequency references in copied plot profiles for each project."""
    global FREQ_REG_DISPLAY_LIMITS_HZ, NOMINAL_FREQUENCY_HZ
    previous_hz = NOMINAL_FREQUENCY_HZ
    NOMINAL_FREQUENCY_HZ = nominal_hz
    FREQ_REG_DISPLAY_LIMITS_HZ = (nominal_hz - 0.6, nominal_hz + 0.6)

    def visit(item):
        if isinstance(item, dict):
            for key, value in item.items():
                if key == "value" and value == previous_hz:
                    item[key] = nominal_hz
                else:
                    visit(value)
        elif isinstance(item, list):
            for value in item:
                visit(value)

    visit(PLOT_PROFILES)

PLOT_PROFILES = {'frt': {'fig_layout': (3, 1),
         'fig_figsize': (12, 9),
         'fig_panels': [{'type': 'Vp',
                         'ylabel': 'V+ (pu)',
                         'ylim': (-0.1, 1.3),
                         'ref_lines': [{'field': 'V_fault_pu', 'ls': ':', 'color': 'gray'}]},
                        {'type': 'Vabc_pu', 'ylabel': 'Vabc (pu)', 'legend': True},
                        {'type': 'Iabc_pu',
                         'ylabel': 'Iabc (pu)',
                         'xlabel': True,
                         'legend': True}]},
 'pwr': {'fig_layout': (4, 1),
         'fig_figsize': (12, 10),
         'fig_panels': [{'type': 'Vabc_pu', 'ylabel': 'Vabc (pu)', 'legend': True},
                        {'type': 'Iabc_pu', 'ylabel': 'Iabc (pu)', 'legend': True},
                        {'type': 'P_pu',
                         'ylabel': 'P (pu)',
                         'power_axis': 'P',
                         'ref_lines': [{'field': 'P_base_pu', 'ls': '--', 'color': 'gray'},
                                       {'field': 'P_step_pu', 'ls': ':', 'color': 'gray'}]},
                        {'type': 'Q_pu',
                         'ylabel': 'Q (pu)',
                         'xlabel': True,
                         'power_axis': 'Q',
                         'ref_lines': [{'field': 'Q_base_pu', 'ls': '--', 'color': 'gray'},
                                       {'field': 'Q_step_pu', 'ls': ':', 'color': 'gray'}]}]},
 'freq_reg': {'fig_layout': (3, 1),
              'fig_figsize': (12, 9),
              'fig_panels': [{'type': 'f',
                              'ylabel': 'f (Hz)',
                              'ref_lines': [{'value': 50.0, 'ls': '--', 'color': 'gray'}]},
                             {'type': 'Vabc_pu', 'ylabel': 'Vabc (pu)', 'legend': True},
                             {'type': 'Iabc_pu',
                              'ylabel': 'Iabc (pu)',
                              'xlabel': True,
                              'legend': True}]},
 'inertia': {'fig_layout': (4, 1),
             'fig_figsize': (12, 10),
             'fig_panels': [{'type': 'Vp', 'ylabel': 'V+ (pu)', 'ylim': (0.8, 1.2)},
                            {'type': 'P_pu',
                             'ylabel': 'P (pu)',
                             'ref_lines': [{'field': 'P_ref_pu',
                                            'ls': '--',
                                            'color': 'gray'}]},
                            {'type': 'Q_pu', 'ylabel': 'Q (pu)'},
                            {'type': 'f',
                             'ylabel': 'f (Hz)',
                             'xlabel': True,
                             'ref_lines': [{'value': 50.0, 'ls': '--', 'color': 'gray'}]}]},
 'volt_prot': {'fig_layout': (3, 1),
               'fig_figsize': (12, 9),
               'fig_panels': [{'type': 'Vp',
                               'ylabel': 'V+ (pu)',
                               'include_ylim': (0.0, 1.0),
                               'ref_lines': [{'field': 'V_target_pu',
                                              'ls': ':',
                                              'color': 'gray'}]},
                              {'type': 'P_pu',
                               'ylabel': 'P (pu)',
                               'ref_lines': [{'field': 'P_ref_pu',
                                              'ls': '--',
                                              'color': 'gray'}]},
                              {'type': 'Q_pu', 'ylabel': 'Q (pu)', 'xlabel': True}]},
 'freq_prot': {'fig_layout': (3, 1),
               'fig_figsize': (12, 9),
               'fig_panels': [{'type': 'f',
                               'ylabel': 'f (Hz)',
                               'ref_lines': [{'value': 50.0, 'ls': '--', 'color': 'gray'},
                                             {'field': 'f_target_hz',
                                              'ls': ':',
                                              'color': 'gray'}]},
                              {'type': 'P_pu',
                               'ylabel': 'P (pu)',
                               'ref_lines': [{'field': 'P_ref_pu',
                                              'ls': '--',
                                              'color': 'gray'}]},
                              {'type': 'Q_pu', 'ylabel': 'Q (pu)', 'xlabel': True}]},
 'scr': {'fig_layout': (4, 1),
         'fig_figsize': (12, 10),
         'fig_panels': [{'type': 'Vp', 'ylabel': 'V+ (pu)'},
                        {'type': 'P_pu',
                         'ylabel': 'P (pu)',
                         'ref_lines': [{'field': 'P_ref_pu', 'ls': '--', 'color': 'gray'}]},
                        {'type': 'Q_pu', 'ylabel': 'Q (pu)'},
                        {'type': 'f',
                         'ylabel': 'f (Hz)',
                         'xlabel': True,
                         'ref_lines': [{'value': 50.0, 'ls': '--', 'color': 'gray'}]}]}}

def freq_reg_display_limits(frequency):
    """Broad common limits, expanding only when finite samples exceed them."""
    lower, upper = FREQ_REG_DISPLAY_LIMITS_HZ
    values = np.asarray(frequency, dtype=float)
    values = values[np.isfinite(values)]
    if values.size:
        if values.min() < lower:
            lower = math.floor((values.min() - 0.05) / FREQ_REG_DISPLAY_TICK_HZ) * FREQ_REG_DISPLAY_TICK_HZ
        if values.max() > upper:
            upper = math.ceil((values.max() + 0.05) / FREQ_REG_DISPLAY_TICK_HZ) * FREQ_REG_DISPLAY_TICK_HZ
    return lower, upper


def pwr_power_axis_limits(signal_pu, meta, channel):
    """Return a centered PWR power-axis range with 0.2 pu grid spacing.

    The center follows all configured steady-state setpoints for the selected
    channel.  The range is at least 0.8 pu wide and expands symmetrically when
    measured transients extend beyond that window.
    """
    channel = str(channel).upper()
    if channel not in ("P", "Q"):
        raise ValueError("channel must be 'P' or 'Q'")

    setpoints = []
    for suffix in ("base", "step", "step2"):
        value = safe_num(meta.get("%s_%s_pu" % (channel, suffix)), float("nan"))
        if math.isfinite(value):
            setpoints.append(value)

    data = np.asarray(signal_pu, dtype=float)
    finite_data = data[np.isfinite(data)]
    if setpoints:
        center = 0.5 * (min(setpoints) + max(setpoints))
    elif finite_data.size:
        center = 0.5 * (float(np.min(finite_data)) + float(np.max(finite_data)))
    else:
        center = 0.0

    required_half_span = PWR_POWER_MIN_SPAN_PU / 2.0
    if setpoints:
        required_half_span = max(
            required_half_span,
            max(abs(value - center) for value in setpoints),
        )
    if finite_data.size:
        required_half_span = max(
            required_half_span,
            abs(float(np.min(finite_data)) - center),
            abs(float(np.max(finite_data)) - center),
        )

    half_span = (
        math.ceil((required_half_span - 1e-12) / PWR_POWER_TICK_STEP_PU)
        * PWR_POWER_TICK_STEP_PU
    )
    return center - half_span, center + half_span


def _draw_figure(prof, tw, sig_wm, meta, P_pu_wm, Q_pu_wm, Vp_pu_wm, f_sig_wm,
                 ev, stem, offset, pre_event_ok, fig_dir, suffix=""):
    """根据 profile 的 fig_panels 配置绘制多子图。

    每个 panel 的 type 决定绘制哪种信号及其样式。
    """
    import matplotlib.pyplot as plt

    nrows, ncols = prof["fig_layout"]
    figsize = prof.get("fig_figsize", (12, 3 * nrows))
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True)
    if nrows == 1 and ncols == 1:
        axes = [axes]
    axes = np.atleast_1d(axes).flatten()

    for i, panel in enumerate(prof["fig_panels"]):
        ax = axes[i]
        ptype = panel["type"]

        # ---- 绘制信号 ----
        if ptype == "Vp":
            ax.plot(tw, Vp_pu_wm, lw=0.8, color="tab:green")
        elif ptype == "Vabc_pu":
            colors = ["tab:blue", "tab:orange", "tab:green"]
            for j, lab in enumerate(["Va", "Vb", "Vc"]):
                ax.plot(tw, sig_wm["Vabc_pu"][:, j], lw=0.4, color=colors[j], label=lab)
            if panel.get("legend"):
                ax.legend(loc="upper right", fontsize=8, ncol=3)
        elif ptype == "Iabc_pu":
            colors = ["tab:blue", "tab:orange", "tab:green"]
            for j, lab in enumerate(["Ia", "Ib", "Ic"]):
                ax.plot(tw, sig_wm["Iabc_pu"][:, j], lw=0.4, color=colors[j], label=lab)
            if panel.get("legend"):
                ax.legend(loc="upper right", fontsize=8, ncol=3)
        elif ptype == "P_pu":
            ax.plot(tw, P_pu_wm, lw=0.8, color="tab:blue")
        elif ptype == "Q_pu":
            ax.plot(tw, Q_pu_wm, lw=0.8, color="tab:orange")
        elif ptype in ("Ip_pu", "Iq_pu"):
            color = "tab:blue" if ptype == "Ip_pu" else "tab:orange"
            ax.plot(tw, sig_wm[ptype], lw=0.8, color=color)
            invalid = ~sig_wm["current_valid"].astype(bool)
            if invalid.any():
                ax.fill_between(tw, 0, 1, where=invalid, step="mid",
                                transform=ax.get_xaxis_transform(), color="gray", alpha=0.15)
                ax.text(0.01, 0.04, "Invalid direction: |V| < 0.05 pu or non-finite input",
                        transform=ax.transAxes, fontsize=8, color="gray")
        elif ptype == "f":
            if f_sig_wm is not None:
                ax.plot(tw, f_sig_wm, lw=0.8, color="tab:red")
            else:
                ax.text(0.5, 0.5, "f signal unavailable", transform=ax.transAxes,
                        ha="center", va="center", fontsize=10, color="gray")

        # ---- 水平参考线 ----
        for rl in panel.get("ref_lines", []):
            if "value" in rl:
                y = rl["value"]
            elif "field" in rl:
                y = safe_num(meta.get(rl["field"]))
                if math.isnan(y):
                    continue
            else:
                continue
            ax.axhline(y, ls=rl.get("ls", "--"), color=rl.get("color", "gray"), lw=0.8)

        # ---- 标签 ----
        ax.set_ylabel(panel.get("ylabel", ""))
        if panel.get("response_active_axis"):
            values = P_pu_wm if ptype == "P_pu" else sig_wm["Ip_pu"]
            # Response curves have no PWR setpoints: center on their data range.
            ax.set_ylim(pwr_power_axis_limits(values, {}, "P"))
            from matplotlib.ticker import MultipleLocator
            ax.yaxis.set_major_locator(MultipleLocator(PWR_POWER_TICK_STEP_PU))
        elif panel.get("power_axis"):
            power_signal = P_pu_wm if panel["power_axis"] == "P" else Q_pu_wm
            ax.set_ylim(pwr_power_axis_limits(power_signal, meta, panel["power_axis"]))
            from matplotlib.ticker import MultipleLocator
            ax.yaxis.set_major_locator(MultipleLocator(PWR_POWER_TICK_STEP_PU))
        elif "ylim" in panel:
            ax.set_ylim(panel["ylim"])
        elif "include_ylim" in panel:
            ymin, ymax = ax.get_ylim()
            include_min, include_max = panel["include_ylim"]
            ax.set_ylim(min(ymin, include_min), max(ymax, include_max))
        if panel.get("ytick_step"):
            from matplotlib.ticker import MultipleLocator
            step = panel["ytick_step"]
            ymin, ymax = ax.get_ylim()
            ax.set_ylim(math.floor(ymin / step) * step,
                        math.ceil(ymax / step) * step)
            ax.yaxis.set_major_locator(MultipleLocator(step))
        if (ptype == "f" and ev.get("profile_name") == "freq_reg"
                and supports_power_current(meta, "freq_reg")):
            from matplotlib.ticker import MultipleLocator
            ax.set_ylim(freq_reg_display_limits(f_sig_wm))
            ax.yaxis.set_major_locator(MultipleLocator(FREQ_REG_DISPLAY_TICK_HZ))
        if panel.get("xlabel"):
            ax.set_xlabel("sim time (s)")

    # ---- 事件标记（所有 panel 共用） ----
    marks = ev.get("event_markers")
    if marks is None:
        marks = [ev["event1_start"], ev["event1_end"]]
        if ev["event2_start"] is not None:
            marks += [ev["event2_start"], ev["event2_end"]]
    for ax in axes:
        for tm in marks:
            ax.axvline(tm, color="red", ls="--", lw=0.8)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(left=0.0)

    axes[0].set_title(ev.get("plot_title", "%s  (offset=%.3fs, pre_ok=%s)" %
                             (stem, offset, pre_event_ok)))
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, stem + suffix + ".png"), dpi=150)
    plt.close(fig)


def _draw_power_current_figure(profile_name, tw, bridged, meta, ev, stem,
                               offset, pre_event_ok, fig_dir):
    disturbance = ({"type": "Vp", "ylabel": "V+ (pu)"} if profile_name == "frt"
                   else {"type": "f", "ylabel": "f (Hz)",
                         "ref_lines": [{"value": NOMINAL_FREQUENCY_HZ}]})
    response_profile = {
        "fig_layout": (5, 1),
        "fig_figsize": (12, 9),
        "fig_panels": [
            disturbance,
            {"type": "P_pu", "ylabel": "P (pu)", "response_active_axis": True,
             "ref_lines": [{"value": 0.0}]},
            {"type": "Ip_pu", "ylabel": "Ip active (pu)",
             "response_active_axis": True,
             "ref_lines": [{"value": 0.0}]},
            {"type": "Q_pu", "ylabel": "Q (pu)",
             "ref_lines": [{"value": 0.0}]},
            {"type": "Iq_pu", "ylabel": "Iq reactive (pu)", "xlabel": True,
             "ref_lines": [{"value": 0.0}]},
        ],
    }
    if profile_name == "freq_reg":
        # Keep a readable baseline range without clipping larger responses.
        for index in (3, 4):
            response_profile["fig_panels"][index]["include_ylim"] = (-0.4, 0.4)
    display_signals = dict(bridged, Ip_pu=bridged["Ip_filtered_pu"],
                           Iq_pu=bridged["Iq_filtered_pu"])
    _draw_figure(response_profile, tw, display_signals, meta, bridged["P_pu"],
                 bridged.get("Q_pu", np.full_like(tw, np.nan)), bridged["Vp_pu"],
                 bridged.get("f_hz"), ev, stem, offset, pre_event_ok, fig_dir,
                 suffix="_power_current")
