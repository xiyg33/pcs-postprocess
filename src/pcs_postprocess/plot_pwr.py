import argparse
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np

from .processing import safe_num, is_freq_reg_deadband, supports_power_current
from .plot_common import (
    _draw_power_current_figure, pwr_power_axis_limits,
    PWR_POWER_TICK_STEP_PU, freq_reg_display_limits,
    FREQ_REG_DISPLAY_TICK_HZ,
)

FIG_SIZE = (12, 10)
DPI = 150
COLORS_ABC = ["tab:blue", "tab:orange", "tab:green"]
F_NOMINAL = 50.0


def render_case(case, fig_dir):
    """Draw the existing mode layout from prepared arrays."""
    stem, meta, ev = case.stem, case.meta, case.events
    data = case.bridged
    tw = data["t_s"]
    power = case.pwr_downsampled
    power_t = power["t_s"]
    Va_pu, Vb_pu, Vc_pu = data["Vabc_pu"].T
    Ia_pu, Ib_pu, Ic_pu = data["Iabc_pu"].T
    P_pu = power["P_pu"]
    Q_pu = power["Q_pu"]
    f_wm = data.get("f_hz")
    Vp_pu_wm = data["Vp_pu"]
    P_base_pu = safe_num(meta.get("P_base_pu"))
    P_step_pu = safe_num(meta.get("P_step_pu"))
    Q_base_pu = safe_num(meta.get("Q_base_pu"))
    Q_step_pu = safe_num(meta.get("Q_step_pu"))

    fig, axes = plt.subplots(4, 1, figsize=FIG_SIZE, sharex=True)

    # ---- Panel 1: 三相电压 ----
    ax = axes[0]
    for j, (y, lab) in enumerate(zip(
        [Va_pu, Vb_pu, Vc_pu], ["Va", "Vb", "Vc"]
    )):
        ax.plot(tw, y, lw=0.4, color=COLORS_ABC[j], label=lab)
    ax.set_ylabel("Vabc (pu)")
    ax.legend(loc="upper right", fontsize=8, ncol=3)

    # ---- Panel 2: 三相电流 ----
    ax = axes[1]
    for j, (y, lab) in enumerate(zip(
        [Ia_pu, Ib_pu, Ic_pu], ["Ia", "Ib", "Ic"]
    )):
        ax.plot(tw, y, lw=0.4, color=COLORS_ABC[j], label=lab)
    ax.set_ylabel("Iabc (pu)")
    ax.legend(loc="upper right", fontsize=8, ncol=3)

    # ---- Panel 3: 有功功率 ----
    ax = axes[2]
    ax.plot(power_t, P_pu, lw=0.8, color="tab:blue")
    if not math.isnan(P_base_pu):
        ax.axhline(P_base_pu, ls="--", color="gray", lw=0.8)
    if not math.isnan(P_step_pu):
        ax.axhline(P_step_pu, ls=":", color="gray", lw=0.8)
    ax.set_ylabel("P (pu)")
    ax.set_ylim(pwr_power_axis_limits(P_pu, meta, "P"))
    ax.yaxis.set_major_locator(MultipleLocator(PWR_POWER_TICK_STEP_PU))

    # ---- Panel 4: 无功功率 ----
    ax = axes[3]
    ax.plot(power_t, Q_pu, lw=0.8, color="tab:orange")
    if not math.isnan(Q_base_pu):
        ax.axhline(Q_base_pu, ls="--", color="gray", lw=0.8)
    if not math.isnan(Q_step_pu):
        ax.axhline(Q_step_pu, ls=":", color="gray", lw=0.8)
    ax.set_ylabel("Q (pu)")
    ax.set_xlabel("sim time (s)")
    ax.set_ylim(pwr_power_axis_limits(Q_pu, meta, "Q"))
    ax.yaxis.set_major_locator(MultipleLocator(PWR_POWER_TICK_STEP_PU))

    # ---- 公共: 事件标记 / 网格 / 标题 ----
    ev_marks = [ev["event1_start"], ev["event1_end"]]
    if ev["event2_start"] is not None:
        ev_marks += [ev["event2_start"], ev["event2_end"]]
    for ax in axes:
        for tm in ev_marks:
            ax.axvline(tm, color="red", ls="--", lw=0.8)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(left=0.0)

    ci = int(safe_num(meta.get("case_index")))
    cn = str(meta.get("case_name", stem))
    fig.suptitle(
        "Case %04d: %s  P_base=%.2f pu  P_step=%.2f pu  Q_base=%.2f pu  Q_step=%.2f pu"
        % (ci, cn, P_base_pu, P_step_pu, Q_base_pu, Q_step_pu),
        fontsize=11,
    )

    fig.tight_layout()
    out_path = os.path.join(fig_dir, stem + ".png")
    fig.savefig(out_path, dpi=DPI)
    plt.close(fig)
    return stem
