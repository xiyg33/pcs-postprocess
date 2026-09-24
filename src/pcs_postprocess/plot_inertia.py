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

FIG_SIZE = (12, 9)
DPI = 150
COLORS_ABC = ["tab:blue", "tab:orange", "tab:green"]
F_NOMINAL = 50.0


def render_case(case, fig_dir):
    """Draw the existing mode layout from prepared arrays."""
    stem, meta, ev = case.stem, case.meta, case.events
    data = case.bridged
    tw = data["t_s"]
    Va_pu, Vb_pu, Vc_pu = data["Vabc_pu"].T
    Ia_pu, Ib_pu, Ic_pu = data["Iabc_pu"].T
    P_pu = data["P_pu"]
    Q_pu = data.get("Q_pu")
    f_wm = data.get("f_hz")
    Vp_pu_wm = data["Vp_pu"]
    dev = safe_num(meta.get("freq_deviation_hz"))
    rate = safe_num(meta.get("freq_ramp_rate_hz_per_s"))
    ramp_dur = abs(dev / rate) if abs(rate) > 1e-12 else 0.0

    fig, axes = plt.subplots(3, 1, figsize=FIG_SIZE, sharex=True)

    # ---- Panel 1: 频率 ----
    ax = axes[0]
    if f_wm is not None:
        ax.plot(tw, f_wm, lw=0.8, color="tab:red")
    else:
        ax.text(
            0.5, 0.5, "f signal unavailable",
            transform=ax.transAxes, ha="center", va="center",
            fontsize=10, color="gray",
        )
    ax.axhline(F_NOMINAL, ls="--", color="gray", lw=0.8)
    # 预期频率斜坡参考线 (灰色虚线)
    if ramp_dur > 0:
        ramp_t = np.array([
            ev["event1_start"],
            ev["event1_start"] + ramp_dur,
        ])
        ramp_f = F_NOMINAL + np.array([0.0, dev])
        ax.plot(ramp_t, ramp_f, "--", color="gray", lw=0.6, alpha=0.7, label="expected ramp")
    ax.set_ylabel("f (Hz)")

    # ---- Panel 2: 三相电压 ----
    ax = axes[1]
    for j, (y, lab) in enumerate(zip(
        [Va_pu, Vb_pu, Vc_pu], ["Va", "Vb", "Vc"]
    )):
        ax.plot(tw, y, lw=0.4, color=COLORS_ABC[j], label=lab)
    ax.set_ylabel("Vabc (pu)")
    ax.legend(loc="upper right", fontsize=8, ncol=3)

    # ---- Panel 3: 三相电流 ----
    ax = axes[2]
    for j, (y, lab) in enumerate(zip(
        [Ia_pu, Ib_pu, Ic_pu], ["Ia", "Ib", "Ic"]
    )):
        ax.plot(tw, y, lw=0.4, color=COLORS_ABC[j], label=lab)
    ax.set_ylabel("Iabc (pu)")
    ax.set_xlabel("normalized time (s)")
    ax.legend(loc="upper right", fontsize=8, ncol=3)

    # ---- 公共: 事件标记 / 网格 / 标题 ----
    for ax in axes:
        ax.axvline(ev["event1_start"], color="red", ls="--", lw=0.8)
        ax.axvline(ev["event1_end"], color="red", ls="--", lw=0.8)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(left=0.0)

    ci = int(safe_num(meta.get("case_index")))
    cn = str(meta.get("case_name", stem))
    P_ref = safe_num(meta.get("P_ref_pu"))
    ramp_rate = safe_num(meta.get("freq_ramp_rate_hz_per_s"))
    fig.suptitle(
        "Case %04d: %s  P=%.2f pu  df/dt=%+.3f Hz/s" % (ci, cn, P_ref, ramp_rate),
        fontsize=11,
    )

    fig.tight_layout()
    out_path = os.path.join(fig_dir, stem + ".png")
    fig.savefig(out_path, dpi=DPI)
    plt.close(fig)
    _draw_power_current_figure("inertia", tw, data, meta, ev, stem,
                               data["rec_offset_s"], case.summary["pre_event_ok"],
                               fig_dir)
    return stem
