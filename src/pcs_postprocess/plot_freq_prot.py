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
    t1 = data["t1_s"]
    t2 = data["t2_s"]
    delta_t = data["delta_t_s"]

    fig, axes = plt.subplots(3, 1, figsize=FIG_SIZE, sharex=True)

    ax = axes[0]
    if f_wm is not None:
        ax.plot(tw, f_wm, lw=0.8, color="tab:red")
    else:
        ax.text(0.5, 0.5, "f signal unavailable",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=10, color="gray")
    ax.axhline(F_NOMINAL, ls="--", color="gray", lw=0.8)
    f_target = safe_num(meta.get("f_target_hz"))
    if not np.isnan(f_target):
        ax.axhline(f_target, ls=":", color="gray", lw=0.8)
    ax.set_ylabel("f (Hz)")

    ax = axes[1]
    ax.plot(tw, P_pu, lw=0.8, color="tab:blue")
    P_ref = safe_num(meta.get("P_ref_pu"))
    if not np.isnan(P_ref):
        ax.axhline(P_ref, ls="--", color="gray", lw=0.8)
    ax.set_ylabel("P (pu)")

    ax = axes[2]
    if Q_pu is not None:
        ax.plot(tw, Q_pu, lw=0.8, color="tab:orange")
    else:
        ax.text(0.5, 0.5, "Q signal unavailable",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=10, color="gray")
    ax.set_ylabel("Q (pu)")
    ax.set_xlabel("sim time (s)")

    for ax in axes:
        ax.axvline(ev["event1_start"], color="red", ls="--", lw=0.8, alpha=0.5)
        ax.axvline(ev["event1_end"], color="red", ls="--", lw=0.8, alpha=0.5)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(left=0.0)

    _add_prot_markers(axes, tw, t1, t2)

    ci = int(safe_num(meta.get("case_index")))
    cn = str(meta.get("case_name", stem))
    ft = safe_num(meta.get("f_target_hz"))

    if not np.isnan(delta_t):
        suptitle = "Case %04d: %s  f_target=%.1f Hz  Δt=%.4f s" % (ci, cn, ft, delta_t)
    elif np.isnan(t2):
        suptitle = "Case %04d: %s  f_target=%.1f Hz  Δt=N/A (no_trip)" % (ci, cn, ft)
    else:
        suptitle = "Case %04d: %s  f_target=%.1f Hz  Δt=N/A (detection_error)" % (ci, cn, ft)

    fig.suptitle(suptitle, fontsize=11)
    fig.tight_layout()
    out_path = os.path.join(fig_dir, stem + ".png")
    fig.savefig(out_path, dpi=DPI)
    plt.close(fig)
    return stem


def _add_prot_markers(axes, tw, t1, t2):
    t_range = tw[-1] - tw[0]
    x_off = t_range * 0.008
    bbox = dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.7)

    if not np.isnan(t1):
        for ax in axes:
            ax.axvline(t1, color="red", ls="--", lw=1.5)
            ylim = ax.get_ylim()
            ax.text(t1 + x_off, ylim[1] - 0.05 * (ylim[1] - ylim[0]),
                    "t1=%.4f" % t1, color="red", fontsize=7,
                    ha="left", va="top", bbox=bbox)

    if not np.isnan(t2):
        for ax in axes:
            ax.axvline(t2, color="green", ls="--", lw=1.5)
            ylim = ax.get_ylim()
            ax.text(t2 + x_off, ylim[0] + 0.05 * (ylim[1] - ylim[0]),
                    "t2=%.4f" % t2, color="green", fontsize=7,
                    ha="left", va="bottom", bbox=bbox)
