"""SCR plot for a prepared canonical case."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from .processing import safe_num

FIG_SIZE = (12, 9)
DPI = 150
COLORS_ABC = ["tab:blue", "tab:orange", "tab:green"]

def render_case(case, fig_dir):
    """Render SCR annotations using aligned prepared arrays."""
    stem, meta, ev = case.stem, case.meta, case.events
    data = case.bridged
    tw = data["t_s"]
    Va_pu, Vb_pu, Vc_pu = data["Vabc_pu"].T
    Ia_pu, Ib_pu, Ic_pu = data["Iabc_pu"].T
    P_pu = data["P_pu"]

    # 5. Breaker 信息
    breaker = data.get("breaker")
    has_breaker_trip = False
    breaker_trip_t = None
    if breaker is not None:
        brk_w = breaker
        brk_chg = np.where(np.diff(brk_w) != 0)[0]
        if len(brk_chg):
            has_breaker_trip = True
            breaker_trip_t = tw[brk_chg[0]]

    # 6. 创建图
    fig, axes = plt.subplots(3, 1, figsize=FIG_SIZE, sharex=True)

    # ---- Panel 1: 三相电压 ----
    ax = axes[0]
    for j, (y, lab) in enumerate(zip(
        [Va_pu, Vb_pu, Vc_pu], ["Va", "Vb", "Vc"]
    )):
        ax.plot(tw, y, lw=0.4, color=COLORS_ABC[j], label=lab)
    ax.set_ylabel("Vabc (pu)")
    ax.set_ylim(-1.5, 1.5)
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
    ax.plot(tw, P_pu, lw=0.8, color="tab:blue")
    P_ref = safe_num(meta.get("P_ref_pu"))
    ax.axhline(P_ref, ls="--", color="gray", lw=0.8)
    ax.set_ylabel("P (pu)")
    ax.set_xlabel("normalized time (s)")

    # ---- 事件标记 ----
    for ax in axes:
        ax.axvline(ev["event1_start"], color="red", ls="--", lw=0.8)
        ax.axvline(ev["event1_end"], color="red", ls="--", lw=0.8)
        if has_breaker_trip and breaker_trip_t is not None:
            ax.axvline(breaker_trip_t, color="red", ls="-", lw=1.2, alpha=0.7)
        ax.axvspan(ev["event1_start"], ev["event1_end"], alpha=0.08, color="red")
        ax.grid(True, alpha=0.3)
        ax.set_xlim(left=0.0)

    # ---- Breaker trip 标注 ----
    if has_breaker_trip and breaker_trip_t is not None:
        axes[0].text(breaker_trip_t, axes[0].get_ylim()[1] * 0.9,
                     "breaker trip", fontsize=8, color="red",
                     rotation=90, ha="right", va="top")

    # ---- 标题 ----
    ci_value = meta.get("case_index")
    cn = str(meta.get("case_name", stem))
    scr = safe_num(meta.get("SCR_value"))
    method = str(meta.get("hil_rec_start_method", "?"))
    stable_str = "UNSTABLE" if has_breaker_trip else "stable"
    if ci_value is None:
        fault_t = safe_num(meta.get("manual_fault_time_s"), float("nan"))
        title = "%s  SCR=%.1f  [%s, fault=%.3f, event_x=%.1f, %s]" % (
            cn, scr, method, fault_t, ev["event1_start"], stable_str,
        )
    else:
        title = "Case %04d: %s  SCR=%.1f  [%s, event_x=%.1f, %s]" % (
            int(safe_num(ci_value)), cn, scr, method,
            ev["event1_start"], stable_str,
        )
    fig.suptitle(title, fontsize=10)

    fig.tight_layout()
    out_path = os.path.join(fig_dir, stem + ".png")
    fig.savefig(out_path, dpi=DPI)
    plt.close(fig)
    return stem




_PROFILE_NAME = "scr"


def _block_voltage_envelope(t_rel, sig):
    """Return 50 ms block positive-sequence voltage envelope in pu."""
    phases = np.column_stack([sig["Va"], sig["Vb"], sig["Vc"]])
    dt = float(np.median(np.diff(t_rel)))
    block = max(8, int(round(0.05 / dt)))
    n = phases.shape[0] // block
    if n < 10:
        return np.array([]), np.array([])
    x = phases[:n * block].reshape(n, block, 3)
    env = np.sqrt(2.0 / 3.0) * np.sqrt(
        np.mean(np.sum(x * x, axis=2), axis=1)
    ) / HIL_V_BASE_PEAK_V
    return t_rel[:n * block:block], env


def detect_manual_fault_time(t_rel, P, breaker, sig=None):
    """Find a manual SCR fault edge using P first and breaker as fallback."""
    t_rel = np.asarray(t_rel, dtype=float)
    P = np.asarray(P, dtype=float)
    if t_rel.size < 100 or P.size != t_rel.size:
        return None, "too_short"
    dt = float(np.median(np.diff(t_rel)))
    if not np.isfinite(dt) or dt <= 0.0:
        return None, "invalid_time"

    win = max(5, int(round(0.048 / dt)))
    P_sm = np.convolve(P, np.ones(win) / win, mode="valid")
    t_sm = t_rel[win // 2:win // 2 + len(P_sm)]
    # 手动录波可能含有较长的控制器启动过程；功率基线取初始稳定段的后半段。
    base_m = (t_sm >= t_rel[0] + 1.0) & (t_sm <= t_rel[0] + 1.8)
    if np.count_nonzero(base_m) < 20:
        base_m = (t_sm >= t_rel[0] + 0.7) & (t_sm <= t_rel[0] + 1.5)
    if not np.any(base_m):
        return None, "no_baseline"
    base_P = float(np.median(P_sm[base_m]))
    if abs(base_P) < 100.0:
        return None, "low_power"

    dev = np.abs((P_sm - base_P) / abs(base_P))
    abnormal = dev > 0.12
    search_start = t_rel[0] + 1.5
    hold = max(2, int(round(0.10 / dt)))
    for i in np.flatnonzero((t_sm >= search_start) & abnormal):
        if i + hold <= len(abnormal) and np.all(abnormal[i:i + hold]):
            return float(t_sm[i]), "P_drop"

    if breaker is not None:
        brk = np.asarray(breaker, dtype=float)
        if brk.size == t_rel.size:
            candidates = np.flatnonzero(np.diff(brk) != 0) + 1
            candidates = candidates[t_rel[candidates] > search_start]
            if candidates.size:
                return float(t_rel[candidates[0]] - 0.04), "breaker"

    # 某些手动录波没有断路器动作，也没有清晰功率沿；只在启动后的预期事件窗搜索，
    # 并标记这是回退判断，避免误认为检测到了可靠事件沿。
    if sig is not None:
        block_t, envelope = _block_voltage_envelope(t_rel, sig)
        if block_t.size:
            base_m = (block_t >= t_rel[0] + 5.0) & (block_t <= t_rel[0] + 8.0)
            search_m = block_t >= t_rel[0] + 8.0
            if np.any(base_m) and np.any(search_m):
                base_v = float(np.median(envelope[base_m]))
                noise_v = float(np.median(np.abs(envelope[base_m] - base_v)))
                threshold = max(0.05, 8.0 * noise_v)
                hit = search_m & (envelope < base_v - threshold)
                hold = max(2, int(round(0.10 / 0.05)))
                for i in np.flatnonzero(hit):
                    if i + hold <= len(hit) and np.all(hit[i:i + hold]):
                        return float(block_t[i]), "voltage_dip_fallback"
    return None, "none"


def parse_manual_name(fname):
    """Parse the filename labels used by the current manual recordings."""
    match = re.search(r"SCR(\d+(?:dot\d+)?)", fname, flags=re.IGNORECASE)
    if not match:
        return None
    scr = float(match.group(1).replace("dot", "."))
    lower = fname.lower()
    if "pm100" in lower or "pm090" in lower:
        p_ref = -1.0 if "pm100" in lower else -0.9
    elif re.search(r"(?:_|^)p100(?:_|\.)", lower):
        p_ref = 1.0
    else:
        return None
    return p_ref, scr


def build_manual_raw_meta(fname, p_ref, scr, csv_row=None):
    """Build plotting metadata when a raw manual file has no case JSON."""
    stem = os.path.splitext(os.path.basename(fname))[0]
    meta = dict(csv_row) if csv_row is not None else {
        "case_index": None,
        "case_name": stem,
        "test_type": "SCR_MANUAL",
        "test_type_code": 14,
        "P_ref_pu": p_ref,
        "Q_ref_pu": 0.0,
        "SCR_value": scr,
        "disturbance_start_s": FAULT_SIM_TIME_S,
        "disturbance_duration_s": FAULT_DUR_S,
        "record_pre_s": 2.0,
        "record_post_s": 8.0,
        "hil_event1_start_eff_s": FAULT_SIM_TIME_S,
    }
    meta["P_ref_pu"] = p_ref
    meta["SCR_value"] = scr
    meta["hil_event1_start_eff_s"] = FAULT_SIM_TIME_S
    meta["record_pre_s"] = safe_num(meta.get("record_pre_s"), 2.0)
    meta["record_post_s"] = safe_num(meta.get("record_post_s"), 8.0)
    meta["disturbance_duration_s"] = safe_num(
        meta.get("disturbance_duration_s"), FAULT_DUR_S,
    )
    return meta


def load_scr_case_rows():
    """Load the current SCR table for raw manual filename matching."""
    csv_path = os.path.join(_PARENT_DIR, "csv_file", "scr_cases.csv")
    if not os.path.exists(csv_path):
        return []
    with open(csv_path, "r", encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def match_scr_case_row(rows, p_ref, scr):
    for row in rows:
        if (abs(safe_num(row.get("P_ref_pu")) - p_ref) <= 1e-9
                and abs(safe_num(row.get("SCR_value")) - scr) <= 1e-9):
            return row
    return None
