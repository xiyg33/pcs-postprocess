"""生成七类完全合成、符合公开契约的小型录波。"""

import json
from pathlib import Path

import numpy as np
from scipy.io import savemat

NAMES = ("t_s", "Va_V", "Vb_V", "Vc_V", "Ia_A", "Ib_A", "Ic_A",
         "P_W", "Q_var", "f_hz", "Vpos_pu", "breaker")
MODES = ("frt", "pwr", "freq_reg", "inertia", "volt_prot", "freq_prot", "scr")
FIELDS = {
    "frt": {"fault_start_s": 2.0, "fault_duration_s": 0.2,
            "V_fault_pu": 0.5, "P_ref_pu": 0.3, "Q_ref_pu": 0},
    "pwr": {"step_start_s": 2.0, "step_dur_s": 1.0,
            "P_base_pu": 0.3, "P_step_pu": 0.5,
            "Q_base_pu": 0.0, "Q_step_pu": 0.0},
    "freq_reg": {"event_start_s": 2.0, "event_hold_s": 0.5,
                 "P_ref_pu": 0.3, "freq_deviation_hz": -0.2},
    "inertia": {"event_start_s": 2.0, "event_hold_s": 0.5,
                "P_ref_pu": 0.3, "freq_deviation_hz": -0.2,
                "freq_ramp_rate_hz_per_s": -2.0},
    "volt_prot": {"event_start_s": 2.0, "event_hold_s": 0.5,
                  "P_ref_pu": 0.3, "V_target_pu": 0.7},
    "freq_prot": {"event_start_s": 2.0, "event_hold_s": 0.5,
                  "P_ref_pu": 0.3, "f_target_hz": 49.5},
    "scr": {"disturbance_start_s": 2.0, "disturbance_duration_s": 0.5,
            "P_ref_pu": 0.3, "SCR_value": 3.0},
}


def create(output: Path, mode: str, time_basis: str = "absolute") -> Path:
    output.mkdir(parents=True, exist_ok=True)
    index = 1000 + MODES.index(mode)
    t = np.arange(0, 5, 0.002)
    event = (t >= 2) & (t < 2.5)
    voltage_pu = np.where(event, 0.5 if mode == "frt" else 0.7, 1.0) if mode in ("frt", "volt_prot") else np.ones_like(t)
    angles = 2 * np.pi * 50 * t[:, None] + np.array([0, -2 * np.pi / 3, 2 * np.pi / 3])
    voltage = 310.269 * voltage_pu[:, None] * np.cos(angles)
    current = 51.57 * 0.3 * np.cos(angles)
    power = np.where(t >= 2, 12000.0, 7200.0) if mode == "pwr" else np.full_like(t, 7200.0)
    freq = np.where(event, 49.5 if mode == "freq_prot" else 49.8, 50.0) if mode in ("freq_reg", "inertia", "freq_prot") else np.full_like(t, 50.0)
    breaker = np.where(t >= 2.3, 0.0, 1.0) if mode in ("volt_prot", "freq_prot") else np.ones_like(t)
    data = np.column_stack((t, voltage, current, power, np.zeros_like(t), freq, voltage_pu, breaker))
    header = np.empty((len(NAMES), 2), dtype=object)
    for row, name in enumerate(NAMES):
        header[row] = (name, "")
    stem = f"case_{index:04d}_{mode}"
    path = output / f"{stem}.mat"
    savemat(path, {"data": data, "header": header}, do_compression=True)
    meta = {"mode": mode, "case_index": index, "case_name": f"synthetic_{mode}",
            "time_basis": time_basis, "record_pre_s": 2.0, "record_post_s": 2.0,
            "sim_time_s": 5.0, **FIELDS[mode]}
    (output / f"{stem}_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    return path


if __name__ == "__main__":
    root = Path(__file__).resolve().parent / "generated"
    for test_mode in MODES:
        print(create(root / test_mode, test_mode))
