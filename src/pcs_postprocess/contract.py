"""校验伙伴提交的 MAT/JSON 录波是否符合公开输入契约。"""

import json
import math
from pathlib import Path
import re

import numpy as np
from scipy.io import loadmat

from .processing import _PROFILES

REQUIRED_SIGNALS = ("t_s", "Va_V", "Vb_V", "Vc_V", "Ia_A", "Ib_A", "Ic_A", "P_W")
FREQUENCY_MODES = frozenset(("freq_reg", "inertia", "freq_prot"))
MODE_FIELDS = {
    "frt": ("V_fault_pu", "P_ref_pu"),
    "pwr": ("P_base_pu", "P_step_pu", "Q_base_pu", "Q_step_pu"),
    "freq_reg": ("P_ref_pu", "freq_deviation_hz"),
    "inertia": ("P_ref_pu", "freq_deviation_hz", "freq_ramp_rate_hz_per_s"),
    "volt_prot": ("P_ref_pu", "V_target_pu"),
    "freq_prot": ("P_ref_pu", "f_target_hz"),
    "scr": ("P_ref_pu", "SCR_value"),
}


def _header_names(header):
    names = []
    for row in header:
        cell = row[0]
        names.append(str(cell[0]) if getattr(cell, "size", 0) else "")
    return names


def validate_case(path: Path, expected_mode: str | None = None) -> dict:
    """先检查事件元数据，再检查 MAT 列名、长度和时间轴。"""
    meta_path = path.with_name(path.stem + "_meta.json")
    if not meta_path.is_file():
        raise ValueError(f"Missing metadata: {meta_path}")
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    for field in ("mode", "case_index", "case_name", "time_basis"):
        if field not in meta:
            raise ValueError(f"{meta_path.name}: missing {field}")
    mode = meta["mode"]
    if mode not in _PROFILES or expected_mode and mode != expected_mode:
        raise ValueError(f"{meta_path.name}: unexpected mode {mode!r}")
    if meta["time_basis"] not in ("recorded", "absolute"):
        raise ValueError(f"{meta_path.name}: time_basis must be recorded or absolute")
    if not isinstance(meta["case_index"], int) or meta["case_index"] < 0:
        raise ValueError(f"{meta_path.name}: case_index must be a nonnegative integer")
    if not isinstance(meta["case_name"], str) or not meta["case_name"]:
        raise ValueError(f"{meta_path.name}: case_name must be nonempty")
    if any(mark in meta["case_name"] for mark in ("/", "\\")):
        raise ValueError(f"{meta_path.name}: case_name cannot contain path separators")
    match = re.fullmatch(r"case_(\d+)_(.+)", path.stem)
    if not match or int(match.group(1)) != meta["case_index"]:
        raise ValueError(f"{path.name}: filename and case_index disagree")
    # 各模式的事件字段由处理器 profile 定义，避免校验与计算要求分叉。
    profile = _PROFILES[mode]
    for field in (profile["event1_field"], profile["event1_dur_field"],
                  *MODE_FIELDS[mode]):
        value = meta.get(field)
        if not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError(f"{meta_path.name}: missing or invalid {field}")
        if field in (profile["event1_field"], profile["event1_dur_field"]) and value < 0:
            raise ValueError(f"{meta_path.name}: {field} must be nonnegative")
    with np.errstate(all="ignore"):
        mat = loadmat(path)
    if "data" not in mat or "header" not in mat:
        raise ValueError(f"{path.name}: expected data and header arrays")
    data = mat["data"]
    names = _header_names(mat["header"])
    if data.ndim != 2 or data.shape[0] < 2 or data.shape[1] != len(names):
        raise ValueError(f"{path.name}: data must be N x len(header), N >= 2")
    # 仅强制要求计算该模式所必需的通道；其他通道可作为可选绘图信息。
    required = set(REQUIRED_SIGNALS)
    if mode in FREQUENCY_MODES:
        required.add("f_hz")
    if mode == "pwr":
        required.add("Q_var")
    missing = required - set(names)
    if missing:
        raise ValueError(f"{path.name}: missing signals {sorted(missing)}")
    time = data[:, names.index("t_s")]
    if not np.all(np.isfinite(time)) or np.any(np.diff(time) <= 0):
        raise ValueError(f"{path.name}: t_s must be finite and strictly increasing")
    if not np.all(np.isfinite(data[:, [names.index(n) for n in required]])):
        raise ValueError(f"{path.name}: required signals must be finite")
    return meta
