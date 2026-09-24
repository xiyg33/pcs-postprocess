"""Project ratings used for per-unit conversion and time alignment."""

from dataclasses import dataclass
import json
import math
from pathlib import Path


@dataclass(frozen=True)
class ProjectConfig:
    power_base_w: float
    voltage_base_peak_v: float
    current_base_peak_a: float
    nominal_frequency_hz: float
    record_start_nominal_s: float = 0.0

    @classmethod
    def load(cls, path: Path) -> "ProjectConfig":
        raw = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(raw, dict):
            raise ValueError("Project config must be a JSON object")
        allowed = set(cls.__dataclass_fields__)
        unexpected = set(raw) - allowed
        if unexpected:
            raise ValueError(f"Unknown project config keys: {sorted(unexpected)}")
        result = cls(**raw)
        for key in allowed - {"record_start_nominal_s"}:
            value = getattr(result, key)
            if not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
                raise ValueError(f"{key} must be a positive finite number")
        if not math.isfinite(result.record_start_nominal_s):
            raise ValueError("record_start_nominal_s must be finite")
        return result
