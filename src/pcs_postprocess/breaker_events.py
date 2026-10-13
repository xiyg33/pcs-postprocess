"""Native-clock breaker evidence, without smoothing binary states."""

import numpy as np


def evidence(t, breaker):
    t = np.asarray(t, float)
    result = {"status": "valid", "reason": "", "falls": [], "incomplete_falls": []}
    if breaker is None:
        return dict(result, status="missing_breaker", reason="缺少 breaker 通道")
    b = np.asarray(breaker, float)
    if (b.shape != t.shape or not np.isfinite(b).all()
            or not np.all(np.isclose(b, 0.) | np.isclose(b, 1.))):
        return dict(result, status="invalid_breaker", reason="breaker 必须是有限的 0/1 状态")
    dt = float(np.median(np.diff(t)))
    for i in np.flatnonzero((b[:-1] >= .5) & (b[1:] < .5)) + 1:
        high_start = max(0, i - 10)
        # A short initial high pulse does not establish a connected condition.
        if i - high_start < 10 or not np.all(b[high_start:i] >= .5):
            continue
        if (i + 10 <= len(t) and np.all(b[i:i + 10] < .5)
                and np.all(np.diff(t[high_start:i + 10]) <= 3 * dt)):
            result["falls"].append(float(t[i]))
        elif np.all(b[i:min(i + 10, len(b))] < .5):
            result["incomplete_falls"].append(float(t[i]))
    return result
