"""Optional case selection and display labels, independent of numerical inputs."""

import math
import re
from pathlib import Path

import numpy as np
from matplotlib.ticker import MultipleLocator


PRESENTATION_FIELDS = ["display_case_name", "reported_P_ref_pu", "source_P_ref_pu", "approximation"]


def select_case_ids(case_id, case_ids, default):
    if case_ids is None:
        return [case_id] if case_id is not None else sorted(default)
    if case_id is not None:
        raise ValueError("case_id and case_ids are mutually exclusive")
    ids = list(case_ids)
    if (not ids or any(type(value) is not int or value < 0 for value in ids)
            or len(ids) != len(set(ids))):
        raise ValueError("case_ids must contain unique nonnegative integers")
    return sorted(ids)


def figure_path(directory, case_id, suffix, presentation=None):
    name = (presentation or {}).get(case_id, {}).get("display_case_name")
    if name is not None and (not name or not re.fullmatch(r"[A-Za-z0-9_.-]+", name)):
        raise ValueError("display_case_name must be a safe filename component")
    prefix = f"case_{case_id:04d}" + (f"_{name}" if name else "")
    return Path(directory) / f"{prefix}_sim_hil_{suffix}.png"


def presentation_record(meta, case_id, presentation):
    if presentation is None:
        return {}
    entry = presentation.get(case_id, {})
    source = meta.get("P_ref_pu", meta.get("P_base_pu"))
    reported = entry.get("reported_P_ref_pu", source)
    return dict(display_case_name=entry.get("display_case_name", meta["case_name"]),
                reported_P_ref_pu=reported, source_P_ref_pu=source,
                approximation=source is not None and reported is not None
                and not math.isclose(float(source), float(reported), abs_tol=1e-8, rel_tol=0.))


def render_plot(renderer, *args, presentation=None):
    if presentation is None:
        return renderer(*args)
    return renderer(*args, presentation=presentation)


def apply_presentation(fig, case_id, presentation):
    """Change only the title and axes; never replace samples or reference values."""
    entry = (presentation or {}).get(case_id)
    if entry is None:
        return
    if entry.get("title"):
        previous = fig._suptitle.get_text() if fig._suptitle is not None else ""
        notes = previous.partition("\n")[2]
        fig.suptitle(entry["title"] + ("\n" + notes if notes else ""), fontsize=10)
    limits = entry.get("power_axis_limits")
    for ax in fig.axes:
        label = ax.get_ylabel().strip()
        is_p = label.startswith(("有功功率", "P (", "P /", "P/"))
        is_q = label.startswith(("无功功率", "Q (", "Q /", "Q/"))
        if not is_p and not is_q:
            continue
        samples = [np.asarray(line.get_ydata(), dtype=float).ravel() for line in ax.lines
                   if line.get_transform() != ax.get_xaxis_transform()]
        values = np.concatenate(samples) if samples else np.array([])
        values = values[np.isfinite(values)]
        if is_p and limits is not None:
            lower, upper = map(float, limits)
            tick = float(entry.get("power_tick_pu", .2))
            if not (math.isfinite(lower) and math.isfinite(upper) and lower < upper
                    and math.isfinite(tick) and tick > 0):
                raise ValueError("Invalid power axis limits/tick")
            if values.size:
                if float(values.min()) < lower:
                    lower = math.floor((float(values.min()) - .02) / tick) * tick
                if float(values.max()) > upper:
                    upper = math.ceil((float(values.max()) + .02) / tick) * tick
            ax.set_ylim(lower, upper)
            ax.yaxis.set_major_locator(MultipleLocator(tick))
        elif is_q and entry.get("broad_q_axis"):
            from .plot_common import pwr_power_axis_limits
            ax.set_ylim(pwr_power_axis_limits(values, {}, "Q"))
            ax.yaxis.set_major_locator(MultipleLocator(float(entry.get("power_tick_pu", .2))
                                                      if limits is not None else .2))
