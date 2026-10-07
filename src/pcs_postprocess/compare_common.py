"""Shared inventory, provenance and display preparation for scalar comparisons."""

import csv
import json
from pathlib import Path

import numpy as np

from . import processing as p
from .compare_pwr import finite, one, recordings, relative
from .config import ProjectConfig
from .contract import validate_case
from .paths import load_paths_config
from .plot_common import configure_nominal_frequency


def read_summary(directory, side, mode, required):
    preferred = directory / f"{side}_summary_{mode}.csv"
    path = preferred if preferred.is_file() else directory / f"summary_{mode}.csv"
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = set(required) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}; "
                             f"rerun pcs-postprocess run --mode {mode} --source {side}")
        rows = {}
        for row in reader:
            rows.setdefault(int(row["case_index"]), []).append(row)
    return path, rows


def validate_pair(metas, numeric_keys):
    for key in ("case_index", "case_name", "test_type", "test_type_code"):
        if metas["sim"].get(key) != metas["hil"].get(key):
            raise ValueError(f"SIM/HIL parameter mismatch: {key}")
    for key in numeric_keys:
        if not np.isclose(finite(metas["sim"].get(key), key), finite(metas["hil"].get(key), key),
                          rtol=0., atol=1e-8):
            raise ValueError(f"SIM/HIL parameter mismatch: {key}")


def signals(path, meta, start, end, power_start=None, power_end=None, currents=False):
    """Apply the saved source clock without detecting events or measuring metrics."""
    raw_t, sig = p.load_case_file(str(path), absolute=meta["time_basis"] == "absolute")
    if not raw_t[0] <= start <= raw_t[-1] or end < start:
        raise ValueError("Saved display timing is outside the recording; rerun single-source analysis")
    t = raw_t - start + 1.
    ev = {"event1_start": 1., "event1_end": end - start + 1.}
    keep = (t >= 0.) & (t <= ev["event1_end"] + float(meta.get("record_post_s", 2.)))
    data = {"t_s": t, "Vp_pu": sig["Vp_pu"], "P_pu": sig["P"] / p.HIL_P_BASE_W,
            "Q_pu": sig["Q"] / p.HIL_P_BASE_W if sig["Q"] is not None else np.full(t.shape, np.nan)}
    if sig["f"] is not None:
        data["f_hz"] = sig["f"]
    vabc = np.column_stack([sig[f"V{phase}"] for phase in "abc"])
    iabc = np.column_stack([sig[f"I{phase}"] for phase in "abc"])
    if currents:
        ip, iq, valid = p.calculate_instantaneous_current(vabc, iabc)
        data["Ip_pu"] = p.filter_current_for_display(raw_t, ip, valid) / p.HIL_I_BASE_PEAK_A
        data["Iq_pu"] = p.filter_current_for_display(raw_t, iq, valid) / p.HIL_I_BASE_PEAK_A
    scalar = p.downsample_scalar_signals(data, tuple(k for k in data if k != "t_s"))
    if power_start is not None and power_end is not None:
        power = p.pwr_downsampled_signals(data, power_start, power_end)
        scalar.update(P_pu=power["P_pu"], Q_pu=power["Q_pu"])
    scalar_keep = (scalar["t_s"] >= 0.) & (scalar["t_s"] <= ev["event1_end"] + float(meta.get("record_post_s", 2.)))
    return {"t_s": t[keep], "Vabc_pu": vabc[keep] / p.HIL_V_BASE_PEAK_V,
            "Iabc_pu": iabc[keep] / p.HIL_I_BASE_PEAK_A,
            "breaker": None if sig["breaker"] is None else sig["breaker"][keep],
            "display": {k: v[scalar_keep] for k, v in scalar.items()}, "events": ev,
            "timing": {"event_start_input_s": start, "event_end_input_s": end,
                       "input_to_normalized_shift_s": 1. - start}}


def run_comparison(mode, paths_config, output, case_id, required, numeric_keys,
                   prepare, plot, make_rows, fields, version_key, version, images_per_case,
                   source_errors=None):
    paths = load_paths_config(paths_config)
    config = ProjectConfig.load(paths.project_config)
    p.configure_project(config)
    configure_nominal_frequency(config.nominal_frequency_hz)
    inventories, summaries, references = {}, {}, {}
    for side in ("sim", "hil"):
        source = paths.source(side)
        inventories[side] = recordings(source.input_dir(mode))
        summary, summaries[side] = read_summary(source.output_dir(mode), side, mode,
                                               ["case_index", "case_name", "test_type", *required])
        references[side] = relative(summary, output)
    ids = ([case_id] if case_id is not None else
           sorted(set(inventories["sim"]) | set(inventories["hil"]) |
                  set(summaries["sim"]) | set(summaries["hil"])))
    if not ids:
        raise ValueError(f"No {mode} cases")
    output.mkdir(parents=True, exist_ok=False)
    (output / "figures").mkdir()
    rows, cases, errors = [], [], []
    for ci in ids:
        try:
            files = {s: one(inventories[s], ci, s + " canonical") for s in inventories}
            sources = {s: one(summaries[s], ci, s + " summary") for s in summaries}
            metas = {s: validate_case(files[s], mode) for s in files}
            for side in files:
                for key in ("case_index", "case_name", "test_type"):
                    if str(sources[side][key]) != str(metas[side][key]):
                        raise ValueError(f"{side}: summary / canonical mismatch: {key}")
                if sources[side][version_key] != version:
                    raise ValueError(f"{side}: outdated {mode} summary; rerun single-source analysis")
                for key in numeric_keys:
                    if key in sources[side] and not np.isclose(
                            finite(sources[side][key], key), finite(metas[side][key], key),
                            rtol=0., atol=1e-8):
                        raise ValueError(f"{side}: summary / canonical parameter mismatch: {key}; "
                                         "rerun single-source analysis")
            validate_pair(metas, numeric_keys)
            image_paths = [output / "figures" / f"case_{ci:04d}_sim_hil_{suffix}.png"
                           for suffix in (["overlay", "power_current"] if images_per_case == 2 else ["overlay"])]
            image_refs = [relative(path, output) for path in image_paths]
            pair = {s: prepare(files[s], metas[s], sources[s]) for s in files}
            case_rows = make_rows(ci, metas["hil"], sources, image_refs)
            plot(pair, metas, image_paths, ci)
            rows.extend(case_rows)
            cases.append({"case_index": ci, "case_name": metas["hil"]["case_name"],
                          "figures": image_refs,
                          "canonical": {s: relative(path, output) for s, path in files.items()},
                          "metadata": {s: relative(path.with_name(path.stem + "_meta.json"), output)
                                       for s, path in files.items()},
                          "alignment": {s: pair[s]["timing"] for s in pair}})
            if source_errors is not None:
                for side, source in sources.items():
                    issue = source_errors(source)
                    if issue:
                        errors.append({"case_index": ci, "source": side, "reason": issue})
            print(f"case {ci}: compared")
        except (ValueError, KeyError, OSError, TypeError) as exc:
            errors.append({"case_index": ci, "reason": str(exc)})
            print(f"case {ci}: error: {exc}")
    csv_name = f"compare_summary_{mode}.csv"
    with (output / csv_name).open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    manifest = {"schema_version": 1, "mode": mode, "metric_version": version,
                "metric_event": "event1", "summaries": references,
                "project_config": relative(paths.project_config, output), "compare_csv": csv_name,
                "completed": len(cases), "exceptions": len(errors), "cases": cases, "errors": errors}
    (output / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    return manifest
