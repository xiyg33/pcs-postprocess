"""Assess saved paired metrics against versioned limits without rereading recordings."""

from collections import Counter, defaultdict
import csv
from datetime import datetime
import hashlib
import io
import json
import math
import os
from pathlib import Path


MODES = ("frt", "conti_frt", "volt_prot", "freq_prot")
DEFAULT_STANDARDS = Path(__file__).with_name("standards") / "reference_limits_v1.json"
PASS, FAIL, UNKNOWN = "达标", "不达标", "无法判定"
FLOAT_EPS = 1e-12  # Numerical equality only, not an allowance for mechanical delay.
GROUPS = ("F1", "F2", "F3", "FG")
METRICS = {
    "f1_a_pre_steady_mean": "F1", "f1_b2_disturbed_steady_mean": "F1",
    "f1_c2_post_steady_mean": "F1", "f2_b1_rise_transient_mae": "F2",
    "f2_c1_recovery_transient_mae": "F2", "f3_a_pre_steady_max": "F3",
    "f3_b2_disturbed_steady_max": "F3", "f3_c2_post_steady_max": "F3",
    "fg_weighted_overall_mae": "FG",
}
SIGNALS = {"Vpos_pu", "P_pu", "Q_pu", "Ip_pu", "Iq_pu"}
SUMMARY_FIELDS = ["mode", "continuous_type", "case_index", "record_id", "case_name", "display_case_name",
                  "rows", "check_count", "pass_count", "fail_count", "undetermined_count",
                  "estimated_row_count", "assessment_verdict", "assessment_reason"]


def number(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def sha256(path):
    # Streaming works on Python 3.10 as well as newer interpreters.
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relative(path, output):
    try:
        return Path(os.path.relpath(Path(path).resolve(), Path(output).resolve())).as_posix()
    except ValueError:
        return Path(path).resolve().as_posix()


def load_standards(path=DEFAULT_STANDARDS):
    data = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    if data.get("schema_version") != 1 or not data.get("version"):
        raise ValueError("Standards require schema_version=1 and a version")
    frt = data["frt"]
    if frt["unit"] != "pu" or frt["metric_groups"] != METRICS or set(frt["limits"]) != SIGNALS:
        raise ValueError("Standards must define the nine FRT metrics and five pu signals")
    for limits in frt["limits"].values():
        if set(limits) != set(GROUPS) or any(number(v) is None or float(v) < 0 for v in limits.values()):
            raise ValueError("FRT limits must contain finite, nonnegative F1/F2/F3/FG values")
    for mode in ("volt_prot", "freq_prot"):
        spec = data["protection"][mode]
        if spec["target_column"] != ("V_target_pu" if mode == "volt_prot" else "f_target_hz"):
            raise ValueError(f"Unexpected target column for {mode}")
        ids = set()
        if not spec["rules"]:
            raise ValueError(f"No protection rules for {mode}")
        for rule in spec["rules"]:
            if rule["id"] in ids or rule["action"] not in ("open", "closed"):
                raise ValueError(f"Invalid/duplicate protection rule in {mode}")
            ids.add(rule["id"])
            for key in ("lower", "upper", "confirmation_time_s", "time_limit_s"):
                if rule[key] is not None and (number(rule[key]) is None or float(rule[key]) < 0):
                    raise ValueError(f"Invalid {key} in {rule['id']}")
            if any(not isinstance(rule[k], bool) for k in ("lower_inclusive", "upper_inclusive")):
                raise ValueError(f"Interval inclusivity must be boolean: {rule['id']}")
            if rule["lower"] is not None and rule["upper"] is not None and rule["lower"] >= rule["upper"]:
                raise ValueError(f"Empty protection interval: {rule['id']}")
            if rule["action"] == "open" and rule["time_limit_s"] is None:
                raise ValueError(f"Missing opening deadline: {rule['id']}")
    return data


def combine(verdicts):
    values = list(verdicts)
    return FAIL if FAIL in values else UNKNOWN if not values or UNKNOWN in values else PASS


def continuous_type(row):
    explicit = row.get("continuous_type")
    if explicit in ("conti_frt", "conti_lvrt"):
        return explicit
    name = str(row.get("record_id", "")).upper()
    if name.startswith("CONTI_LVRT_"):
        return "conti_lvrt"
    if name.startswith("CONTI_FRT_"):
        return "conti_frt"
    raise ValueError("Continuous performance requires continuous_type or a CONTI_FRT/LVRT record_id")


def assess_frt(row, standards):
    if row["signal"] not in SIGNALS or row["unit"] != "pu":
        raise ValueError(f"Unknown FRT signal/unit: {row['signal']} / {row['unit']}")
    result = dict(row)
    limits = standards["frt"]["limits"][row["signal"]]
    result.update({group.lower() + "_limit_pu": limits[group] for group in GROUPS})
    verdicts, reasons = [], []
    for metric, group in METRICS.items():
        value = number(row[metric])
        verdict = UNKNOWN if value is None else PASS if abs(value) <= limits[group] + FLOAT_EPS else FAIL
        result[metric + "_verdict"] = verdict
        verdicts.append(verdict)
        if verdict != PASS:
            reasons.append(f"{metric}: " + ("缺少有效数值" if value is None else f"超出 {limits[group]} pu"))
    result.update(assessment_verdict=combine(verdicts), assessment_reason="; ".join(reasons),
                  assessment_basis="estimated" if row.get("status") == "estimated" else
                                   "partial_metrics" if row.get("status") == "unavailable" else "saved_metrics")
    return result


def select_rule(target, rules):
    if target is None or target < 0:
        return None
    matches = []
    for rule in rules:
        lower, upper = rule["lower"], rule["upper"]
        if lower is not None and (target < lower or (target == lower and not rule["lower_inclusive"])):
            continue
        if upper is not None and (target > upper or (target == upper and not rule["upper_inclusive"])):
            continue
        matches.append(rule)
    if len(matches) > 1:
        raise ValueError(f"Overlapping protection rules at target {target}")
    return matches[0] if matches else None


def protection_verdict(row, side, rule):
    state = row[side + "_protection_status"]
    delay = number(row[side + "_value"])
    observed = number(row[side + "_protection_observed_duration_s"])
    if rule is None:
        return UNKNOWN, "目标值缺失、无效或未匹配标准规则"
    if state not in ("tripped", "no_trip"):
        return UNKNOWN, row.get(side + "_protection_reason") or row.get(side + "_reason") or state or "缺少保护状态"
    if observed is None or observed <= 0:
        return UNKNOWN, "缺少有效观察时长"
    if state == "tripped":
        start, end = number(row[side + "_t1_s"]), number(row[side + "_t2_s"])
        if (delay is None or delay < 0 or start is None or end is None
                or abs(end - start - delay) > 1e-8 or delay > observed + FLOAT_EPS):
            return UNKNOWN, "跳闸时间不完整或与 t2−t1/观察时长不一致"
        if rule["action"] == "closed":
            return FAIL, "要求 breaker 保持闭合，但检测到跳闸"
        if delay <= rule["time_limit_s"] + FLOAT_EPS:
            return PASS, ""
        return FAIL, f"响应时间超出 {rule['time_limit_s']} s 上限"
    if row[side + "_value"] not in (None, "") or row[side + "_t2_s"] not in (None, ""):
        return UNKNOWN, "未跳闸状态包含跳闸时间，需重新生成性能表"
    minimum = rule["time_limit_s"] if rule["action"] == "open" else rule["confirmation_time_s"]
    if minimum is not None and observed + FLOAT_EPS < minimum:
        return UNKNOWN, "观察时长不足以评价该保护要求"
    if rule["action"] == "closed":
        return PASS, "完整观察窗口内 breaker 保持闭合" + (
            "；不评价Ⅰ段功率转换" if rule["confirmation_time_s"] is not None else "")
    return FAIL, "完整观察窗口内未检测到要求的 breaker 断开"


def assess_protection(row, mode, standards):
    result = dict(row)
    spec = standards["protection"][mode]
    rule = select_rule(number(row[spec["target_column"]]), spec["rules"])
    result.update(rule_id=rule["id"] if rule else "", protection_stage=rule["stage"] if rule else "",
                  expected_breaker_action=rule["action"] if rule else "",
                  confirmation_time_s=rule["confirmation_time_s"] if rule else None,
                  breaker_time_limit_s=rule["time_limit_s"] if rule else None)
    verdicts = []
    for side in ("sim", "hil"):
        verdict, reason = protection_verdict(row, side, rule)
        result[side + "_verdict"] = verdict
        result[side + "_assessment_reason"] = reason
        verdicts.append(verdict)
    result.update(assessment_verdict=combine(verdicts), assessment_basis="estimated" if any(
        row.get(s + "_status") == "estimated" for s in ("sim", "hil")) else "saved_breaker_evidence",
        assessment_reason="; ".join(f"{s}: {result[s + '_assessment_reason']}" for s in ("sim", "hil")
                                    if result[s + "_verdict"] != PASS))
    return result


def read_csv(path, mode, standards):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        if len(fields) != len(set(fields)):
            raise ValueError(f"Duplicate CSV columns: {path}")
        rows = list(reader)
    validate_rows(fields, rows, path, mode, standards)
    return fields, rows


def validate_rows(fields, rows, path, mode, standards):
    required = {"case_index"} if mode != "conti_frt" else {"record_id", "event_index"}
    if mode in ("frt", "conti_frt"):
        required |= {"signal", "unit", "status", *METRICS}
    else:
        required.add("unit")
        required.add(standards["protection"][mode]["target_column"])
        required |= {s + "_" + k for s in ("sim", "hil") for k in (
            "value", "t1_s", "t2_s", "protection_status", "protection_observed_duration_s")}
    if required - set(fields):
        raise ValueError(f"{path}: missing columns {sorted(required - set(fields))}; regenerate performance metrics")
    if not rows:
        raise ValueError(f"Empty performance CSV: {path}")
    seen, group_signals = set(), defaultdict(set)
    for row in rows:
        if None in row or any(v is None for v in row.values()):
            raise ValueError(f"Malformed CSV row: {path}")
        group = (row["record_id"], row["event_index"]) if mode == "conti_frt" else (row["case_index"],)
        if not all(group):
            raise ValueError(f"Missing case/segment identity: {path}")
        key = (*group, row["signal"]) if mode in ("frt", "conti_frt") else group
        if key in seen:
            raise ValueError(f"Duplicate performance row {key}: {path}")
        seen.add(key)
        if mode in ("frt", "conti_frt"):
            if row["status"] not in ("valid", "estimated", "unavailable"):
                raise ValueError(f"Unknown performance status: {path}")
            group_signals[group].add(row["signal"])
        elif row["unit"] != "s":
            raise ValueError(f"Protection response must use seconds: {path}")
    if any(signals != SIGNALS for signals in group_signals.values()):
        raise ValueError(f"{path}: each FRT case/segment requires all five signals; regenerate performance metrics")


def summarize(mode, rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row["record_id"] if mode == "conti_frt" else row["case_index"]].append(row)
    result = []
    for group in grouped.values():
        first = group[0]
        verdicts = [r[k + "_verdict"] for r in group for k in (
            METRICS if mode in ("frt", "conti_frt") else ("sim", "hil"))]
        counts = Counter(verdicts)
        record = {key: first.get(key, "") for key in SUMMARY_FIELDS[:6]}
        record.update(mode=mode, rows=len(group), check_count=len(verdicts), pass_count=counts[PASS],
                      fail_count=counts[FAIL], undetermined_count=counts[UNKNOWN],
                      estimated_row_count=sum(r["assessment_basis"] == "estimated" for r in group),
                      assessment_verdict=combine(verdicts), assessment_reason="; ".join(dict.fromkeys(
                          r["assessment_reason"] for r in group if r["assessment_reason"])))
        result.append(record)
    return result


def write_csv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def default_output(compare_dir):
    return Path(compare_dir) / "assessment" / datetime.now().strftime("%Y-%m-%d_%H%M%S_%f")


def run(compare_dir, output=None, *, standards_path=None, modes=None, continuous_performance=None,
        prepared_inputs=None, preparation=None):
    """Evaluate CSV tables; prepared_inputs is an in-memory adapter contract, not a recording reader."""
    compare_dir = Path(compare_dir).resolve()
    output = Path(output or default_output(compare_dir)).resolve()
    if output.exists():
        raise FileExistsError(f"Assessment output already exists: {output}")
    modes = list(MODES if modes is None else modes)
    if not modes or len(modes) != len(set(modes)) or set(modes) - set(MODES):
        raise ValueError("Assessment modes must be a nonempty, unique subset of " + ", ".join(MODES))
    standards_path = Path(standards_path or DEFAULT_STANDARDS).resolve()
    standards = load_standards(standards_path)
    inputs, tables, summary, provenance = {}, {}, [], []
    prepared_inputs = prepared_inputs or {}
    for mode in modes:
        if mode in prepared_inputs:
            spec = prepared_inputs[mode]
            # Validate prepared tables with the same CSV contract as saved tables.
            buffer = io.StringIO(newline="")
            writer = csv.DictWriter(buffer, fieldnames=spec["fields"])
            writer.writeheader()
            writer.writerows(spec["rows"])
            content = buffer.getvalue().encode("utf-8-sig")
            source = output / "inputs" / spec["filename"]
            if Path(spec["filename"]).name != spec["filename"]:
                raise ValueError("Prepared input filename must be a basename")
            parsed = csv.DictReader(io.StringIO(buffer.getvalue()))
            fields, raw_rows = list(parsed.fieldnames), list(parsed)
            validate_rows(fields, raw_rows, source, mode, standards)
            inputs[mode] = (source, content)
            provenance.append(dict(mode=mode, csv=relative(source, output), sha256=hashlib.sha256(content).hexdigest(),
                                   prepared=True))
        else:
            name = "compare_performance_conti_frt.csv" if mode == "conti_frt" else f"compare_summary_{mode}.csv"
            source = Path(continuous_performance).resolve() if mode == "conti_frt" and continuous_performance else compare_dir / mode / name
            fields, raw_rows = read_csv(source, mode, standards)
            provenance.append(dict(mode=mode, csv=relative(source, output), sha256=sha256(source), prepared=False))
        evaluated = []
        for row in raw_rows:
            checked = assess_frt(row, standards) if mode in ("frt", "conti_frt") else assess_protection(row, mode, standards)
            if mode == "conti_frt":
                checked["continuous_type"] = continuous_type(row)
            checked.update(standards_version=standards["version"], source_csv=relative(source, output))
            evaluated.append(checked)
        tables[mode] = (list(dict.fromkeys([*fields, *(k for row in evaluated for k in row)])), evaluated)
        summary.extend(summarize(mode, evaluated))
    # Finish validation before creating any assessment artifacts.
    if any(not table[1] for table in tables.values()):
        raise ValueError("No prepared performance rows")
    output.mkdir(parents=True, exist_ok=False)
    if inputs:
        (output / "inputs").mkdir()
        for source, content in inputs.values():
            source.write_bytes(content)
    for mode, (fields, rows) in tables.items():
        write_csv(output / f"compliance_{mode}.csv", fields, rows)
    write_csv(output / "compliance_summary.csv", SUMMARY_FIELDS, summary)
    (output / "standards_snapshot.json").write_text(json.dumps(standards, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    manifest = dict(schema_version=1, assessment_version="paired_compliance_v1", standards_version=standards["version"],
                    standards_source=relative(standards_path, output), standards_sha256=sha256(standards_path),
                    standards_snapshot="standards_snapshot.json", compare_dir=relative(compare_dir, output),
                    inputs=provenance, modes={m: dict(csv=f"compliance_{m}.csv", rows=len(rows),
                        verdict_counts=dict(Counter(row["assessment_verdict"] for row in rows))) for m, (_, rows) in tables.items()},
                    summary_csv="compliance_summary.csv", summary_cases=len(summary),
                    verdict_counts=dict(Counter(row["assessment_verdict"] for row in summary)),
                    interpretation=dict(frt="Absolute paired deviation <= group limit; nine checks per signal.",
                        protection=standards["protection"]["timing_interpretation"],
                        frequency_stage_one=standards["protection"]["frequency_stage_one_scope"],
                        estimated="Finite estimated metrics are evaluated and their evidence status is retained.",
                        aggregation="fail first, then undetermined, otherwise pass", engineering_tolerance_s=0),
                    run_status="complete")
    if preparation is not None:
        (output / "inputs").mkdir(exist_ok=True)
        (output / "inputs" / "continuous_preparation.json").write_text(
            json.dumps(preparation, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
        manifest["continuous_preparation"] = "inputs/continuous_preparation.json"
    (output / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    return manifest


def add_arguments(parser):
    parser.add_argument("--compare-dir", type=Path, required=True)
    parser.add_argument("--standards", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--mode", nargs="+", choices=MODES)
    parser.add_argument("--continuous-performance", type=Path,
                        help="逐故障段 F1/F2/F3/FG 表；不能使用旧整段统计表")


def execute(args, parser):
    try:
        output = args.output or default_output(args.compare_dir)
        manifest = run(args.compare_dir, output, standards_path=args.standards, modes=args.mode,
                       continuous_performance=args.continuous_performance)
        print(json.dumps(dict(output=str(output.resolve()), verdict_counts=manifest["verdict_counts"]), ensure_ascii=False))
        # A failed check is a completed evaluation, not a processing exception.
        return 0
    except (ValueError, KeyError, OSError, TypeError) as exc:
        parser.exit(2, f"error: {exc}\n")
