"""Reference limits, protection endpoints and saved-table assessment contracts."""

from copy import deepcopy
import csv
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from pcs_postprocess import assessment as a
from pcs_postprocess.cli import main


def frt_rows(case=1, record=None, event=1):
    rows = []
    for signal in sorted(a.SIGNALS):
        row = dict(case_index=str(case), case_name="source_name", display_case_name="display_name",
                   signal=signal, unit="pu", status="estimated", reason="timing_estimated",
                   **{metric: "0.01" for metric in a.METRICS})
        if record:
            row.update(record_id=record, event_index=str(event), target_voltage_pu="0.2")
        rows.append(row)
    return rows


def protection_row(target=1.15, mode="volt_prot", delay=.1, state="tripped", observed=60.):
    row = dict(case_index="1", case_name="case", unit="s")
    row["V_target_pu" if mode == "volt_prot" else "f_target_hz"] = target
    for side in ("sim", "hil"):
        row.update({side + "_value": delay if state == "tripped" else "",
                    side + "_t1_s": 1., side + "_t2_s": 1. + delay if state == "tripped" else "",
                    side + "_protection_status": state, side + "_protection_observed_duration_s": observed,
                    side + "_protection_reason": "", side + "_status": "estimated"})
    return row


class LimitTests(unittest.TestCase):
    def setUp(self):
        self.standard = a.load_standards()

    def test_reference_transcription_and_all_nine_metric_boundaries(self):
        expected = {s: (.02, .05, .05, .05) if s == "Vpos_pu" else (.1, .2, .15, .15) for s in a.SIGNALS}
        for row in frt_rows():
            limits = self.standard["frt"]["limits"][row["signal"]]
            self.assertEqual(tuple(limits[k] for k in a.GROUPS), expected[row["signal"]])
            for metric, group in a.METRICS.items():
                for delta, verdict in ((-1e-6, a.PASS), (0., a.PASS), (1e-6, a.FAIL)):
                    with self.subTest(signal=row["signal"], metric=metric, delta=delta):
                        copy = dict(row)
                        copy[metric] = limits[group] + delta
                        result = a.assess_frt(copy, self.standard)
                        self.assertEqual(result[metric + "_verdict"], verdict)
                        self.assertEqual(result["assessment_basis"], "estimated")
                        self.assertEqual(result[metric], copy[metric])

    def test_missing_one_window_leaves_other_checks_and_fail_has_priority(self):
        row = frt_rows()[0]
        row.update(status="unavailable", f1_b2_disturbed_steady_mean="", fg_weighted_overall_mae="nan")
        result = a.assess_frt(row, self.standard)
        self.assertEqual(result["f1_b2_disturbed_steady_mean_verdict"], a.UNKNOWN)
        self.assertEqual(result["f1_a_pre_steady_mean_verdict"], a.PASS)
        self.assertEqual(result["assessment_verdict"], a.UNKNOWN)
        row["f3_c2_post_steady_max"] = .3
        self.assertEqual(a.assess_frt(row, self.standard)["assessment_verdict"], a.FAIL)
        for invalid in ("", "NaN", "inf", "-inf", "bad", None):
            self.assertIsNone(a.number(invalid))

    def test_exact_protection_interval_endpoints(self):
        cases = {
            "volt_prot": [(0., "voltage_under_II"), (.499999, "voltage_under_II"), (.5, "voltage_under_I"),
                          (.849999, "voltage_under_I"), (.85, "voltage_untriggered"), (1.1, "voltage_untriggered"),
                          (1.100001, "voltage_over_I"), (1.2, "voltage_over_I"), (1.200001, "voltage_over_II")],
            "freq_prot": [(46.499999, "frequency_under_II"), (46.5, "frequency_under_I"),
                          (48.499999, "frequency_under_I"), (48.5, "frequency_untriggered"),
                          (50.5, "frequency_untriggered"), (50.500001, "frequency_over_I"),
                          (51.5, "frequency_over_I"), (51.500001, "frequency_over_II")]}
        for mode, targets in cases.items():
            rules = self.standard["protection"][mode]["rules"]
            for target, expected in targets:
                with self.subTest(mode=mode, target=target):
                    self.assertEqual(a.select_rule(target, rules)["id"], expected)

    def test_all_opening_deadlines_and_no_engineering_tolerance(self):
        for mode, target, limit in (("volt_prot", 1.15, 10.1), ("volt_prot", 1.3, 1.1),
                                    ("volt_prot", .7, 2.1), ("volt_prot", .3, 1.5),
                                    ("freq_prot", 52., .1), ("freq_prot", 46., .1)):
            for delta, expected in ((-1e-6, a.PASS), (0., a.PASS), (1e-6, a.FAIL)):
                row = protection_row(target, mode, limit + delta)
                row["deviation"] = 0  # Matching SIM/HIL delays do not imply compliance.
                result = a.assess_protection(row, mode, self.standard)
                self.assertEqual(result["breaker_time_limit_s"], limit)
                self.assertEqual(result["sim_verdict"], expected)
                self.assertEqual(result["assessment_verdict"], expected)

    def test_frequency_stage_one_and_normal_conditions_require_closed_breaker(self):
        for mode, target in (("freq_prot", 51.), ("freq_prot", 47.), ("freq_prot", 50.),
                             ("volt_prot", 1.0), ("volt_prot", 1.1)):
            row = protection_row(target, mode, state="no_trip")
            self.assertEqual(a.assess_protection(row, mode, self.standard)["assessment_verdict"], a.PASS)
            row = protection_row(target, mode)
            result = a.assess_protection(row, mode, self.standard)
            self.assertEqual(result["sim_verdict"], a.FAIL)
            self.assertIsNone(result["breaker_time_limit_s"])

    def test_unavailable_evidence_no_trip_short_observation_and_inconsistent_times(self):
        for state in ("disturbance_unconfirmed", "missing_breaker", "invalid_breaker", "already_open", "insufficient_observation"):
            row = protection_row(state=state)
            self.assertEqual(a.assess_protection(row, "volt_prot", self.standard)["sim_verdict"], a.UNKNOWN)
        row = protection_row(state="no_trip", observed=60.)
        self.assertEqual(a.assess_protection(row, "volt_prot", self.standard)["sim_verdict"], a.FAIL)
        row["sim_protection_observed_duration_s"] = .05
        result = a.assess_protection(row, "volt_prot", self.standard)
        self.assertEqual(result["sim_verdict"], a.UNKNOWN)
        self.assertEqual(result["assessment_verdict"], a.FAIL)
        for key, value in (("sim_value", ""), ("sim_value", -.1), ("sim_t2_s", 2.),
                           ("sim_protection_observed_duration_s", .01)):
            row = protection_row()
            row[key] = value
            self.assertEqual(a.assess_protection(row, "volt_prot", self.standard)["sim_verdict"], a.UNKNOWN)


class AssessmentWorkflowTests(unittest.TestCase):
    def make_delivery(self, base):
        for mode in a.MODES:
            if mode == "frt":
                rows = frt_rows()
            elif mode == "conti_frt":
                rows = [row for record, events in (("CONTI_FRT_p020", range(1, 7)), ("CONTI_LVRT_p020", (1, 3)))
                        for event in events for row in frt_rows(record=record, event=event)]
            else:
                rows = [protection_row(52. if mode == "freq_prot" else 1.15, mode)]
            name = "compare_performance_conti_frt.csv" if mode == "conti_frt" else f"compare_summary_{mode}.csv"
            path = base / mode / name
            path.parent.mkdir(parents=True)
            a.write_csv(path, list(rows[0]), rows)

    def test_complete_output_combined_continuous_summary_snapshot_and_input_preservation(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp) / "delivery"
            self.make_delivery(base)
            original = {p: p.read_bytes() for p in base.glob("*/*.csv")}
            output = base / "assessment/result"
            result = a.run(base, output)
            self.assertEqual(result["modes"]["conti_frt"]["rows"], 40)
            with (output / "compliance_conti_frt.csv").open(encoding="utf-8-sig") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len({(r["record_id"], r["event_index"]) for r in rows}), 8)
            self.assertEqual({r["continuous_type"] for r in rows}, {"conti_frt", "conti_lvrt"})
            self.assertTrue(all(r["assessment_verdict"] == a.PASS for r in rows))
            self.assertEqual(result["summary_cases"], 5)
            self.assertEqual(a.load_standards(output / "standards_snapshot.json"), a.load_standards())
            self.assertEqual({p: p.read_bytes() for p in original}, original)
            self.assertTrue((output / "compliance_summary.csv").read_bytes().startswith(b"\xef\xbb\xbf"))
            with self.assertRaises(FileExistsError):
                a.run(base, output)

    def test_bad_missing_duplicate_and_old_tables_fail_before_output_creation(self):
        for kind in ("missing_column", "duplicate", "missing_signal", "old_continuous", "bad_units", "bad_status"):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory() as tmp:
                base = Path(tmp)
                self.make_delivery(base)
                path = base / "frt/compare_summary_frt.csv"
                fields, rows = a.read_csv(path, "frt", a.load_standards())
                if kind == "missing_column":
                    fields.remove(next(iter(a.METRICS)))
                    for row in rows:
                        row.pop(next(iter(a.METRICS)))
                elif kind == "duplicate":
                    rows.append(rows[0])
                elif kind == "missing_signal":
                    rows.pop()
                elif kind == "bad_units":
                    rows[0]["unit"] = "W"
                elif kind == "bad_status":
                    rows[0]["status"] = "bad"
                else:
                    path = base / "conti_frt/compare_performance_conti_frt.csv"
                    fields, rows = ["record_id", "mae_pu"], [{"record_id": "old", "mae_pu": "0"}]
                a.write_csv(path, fields, rows)
                with self.assertRaises(ValueError):
                    a.run(base, base / "output")
                self.assertFalse((base / "output").exists())

    def test_custom_standards_cli_and_subset(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            self.make_delivery(base)
            data = a.load_standards()
            data["version"] = "custom"
            data["frt"]["limits"]["P_pu"]["F1"] = .001
            standards = base / "custom.json"
            standards.write_text(json.dumps(data), encoding="utf-8")
            output = base / "custom_assessment"
            self.assertEqual(main(["assess", "--compare-dir", str(base), "--mode", "frt", "--standards",
                                   str(standards), "--output", str(output)]), 0)
            manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["standards_version"], "custom")
            self.assertEqual(manifest["verdict_counts"], {a.FAIL: 1})
            self.assertFalse((output / "compliance_freq_prot.csv").exists())

    def test_prepared_continuous_uses_same_contract_and_immutable_input_snapshot(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            rows = frt_rows(record="CONTI_FRT_p020")
            prepared = {"conti_frt": dict(filename="performance.csv", fields=list(rows[0]), rows=rows)}
            result = a.run(base, base / "out", modes=["conti_frt"], prepared_inputs=prepared)
            self.assertTrue(result["inputs"][0]["prepared"])
            self.assertTrue((base / "out/inputs/performance.csv").exists())
            prepared["conti_frt"]["rows"] = []
            with self.assertRaises(ValueError):
                a.run(base, base / "bad", modes=["conti_frt"], prepared_inputs=prepared)
            self.assertFalse((base / "bad").exists())


if __name__ == "__main__":
    unittest.main()
