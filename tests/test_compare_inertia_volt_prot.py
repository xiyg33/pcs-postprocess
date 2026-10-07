"""Measured event clocks, native breaker edges and referenced comparison outputs."""

import csv
import json
import shutil
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from scipy.io import savemat

from pcs_postprocess import inertia_response as inertia, protection_timing as protection, processing as p
from pcs_postprocess import compare_inertia, compare_volt_prot
from pcs_postprocess.cli import main
from pcs_postprocess.config import ProjectConfig


PACKAGE = Path(__file__).resolve().parents[1]
TMP_PARENT = PACKAGE / "outputs" / "test_tmp"
TMP_PARENT.mkdir(parents=True, exist_ok=True)


class InertiaResponseTests(unittest.TestCase):
    def setUp(self):
        p.configure_project(ProjectConfig.load(PACKAGE / "examples/project.json"))
        self.t = np.arange(0., 6.001, .001)
        self.clock = dict(inertia_metrics_version=inertia.VERSION, inertia_event_start_input_s=2.,
                          inertia_metric_end_input_s=5., inertia_timing_method="frequency_ramp_fit",
                          inertia_end_method="measured_frequency_recovery")

    def measure(self, y, deviation=-.5, clock=None):
        return inertia.summary(self.t, y, np.zeros_like(y), {"freq_deviation_hz": deviation}, clock or self.clock)

    def test_positive_negative_ramps_charge_discharge_and_cross_zero(self):
        for base in (.35, -.35):
            for delta in (.6, -.6):
                with self.subTest(base=base, delta=delta):
                    y = base + delta * np.clip(self.t - 2., 0., 1.)
                    result = self.measure(y, -np.sign(delta) * .5)
                    self.assertAlmostEqual(result["T_up_s"], .9, delta=.015)
                    self.assertLess(result["M_p_pct"], .01)
                    self.assertEqual(result["step_response_start_s"], 1.)
                    self.assertEqual(result["T_s_method"], "continuous_5pct")
                    self.assertAlmostEqual(result["inertia_P_change_pu"], delta, places=5)

    def test_unstable_tail_keeps_peak_time_but_not_final_value_metrics(self):
        y = .35 + .2 * np.clip(self.t - 2., 0., 1.) + .04 * np.maximum(self.t - 3., 0.)
        result = self.measure(y)
        self.assertIsNotNone(result["T_p_s"])
        self.assertIsNone(result["M_p_pct"])
        self.assertIsNone(result["T_s_s"])
        self.assertIn("末段", result["T_s_reason"])

    def test_restoration_samples_cannot_change_first_response(self):
        y = .35 + .2 * np.clip(self.t - 2., 0., 1.)
        original = self.measure(y)
        y[self.t >= 5.] = -100.
        changed = self.measure(y)
        for key in ("M_p_pct", "T_up_s", "T_p_s", "T_s_s", "inertia_final_P_pu"):
            self.assertAlmostEqual(original[key], changed[key], places=5)

    def test_flat_wrong_direction_missing_and_unconfirmed_timing(self):
        for kind in ("flat", "wrong", "missing", "clock", "short"):
            with self.subTest(kind=kind):
                y = .35 + .2 * np.clip(self.t - 2., 0., 1.)
                clock = dict(self.clock)
                if kind == "flat":
                    y[:] = .35
                if kind == "wrong":
                    y = .7 - y
                if kind == "missing":
                    y[(self.t >= 3.5) & (self.t < 3.6)] = np.nan
                if kind == "clock":
                    clock["inertia_timing_method"] = "timing_unavailable"
                if kind == "short":
                    clock["inertia_end_method"] = "recording_truncated"
                result = self.measure(y, clock=clock)
                self.assertIsNone(result["M_p_pct"])
                self.assertTrue(result["M_p_reason"])

    def test_measured_ramp_clock_and_recovery(self):
        start = 2.13
        f = 50. + .5 * np.clip(self.t - start, 0., 1.)
        f[self.t >= start + 2.] = 50.
        meta = {"mode": "inertia", "time_basis": "absolute", "event_start_s": 2.,
                "event_hold_s": 1., "freq_deviation_hz": .5, "freq_ramp_rate_hz_per_s": .5}
        ev = p.event_times_from_meta(meta)
        clock = inertia.timing(self.t, {"f": f}, meta, ev, 0., "absolute_time")
        self.assertAlmostEqual(clock["inertia_event_start_input_s"], start, delta=.005)
        self.assertEqual(clock["inertia_timing_method"], "frequency_ramp_fit")
        self.assertLessEqual(clock["inertia_metric_end_input_s"], start + 2.)


class VoltageTimingTests(unittest.TestCase):
    def setUp(self):
        p.configure_project(ProjectConfig.load(PACKAGE / "examples/project.json"))
        self.t = np.arange(0., 5.001, .001)
        self.voltage = np.where(self.t >= 2.15, .7, 1.)
        self.meta = {"mode": "volt_prot", "time_basis": "absolute", "event_start_s": 2.,
                     "event_hold_s": 1., "V_target_pu": .7}

    def measure(self, breaker, t=None, voltage=None, meta=None):
        meta = meta or self.meta
        return protection.measure(self.t if t is None else t, self.voltage if voltage is None else voltage,
                                  breaker, meta, p.event_times_from_meta(meta), 0.)

    def test_actual_voltage_anchor_native_trip_and_short_glitch(self):
        b = np.where(self.t >= 2.45, 0., 1.)
        b[(self.t >= 2.2) & (self.t < 2.205)] = 0.
        result = self.measure(b)
        self.assertEqual(result["protection_status"], "tripped")
        self.assertAlmostEqual(result["protection_event_start_input_s"], 2.15, delta=.015)
        self.assertAlmostEqual(result["protection_trip_input_s"], 2.45, places=5)
        self.assertEqual(result["t1_s"], 1.)
        self.assertAlmostEqual(result["delta_t_s"], result["t2_s"] - 1., places=8)
        self.assertLess(result["delta_t_s"], .32)

    def test_time_translation_invariance(self):
        b = np.where(self.t >= 2.45, 0., 1.)
        before = self.measure(b)
        after = self.measure(b, t=self.t + 17., meta=dict(self.meta, event_start_s=19.))
        self.assertAlmostEqual(before["delta_t_s"], after["delta_t_s"], places=8)

    def test_no_trip_missing_already_open_invalid_and_unconfirmed(self):
        for b, v, expected in ((np.ones_like(self.t), self.voltage, "no_trip"),
                               (None, self.voltage, "missing_breaker"),
                               (np.zeros_like(self.t), self.voltage, "already_open"),
                               (np.full_like(self.t, 2.), self.voltage, "invalid_breaker"),
                               (np.ones_like(self.t), np.ones_like(self.t), "disturbance_unconfirmed")):
            with self.subTest(expected=expected):
                result = self.measure(b, voltage=v)
                self.assertEqual(result["protection_status"], expected)
                self.assertIsNone(result["delta_t_s"])
                self.assertTrue(result["protection_reason"])

    def test_short_recording_or_unconfirmed_short_tail_is_not_trip(self):
        keep = self.t <= 2.7
        b = np.ones(int(keep.sum()))
        b[-4:] = 0.
        result = self.measure(b, t=self.t[keep], voltage=self.voltage[keep])
        self.assertEqual(result["protection_status"], "insufficient_observation")
        self.assertIsNone(result["delta_t_s"])

    def test_breaker_transition_across_missing_samples_is_unavailable(self):
        keep = ~((self.t >= 2.3) & (self.t < 2.5))
        result = self.measure(np.where(self.t[keep] >= 2.4, 0., 1.),
                              t=self.t[keep], voltage=self.voltage[keep])
        self.assertEqual(result["protection_status"], "insufficient_observation")
        self.assertIsNone(result["delta_t_s"])


class ComparisonWorkflowTests(unittest.TestCase):
    def project(self, root, mode, missing_breaker=False):
        config = root / "paths.json"
        config.write_text(json.dumps({"project_config": str(PACKAGE / "examples/project.json"),
                                      "sources": {s: {"results_root": s} for s in ("sim", "hil")}}), encoding="utf-8")
        for side in ("sim", "hil"):
            directory = root / side / ("results_" + mode) / "canonical"
            directory.mkdir(parents=True)
            t = np.arange(0., 7., .002)
            start = 2. + (.12 if side == "hil" else 0.)
            fields = dict(event_start_s=2., event_hold_s=2., P_ref_pu=.35,
                          hil_rec_start_target_s=0., test_type_code=8 if mode == "inertia" else 11)
            if mode == "inertia":
                fields.update(freq_deviation_hz=-.5, freq_ramp_rate_hz_per_s=-.5)
                f = 50. - .5 * np.clip(t - start, 0., 1.)
                f[t >= start + 3.] = 50.
                power = .35 + .2 * np.clip(t - start, 0., 1.)
                v = np.ones_like(t)
            else:
                fields.update(V_target_pu=.7)
                f = np.full_like(t, 50.)
                power = np.full_like(t, .35)
                v = np.where((t >= start) & (t < start + 2.), .7, 1.)
            angles = 2 * np.pi * 50 * t[:, None] + np.array([0., -2 * np.pi / 3, 2 * np.pi / 3])
            b = np.where(t >= start + .4 + (.02 if side == "sim" else 0.), 0., 1.)
            names = ["t_s", "Va_V", "Vb_V", "Vc_V", "Ia_A", "Ib_A", "Ic_A", "P_W", "Q_var", "f_hz", "Vpos_pu"]
            columns = [t, *(310.269 * v[:, None] * np.cos(angles)).T,
                       *(51.57 * .35 * np.cos(angles)).T, power * 24000., np.zeros_like(t), f, v]
            if not (missing_breaker and side == "hil"):
                names.append("breaker")
                columns.append(b)
            header = np.empty((len(names), 2), object)
            for i, name in enumerate(names):
                header[i] = (name, "")
            path = directory / f"case_1_{mode}.mat"
            savemat(path, {"data": np.column_stack(columns), "header": header})
            meta = dict(mode=mode, case_index=1, case_name=mode, test_type=mode.upper(),
                        time_basis="absolute" if side == "sim" else "recorded", **fields)
            path.with_name(path.stem + "_meta.json").write_text(json.dumps(meta), encoding="utf-8")
            self.assertEqual(main(["run", "--mode", mode, "--source", side,
                                   "--paths-config", str(config), "--no-fig"]), 0)
        return config

    def test_public_commands_reference_metrics_and_produce_expected_figures(self):
        for mode, count, rows in (("inertia", 2, 4), ("volt_prot", 1, 1)):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory(dir=TMP_PARENT) as directory:
                root = Path(directory)
                config = self.project(root, mode)
                output = root / "compare"
                # Comparing must never invoke either source metric calculator.
                with patch.object(inertia, "summary", side_effect=AssertionError("remeasured inertia")), \
                     patch.object(protection, "measure", side_effect=AssertionError("remeasured protection")):
                    self.assertEqual(main(["compare", "--mode", mode, "--all",
                                           "--paths-config", str(config), "--output", str(output)]), 0)
                self.assertEqual(len(list((output / "figures").glob("*.png"))), count)
                with (output / f"compare_summary_{mode}.csv").open(encoding="utf-8-sig") as handle:
                    records = list(csv.DictReader(handle))
                self.assertEqual(len(records), rows)
                source_path = root / "hil" / ("results_" + mode) / "postprocess" / f"hil_summary_{mode}.csv"
                with source_path.open(encoding="utf-8-sig") as handle:
                    source = next(csv.DictReader(handle))
                self.assertEqual(float(records[0]["hil_value"]), float(source["M_p_pct" if mode == "inertia" else "delta_t_s"]))
                self.assertEqual(json.loads((output / "manifest.json").read_text(encoding="utf-8"))["exceptions"], 0)

    def test_missing_breaker_is_reported_in_manifest_with_nonzero_exit(self):
        with tempfile.TemporaryDirectory(dir=TMP_PARENT) as directory:
            root = Path(directory)
            config = self.project(root, "volt_prot", missing_breaker=True)
            output = root / "compare"
            self.assertEqual(main(["compare", "--mode", "volt_prot", "--case", "1", "--paths-config", str(config),
                                   "--output", str(output)]), 1)
            manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["exceptions"], 1)
            self.assertIn("breaker", manifest["errors"][0]["reason"])

    def test_old_summary_columns_rejected_and_parameter_mismatch_recorded(self):
        with tempfile.TemporaryDirectory(dir=TMP_PARENT) as directory:
            root = Path(directory)
            config = self.project(root, "inertia")
            meta_path = next((root / "hil/results_inertia/canonical").glob("*_meta.json"))
            meta = json.loads(meta_path.read_text())
            meta["freq_ramp_rate_hz_per_s"] = -.1
            meta_path.write_text(json.dumps(meta))
            result = compare_inertia.run(config, root / "mismatch")
            self.assertEqual(result["exceptions"], 1)
            self.assertIn("parameter mismatch", result["errors"][0]["reason"])
            summary = root / "sim/results_inertia/postprocess/sim_summary_inertia.csv"
            summary.write_text("case_index,case_name,test_type\n1,inertia,INERTIA\n")
            with self.assertRaisesRegex(ValueError, "rerun pcs-postprocess run"):
                compare_inertia.run(config, root / "old")

    def test_no_trip_and_one_sided_trip_are_not_processing_errors(self):
        common = {"t1_s": "1.0", "t2_s": "", "delta_t_s": "", **{key: "" for key in protection.SUMMARY_FIELDS}}
        common.update(protection_status="no_trip", protection_reason="完整观察窗口内未跳闸")
        meta = {"case_name": "sample", "test_type": "VOLT_PROT_UP", "P_ref_pu": .9, "V_target_pu": 1.1}
        both = compare_volt_prot.comparison_rows(1, meta, {"sim": common, "hil": common}, ["figure.png"])[0]
        self.assertIsNone(both["deviation"])
        self.assertEqual(both["response_match"], "same")
        self.assertIsNone(compare_volt_prot.source_error(common))
        trip = dict(common, protection_status="tripped", protection_reason="", t2_s="1.3", delta_t_s=".3")
        mixed = compare_volt_prot.comparison_rows(1, meta, {"sim": trip, "hil": common}, ["figure.png"])[0]
        self.assertEqual(mixed["response_match"], "trip_mismatch")
        self.assertIsNone(mixed["deviation"])

    def test_batch_continues_after_duplicate_or_missing_source(self):
        for problem in ("duplicate", "missing_source", "missing_summary", "old_version"):
            with self.subTest(problem=problem), tempfile.TemporaryDirectory(dir=TMP_PARENT) as directory:
                root = Path(directory)
                config = self.project(root, "inertia")
                for side in ("sim", "hil"):
                    folder = root / side / "results_inertia/canonical"
                    original = folder / "case_1_inertia.mat"
                    second = folder / "case_2_inertia.mat"
                    shutil.copyfile(original, second)
                    meta = json.loads(original.with_name(original.stem + "_meta.json").read_text())
                    meta["case_index"] = 2
                    second.with_name(second.stem + "_meta.json").write_text(json.dumps(meta), encoding="utf-8")
                    self.assertEqual(main(["run", "--mode", "inertia", "--source", side,
                                           "--paths-config", str(config), "--no-fig"]), 0)
                folder = root / "hil/results_inertia/canonical"
                if problem == "duplicate":
                    shutil.copyfile(folder / "case_1_inertia.mat", folder / "case_1_variant.mat")
                elif problem == "missing_source":
                    (folder / "case_1_inertia.mat").unlink()
                else:
                    summary = root / "hil/results_inertia/postprocess/hil_summary_inertia.csv"
                    with summary.open(encoding="utf-8-sig", newline="") as handle:
                        reader = csv.DictReader(handle)
                        columns, rows = reader.fieldnames, list(reader)
                    if problem == "missing_summary":
                        rows = rows[1:]
                    else:
                        rows[0]["inertia_metrics_version"] = "old"
                    with summary.open("w", encoding="utf-8-sig", newline="") as handle:
                        writer = csv.DictWriter(handle, fieldnames=columns)
                        writer.writeheader()
                        writer.writerows(rows)
                result = compare_inertia.run(config, root / "comparison")
                self.assertEqual(result["exceptions"], 1)
                self.assertEqual(result["completed"], 1)
                self.assertEqual(result["cases"][0]["case_index"], 2)


if __name__ == "__main__":
    unittest.main()
