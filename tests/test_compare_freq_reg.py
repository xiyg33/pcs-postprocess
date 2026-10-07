"""Regulation response metrics, referenced comparisons and deadband records."""

import csv
import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from scipy.io import loadmat, savemat

from pcs_postprocess import compare_freq_reg as compare, processing
from pcs_postprocess.cli import main
from pcs_postprocess.config import ProjectConfig

PACKAGE = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("freq_synthetic", PACKAGE / "examples/make_synthetic.py")
synthetic = importlib.util.module_from_spec(spec)
spec.loader.exec_module(synthetic)


class RegulationMetricsTests(unittest.TestCase):
    def setUp(self):
        processing.configure_project(ProjectConfig.load(PACKAGE / "examples/project.json"))

    def response(self, baseline=.35, delta=.2, delay=.08):
        t = np.arange(0, 5.001, .001)
        u = np.maximum(t - 1. - delay, 0.)
        y = baseline + delta * (1. - np.exp(-u / .08))
        y[t >= 4.] = baseline  # restoration must never enter first-step metrics
        return {"t_s": t, "P_pu": y, "Q_pu": np.zeros_like(t),
                "event1_start_s": 1., "event1_end_s": 4.,
                "response_edge_detected": 1, "response_edge_normalized_s": 1.,
                "command_time_reliable": 1, "command_time_source": "synchronized",
                "command_event1_normalized_s": 1.}

    def test_both_power_directions_and_charge_discharge(self):
        for baseline in (.35, -.65):
            for delta in (.2, -.2):
                with self.subTest(baseline=baseline, delta=delta):
                    data = self.response(baseline, delta)
                    original = data["P_pu"].copy()
                    result = processing.freq_reg_step_summary(data, {"freq_deviation_hz": -np.sign(delta) * .2})
                    self.assertGreater(result["step_response_start_s"], 1.08)
                    self.assertEqual(result["step_timing_basis"], "five_point_response_estimate")
                    self.assertAlmostEqual(result["T_up_s"], .18, delta=.025)
                    self.assertLess(result["M_p_pct"], .01)
                    self.assertEqual(result["metric_sampling_ms"], 5)
                    self.assertEqual(result["T_s_method"], "continuous_5pct")
                    np.testing.assert_array_equal(data["P_pu"], original)

    def test_fourth_entry_and_drifting_tail_keep_pwr_method(self):
        data = self.response()
        t = data["t_s"]
        u = np.maximum(t - 1.08, 0.)
        data["P_pu"] = .35 + .2 * (1. - np.exp(-u / .08)) + .02 * np.sin(2 * np.pi * 5 * u)
        data["P_pu"][t <= 1.08] = .35
        result = processing.freq_reg_step_summary(data, {"freq_deviation_hz": -.2})
        self.assertEqual(result["T_s_method"], "fourth_band_entry_estimate")
        self.assertIn("第4次", result["T_s_reason"])
        drifting = self.response()
        drifting["P_pu"] += np.maximum(t - 2., 0.) * .03
        result = processing.freq_reg_step_summary(drifting, {"freq_deviation_hz": -.2})
        self.assertIsNotNone(result["M_p_pct"])
        self.assertIn("数值仅供参考", result["M_p_reason"])

    def test_no_response_wrong_direction_and_short_recording(self):
        for kind, reason in (("flat", "幅值过小"), ("wrong", "方向与指令不符"),
                             ("short", "观察时长不足"), ("gap", "缺失采样")):
            with self.subTest(kind=kind):
                data = self.response()
                if kind == "flat":
                    data["P_pu"][:] = .35
                elif kind == "wrong":
                    data["P_pu"] = .7 - data["P_pu"]
                elif kind == "short":
                    data["event1_end_s"] = 1.2
                else:
                    data["P_pu"][(data["t_s"] >= 3.8) & (data["t_s"] < 3.9)] = np.nan
                result = processing.freq_reg_step_summary(data, {"freq_deviation_hz": -.2})
                self.assertEqual(result["M_p_pct"], "")
                self.assertIn(reason, result["M_p_reason"])

    def test_missing_samples_inside_response_are_not_interpolated(self):
        data = self.response()
        data["P_pu"][(data["t_s"] > 2.) & (data["t_s"] < 2.1)] = np.nan
        result = processing.freq_reg_step_summary(data, {"freq_deviation_hz": -.2})
        self.assertIsNone(result["M_p_pct"])
        self.assertIn("缺失采样", result["M_p_reason"])


class RegulationComparisonTests(unittest.TestCase):
    def make_project(self, base, deadband=False):
        paths = base / "paths.json"
        paths.write_text(json.dumps({"project_config": str(PACKAGE / "examples/project.json"),
                                     "sources": {s: {"results_root": s} for s in ("sim", "hil")}}), encoding="utf-8")
        files, rows = {}, {}
        for side in ("sim", "hil"):
            root = base / side / "results_freq_reg"
            files[side] = synthetic.create(root / "canonical", "freq_reg", "absolute" if side == "sim" else "recorded")
            meta_path = files[side].with_name(files[side].stem + "_meta.json")
            meta = json.loads(meta_path.read_text())
            meta.update(test_type="FREQ_REG_DEADBAND" if deadband else "FREQ_REG", test_type_code=6,
                        event_hold_s=2., freq_deviation_hz=-.02 if deadband else -.2)
            if side == "sim":
                meta["data_source"] = "sim"
            meta_path.write_text(json.dumps(meta), encoding="utf-8")
            mat = loadmat(files[side])
            t = mat["data"][:, 0]
            u = np.maximum(t - 2.08, 0.)
            mat["data"][:, 7] = 24000 * (.3 + .2 * (1. - np.exp(-u / .08)))
            savemat(files[side], {"data": mat["data"], "header": mat["header"]})
            rows[side] = {"case_index": 1002, "case_name": meta["case_name"], "test_type": meta["test_type"],
                          "rec_offset_s": 0, "offset_source": "absolute_time" if side == "sim" else "edge",
                          "P_change_abs_pu": .05, "P_change_5pct_pass": 1,
                          "M_p_pct": 12.25 if side == "sim" else 10., "T_up_s": .1,
                          "T_p_s": .2, "T_s_s": .3, "M_p_reason": "", "T_up_reason": "",
                          "T_p_reason": "", "T_s_reason": "", "step_timing_basis": "five_point_response_estimate",
                          "step_response_start_s": 1.08, "T_s_method": "fourth_band_entry_estimate",
                          "metric_sampling_ms": 5, "人工确认Ts": ""}
            rows[side].update({key: "" for key in compare.TIMING_FIELDS})
            rows[side].update(freq_event_start_recorder_s=2., freq_event_end_recorder_s=4.,
                              freq_timing_method="absolute_time" if side == "sim" else "frequency_piecewise_fit",
                              freq_timing_confidence="absolute" if side == "sim" else "accepted_estimate",
                              freq_timing_review_required=0, freq_metric_end_recorder_s=4.)
            self.write_summary(base, side, rows[side])
        return paths, files, rows

    def write_summary(self, base, side, row, filename=None):
        path = base / side / "results_freq_reg/postprocess" / (filename or f"{side}_summary_freq_reg.csv")
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(row))
            writer.writeheader()
            writer.writerow(row)
        return path

    def test_referenced_metrics_two_figures_and_tracing(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            paths, files, _ = self.make_project(base)
            originals = [compare.summary_path(base / s / "results_freq_reg/postprocess", s).read_bytes() for s in files]
            out = base / "compare"
            captured = []
            close = compare.plt.close
            with patch.object(processing, "prepare_case", side_effect=AssertionError("must not calculate")), \
                 patch.object(processing, "freq_reg_step_summary", side_effect=AssertionError("must not calculate")), \
                 patch.object(processing, "detect_pwr_edge_offset", side_effect=AssertionError("must not detect")), \
                 patch.object(compare.plt, "close", side_effect=captured.append):
                manifest = compare.run(paths, out)
            self.assertEqual(manifest["completed"], 1)
            self.assertEqual(manifest["exceptions"], 0)
            self.assertEqual(len(list((out / "figures").glob("*.png"))), 2)
            figs = [f for f in captured if hasattr(f, "axes")]
            self.assertEqual([len(f.axes) for f in figs], [3, 5])
            for fig in figs:
                for ax in fig.axes:
                    self.assertEqual(ax.get_xlim()[0], 0.)
                    self.assertEqual(ax.get_legend_handles_labels()[1], ["HIL", "SIM"])
                close(fig)
            for side in files:
                self.assertEqual((out / manifest["cases"][0]["canonical"][side]).resolve(), files[side].resolve())
                self.assertEqual(manifest["cases"][0]["alignment"][side]["step_response_start_s"], 1.08)
            self.assertEqual([compare.summary_path(base / s / "results_freq_reg/postprocess", s).read_bytes() for s in files], originals)
            with (out / "compare_summary_freq_reg.csv").open(encoding="utf-8-sig") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 4)
            self.assertEqual(float(rows[0]["deviation"]), 2.25)
            self.assertEqual(rows[1]["status"], "estimated")
            self.assertTrue((out / rows[0]["power_current_figure_png"]).exists())
            self.assertTrue((out / "compare_summary_freq_reg.csv").read_bytes().startswith(b"\xef\xbb\xbf"))
            with self.assertRaises(FileExistsError):
                compare.run(paths, out)

    def test_deadband_threshold_no_figures_and_missing_hil(self):
        for missing in (False, True):
            with self.subTest(missing=missing), tempfile.TemporaryDirectory() as tmp:
                base = Path(tmp)
                paths, files, sources = self.make_project(base, True)
                sources["sim"].update(P_change_abs_pu=.0501, P_change_5pct_pass=0)
                self.write_summary(base, "sim", sources["sim"])
                if missing:
                    files["hil"].unlink()
                out = base / "compare"
                with patch.object(compare, "plot_pair", side_effect=AssertionError("deadband must not draw")):
                    result = compare.run(paths, out)
                with (out / "compare_summary_freq_reg.csv").open(encoding="utf-8-sig") as handle:
                    rows = list(csv.DictReader(handle))
                self.assertEqual(len(rows), 1)
                self.assertEqual(rows[0]["sim_5pct_pass"], "0.0")
                self.assertEqual(rows[0]["hil_5pct_pass"], "" if missing else "1.0")
                self.assertEqual(result["exceptions"], int(missing))
                self.assertEqual(result["partial"], int(missing))
                self.assertFalse(list((out / "figures").iterdir()))

    def test_old_summary_requires_regeneration_and_preferred_name(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            _, _, rows = self.make_project(base)
            row = rows["sim"]
            row.pop("T_up_s")
            path = self.write_summary(base, "sim", row)
            with self.assertRaisesRegex(ValueError, "rerun pcs-postprocess"):
                compare.read_summary(path, "sim")
            fallback = self.write_summary(base, "sim", row, "summary_freq_reg.csv")
            self.assertEqual(compare.summary_path(path.parent, "sim"), path)
            path.unlink()
            self.assertEqual(compare.summary_path(path.parent, "sim"), fallback)

    def test_duplicate_and_parameter_mismatch(self):
        for kind in ("duplicate", "mismatch", "missing_summary"):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory() as tmp:
                base = Path(tmp)
                paths, files, _ = self.make_project(base)
                if kind == "duplicate":
                    shutil.copyfile(files["sim"], files["sim"].with_name("case_1002_duplicate.mat"))
                elif kind == "mismatch":
                    path = files["sim"].with_name(files["sim"].stem + "_meta.json")
                    meta = json.loads(path.read_text())
                    meta["freq_deviation_hz"] = -.3
                    path.write_text(json.dumps(meta))
                else:
                    path = compare.summary_path(base / "sim/results_freq_reg/postprocess", "sim")
                    with path.open(encoding="utf-8-sig") as handle:
                        header = handle.readline()
                    path.write_text(header, encoding="utf-8-sig")
                result = compare.run(paths, base / "compare")
                self.assertEqual(result["exceptions"], 1)
                self.assertEqual(result["completed"], 0)

    def test_cli_dispatch(self):
        with patch.object(compare, "run", return_value={"completed": 1, "exceptions": 0}) as runner:
            self.assertEqual(main(["compare", "--mode", "freq_reg", "--case", "3019", "--output", "freq-test"]), 0)
        runner.assert_called_once_with(None, Path("freq-test"), 3019)


if __name__ == "__main__":
    unittest.main()
