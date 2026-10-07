"""Regression tests for referenced PWR summaries and complete recordings."""

import csv
import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from pcs_postprocess import compare_pwr as compare
from pcs_postprocess.cli import main as cli_main


PACKAGE = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("synthetic_compare", PACKAGE / "examples/make_synthetic.py")
synthetic = importlib.util.module_from_spec(spec)
spec.loader.exec_module(synthetic)


class ComparePwrTests(unittest.TestCase):
    def make_project(self, base):
        config = base / "paths.json"
        config.write_text(json.dumps({"project_config": str(PACKAGE / "examples/project.json"),
                                     "sources": {s: {"results_root": s} for s in ("sim", "hil")}}), encoding="utf-8")
        rows, files = {}, {}
        for side in ("sim", "hil"):
            mode_dir = base / side / "results_pwr"
            files[side] = synthetic.create(mode_dir / "canonical", "pwr", "absolute" if side == "sim" else "recorded")
            meta_path = files[side].with_name(files[side].stem + "_meta.json")
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
            meta.update(test_type="PWR_P_SINGLE", test_type_code=3)
            if side == "sim":
                meta["data_source"] = "sim"
            meta_path.write_text(json.dumps(meta), encoding="utf-8")
            rows[side] = {"case_index": 1001, "case_name": "synthetic_pwr", "test_type": "PWR_P_SINGLE",
                          "rec_offset_s": 0, "offset_source": "absolute_time" if side == "sim" else "edge",
                          "M_p_pct": 12.25 if side == "sim" else 10.0, "T_up_s": .1,
                          "T_p_s": .2, "T_s_s": .3, "M_p_reason": "", "T_up_reason": "",
                          "T_p_reason": "", "T_s_reason": "",
                          "step_timing_basis": "sim_command" if side == "sim" else "response_edge_estimate",
                          "step_response_start_s": 1, "T_s_method": "continuous_5pct",
                          "metric_sampling_ms": 5, "人工确认Ts": ""}
            self.write_summary(base, side, [rows[side]])
        return config, rows, files

    def write_summary(self, base, side, rows, filename=None):
        directory = base / side / "results_pwr/postprocess"
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / (filename or f"{side}_summary_pwr.csv")
        with path.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        return path

    def test_references_metrics_without_calculating_or_copying(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            config, _, files = self.make_project(base)
            summaries = [base / s / f"results_pwr/postprocess/{s}_summary_pwr.csv" for s in ("sim", "hil")]
            originals = [p.read_bytes() for p in summaries]
            output = base / "compare"
            with patch.object(compare.processing, "prepare_case", side_effect=AssertionError("must not calculate")), \
                 patch.object(compare.processing, "pwr_step_summary", side_effect=AssertionError("must not calculate")):
                manifest = compare.run(config, output)
            self.assertEqual(manifest["completed"], 1)
            self.assertEqual(manifest["exceptions"], 0)
            self.assertEqual({p.name for p in output.iterdir()}, {"figures", "compare_summary_pwr.csv", "manifest.json"})
            self.assertEqual([p.read_bytes() for p in summaries], originals)
            for side, path in zip(("sim", "hil"), summaries):
                self.assertEqual((output / manifest["summaries"][side]).resolve(), path.resolve())
                self.assertEqual((output / manifest["cases"][0]["canonical"][side]).resolve(), files[side].resolve())
                self.assertEqual(manifest["cases"][0]["recording_range_s"][side], [0., 3.998])
            with (output / "compare_summary_pwr.csv").open(encoding="utf-8-sig", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 4)
            self.assertEqual(float(rows[0]["deviation"]), 2.25)
            self.assertEqual(rows[0]["deviation_unit"], "percentage_points")
            self.assertEqual(rows[0]["status"], "valid")
            self.assertEqual(rows[1]["status"], "estimated")
            self.assertTrue((output / rows[0]["figure_png"]).is_file())
            self.assertTrue((output / "compare_summary_pwr.csv").read_bytes().startswith(b"\xef\xbb\xbf"))
            with self.assertRaises(FileExistsError):
                compare.run(config, output)

    def test_preferred_summary_and_legacy_column_aliases(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            _, rows, _ = self.make_project(base)
            legacy = dict(rows["sim"], M_p_pct=999)
            self.write_summary(base, "sim", [legacy], "summary_pwr.csv")
            selected = compare.summary_path(base / "sim/results_pwr/postprocess", "sim")
            self.assertEqual(compare.read_summary(selected, "sim")[1001][0]["M_p_pct"], "12.25")
            selected.unlink()
            self.assertEqual(compare.summary_path(selected.parent, "sim").name, "summary_pwr.csv")
            aliased = dict(rows["hil"])
            for metric, _, column, _, _ in compare.METRICS:
                aliased[column.replace(metric + "_", metric + "_hil_", 1)] = aliased.pop(column)
            path = self.write_summary(base, "hil", [aliased])
            self.assertEqual(compare.read_summary(path, "hil")[1001][0]["T_s_s"], "0.3")

    def test_empty_values_quality_and_q_signal(self):
        with tempfile.TemporaryDirectory() as temp:
            _, rows, _ = self.make_project(Path(temp))
            rows["sim"].update(M_p_pct="", M_p_reason="baseline_unstable")
            rows["hil"].update(T_s_method="fourth_band_entry_estimate", 人工确认Ts="待确认")
            output = compare.comparison_rows(1001, {"case_name": "q", "test_type": "PWR_Q_SINGLE", "test_type_code": 4}, rows, "figures/q.png")
            self.assertIsNone(output[0]["deviation"])
            self.assertEqual(output[0]["status"], "unavailable")
            self.assertEqual(output[3]["hil_status"], "estimated")
            self.assertEqual(output[3]["hil_manual_T_s"], "待确认")
            self.assertEqual(output[0]["signal"], "Q_pu")

    def test_display_starts_one_second_before_event_and_keeps_second_event(self):
        with tempfile.TemporaryDirectory() as temp:
            config, rows, files = self.make_project(Path(temp))
            compare.processing.configure_project(compare.ProjectConfig.load(PACKAGE / "examples/project.json"))
            meta = json.loads(files["sim"].with_name(files["sim"].stem + "_meta.json").read_text(encoding="utf-8"))
            meta.update(test_type_code=5, step2_start_s=3.5, step2_dur_s=.5, P_step2_pu=.3, Q_step2_pu=0)
            data = compare.full_signals(files["sim"], meta, rows["sim"])
            self.assertAlmostEqual(data["power"]["t_s"][0], 0.)
            self.assertAlmostEqual(data["t_s"][0], 0.)
            self.assertAlmostEqual(data["power"]["t_s"][-1], 3.995)
            np.testing.assert_allclose(np.diff(data["power"]["t_s"]), .005)
            self.assertEqual(data["events"]["event2_start"], 2.5)
            self.assertEqual(data["events"]["event2_end"], 3.)

    def test_plot_uses_two_source_colors_and_nonnegative_axis(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            _, rows, files = self.make_project(base)
            compare.processing.configure_project(compare.ProjectConfig.load(PACKAGE / "examples/project.json"))
            metas = {side: compare.validate_case(files[side], "pwr") for side in files}
            pair = {side: compare.full_signals(files[side], metas[side], rows[side]) for side in files}
            figures = []
            real_close = compare.plt.close
            try:
                with patch.object(compare.plt, "close", side_effect=figures.append):
                    compare.plot_pair(pair, metas, base / "plot.png", 1001)
                fig = figures[-1]
                for ax in fig.axes:
                    self.assertEqual(ax.get_xlim()[0], 0.)
                    self.assertEqual(ax.get_legend_handles_labels()[1], ["HIL", "SIM"])
                for ax in fig.axes[:2]:
                    self.assertEqual([line.get_color() for line in ax.lines[:6]], ["tab:red"] * 3 + ["tab:blue"] * 3)
                    self.assertTrue(all(np.min(line.get_xdata()) >= 0 for line in ax.lines[:6]))
            finally:
                for fig in figures:
                    if hasattr(fig, "axes"):
                        real_close(fig)

    def test_package_compare_command(self):
        result = {"completed": 1, "exceptions": 0}
        with patch.object(compare, "run", return_value=result) as runner:
            self.assertEqual(cli_main(["compare", "--mode", "pwr", "--case", "2001", "--output", "compare-test"]), 0)
        runner.assert_called_once_with(None, Path("compare-test"), 2001)

    def test_duplicate_missing_summary_and_mismatch_are_errors(self):
        for problem in ("duplicate", "missing_summary", "mismatch", "missing_source", "duplicate_input"):
            with self.subTest(problem=problem), tempfile.TemporaryDirectory() as temp:
                base = Path(temp)
                config, rows, files = self.make_project(base)
                if problem == "duplicate":
                    self.write_summary(base, "sim", [rows["sim"], rows["sim"]])
                elif problem == "missing_summary":
                    self.write_summary(base, "sim", [dict(rows["sim"], case_index=9999)])
                elif problem == "mismatch":
                    path = files["sim"].with_name(files["sim"].stem + "_meta.json")
                    meta = json.loads(path.read_text(encoding="utf-8"))
                    meta["P_step_pu"] = .9
                    path.write_text(json.dumps(meta), encoding="utf-8")
                elif problem == "missing_source":
                    files["sim"].unlink()
                else:
                    shutil.copyfile(files["sim"], files["sim"].with_name("case_1001_duplicate.mat"))
                output = base / "compare"
                self.assertEqual(compare.main(["--all", "--paths-config", str(config), "--output", str(output)]), 1)
                manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
                self.assertEqual(manifest["completed"], 0)
                self.assertEqual(manifest["exceptions"], 1)
                self.assertEqual(manifest["errors"][0]["case_index"], 1001)

    def test_batch_continues_after_an_unpaired_case(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            config, _, files = self.make_project(base)
            shutil.copyfile(files["sim"], files["sim"].with_name("case_1000_unpaired.mat"))
            output = base / "compare"
            with patch.object(compare, "plot_pair"):
                manifest = compare.run(config, output)
            self.assertEqual(manifest["completed"], 1)
            self.assertEqual(manifest["exceptions"], 1)
            self.assertEqual(manifest["cases"][0]["case_index"], 1001)
            self.assertEqual(manifest["errors"][0]["case_index"], 1000)


if __name__ == "__main__":
    unittest.main()
