"""FRT segmentation, timing, skip policy, and public comparison regression."""

import csv
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from pcs_postprocess import compare_frt as compare
from pcs_postprocess.cli import main
from test_workflow import ROOT, synthetic


def pair_data(recovery=1.15, bias=0.):
    t = np.arange(0., 4., .005)
    data = {}
    for side in ("sim", "hil"):
        values = np.ones(t.size) + (bias if side == "sim" else 0.)
        data[side] = {"t_s": t, "display": {"t_s": t, **{s: values.copy() for s in compare.SIGNALS}},
                      "recovery_s": recovery, "recovery_basis": "detected",
                      "timing": {"offset_source": "absolute_time" if side == "sim" else "edge"}}
    return data


class FrtComparisonTests(unittest.TestCase):
    def make_project(self, base):
        config = base / "paths.json"
        config.write_text(json.dumps({"project_config": str(ROOT / "examples/project.json"),
                                     "sources": {side: {"results_root": side} for side in ("sim", "hil")}}), encoding="utf-8")
        files = {}
        for side in ("sim", "hil"):
            directory = base / side / "results_frt"
            path = synthetic.create(directory / "canonical", "frt", "absolute" if side == "sim" else "recorded")
            meta_path = path.with_name(path.stem + "_meta.json")
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
            meta.update(test_type="LVRT_SINGLE", fault_duration_s=.5)
            if side == "sim":
                meta["data_source"] = "sim"
            meta_path.write_text(json.dumps(meta), encoding="utf-8")
            out = directory / "postprocess"
            out.mkdir()
            summary = out / ("summary_frt.csv" if side == "sim" else "hil_summary_frt.csv")
            with summary.open("w", encoding="utf-8-sig", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["case_index", "case_name", "test_type", "rec_offset_s", "offset_source"])
                writer.writeheader()
                writer.writerow(dict(case_index=1000, case_name="synthetic_frt", test_type="LVRT_SINGLE",
                                     rec_offset_s=0, offset_source="absolute_time" if side == "sim" else "edge"))
            files[side] = path
        return config, files

    def test_short_fault_has_positive_nonoverlapping_windows(self):
        pair = pair_data()
        definitions = compare.windows(pair)
        indexed = {w["code"]: w for w in definitions}
        self.assertAlmostEqual(indexed["B1"]["end_s"], 1.0375)
        self.assertAlmostEqual(indexed["B2"]["end_s"], 1.135)
        self.assertEqual(indexed["B1"]["end_s"], indexed["B2"]["start_s"])
        self.assertEqual(indexed["C1"]["end_s"], indexed["C2"]["start_s"])
        rows, details = compare.deviations(pair)
        self.assertEqual(len(rows), 3)
        self.assertTrue(all(d["sample_count"] > 0 for d in details["Vpos_pu"].values()))
        self.assertTrue(all(row[c] == 0. for row in rows for c in compare.ALL_METRICS))

    def test_known_bias_and_weighted_overall_error(self):
        rows, _ = compare.deviations(pair_data(bias=.2))
        for row in rows:
            for column in compare.ALL_METRICS:
                self.assertAlmostEqual(row[column], .2)
            self.assertEqual(row["status"], "valid")
        pair = pair_data()
        for signal in compare.SIGNALS:
            pair["sim"]["display"][signal][:] = np.where(pair["sim"]["display"]["t_s"] < 1., 3., 1.)
        rows, details = compare.deviations(pair)
        selected = details["P_pu"].values()
        expected = sum(d["mae"] * d["sample_count"] for d in selected) / sum(d["sample_count"] for d in selected)
        self.assertAlmostEqual(rows[1]["fg_weighted_overall_mae"], expected)

    def test_recovery_alignment_removes_artificial_recovery_error(self):
        pair = pair_data(recovery=1.15)
        pair["sim"]["recovery_s"] = 1.17
        for side, data in pair.items():
            t = data["display"]["t_s"]
            for signal in compare.SIGNALS:
                # Linear recovery gives an exact interpolation reference on either time axis.
                data["display"][signal] = np.maximum(0., 1. - 10. * (t - data["recovery_s"]))
        aligned, details = compare.deviations(pair)
        onset, _ = compare.deviations(pair, recovery_aligned=False)
        column = "f2_c1_recovery_transient_mae"
        self.assertAlmostEqual(aligned[0][column], 0., places=12)
        self.assertGreater(onset[0][column], .1)
        assessment = compare.duration_assessment(pair, {"sim": {"fault_duration_s": .17}, "hil": {"fault_duration_s": .15}}, aligned, details)
        self.assertEqual(assessment["effect"], "affected")
        self.assertTrue(assessment["onset_only_metric_changes"]["Vpos_pu"][column]["changed"])

    def test_missing_samples_and_empty_windows_leave_metrics_blank(self):
        pair = pair_data()
        pair["sim"]["display"]["P_pu"][216] = np.nan
        rows, details = compare.deviations(pair)
        p = rows[1]
        self.assertIsNone(p["f1_b2_disturbed_steady_mean"])
        self.assertIsNone(p["fg_weighted_overall_mae"])
        self.assertEqual(p["f2_b1_rise_transient_mae"], 0.)
        self.assertEqual(details["P_pu"]["B2"]["reason"], "missing_samples")
        empty = {"code": "B2", "start_s": 1.2, "end_s": 1.1, "origins_s": {"sim": 0., "hil": 0.}}
        self.assertEqual(compare.segment_metrics(pair, "P_pu", empty)["reason"], "empty_common_window")

    def test_metadata_recovery_is_estimated(self):
        pair = pair_data()
        pair["hil"]["recovery_basis"] = "metadata_estimate"
        rows, _ = compare.deviations(pair)
        self.assertTrue(all(row["status"] == "estimated" for row in rows))

    def test_run_references_sources_and_produces_two_figures(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            config, files = self.make_project(base)
            originals = {s: compare.summary_path(p.parent.parent / "postprocess", s).read_bytes() for s, p in files.items()}
            output = base / "compare"
            with patch.object(compare.processing, "prepare_case", side_effect=AssertionError("must not rewrite summary")):
                self.assertEqual(main(["compare", "--mode", "frt", "--case", "1000", "--paths-config", str(config), "--output", str(output)]), 0)
            manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["completed"], 1)
            self.assertEqual(manifest["skipped"], 0)
            self.assertEqual(manifest["exceptions"], 0)
            self.assertEqual(len(list((output / "figures").glob("*.png"))), 2)
            with (output / "compare_summary_frt.csv").open(encoding="utf-8-sig", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 3)
            self.assertEqual(set(row["signal"] for row in rows), set(compare.SIGNALS))
            for s in originals:
                self.assertEqual((output / manifest["summaries"][s]).read_bytes(), originals[s])
            self.assertEqual(set(p.name for p in output.iterdir()), {"figures", "compare_summary_frt.csv", "manifest.json"})
            with self.assertRaises(FileExistsError):
                compare.run(config, output)

    def test_scan_and_missing_sim_are_successful_skips(self):
        for reason in ("boundary_scan_excluded", "missing_sim"):
            with self.subTest(reason=reason), tempfile.TemporaryDirectory() as temp:
                base = Path(temp)
                config, files = self.make_project(base)
                if reason == "missing_sim":
                    files["sim"].unlink()
                else:
                    for path in files.values():
                        meta_path = path.with_name(path.stem + "_meta.json")
                        meta = json.loads(meta_path.read_text(encoding="utf-8"))
                        meta["test_type"] = "LVRT_BND"
                        meta_path.write_text(json.dumps(meta), encoding="utf-8")
                output = base / "compare"
                self.assertEqual(main(["compare", "--mode", "frt", "--all", "--paths-config", str(config), "--output", str(output)]), 0)
                manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
                self.assertEqual(manifest["completed"], 0)
                self.assertEqual(manifest["skipped"], 1)
                self.assertEqual(manifest["skips"][0]["reason"], reason)

    def test_parameter_mismatch_and_duplicate_are_errors(self):
        for problem in ("parameters", "duplicate"):
            with self.subTest(problem=problem), tempfile.TemporaryDirectory() as temp:
                base = Path(temp)
                config, files = self.make_project(base)
                if problem == "duplicate":
                    shutil.copyfile(files["sim"], files["sim"].with_name("case_1000_duplicate.mat"))
                else:
                    path = files["sim"].with_name(files["sim"].stem + "_meta.json")
                    meta = json.loads(path.read_text(encoding="utf-8"))
                    meta["V_fault_pu"] = .7
                    path.write_text(json.dumps(meta), encoding="utf-8")
                output = base / "compare"
                self.assertEqual(main(["compare", "--mode", "frt", "--all", "--paths-config", str(config), "--output", str(output)]), 1)
                manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
                self.assertEqual(manifest["exceptions"], 1)

    def test_plot_layout_colors_and_invalid_current(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            config, files = self.make_project(base)
            compare.processing.configure_project(compare.ProjectConfig.load(ROOT / "examples/project.json"))
            pair = {side: compare.prepare_signals(path, compare.validate_case(path, "frt"),
                    {"rec_offset_s": "0", "offset_source": "edge" if side == "hil" else "absolute_time"}) for side, path in files.items()}
            # Exercise the low-voltage mask independently of the synthetic 0.5 pu fault.
            voltage = np.zeros((10, 3))
            ip, iq, valid = compare.processing.calculate_instantaneous_current(voltage, voltage)
            self.assertFalse(valid.any())
            self.assertTrue(np.isnan(ip).all() and np.isnan(iq).all())
            figures = []
            real_close = compare.plt.close
            try:
                with patch.object(compare.plt, "close", side_effect=figures.append):
                    compare.plot_pair(pair, (base / "main.png", base / "power.png"), 1000, "synthetic_frt")
                rendered = [fig for fig in figures if hasattr(fig, "axes")]
                self.assertEqual([len(fig.axes) for fig in rendered], [3, 5])
                for fig in rendered:
                    for ax in fig.axes:
                        self.assertEqual(ax.get_xlim()[0], 0.)
                        handles, labels = ax.get_legend_handles_labels()
                        self.assertEqual(labels, ["HIL", "SIM"])
                        self.assertEqual([line.get_color() for line in handles], ["tab:red", "tab:blue"])
            finally:
                for fig in figures:
                    if hasattr(fig, "axes"):
                        real_close(fig)


if __name__ == "__main__":
    unittest.main()
