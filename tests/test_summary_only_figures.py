"""Summary-only plotting policy for both recording time bases."""

import csv
import io
import json
from contextlib import redirect_stdout
from pathlib import Path
import shutil
import tempfile
import unittest

from pcs_postprocess.cli import _is_summary_only_case, main, run

from test_workflow import ROOT, TMP_PARENT, synthetic


class SummaryOnlyFigureTests(unittest.TestCase):
    def make_case(self, raw, mode, kind, time_basis, name_only=False):
        path = synthetic.create(raw, mode, time_basis)
        meta_path = path.with_name(path.stem + "_meta.json")
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        meta["case_name"] = kind + "_synthetic"
        if not name_only:
            meta["test_type"] = kind
        meta_path.write_text(json.dumps(meta), encoding="utf-8")
        return path, meta

    def read_rows(self, output, mode, source=None):
        prefix = source + "_" if source else ""
        with (output / f"{prefix}summary_{mode}.csv").open(encoding="utf-8-sig") as file:
            return list(csv.DictReader(file))

    def test_detection_is_mode_and_name_based(self):
        for mode, kind in (("freq_reg", "FREQ_REG_DEADBAND"),
                           ("frt", "LVRT_BND"), ("frt", "HVRT_BND")):
            with self.subTest(mode=mode, kind=kind):
                self.assertTrue(_is_summary_only_case(mode, {"test_type": kind.lower()}))
                self.assertTrue(_is_summary_only_case(mode, {"case_name": " " + kind.lower() + "_old "}))
                self.assertTrue(_is_summary_only_case(mode, {"test_type": "ordinary", "case_name": kind + "_old"}))
                for other in synthetic.MODES:
                    if other != mode:
                        self.assertFalse(_is_summary_only_case(other, {"test_type": kind, "case_name": kind + "_old"}))
        self.assertFalse(_is_summary_only_case("freq_reg", {
            "case_index": 3001, "test_type": "FREQ_REG", "case_name": "FREQ_REG_normal",
            "freq_deviation_hz": 0.005,
        }))
        self.assertFalse(_is_summary_only_case("frt", {"case_index": 1097, "test_type": "LVRT"}))

    def test_default_preserves_summary_and_processed_for_sim_and_hil(self):
        for source, time_basis in (("sim", "absolute"), ("hil", "recorded")):
            for mode, kind in (("freq_reg", "FREQ_REG_DEADBAND"),
                               ("frt", "LVRT_BND"), ("frt", "HVRT_BND")):
                for name_only in (False, True):
                    with self.subTest(source=source, kind=kind, name_only=name_only):
                        with tempfile.TemporaryDirectory(dir=TMP_PARENT) as tmp:
                            base = Path(tmp)
                            raw, out, baseline = base / "raw", base / "out", base / "baseline"
                            _, meta = self.make_case(raw, mode, kind, time_basis, name_only)
                            log = io.StringIO()
                            with redirect_stdout(log):
                                self.assertEqual(run(mode, raw, out, ROOT / "examples/project.json",
                                                     {meta["case_index"]}, save_processed=True, source=source), 0)
                                self.assertEqual(run(mode, raw, baseline, ROOT / "examples/project.json",
                                                     no_fig=True, source=source), 0)
                            self.assertIn("processed, summary only", log.getvalue())
                            self.assertFalse((out / "figures").exists())
                            self.assertEqual(len(list((out / "processed").glob("*.mat"))), 1)
                            rows = self.read_rows(out, mode, source)
                            self.assertEqual(len(rows), 1)
                            self.assertEqual(rows[0]["case_name"], meta["case_name"])
                            self.assertEqual(rows, self.read_rows(baseline, mode, source))
                            if mode == "freq_reg":
                                self.assertNotEqual(rows[0]["event_end_P_mean_pu"], "")

    def test_cli_override_and_no_fig_precedence(self):
        for source, time_basis in (("sim", "absolute"), ("hil", "recorded")):
            for mode, kind in (("freq_reg", "FREQ_REG_DEADBAND"),
                               ("frt", "LVRT_BND"), ("frt", "HVRT_BND")):
                with self.subTest(source=source, kind=kind):
                    with tempfile.TemporaryDirectory(dir=TMP_PARENT) as tmp:
                        base = Path(tmp)
                        raw = base / "raw"
                        _, meta = self.make_case(raw, mode, kind, time_basis)
                        common = ["run", "--mode", mode, "--input", str(raw), "--config",
                                  str(ROOT / "examples/project.json"), "--case", str(meta["case_index"])]
                        outputs = []
                        for label, flags, expected in (
                                ("default", [], 0),
                                ("included", ["--include-summary-only-figures"], 1),
                                ("disabled", ["--include-summary-only-figures", "--no-fig"], 0)):
                            out = base / label
                            log = io.StringIO()
                            with redirect_stdout(log):
                                self.assertEqual(main([*common, "--output", str(out), *flags]), 0)
                            self.assertEqual(len(list((out / "figures").glob("*.png"))), expected)
                            self.assertIn("processed with figures" if expected else "processed, summary only",
                                          log.getvalue())
                            outputs.append(self.read_rows(out, mode))
                        self.assertEqual(outputs[0], outputs[1])
                        self.assertEqual(outputs[0], outputs[2])

    def test_mixed_batch_keeps_ordinary_figures_and_all_rows(self):
        for mode, kind in (("freq_reg", "FREQ_REG_DEADBAND"), ("frt", "LVRT_BND")):
            with self.subTest(mode=mode):
                with tempfile.TemporaryDirectory(dir=TMP_PARENT) as tmp:
                    base = Path(tmp)
                    raw, out = base / "raw", base / "out"
                    ordinary = synthetic.create(raw, mode)
                    meta_path = ordinary.with_name(ordinary.stem + "_meta.json")
                    meta = json.loads(meta_path.read_text(encoding="utf-8"))
                    meta["case_index"] += 100
                    meta["case_name"] = kind + "_synthetic"
                    target = raw / f"case_{meta['case_index']}_{meta['case_name']}.mat"
                    shutil.copyfile(ordinary, target)
                    target.with_name(target.stem + "_meta.json").write_text(json.dumps(meta), encoding="utf-8")
                    self.assertEqual(run(mode, raw, out, ROOT / "examples/project.json"), 0)
                    figures = list((out / "figures").glob("*.png"))
                    self.assertEqual(len(figures), 2)
                    self.assertTrue(all(image.stem.startswith(ordinary.stem) for image in figures))
                    self.assertEqual(len(self.read_rows(out, mode)), 2)


if __name__ == "__main__":
    unittest.main()
