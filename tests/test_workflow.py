import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

from scipy.io import loadmat, savemat

from pcs_postprocess.cli import main
from pcs_postprocess.contract import validate_case

ROOT = Path(__file__).resolve().parents[1]
TMP_PARENT = ROOT / "outputs" / "test_tmp"
TMP_PARENT.mkdir(parents=True, exist_ok=True)
spec = importlib.util.spec_from_file_location("make_synthetic", ROOT / "examples/make_synthetic.py")
synthetic = importlib.util.module_from_spec(spec)
spec.loader.exec_module(synthetic)


class WorkflowTests(unittest.TestCase):
    def test_all_modes_produce_processed_data_and_figures(self):
        with tempfile.TemporaryDirectory(dir=TMP_PARENT) as tmp:
            base = Path(tmp)
            for mode in synthetic.MODES:
                with self.subTest(mode=mode):
                    raw = base / mode / "raw"
                    synthetic.create(raw, mode)
                    out = base / mode / "out"
                    self.assertEqual(main(["validate", "--mode", mode, "--input", str(raw)]), 0)
                    self.assertEqual(main(["run", "--mode", mode, "--input", str(raw),
                                           "--output", str(out), "--config",
                                           str(ROOT / "examples/project.json")]), 0)
                    self.assertEqual(len(list((out / "processed").glob("*.mat"))), 1)
                    self.assertEqual(len(list((out / "figures").glob("*.png"))),
                                     2 if mode in ("frt", "freq_reg", "inertia") else 1)
                    with (out / f"summary_{mode}.csv").open(encoding="utf-8-sig") as file:
                        self.assertEqual(len(list(csv.DictReader(file))), 1)
                    if mode == "pwr":
                        self.assertEqual(len(list((out / "downsampled").glob("*.csv"))), 1)

    def test_recorded_time_and_no_fig(self):
        with tempfile.TemporaryDirectory(dir=TMP_PARENT) as tmp:
            base = Path(tmp)
            raw = base / "raw"
            synthetic.create(raw, "frt", "recorded")
            out = base / "out"
            self.assertEqual(main(["run", "--mode", "frt", "--input", str(raw),
                                   "--output", str(out), "--config",
                                   str(ROOT / "examples/project.json"), "--case", "1000",
                                   "--no-fig"]), 0)
            self.assertFalse((out / "figures").exists())
            mat = loadmat(next((out / "processed").glob("*.mat")))
            self.assertGreater(mat["t_s"].size, 0)

    def test_absolute_time_preserves_nonzero_origin(self):
        with tempfile.TemporaryDirectory(dir=TMP_PARENT) as tmp:
            base = Path(tmp)
            raw = base / "raw"
            path = synthetic.create(raw, "frt", "absolute")
            mat = loadmat(path)
            mat["data"][:, 0] += 10.0
            savemat(path, {"data": mat["data"], "header": mat["header"]})
            meta_path = path.with_name(path.stem + "_meta.json")
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
            meta["fault_start_s"] += 10.0
            meta["sim_time_s"] += 10.0
            meta_path.write_text(json.dumps(meta), encoding="utf-8")
            out = base / "out"
            self.assertEqual(main(["run", "--mode", "frt", "--input", str(raw),
                                   "--output", str(out), "--config",
                                   str(ROOT / "examples/project.json"), "--no-fig"]), 0)
            result = loadmat(next((out / "processed").glob("*.mat")))
            # 绘图窗口把事件统一移到 1 s；绝对时钟平移后仍应找到该事件。
            self.assertAlmostEqual(float(result["event1_start_s"].squeeze()), 1.0)

    def test_missing_required_metadata_fails(self):
        with tempfile.TemporaryDirectory(dir=TMP_PARENT) as tmp:
            path = synthetic.create(Path(tmp), "frt")
            meta_path = path.with_name(path.stem + "_meta.json")
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
            del meta["fault_start_s"]
            meta_path.write_text(json.dumps(meta), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "fault_start_s"):
                validate_case(path)

    def test_csv_adapter(self):
        with tempfile.TemporaryDirectory(dir=TMP_PARENT) as tmp:
            base = Path(tmp)
            source = synthetic.create(base / "source", "frt")
            from scipy.io import loadmat
            mat = loadmat(source)
            names = synthetic.NAMES
            csv_path = base / "wave.csv"
            with csv_path.open("w", newline="", encoding="utf-8") as file:
                writer = csv.writer(file)
                writer.writerow(names)
                writer.writerows(mat["data"])
            spec = importlib.util.spec_from_file_location("adapter", ROOT / "adapters/csv_to_canonical.py")
            adapter = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(adapter)
            result = adapter.convert(csv_path, source.with_name(source.stem + "_meta.json"), base / "converted")
            self.assertEqual(validate_case(result)["mode"], "frt")


if __name__ == "__main__":
    unittest.main()
