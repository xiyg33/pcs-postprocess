"""个人目录配置与旧命令并存的回归检查。"""

import json
from pathlib import Path
import tempfile
import unittest

from pcs_postprocess.cli import main
from pcs_postprocess.paths import PathsConfig

from test_workflow import ROOT, TMP_PARENT, synthetic


class PathsTests(unittest.TestCase):
    def make_config(self, base: Path) -> Path:
        config = base / "postprocess_paths.local.json"
        config.write_text(json.dumps({
            "project_config": str(ROOT / "examples/project.json"),
            "sources": {
                "hil": {"results_root": "hil_results", "input_subdir": "canonical",
                        "output_subdir": "postprocess", "raw_subdir": "hil_raw"},
                "sim": {"results_root": "sim_results"},
            },
        }), encoding="utf-8")
        return config

    def test_profiles_find_all_seven_mode_directories(self):
        with tempfile.TemporaryDirectory(dir=TMP_PARENT) as tmp:
            base = Path(tmp)
            paths = PathsConfig.load(self.make_config(base))
            for name, root in (("hil", "hil_results"), ("sim", "sim_results")):
                for mode in synthetic.MODES:
                    profile = paths.source(name)
                    self.assertEqual(profile.input_dir(mode), base / root / f"results_{mode}" / "canonical")
                    self.assertEqual(profile.output_dir(mode), base / root / f"results_{mode}" / "postprocess")

    def test_profile_validate_run_and_explicit_overrides(self):
        with tempfile.TemporaryDirectory(dir=TMP_PARENT) as tmp:
            base = Path(tmp)
            paths_file = self.make_config(base)
            profile = PathsConfig.load(paths_file).source("hil")
            synthetic.create(profile.input_dir("frt"), "frt")
            # 目录可能含尚未迁移的其他编号，--case 只校验指定工况。
            (profile.input_dir("frt") / "case_9999_old.mat").write_bytes(b"legacy")
            common = ["--mode", "frt", "--source", "hil", "--paths-config", str(paths_file)]
            self.assertEqual(main(["validate", *common, "--case", "1000"]), 0)
            self.assertEqual(main(["run", *common, "--case", "1000", "--no-fig"]), 0)
            self.assertTrue((profile.output_dir("frt") / "summary_frt.csv").is_file())
            override_input = base / "other_input"
            synthetic.create(override_input, "frt")
            override_output = base / "other_output"
            self.assertEqual(main(["run", *common, "--input", str(override_input),
                                   "--output", str(override_output), "--config",
                                   str(ROOT / "examples/project.json"), "--no-fig"]), 0)
            self.assertTrue((override_output / "summary_frt.csv").is_file())

    def test_missing_config_and_unknown_source_fail(self):
        with tempfile.TemporaryDirectory(dir=TMP_PARENT) as tmp:
            base = Path(tmp)
            missing = base / "absent.json"
            with self.assertRaises(SystemExit) as error:
                main(["validate", "--mode", "frt", "--source", "hil",
                      "--paths-config", str(missing)])
            self.assertEqual(error.exception.code, 2)
            paths_file = self.make_config(base)
            with self.assertRaises(SystemExit) as error:
                main(["validate", "--mode", "frt", "--source", "unknown",
                      "--paths-config", str(paths_file)])
            self.assertEqual(error.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
