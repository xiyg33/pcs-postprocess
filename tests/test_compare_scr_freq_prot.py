"""SCR paired horizons and voltage-phase protection timing regressions."""

import csv
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from scipy.io import savemat

from pcs_postprocess import frequency_protection as frequency, scr_stability as scr, processing as p
from pcs_postprocess import compare_scr, compare_freq_prot
from pcs_postprocess.cli import main
from pcs_postprocess.config import ProjectConfig


PACKAGE = Path(__file__).resolve().parents[1]
TMP_PARENT = PACKAGE / "outputs/test_tmp"
TMP_PARENT.mkdir(parents=True, exist_ok=True)


def abc(t, start, change, amplitude=1., harmonic=0.):
    angle = 2 * np.pi * (50 * t + change * np.maximum(t - start, 0.))
    return np.column_stack([p.HIL_V_BASE_PEAK_V * amplitude *
                            (np.sin(angle - phase * 2 * np.pi / 3)
                             + harmonic * np.sin(5 * angle + phase * 2 * np.pi / 3)) for phase in range(3)])


class FrequencyProtectionTests(unittest.TestCase):
    def setUp(self):
        p.configure_project(ProjectConfig.load(PACKAGE / "examples/project.json"))
        self.t = np.arange(0., 5.001, .001)

    def measure(self, change=2., start=2.13, trip=2.4, pll_delay=.2, shift=0., breaker=True, harmonic=0., t=None):
        t = self.t if t is None else t
        voltages = abc(t, start, change, harmonic=harmonic)
        b = np.ones(t.shape)
        if trip is not None:
            b[t >= trip] = 0.
        sig = {"V" + phase: voltages[:, i] for i, phase in enumerate("abc")}
        sig.update(breaker=b if breaker else None, f=50. + np.where(t >= start + pll_delay, change, 0.))
        meta = {"mode": "freq_prot", "time_basis": "absolute", "event_start_s": 2. + shift,
                "event_hold_s": 1., "f_target_hz": 50. + change,
                "frequency_protection_reference": {"over": {"threshold_hz": 51.5, "delay_s": .1},
                                                   "under": {"threshold_hz": 46.5, "delay_s": .1}}}
        return frequency.measure(t + shift, sig, meta, p.event_times_from_meta(meta), 0.)[0]

    def test_up_down_noise_and_disturbance_offset(self):
        for change in (2., -4.):
            for harmonic in (0., .05):
                with self.subTest(change=change, harmonic=harmonic):
                    result = self.measure(change=change, harmonic=harmonic)
                    self.assertEqual(result["protection_status"], "tripped")
                    self.assertAlmostEqual(result["protection_event_start_input_s"], 2.13, delta=.008)
                    self.assertAlmostEqual(result["delta_t_s"], .27, delta=.008)
                    self.assertAlmostEqual(result["freq_voltage_start_input_s"], 2.13, delta=.008)

    def test_pll_changes_do_not_select_disturbance_clock(self):
        a, b = self.measure(pll_delay=.05), self.measure(pll_delay=.8)
        self.assertEqual(a["delta_t_s"], b["delta_t_s"])
        self.assertEqual(a["protection_event_start_input_s"], b["protection_event_start_input_s"])

    def test_pll_crossing_after_trip_keeps_total_response_time(self):
        result = self.measure(pll_delay=.4)
        self.assertEqual(result["protection_status"], "tripped")
        self.assertAlmostEqual(result["delta_t_s"], .27, delta=.008)
        self.assertAlmostEqual(result["freq_breaker_trip_input_s"], 2.4)

    def test_unconfirmed_voltage_never_uses_pll_or_nominal_start(self):
        t = self.t
        v = abc(t, 2., 0.)
        sig = {"V" + phase: v[:, i] for i, phase in enumerate("abc")}
        sig.update(f=np.where(t >= 2.2, 52., 50.), breaker=np.where(t >= 2.4, 0., 1.))
        meta = {"mode": "freq_prot", "time_basis": "absolute", "event_start_s": 2.,
                "event_hold_s": 1., "f_target_hz": 52.}
        result = frequency.measure(t, sig, meta, p.event_times_from_meta(meta), 0.)[0]
        self.assertEqual(result["protection_status"], "disturbance_unconfirmed")
        self.assertIsNone(result["delta_t_s"])
        self.assertIsNone(result["t1_s"])

    def test_time_translation_invariance(self):
        a, b = self.measure(), self.measure(shift=1000.)
        self.assertAlmostEqual(a["delta_t_s"], b["delta_t_s"], places=7)
        self.assertAlmostEqual(b["protection_event_start_input_s"] - a["protection_event_start_input_s"], 1000., places=7)

    def test_no_trip_missing_breaker_and_unconfirmed_voltage(self):
        self.assertEqual(self.measure(trip=None)["protection_status"], "no_trip")
        missing = self.measure(breaker=False)
        self.assertEqual(missing["protection_status"], "missing_breaker")
        self.assertIsNone(missing["delta_t_s"])
        flat = self.measure(change=0.)
        self.assertEqual(flat["protection_status"], "disturbance_unconfirmed")
        self.assertIsNone(flat["t1_s"])

    def test_breaker_glitch_preopen_and_short_recording(self):
        t = self.t
        v = abc(t, 2.13, 2.)
        meta = {"mode": "freq_prot", "time_basis": "absolute", "event_start_s": 2., "event_hold_s": 1., "f_target_hz": 52.}
        meta["frequency_protection_reference"] = {"over": {"threshold_hz": 51.5, "delay_s": .1}}
        sig = {"V" + phase: v[:, i] for i, phase in enumerate("abc")}
        b = np.ones(t.shape)
        b[(t >= 2.2) & (t < 2.204)] = 0.
        b[t >= 2.4] = 0.
        sig.update(f=np.where(t >= 2.15, 52., 50.), breaker=b)
        result, _ = frequency.measure(t, sig, meta, p.event_times_from_meta(meta), 0.)
        self.assertAlmostEqual(result["protection_trip_input_s"], 2.4)
        sig["breaker"] = np.zeros(t.shape)
        self.assertEqual(frequency.measure(t, sig, meta, p.event_times_from_meta(meta), 0.)[0]["protection_status"], "already_open")
        self.assertEqual(self.measure(trip=None, t=t[t < 2.8])["protection_status"], "insufficient_observation")

    def test_phase_fit_bad_direction_boundary_and_gaps(self):
        t = self.t
        br = {"falls": []}
        bad = frequency.fit_start(t, abc(t, 2.13, -2.), 2., 52., br)
        self.assertIsNone(bad["start"])
        boundary = frequency.fit_start(t, abc(t, 2.4, 2.), 2., 52., br)
        self.assertIsNone(boundary["start"])
        gapped = t[(t < 2.04) | (t > 2.06)]
        self.assertIsNone(frequency.fit_start(gapped, abc(gapped, 2.13, 2.), 2., 52., br)["start"])

    def test_protection_settings_are_not_required_for_line_clock(self):
        t = self.t
        v = abc(t, 2.13, 2.)
        sig = {"V" + phase: v[:, i] for i, phase in enumerate("abc")}
        sig.update(f=None, breaker=np.where(t >= 2.4, 0., 1.))
        meta = {"mode": "freq_prot", "time_basis": "absolute", "event_start_s": 2.,
                "event_hold_s": 1., "f_target_hz": 52.}
        result = frequency.measure(t, sig, meta, p.event_times_from_meta(meta), 0.)[0]
        self.assertEqual(result["protection_status"], "tripped")
        self.assertAlmostEqual(result["delta_t_s"], .27, delta=.008)


class ScrStabilityTests(unittest.TestCase):
    def setUp(self):
        p.configure_project(ProjectConfig.load(PACKAGE / "examples/project.json"))
        self.t = np.arange(0., 5.001, .005)
        self.summary = {"scr_trip_times_input_s": "[]", "scr_incomplete_trip_times_input_s": "[]",
                        "scr_event_start_input_s": 2., "scr_recovery_input_s": 2.2,
                        "scr_timing_status": "confirmed", "scr_timing_reason": "", "scr_pre_disconnected": 0,
                        "scr_breaker_status": "valid", "scr_breaker_reason": "", "scr_input_sampling_s": .005}

    def classify(self, power, horizon=2.5, summary=None):
        data = {"display": {"t_s": self.t, "P_pu": power}, "timing": {"event_start_input_s": 2.},
                "input_t_s": self.t + 1.}
        return scr.classify(data, summary or self.summary, horizon)

    def test_growing_decaying_and_noise_only(self):
        elapsed = np.maximum(self.t - 1.5, 0.)
        growing = 1. + (.02 + .08 * elapsed) * np.sin(2 * np.pi * 5 * self.t)
        decaying = 1. + .15 * np.exp(-elapsed) * np.sin(2 * np.pi * 5 * self.t)
        noise = 1. + .0001 * np.sin(2 * np.pi * 5 * self.t) * (1 + elapsed)
        self.assertEqual(self.classify(growing)["state"], "UNSTABLE")
        self.assertEqual(self.classify(decaying)["state"], "STABLE")
        self.assertEqual(self.classify(noise)["state"], "STABLE")

    def test_overlimit_requires_twenty_ms_and_common_horizon(self):
        power = np.ones(self.t.shape)
        power[(self.t >= 3.4) & (self.t <= 3.405)] = 2.2
        self.assertNotEqual(self.classify(power)["reason"], "power_overlimit")
        power[(self.t >= 3.4) & (self.t <= 3.44)] = 2.2
        self.assertEqual(self.classify(power)["reason"], "power_overlimit")
        self.assertEqual(self.classify(power, 1.5)["state"], "STABLE")

    def test_pretrip_survives_unknown_timing_or_short_window(self):
        summary = dict(self.summary, scr_trip_times_input_s="[1.2]", scr_timing_status="unconfirmed")
        result = self.classify(np.ones(self.t.shape), horizon=.5, summary=summary)
        self.assertEqual(result["state"], "TRIP")
        self.assertEqual(result["reason"], "pre_disturbance_trip")
        self.assertEqual(self.classify(np.ones(self.t.shape), horizon=.5)["state"], "UNDETERMINED")

    def test_short_pulses_and_invalid_breaker(self):
        from pcs_postprocess.breaker_events import evidence
        b = np.ones(self.t.shape)
        b[(self.t >= 1.) & (self.t <= 1.01)] = 0.
        self.assertEqual(evidence(self.t, b)["falls"], [])
        self.assertEqual(evidence(self.t, None)["status"], "missing_breaker")

    def test_sparse_single_sided_nonmonotonic_and_unknown_boundaries(self):
        def points(values):
            return [{"scr": s, "state": state} for s, state in values]
        result = scr.boundary(points([(3., "STABLE"), (2.8, "STABLE"), (2.7, "TRIP"), (1.5, "TRIP")]))
        self.assertEqual((result["lower_scr"], result["upper_scr"]), (2.7, 2.8))
        self.assertEqual(scr.boundary(points([(3., "STABLE")]))["status"], "stable_side_only")
        self.assertEqual(scr.boundary(points([(3., "TRIP"), (2.8, "STABLE")]))["status"], "nonmonotonic")
        self.assertEqual(scr.boundary(points([(3., "STABLE"), (2.9, "UNDETERMINED"), (2.8, "TRIP")]))["status"], "incomplete_boundary")

    def test_local_voltage_drift_and_startup_do_not_replace_pulse(self):
        t = np.arange(0., 5.001, .001)
        v = 1. + .01 * t
        v[(t >= 2.13) & (t < 2.33)] -= .1
        sig = {"Vp_pu": v, "breaker": np.ones(t.shape)}
        meta = {"mode": "scr", "time_basis": "recorded", "disturbance_start_s": 2.,
                "disturbance_duration_s": .2, "V_target_pu": .9, "scr_event_search_center_input_s": 2.13}
        timing = scr.timing(t, sig, meta, p.event_times_from_meta(meta), 0.)
        self.assertAlmostEqual(timing["scr_event_start_input_s"], 2.13, delta=.012)
        self.assertAlmostEqual(timing["scr_recovery_input_s"], 2.33, delta=.012)


class ComparisonWorkflowTests(unittest.TestCase):
    def make_project(self, base, mode, ids=(1,)):
        p.configure_project(ProjectConfig.load(PACKAGE / "examples/project.json"))
        config = base / "paths.json"
        config.write_text(json.dumps({"project_config": str(PACKAGE / "examples/project.json"),
                                      "sources": {side: {"results_root": side} for side in ("sim", "hil")}}), encoding="utf-8")
        for side in ("sim", "hil"):
            for ci in ids:
                folder = base / side / f"results_{mode}/canonical"
                folder.mkdir(parents=True, exist_ok=True)
                t = np.arange(0., 5.001, .001)
                v = abc(t, 2.1, 2. if mode == "freq_prot" else 0.)
                b = np.ones(t.shape)
                if mode == "freq_prot":
                    b[t >= (2.4 if side == "sim" else 2.5)] = 0.
                else:
                    v[(t >= 2.1) & (t < 2.3)] *= .9
                names = ["t_s", "Va_V", "Vb_V", "Vc_V", "Ia_A", "Ib_A", "Ic_A", "P_W", "Q_var", "f_hz", "breaker"]
                f = np.where(t >= 2.2, 52., 50.) if mode == "freq_prot" else np.full(t.shape, 50.)
                data = np.column_stack((t, v, v * .05, np.full(t.shape, 24000.), np.zeros(t.shape), f, b))
                header = np.empty((len(names), 2), object)
                for i, name in enumerate(names):
                    header[i] = (name, "")
                path = folder / f"case_{ci:04d}_example.mat"
                savemat(path, {"data": data, "header": header})
                meta = {"mode": mode, "time_basis": "absolute", "case_index": ci, "case_name": "example",
                        "P_ref_pu": 1., "record_post_s": 2.}
                if mode == "freq_prot":
                    meta.update(event_start_s=2., event_hold_s=1., f_target_hz=52., test_type="FREQ_PROT_UP", test_type_code=12)
                    meta["frequency_protection_reference"] = {"over": {"threshold_hz": 51.5, "delay_s": .1},
                                                              "sha256": "fixture"}
                else:
                    meta.update(disturbance_start_s=2.1, disturbance_duration_s=.2, SCR_value=3.1 - ci * .1,
                                V_target_pu=.9, test_type="SCR_DISCHARGE", test_type_code=14, sim_time_s=4.3, record_pre_s=2.1)
                path.with_name(path.stem + "_meta.json").write_text(json.dumps(meta), encoding="utf-8")
        for side in ("sim", "hil"):
            self.assertEqual(main(["run", "--mode", mode, "--source", side, "--paths-config", str(config), "--no-fig"]), 0)
        return config

    def test_freq_compare_reads_saved_metrics_and_outputs_diagnostics(self):
        with tempfile.TemporaryDirectory(dir=TMP_PARENT) as tmp:
            base = Path(tmp)
            config = self.make_project(base, "freq_prot")
            with patch.object(frequency, "measure", side_effect=AssertionError("compare must not refit timing")):
                self.assertEqual(main(["compare", "--mode", "freq_prot", "--all", "--paths-config", str(config), "--output", str(base / "compare")]), 0)
            with (base / "compare/compare_summary_freq_prot.csv").open(encoding="utf-8-sig") as handle:
                row = next(csv.DictReader(handle))
            self.assertAlmostEqual(float(row["deviation"]), -.1, places=6)
            self.assertEqual(len(list((base / "compare/figures").glob("*.png"))), 1)
            with (base / "compare/timing_diagnostics_freq_prot.csv").open(encoding="utf-8-sig") as handle:
                self.assertEqual(len(list(csv.DictReader(handle))), 2)

    def test_scr_only_hil_present_cases_and_boundary_output(self):
        with tempfile.TemporaryDirectory(dir=TMP_PARENT) as tmp:
            base = Path(tmp)
            config = self.make_project(base, "scr", ids=(1, 2))
            next((base / "hil/results_scr/canonical").glob("case_0002_*.mat")).unlink()
            self.assertEqual(main(["compare", "--mode", "scr", "--all", "--paths-config", str(config), "--output", str(base / "compare")]), 0)
            manifest = json.loads((base / "compare/manifest.json").read_text())
            self.assertEqual(manifest["completed"], 1)
            self.assertEqual(len(manifest["cases"]), 1)
            self.assertTrue((base / "compare/compare_boundary_scr.csv").is_file())

    def test_old_summary_and_duplicate_cases_report_errors_and_continue(self):
        with tempfile.TemporaryDirectory(dir=TMP_PARENT) as tmp:
            base = Path(tmp)
            config = self.make_project(base, "freq_prot", ids=(1, 2))
            summary = base / "hil/results_freq_prot/postprocess/hil_summary_freq_prot.csv"
            summary.write_text("case_index,case_name,test_type\n1,example,FREQ_PROT_UP\n", encoding="utf-8")
            self.assertEqual(main(["compare", "--mode", "freq_prot", "--all", "--paths-config", str(config), "--output", str(base / "old")]), 1)
            manifest = json.loads((base / "old/manifest.json").read_text())
            self.assertEqual(manifest["exceptions"], 2)
            main(["run", "--mode", "freq_prot", "--source", "hil", "--paths-config", str(config), "--no-fig"])
            file = next((base / "hil/results_freq_prot/canonical").glob("case_0001_*.mat"))
            file.with_name(file.stem + "_duplicate.mat").write_bytes(file.read_bytes())
            self.assertEqual(main(["compare", "--mode", "freq_prot", "--all", "--paths-config", str(config), "--output", str(base / "duplicate")]), 1)
            self.assertEqual(json.loads((base / "duplicate/manifest.json").read_text(encoding="utf-8"))["completed"], 1)

    def test_pll_crossing_after_trip_keeps_comparison_values(self):
        with tempfile.TemporaryDirectory(dir=TMP_PARENT) as tmp:
            base = Path(tmp)
            config = self.make_project(base, "freq_prot")
            from scipy.io import loadmat
            file = next((base / "hil/results_freq_prot/canonical").glob("*.mat"))
            mat = loadmat(file)
            mat["data"][:, -2] = np.where(mat["data"][:, 0] >= 2.6, 52., 50.)
            savemat(file, {"data": mat["data"], "header": mat["header"]})
            main(["run", "--mode", "freq_prot", "--source", "hil", "--paths-config", str(config), "--no-fig"])
            self.assertEqual(main(["compare", "--mode", "freq_prot", "--all", "--paths-config", str(config), "--output", str(base / "compare")]), 0)
            with (base / "compare/compare_summary_freq_prot.csv").open(encoding="utf-8-sig") as handle:
                row = next(csv.DictReader(handle))
            self.assertEqual(row["hil_protection_status"], "tripped")
            self.assertAlmostEqual(float(row["deviation"]), -.1, delta=.005)
            self.assertAlmostEqual(float(row["hil_value"]), .4, delta=.005)

    def test_default_compare_archives_timestamp_runs_and_preserves_sources(self):
        with tempfile.TemporaryDirectory(dir=TMP_PARENT) as tmp:
            base = Path(tmp)
            config = self.make_project(base, "freq_prot")
            source = base / "hil/results_freq_prot/postprocess"
            figures = source / "figures"
            figures.mkdir()
            single_source = figures / "case_0001_example.png"
            single_source.write_bytes(b"single-source picture")
            original = (source / "hil_summary_freq_prot.csv").read_bytes()
            with patch("pcs_postprocess.compare.Path.cwd", return_value=base):
                for _ in range(2):
                    self.assertEqual(main(["compare", "--mode", "freq_prot", "--all", "--paths-config", str(config)]), 0)
            runs = sorted((base / "result_compare/freq_prot").iterdir())
            self.assertEqual(len(runs), 2)
            for output in runs:
                self.assertRegex(output.name, r"^\d{4}-\d{2}-\d{2}_\d{6}_\d+$")
                self.assertEqual([p.name for p in output.glob("*.json")], ["manifest.json"])
                self.assertTrue((output / "figures/case_0001_sim_hil_overlay.png").is_file())
                self.assertNotIn("freq_settings_reference", (output / "timing_notes.md").read_text(encoding="utf-8"))
            self.assertEqual((source / "hil_summary_freq_prot.csv").read_bytes(), original)
            self.assertEqual(single_source.read_bytes(), b"single-source picture")
            self.assertFalse((source / "compare_summary_freq_prot.csv").exists())

    def test_scr_representative_gallery_and_full_capture_trip_evidence(self):
        with tempfile.TemporaryDirectory(dir=TMP_PARENT) as tmp:
            base = Path(tmp)
            config = self.make_project(base, "scr", ids=(1, 2, 3))
            settings = json.loads(config.read_text())
            settings["compare_cases"] = {"scr": [1, 3]}
            config.write_text(json.dumps(settings), encoding="utf-8")
            from scipy.io import loadmat
            file = next((base / "hil/results_scr/canonical").glob("case_0001_*.mat"))
            mat = loadmat(file)
            mat["data"][:, -1] = np.where(mat["data"][:, 0] >= 1.5, 0., 1.)
            savemat(file, {"data": mat["data"], "header": mat["header"]})
            main(["run", "--mode", "scr", "--source", "hil", "--paths-config", str(config), "--no-fig"])
            output = base / "selected_compare"
            figures = output / "figures"
            renderer = compare_scr.plot_pair

            def render(pair, results, meta, path, horizon):
                self.assertEqual(pair["hil"]["recording"]["t_s"][0], 0.)
                if meta["case_index"] == 1:
                    self.assertEqual(results["hil"]["reason"], "pre_disturbance_trip")
                    self.assertGreater(results["hil"]["trip_input_s"], pair["hil"]["recording"]["input_origin_s"])
                return renderer(pair, results, meta, path, horizon)

            with patch.object(compare_scr, "plot_pair", side_effect=render), patch("pcs_postprocess.compare.Path.cwd", return_value=base):
                self.assertEqual(main(["compare", "--mode", "scr", "--all", "--paths-config", str(config), "--output", str(output)]), 0)
            self.assertEqual({path.name for path in figures.glob("*.png")}, {"case_0001_sim_hil_overlay.png", "case_0003_sim_hil_overlay.png"})
            manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["selected_case_ids"], [1, 3])
            self.assertEqual(manifest["completed"], 2)
            with (output / "compare_summary_scr.csv").open(encoding="utf-8-sig") as handle:
                self.assertEqual({int(r["case_index"]) for r in csv.DictReader(handle)}, {1, 3})

    def test_single_source_scr_picture_keeps_early_trip_and_source_only_legend(self):
        with tempfile.TemporaryDirectory(dir=TMP_PARENT) as tmp:
            base = Path(tmp)
            self.make_project(base, "scr")
            from scipy.io import loadmat
            from pcs_postprocess import plot_scr
            file = next((base / "hil/results_scr/canonical").glob("*.mat"))
            mat = loadmat(file)
            mat["data"][:, -1] = np.where(mat["data"][:, 0] >= 1.5, 0., 1.)
            savemat(file, {"data": mat["data"], "header": mat["header"]})
            meta = json.loads(file.with_name(file.stem + "_meta.json").read_text())
            case = p.prepare_case(str(file), meta, 0., "absolute_time")
            output = base / "single_source"
            output.mkdir()
            close = plot_scr.plt.close
            try:
                with patch.object(plot_scr.plt, "close"):
                    plot_scr.render_case(case, str(output))
                    figure = plot_scr.plt.gcf()
                    self.assertEqual(figure.axes[0].lines[0].get_xdata()[0], 0.)
                    self.assertIn("不正常（扰动前脱网）", figure._suptitle.get_text())
                    self.assertTrue(any(np.allclose(line.get_xdata(), [1.5, 1.5])
                                        for line in figure.axes[0].lines if len(line.get_xdata()) == 2))
                    self.assertFalse(any(line.get_label() in ("SIM", "HIL") for ax in figure.axes for line in ax.lines))
            finally:
                close("all")

    def test_scr_retested_timing_cannot_mix_with_old_two_second_source(self):
        with tempfile.TemporaryDirectory(dir=TMP_PARENT) as tmp:
            base = Path(tmp)
            config = self.make_project(base, "scr")
            meta_file = next((base / "hil/results_scr/canonical").glob("*_meta.json"))
            meta = json.loads(meta_file.read_text())
            meta.update(disturbance_start_s=10., sim_time_s=12.2, record_pre_s=10.)
            meta_file.write_text(json.dumps(meta), encoding="utf-8")
            summary = base / "hil/results_scr/postprocess/hil_summary_scr.csv"
            with summary.open(encoding="utf-8-sig", newline="") as handle:
                reader = csv.DictReader(handle)
                fields = reader.fieldnames
                rows = list(reader)
            rows[0].update(disturbance_start_s="10.0", sim_time_s="12.2", record_pre_s="10.0")
            with summary.open("w", encoding="utf-8-sig", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerows(rows)
            output = base / "mismatch"
            self.assertEqual(main(["compare", "--mode", "scr", "--all", "--paths-config", str(config), "--output", str(output)]), 1)
            manifest = json.loads((output / "manifest.json").read_text())
            self.assertEqual(manifest["completed"], 0)
            self.assertIn("disturbance_start_s", manifest["errors"][0]["reason"])

    def test_retested_scr_picture_uses_configured_ten_second_marker_and_end(self):
        p.configure_project(ProjectConfig.load(PACKAGE / "examples/project.json"))
        with tempfile.TemporaryDirectory(dir=TMP_PARENT) as tmp:
            root = Path(tmp)
            t = np.arange(0., 12.201, .001)
            voltages = abc(t, 10., 0.)
            voltages[(t >= 10.) & (t < 10.2)] *= .9
            names = ["t_s", "Va_V", "Vb_V", "Vc_V", "Ia_A", "Ib_A", "Ic_A", "P_W", "Q_var", "f_hz", "breaker"]
            header = np.empty((len(names), 2), object)
            for i, name in enumerate(names):
                header[i] = (name, "")
            file = root / "case_7001_example.mat"
            savemat(file, {"data": np.column_stack((t, voltages, voltages * .05, np.full(t.shape, 24000.),
                                                   np.zeros(t.shape), np.full(t.shape, 50.), np.ones(t.shape))), "header": header})
            meta = {"mode": "scr", "time_basis": "absolute", "case_index": 7001, "case_name": "example",
                    "test_type": "SCR_DISCHARGE", "test_type_code": 14, "disturbance_start_s": 10.,
                    "disturbance_duration_s": .2, "sim_time_s": 12.2, "record_pre_s": 10., "record_post_s": 2.,
                    "P_ref_pu": 1., "SCR_value": 3., "V_target_pu": .9}
            case = p.prepare_case(str(file), meta, 0., "absolute_time")
            from pcs_postprocess import plot_scr
            close = plot_scr.plt.close
            try:
                with patch.object(plot_scr.plt, "close"):
                    plot_scr.render_case(case, str(root))
                    figure = plot_scr.plt.gcf()
                    self.assertAlmostEqual(figure.axes[-1].get_xlim()[1], 12.2)
                    self.assertTrue(any(np.allclose(line.get_xdata(), [10., 10.])
                                        for line in figure.axes[0].lines if len(line.get_xdata()) == 2))
            finally:
                close("all")


if __name__ == "__main__":
    unittest.main()
