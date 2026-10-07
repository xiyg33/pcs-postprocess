"""Deterministic frequency event clocks and independent raw-axis P metrics."""

import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from pcs_postprocess import freq_reg_timing as timing, processing as p
from pcs_postprocess.config import ProjectConfig


class FrequencyTimingTests(unittest.TestCase):
    def setUp(self):
        p.configure_project(ProjectConfig.load(Path(__file__).resolve().parents[1] / "examples/project.json"))

    def case(self, deviation=.2, pll_delay=0., start=2.6, hold=3.):
        t = np.arange(0., 8., .001)
        up = np.clip((t - start - pll_delay) / .08, 0., 1.)
        down = np.clip((t - start - hold - pll_delay) / .08, 0., 1.)
        f = 50. + deviation * (up - down) + .003 * np.sin(2 * np.pi * 60 * t)
        power = .35 - np.sign(deviation) * .2 * (1. - np.exp(-np.maximum(t - start - .04, 0.) / .05))
        power[t >= start + hold + .04] = .35
        meta = {"mode": "freq_reg", "test_type": "FREQ_REG", "test_type_code": 6,
                "time_basis": "recorded", "case_index": 3019, "case_name": "synthetic_timing",
                "event_start_s": 2., "event_hold_s": hold, "freq_deviation_hz": deviation,
                "P_ref_pu": .35, "hil_settle_shift_s": 6., "hil_rec_start_target_s": 5.5}
        return t, {"f": f, "P": power * 24000}, meta, p.event_times_from_meta(meta)

    def test_ripple_spike_small_steps_and_both_frequency_directions(self):
        for dev in (.05, -.05, .2, -.2):
            with self.subTest(deviation=dev):
                t, sig, meta, ev = self.case(dev)
                sig["f"][2300:2303] += 1.
                result = timing.resolve(t, sig, meta, ev, 5.5, 24000.)
                self.assertEqual(result["freq_timing_review_required"], 0, result["freq_timing_reason"])
                self.assertAlmostEqual(result["freq_event_start_recorder_s"], 2.6, delta=.015)
                self.assertAlmostEqual(result["freq_event_end_recorder_s"], 5.6, delta=.015)
                self.assertGreater(result["freq_power_start_recorder_s"], result["freq_event_start_recorder_s"])

    def test_no_frequency_event_and_missing_recovery_require_review(self):
        for problem in ("no_event", "no_recovery", "gap"):
            with self.subTest(problem=problem):
                t, sig, meta, ev = self.case()
                if problem == "no_event":
                    sig["f"] = np.full_like(t, 50.)
                elif problem == "no_recovery":
                    sig["f"][t > 5.6] = 50.2
                else:
                    sig["f"][(t > 2.5) & (t < 2.6)] = np.nan
                result = timing.resolve(t, sig, meta, ev, 5.5, 24000.)
                self.assertEqual(result["freq_timing_review_required"], 1)
                self.assertEqual(result["freq_event_start_recorder_s"], 2.5)
                self.assertTrue(result["freq_timing_reason"])

    def test_late_pll_does_not_replace_power_onset(self):
        t, sig, meta, ev = self.case(pll_delay=.1)
        result = timing.resolve(t, sig, meta, ev, 5.5, 24000.)
        self.assertEqual(result["freq_timing_review_required"], 1)
        self.assertIn("power_response_precedes_frequency_estimate", result["freq_timing_reason"])
        self.assertLess(result["freq_power_start_recorder_s"], result["freq_candidate_start_recorder_s"])
        self.assertEqual(result["freq_event_start_recorder_s"], 2.5)

    def test_manual_confirmation_and_synchronized_marker_precedence(self):
        t, sig, meta, ev = self.case()
        meta.update(freq_timing_override_source="manual_confirmed", freq_event1_start_recorder_s=2.59,
                    freq_event1_end_recorder_s=5.59)
        manual = timing.resolve(t, sig, meta, ev, 5.5, 24000.)
        self.assertEqual(manual["freq_timing_confidence"], "manual_confirmed")
        self.assertEqual(manual["freq_event_start_recorder_s"], 2.59)
        meta.update(hil_command_time_source="synchronized", hil_command_event1_recorder_s=2.6)
        synchronized = timing.resolve(t, sig, meta, ev, 5.5, 24000.)
        self.assertEqual(synchronized["freq_event_start_recorder_s"], 2.6)
        self.assertEqual(synchronized["freq_timing_confidence"], "synchronized")
        meta["hil_command_event1_recorder_s"] = -1.
        with self.assertRaisesRegex(ValueError, "Frequency timing override"):
            timing.resolve(t, sig, meta, ev, 5.5, 24000.)

    def test_single_and_batch_match_without_ordinary_median_fallback(self):
        t, sig, meta, ev = self.case()
        sig["f"][:] = 50.
        other_meta = dict(meta, case_index=3020)
        with patch.object(p, "load_case_file", return_value=(t, sig)):
            single = p.resolve_case_offsets([(3019, "a.mat", dict(meta))])[0]
            batch = p.resolve_case_offsets([(3019, "a.mat", dict(meta)), (3020, "b.mat", other_meta)])[0]
        self.assertEqual(single[3019], batch[3019])
        self.assertEqual(single[3019], (5.5, "frequency_nominal_unconfirmed"))

    def test_full_recording_power_start_can_precede_display_frequency_anchor(self):
        t, sig, meta, ev = self.case()
        result = timing.resolve(t, sig, meta, ev, 5.5, 24000.)
        view = {"t_s": t, "P_pu": sig["P"] / 24000., "Q_pu": np.zeros_like(t),
                "event1_start_s": 2.8, "event1_end_s": 5.6, "power_search_start_s": 2.5,
                "independent_power_start_s": result["freq_power_start_recorder_s"],
                "independent_power_end_s": result["freq_event_end_recorder_s"]}
        summary = p.freq_reg_step_summary(view, meta)
        self.assertLess(summary["step_response_start_s"], 2.8)
        self.assertIsNotNone(summary["T_up_s"])
        self.assertIsNotNone(summary["M_p_pct"])

    def test_power_response_before_nominal_event_uses_clean_early_baseline(self):
        t, sig, meta, ev = self.case(start=2.3)
        result = timing.resolve(t, sig, meta, ev, 5.5, 24000.)
        self.assertIsNotNone(result["freq_power_start_recorder_s"])
        self.assertAlmostEqual(result["freq_power_start_recorder_s"], 2.343, delta=.02)

    def test_absolute_sim_power_clock_never_uses_hil_recording_prior(self):
        t, sig, meta, _ = self.case(start=2.)
        meta.update(time_basis="absolute", data_source="sim")
        meta.pop("hil_settle_shift_s")
        meta.pop("hil_rec_start_target_s")
        ev = p.event_times_from_meta(meta)
        result = timing.resolve(t, sig, meta, ev, .75, 24000.)
        self.assertEqual(result["freq_event_start_recorder_s"], 2.)
        self.assertEqual(result["freq_timing_review_required"], 0)
        self.assertIsNotNone(result["freq_power_start_recorder_s"])
        self.assertAlmostEqual(result["freq_power_start_recorder_s"], 2.043, delta=.02)

    def test_deadband_and_inertia_keep_legacy_detector(self):
        t, sig, meta, ev = self.case()
        for mode, test_type in (("freq_reg", "FREQ_REG_DEADBAND"), ("inertia", "INERTIA")):
            with self.subTest(mode=mode):
                changed = dict(meta, mode=mode, test_type=test_type)
                with patch.object(timing, "resolve", side_effect=AssertionError("legacy only")), \
                     patch.object(p, "_legacy_frequency_edge_offset", return_value=5.4):
                    self.assertEqual(p.detect_frequency_edge_offset(t, sig, changed, ev, 5.5), 5.4)

    def test_boundary_and_window_sensitivity_require_review(self):
        for reason in ("boundary", "window"):
            with self.subTest(reason=reason):
                t, sig, meta, ev = self.case()
                def fit(grid, values, nominal, deviation, noise, padding=.30):
                    start = nominal + .1 + (.03 if reason == "window" and padding == .20 else 0.)
                    return {"start_s": start, "duration_s": .08, "amplitude_hz": deviation,
                            "baseline_hz": 50., "rmse_hz": .001, "outlier_fraction": 0.,
                            "boundary": reason == "boundary"}
                with patch.object(timing, "_fit", side_effect=fit):
                    result = timing.resolve(t, sig, meta, ev, 5.5, 24000.)
                self.assertEqual(result["freq_timing_review_required"], 1)
                self.assertIn("fit_at_search_boundary" if reason == "boundary" else "fit_sensitive_to_window",
                              result["freq_timing_reason"])


if __name__ == "__main__":
    unittest.main()
