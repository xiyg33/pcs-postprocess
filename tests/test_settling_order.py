"""Strict rise/settling order on measured samples and referenced summaries."""

import copy
import unittest
from unittest.mock import patch

import numpy as np

from pcs_postprocess import processing as p
from pcs_postprocess.compare_pwr import comparison_rows
from pcs_postprocess.step_response import measure_step, peak_order_reason, settling_order_reason


class SettlingOrderTests(unittest.TestCase):
    def test_continuous_response_in_both_directions(self):
        t = np.arange(0., 4.001, .005)
        for direction in (-1., 1.):
            y = .35 + direction * np.clip((t - 1.) / .5, 0., 1.)
            measured = measure_step(t, y, 1., 4., direction, smooth_ms=0.)
            self.assertLess(measured["values"]["T_up"], measured["values"]["T_p"])
            self.assertLess(measured["values"]["T_p"], measured["values"]["T_s"])

    def test_instant_response_does_not_invent_a_peak_after_equal_rise(self):
        t = np.arange(0., 3.001, .005)
        measured = measure_step(t, (t >= 1.).astype(float), 1., 3., 1., smooth_ms=0.)
        self.assertEqual(measured["values"]["T_up"], 0.)
        self.assertIsNone(measured["values"]["T_p"])
        self.assertIsNone(measured["values"]["T_s"])
        self.assertEqual(measured["reasons"]["T_p"], "peak_not_after_rise")

    def test_late_in_band_peak_still_requires_observation_after_peak(self):
        t = np.arange(0., 4.001, .005)
        y = np.clip((t - 1.) / .2, 0., 1.)
        y[np.argmin(abs(t - 3.9))] += .001
        measured = measure_step(t, y, 1., 4., 1., smooth_ms=0.)
        self.assertAlmostEqual(measured["values"]["T_p"], 2.9)
        self.assertIsNone(measured["values"]["T_s"])
        self.assertEqual(measured["reasons"]["T_s"], "not_settled_within_observation")

    def test_unconfirmed_rise_or_short_observation_has_no_settling(self):
        t = np.arange(0., 3.001, .005)
        for end, unstable in ((1.2, False), (3., True)):
            y = (t >= 1.).astype(float)
            if unstable:
                y[t < 1.] = .4 * t[t < 1.]
            measured = measure_step(t, y, 1., end, 1., smooth_ms=0.)
            self.assertIsNone(measured["values"]["T_up"])
            self.assertIsNone(measured["values"]["T_s"])

    def test_fourth_entry_must_follow_available_rise(self):
        t = np.arange(0., 2.001, .005)
        y = np.zeros_like(t)
        for start in (1.01, 1.03, 1.05, 1.07):
            y[(t >= start) & (t < start + .006)] = 1.
        measured = dict(initial_pu=0., final_pu=1., values=dict(T_up=.01, T_p=.02))
        candidate = p.pwr_fourth_band_entry(t, y, 1., 2., measured, smooth_ms=0.)
        self.assertIsNotNone(candidate)
        for rise in (candidate, candidate + .1, None):
            measured["values"]["T_up"] = rise
            self.assertIsNone(p.pwr_fourth_band_entry(t, y, 1., 2., measured, smooth_ms=0.))
        measured["values"]["T_up"] = .01
        for peak in (candidate, candidate + .1, None):
            measured["values"]["T_p"] = peak
            self.assertIsNone(p.pwr_fourth_band_entry(t, y, 1., 2., measured, smooth_ms=0.))

    def test_rejected_fourth_entry_falls_back_to_continuous_result(self):
        t = np.arange(0., 3.001, .005)
        data = dict(t_s=t, P_pu=(t >= 1.).astype(float), Q_pu=np.zeros_like(t),
                    event1_start_s=1., command_time_reliable=1, command_time_source="synchronized",
                    command_event1_normalized_s=1.)
        measured = dict(initial_pu=0., final_pu=1., baseline_stable=True, tail_stable=True,
                        values=dict(M_p=0., T_up=.2, T_p=.3, T_s=.4), reasons={})
        with patch("pcs_postprocess.step_response.measure_step", return_value=copy.deepcopy(measured)):
            result = p.pwr_step_summary(data, dict(test_type_code=3, P_base_pu=0., P_step_pu=1.,
                                                   step_start_s=1., step_dur_s=2.))
        self.assertEqual(result["T_s_s"], .4)
        self.assertEqual(result["T_s_method"], "continuous_5pct")

    def test_compare_rejects_old_equal_earlier_or_missing_rise_on_one_side(self):
        source = dict(M_p_pct=0., T_up_s=.1, T_p_s=.2, T_s_s=.3,
                      M_p_reason="", T_up_reason="", T_p_reason="", T_s_reason="original reason",
                      step_timing_basis="synchronized_command", step_response_start_s=1.,
                      T_s_method="continuous_5pct", metric_sampling_ms=5, 人工确认Ts="")
        for rise, peak, settling in ((.1, .3, .3), (.1, .4, .3), ("", .2, .3),
                                     (.1, "", .3), (.2, .2, .3), (.3, .2, .4), (.1, .2, "")):
            with self.subTest(rise=rise, settling=settling):
                summaries = dict(sim=dict(source, T_up_s=rise, T_p_s=peak, T_s_s=settling), hil=source)
                original = copy.deepcopy(summaries)
                rows = comparison_rows(2009, dict(case_name="case", test_type="PWR", test_type_code=3), summaries, "plot.png")
                row = next(r for r in rows if r["metric"] == "T_s")
                self.assertIsNone(row["sim_value"])
                self.assertEqual(row["hil_value"], .3)
                self.assertIsNone(row["deviation"])
                self.assertIsNone(row["absolute_deviation"])
                self.assertEqual(row["status"], "unavailable")
                self.assertEqual(summaries, original)
                if settling != "":
                    self.assertIn("rerun single-source", row["sim_reason"])
                    self.assertIn("original reason", row["sim_reason"])
                    self.assertEqual(row["sim_T_s_method"], "unavailable")
                if rise != "" and peak != "" and peak <= rise:
                    self.assertIsNone(next(r for r in rows if r["metric"] == "T_p")["sim_value"])

    def test_strict_comparison_does_not_add_a_sample_or_round(self):
        self.assertEqual(peak_order_reason(.1, .1), "peak_not_after_rise")
        self.assertEqual(settling_order_reason(.1, .2, .2), "settling_not_after_peak")
        self.assertEqual(settling_order_reason(.1, .2, np.nextafter(.2, 1.)), "")


if __name__ == "__main__":
    unittest.main()
