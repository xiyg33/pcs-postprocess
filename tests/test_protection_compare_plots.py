"""Protection overlays retain the full recording regardless of trip state."""

import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from matplotlib.figure import Figure

from pcs_postprocess import compare_freq_prot, compare_volt_prot, processing
from pcs_postprocess.config import ProjectConfig


class ProtectionComparePlotTests(unittest.TestCase):
    def setUp(self):
        processing.configure_project(ProjectConfig.load(Path(__file__).resolve().parents[1] / "examples/project.json"))

    def pair(self, trips, end=63.):
        t = np.arange(0., end, .005)
        pair = {}
        for side, trip in zip(("hil", "sim"), trips):
            pair[side] = {
                "t_s": t, "trip_s": trip,
                # This deliberately differs from t2 - 1: annotations must use
                # the saved metric, rather than measuring it in the renderer.
                "delta_t_s": .1234 if trip is not None else None,
                "breaker": np.ones(t.shape),
                "timing": {"protection_status": "tripped" if trip is not None else "no_trip"},
                "display": {"t_s": t, "Vp_pu": np.ones(t.shape),
                            "P_pu": np.zeros(t.shape), "Q_pu": np.zeros(t.shape),
                            "voltage_frequency_hz": np.full(t.shape, 52.), "f_hz": np.full(t.shape, 52.)},
            }
        return pair

    def test_full_window_three_panels_and_saved_timing_annotations(self):
        metas = {side: {"case_name": "synthetic", "f_target_hz": 52.,
                        "event_hold_s": 60., "record_post_s": 2.} for side in ("hil", "sim")}
        for module in (compare_freq_prot, compare_volt_prot):
            for trips in ((1.2, 1.3), (1.2, None), (None, None)):
                with self.subTest(mode=module.__name__, trips=trips):
                    def inspect(fig, *args, **kwargs):
                        self.assertEqual(len(fig.axes), 3)
                        for ax in fig.axes:
                            self.assertEqual(ax.get_xlim(), (0., 63.))
                            self.assertNotIn("breaker", ax.get_ylabel())
                            self.assertTrue(all(line.get_drawstyle() == "default" for line in ax.lines))
                            verticals = [float(line.get_xdata()[0]) for line in ax.lines
                                         if len(line.get_xdata()) == 2]
                            self.assertIn(1., verticals)
                            for trip in trips:
                                if trip is not None:
                                    self.assertIn(trip, verticals)
                        texts = fig.axes[0].texts
                        self.assertEqual(len(texts), 2)
                        self.assertNotEqual(texts[0].get_position(), texts[1].get_position())
                        for text, trip in zip(texts, trips):
                            self.assertIs(text.get_transform(), fig.axes[0].transAxes)
                            self.assertIn("Δt=0.1234 s" if trip is not None else "no_trip", text.get_text())
                    with patch.object(Figure, "savefig", autospec=True, side_effect=inspect):
                        module.plot_pair(self.pair(trips), metas, [Path("unused.png")], 6001)

    def test_short_recording_is_not_extrapolated(self):
        pair = self.pair((1.2, 1.3), end=4.)
        metas = {side: {"event_hold_s": 60., "record_post_s": 2.} for side in pair}
        self.assertAlmostEqual(compare_volt_prot.full_recording_stop(pair, metas), 3.995)
        for data in pair.values():
            self.assertLess(data["display"]["t_s"][-1], 4.)


if __name__ == "__main__":
    unittest.main()
