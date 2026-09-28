"""标量曲线的绘图降采样不改变指标输入。"""

import unittest

import numpy as np

from pcs_postprocess.plot_common import scalar_plot_data
from pcs_postprocess.processing import pwr_downsampled_signals


class ScalarPlotTests(unittest.TestCase):
    def test_uses_pwr_filter_and_five_ms_grid_without_mutating_input(self):
        t = np.arange(0.0, 0.101, 0.001)
        p = 0.5 + 0.1 * (-1.0) ** np.arange(t.size)
        original = p.copy()
        bridged = {"t_s": t, "P_pu": p, "Vp_pu": np.ones_like(t),
                   "Ip_pu": np.linspace(0.0, 1.0, t.size)}
        events = {"event1_start": 0.04, "event1_end": 0.08}
        display = scalar_plot_data(bridged, events)
        pwr = pwr_downsampled_signals(bridged, 0.04, 0.08)
        np.testing.assert_allclose(display["t_s"], np.arange(0, 0.101, 0.005))
        np.testing.assert_allclose(display["P_pu"], pwr["P_pu"])
        self.assertLess(np.ptp(display["P_pu"]), np.ptp(p))
        np.testing.assert_array_equal(p, original)

    def test_missing_values_and_time_gaps_are_not_bridged(self):
        t = np.r_[np.arange(0.0, 0.041, 0.001), np.arange(0.07, 0.101, 0.001)]
        p = np.ones_like(t)
        p[10:15] = np.nan
        display = scalar_plot_data({"t_s": t, "P_pu": p, "Vp_pu": np.ones_like(t)},
                                   {"event1_start": 0.04, "event1_end": 0.08})
        self.assertTrue(np.isnan(display["P_pu"][2]))  # t=0.010，原始数据缺失
        self.assertTrue(np.isnan(display["P_pu"][10]))  # t=0.050，原始时间缺口
        self.assertTrue(np.isnan(display["Q_pu"]).all())


if __name__ == "__main__":
    unittest.main()
