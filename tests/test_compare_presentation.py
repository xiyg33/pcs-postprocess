"""Selection must isolate cases; delivery styling must preserve numerical evidence."""

import csv
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from pcs_postprocess import compare_pwr as compare
from pcs_postprocess import compare_scr
from pcs_postprocess.compare_presentation import apply_presentation, figure_path, select_case_ids
import test_compare_pwr as fixtures


class PresentationTests(unittest.TestCase):
    def test_selection_rejects_conflicting_empty_duplicate_or_noninteger_ids(self):
        self.assertEqual(select_case_ids(None, [2, 1], {9}), [1, 2])
        self.assertEqual(select_case_ids(3, None, {9}), [3])
        for ids in ([], [1, 1], [True], [-1], ["1"]):
            with self.assertRaises(ValueError):
                select_case_ids(None, ids, {9})
        with self.assertRaises(ValueError):
            select_case_ids(1, [1], {9})
        with self.assertRaises(ValueError):
            figure_path(Path("figures"), 1, "overlay", {1: {"display_case_name": "../../escape"}})

    def test_axes_include_outliers_ignore_vertical_markers_and_keep_samples(self):
        fig, axes = compare.plt.subplots(2, 1)
        try:
            p_values = np.array([-1., -.9, -1.6])
            q_values = np.array([0., .01, -.01])
            lines = [axes[0].plot([0, 1, 2], p_values)[0], axes[1].plot([0, 1, 2], q_values)[0]]
            axes[0].set_ylabel("有功功率 / pu")
            axes[1].set_ylabel("无功功率 / pu")
            for ax in axes:
                ax.axvline(1.)
            fig.suptitle("original\nstatus retained")
            apply_presentation(fig, 1, {1: dict(title="delivery -0.9Pn", power_axis_limits=[-1.2, .2],
                                               power_tick_pu=.2, broad_q_axis=True)})
            self.assertLess(axes[0].get_ylim()[0], -1.6)
            self.assertEqual(axes[0].get_ylim()[1], .2)
            self.assertLess(axes[1].get_ylim()[1], 1.)
            self.assertEqual(fig._suptitle.get_text(), "delivery -0.9Pn\nstatus retained")
            np.testing.assert_array_equal(lines[0].get_ydata(), p_values)
            np.testing.assert_array_equal(lines[1].get_ydata(), q_values)
        finally:
            compare.plt.close(fig)

    def test_styled_selected_comparison_preserves_metrics_and_input_bytes(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            config, _, files = fixtures.ComparePwrTests().make_project(base)
            paths = list(files.values()) + [p.with_name(p.stem + "_meta.json") for p in files.values()]
            paths += [base / side / f"results_pwr/postprocess/{side}_summary_pwr.csv" for side in files]
            before = [path.read_bytes() for path in paths]
            (base / "sim/results_pwr/canonical/case_9999_invalid.mat").write_bytes(b"unselected invalid data")
            presentation = {1001: dict(display_case_name="reported_P030", reported_P_ref_pu=.3,
                                       title="reported +0.3Pn", power_axis_limits=[-.2, 1.2])}
            with patch.object(compare.plt.Figure, "savefig"):
                plain = compare.run(config, base / "plain", case_ids=[1001])
                styled = compare.run(config, base / "styled", case_ids=[1001], presentation=presentation)
            self.assertEqual(plain["exceptions"], 0)
            self.assertEqual(styled["exceptions"], 0)
            self.assertEqual(styled["selected_case_ids"], [1001])
            self.assertEqual([path.read_bytes() for path in paths], before)
            def rows(folder):
                with (folder / "compare_summary_pwr.csv").open(encoding="utf-8-sig", newline="") as handle:
                    return list(csv.DictReader(handle))
            for first, second in zip(rows(base / "plain"), rows(base / "styled")):
                for key in compare.FIELDS:
                    if key != "figure_png":
                        self.assertEqual(first[key], second[key], key)
                self.assertEqual(second["display_case_name"], "reported_P030")
                self.assertIn("reported_P030", second["figure_png"])

    def test_reported_power_does_not_bypass_actual_pair_validation(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            config, _, files = fixtures.ComparePwrTests().make_project(base)
            meta_path = files["sim"].with_name(files["sim"].stem + "_meta.json")
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
            meta["P_base_pu"] += .05
            meta_path.write_text(json.dumps(meta), encoding="utf-8")
            result = compare.run(config, base / "compare", case_ids=[1001],
                                 presentation={1001: dict(display_case_name="P030", reported_P_ref_pu=.3)})
            self.assertEqual(result["completed"], 0)
            self.assertEqual(result["exceptions"], 1)
            self.assertIn("P_base_pu", result["errors"][0]["reason"])

    def test_scr_boundary_preserves_actual_power_and_records_reported_power(self):
        rows = [dict(case_index=1, P_ref_pu=-1., reported_P_ref_pu=-.9,
                     SCR_value=2.8, sim_state="STABLE", hil_state="STABLE"),
                dict(case_index=2, P_ref_pu=-1., reported_P_ref_pu=-.9,
                     SCR_value=2.7, sim_state="UNSTABLE", hil_state="UNSTABLE")]
        row = compare_scr.boundary_rows(rows)[0]
        self.assertEqual(row["P_ref_pu"], -1.)
        self.assertEqual(row["source_P_ref_pu"], -1.)
        self.assertEqual(row["reported_P_ref_pu"], -.9)
        self.assertTrue(row["approximation"])
        self.assertEqual(row["boundary_match"], "same_tested_states")


if __name__ == "__main__":
    unittest.main()
