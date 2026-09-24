import argparse
import os
import sys

from .processing import supports_power_current
from .plot_common import draw_processed_case


def render_case(case, fig_dir):
    return draw_processed_case(case, fig_dir,
                               include_power_current=supports_power_current(case.meta, "frt"))
