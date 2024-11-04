__doc__ = """\
Package for data analysis of optical trap experiments
"""

import os
import glob



SOURCE_FILES = glob.glob(os.path.dirname(__file__) + "/*.py")
__all__ = [os.path.basename(f)[: -3] for f in SOURCE_FILES]

__version__ = (0, 0, 1)

# tweezer related imports ----

# This has to be declared before any foldometer imports
_ROOT = os.path.abspath(os.path.dirname(__file__))

from .ixo.binary import read_file
from .ixo.lumicks_c_trap import lumicks_file
from .ixo.data_conversion import process_file, analyse_file

from . import physics
from .simulate.trap import simulate_trap
from .tools.plots import force_extension_curve, plot_xcorr
from .tools.maths import cross_correlation

from .analysis.thermal_calibration import calibration_file, calibration_data
from .analysis.region_classification import assign_regions
from .analysis.event_classification import find_unfolding_events
from .analysis.wlc_curve_fit import wlc_fit_data
from .analysis.tweezers_parameters import MIRRORVOLTDISTANCEFACTOR

from .core.main import Folding
from .core.subclasses import Thread, LifeTimeMeasurement
