# coding=utf-8
__doc__ = """\
General physical and chemical functions and constants
"""

import os
import glob

SOURCE_FILES = glob.glob(os.path.dirname(__file__) + "/*.py")
__all__ = [os.path.basename(f)[: -3] for f in SOURCE_FILES]

# imports
from .thermodynamics import thermal_energy

from .hydrodynamics import drag_sphere

from .viscosity import (dynamic_viscosity_of_mixture,
                        water_dynamic_viscosity,
                        water_density,
                        glycerol_dynamic_viscosity,
                        glycerol_density)

from .constants import (kB, h, NA, c, vacuumPermittivity)

from . import polymer #, electrostatics, materials

from .utils import (mass_sphere, volume_sphere, as_Celsius, as_Kelvin)
