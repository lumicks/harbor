#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
General scientific constants
"""
import builtins


class Constant(float):

    def __init__(self, unit=None):
        super().__init__()
        self.unit = unit
        self.description = None


kB = Constant(1.3806488e-23)
kB.unit = 'J/K'
kB.__doc__ = """\
Boltzmann constant [J/K]
"""

h = Constant(6.62606957e-34)
h.unit = 'J * s'
h.__doc__ = """\
Planck's constant [J * s]
"""

NA = Constant(6.02214129e23)
NA.unit = "1/mol"
NA.__doc__ = """\
Avogadro's constant [1 / mol]
"""

vacuumPermittivity = Constant(8.854187817620e-12)
vacuumPermittivity.unit = "F/m"
vacuumPermittivity.__doc__ = """\
Permittivity in free space.
"""

c = Constant(299792458)
c.unit = "m/s"
c.__doc__ = """\
Speed of light
"""
