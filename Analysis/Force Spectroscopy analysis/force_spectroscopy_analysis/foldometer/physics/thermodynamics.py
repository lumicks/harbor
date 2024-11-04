#!/usr/bin/env python
# -*- coding: utf-8 -*-

from foldometer.physics.constants import Constant
from foldometer.physics.constants import kB


def thermal_energy(temperature=298, units='pN nm'):
    """
    Thermal energy in units of [pN nm]

    Args:
        temperature (float): temperature in units of [K]. (optional, default 298)
        units (str): unit of the returned energy value, allowed values: 'pN nm', 'J', \
        or None. (optional, default 'pN nm')

    Returns:
        energy (float or Constant): thermal energy in [units] (default is [pN nm])

    Usage:

        >>> thermal_energy(273.15, 'J')
        >>> 4.1143334240000004e-21

        >>> thermal_energy()
        >>> 4.114333424000001

    .. note::

        The default temperature corresponds to 25°C.

    """
    if units is None:
        energy = kB * temperature * 10**21
    elif 'pN nm' in units:
        energy = Constant(kB * temperature * 10**21)
        energy.unit = 'pN nm'
    elif 'J' in units:
        energy = Constant(kB * temperature)
        energy.unit = 'J'
    else:
        raise BaseException("Can't figure out how to return thermal energy value")

    return energy
