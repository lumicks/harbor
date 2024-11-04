#!/usr/bin/env python
#-*- coding: utf-8 -*-

import numpy as np

from foldometer.physics.thermodynamics import thermal_energy
from foldometer.physics.utils import as_Kelvin


def drag_sphere(radius=1000, dynamicViscosity=0.9e-9, verbose=False):
    """
    Calculates the simple Stokes' drag coefficient of a sphere in a Newtonian fluid at low Reynolds number.

    Args:
        radius (float): Radius of solid sphere in [nm]. Default: 1000 nm
        dynamicViscosity (float): Dynamic viscosity in [pN/nm^2 s]. Default: 0.9e-9 pN/nm^2 s
        verbose (bool): Print parameters and results with units. Default: False

    Returns:
        dragCoefficient (float): Stokes drag coefficient in [pN/nm s]

    """
    if not isinstance(radius, (int, float, np.float)):
        try:
            radius = np.float(radius)
        except ValueError:
            print('Radius must be a number, not a {}'.format(type(radius)))

    if not isinstance(dynamicViscosity, (int, float, np.float)):
        try:
            dynamicViscosity = np.float(dynamicViscosity)
        except ValueError:
            print('Viscosity must be a number, not a {}'.format(type(dynamicViscosity)))

    assert (radius > 0), 'Radius of sphere must be positive'

    dragCoefficient = 6 * np.pi * radius * dynamicViscosity

    if verbose:
        print("In:")
        print("Radius: r = {} nm".format(radius))
        print("Viscosity: eta = {} pN/nm^2 s\n".format(dynamicViscosity))

        print("Out:")
        print("Drag coefficient: gamma = {} pN/nm s".format(round(dragCoefficient, 12)))

    return dragCoefficient


def diffusion_coefficient(radius=1000, temperature=25, dynamicViscosity=1e-9, verbose=False):
    """
    Calculates the diffusion coefficient for a sphere based on Stokes drag and the Stokes-Einstein relation:
    D = kT / gamma

    Args:
        radius (float): Radius of sphere in [nm]. Default: 1000 nm
        temperature (float): Solvent temperature in °C. Default: 25
        dynamicViscosity (float): Dynamic viscosity in [pN/nm^2 s]. Default: 0.9e-9 pN/nm^2 s
        verbose (bool): Print parameters and results with units. Default: False

    Returns:
        diffusionConstant (float): Diffusion constant in [nm^2 / s]
    """
    assert radius > 0
    assert temperature >= -273.15
    assert dynamicViscosity > 0

    kT = thermal_energy(as_Kelvin(temperature))
    drag = drag_sphere(radius=radius, dynamicViscosity=dynamicViscosity)

    diffusionConstant = kT / drag

    if verbose:
        print("In:")
        print("Radius: r = {} nm".format(radius))
        print("Temperature: T = {} °C".format(temperature))
        print("Dynamic viscosity: eta = {} pN/nm^2 s\n".format(dynamicViscosity))

        print("Out:")
        print("Diffusion constant: D = {} nm^2 / s".format(round(diffusionConstant, 12)))

    return diffusionConstant



class StokesDragSphere(object):
    def __init__(self):
        raise NotImplementedError('Nope')

class TwoSphereHydrodynamicInteractions(object):
    def __init__(self):
        raise NotImplementedError('Nope')


class TemperatureDependentViscosity(object):
    def __init__(self):
        raise NotImplementedError('Nope')

