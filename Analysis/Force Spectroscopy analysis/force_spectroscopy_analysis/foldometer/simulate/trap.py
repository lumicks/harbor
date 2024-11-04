# coding=utf-8
from collections import namedtuple

import numpy as np
import pandas as pd

from foldometer.physics import drag_sphere, mass_sphere
from foldometer.physics.hydrodynamics import diffusion_coefficient

__doc__ = """\
Simulation of optical trap time series according to Norrelykke et al.
"""


def eigenvalues(dragCoefficient=drag_sphere(1000),
                massSphere=mass_sphere(1000),
                trapStiffness=0.1):
    """
    Calculates eigenvalues in the simulation of OTs according to Norrelykke et al.

    Args:
        dragCoefficient (float): Stokes drag coefficient in [pN/nm s]
        massSphere (float): Mass in [g]
        trapStiffness (float): Stiffness in [pN/nm] (Default: 0.1 pN/nm)

    Returns:
        eigenvalues (namedtuple): Eigenvalues of Langevin equation; fields are 'plus' and 'minus'
    """
    g = dragCoefficient
    m = massSphere
    k = trapStiffness

    discriminant = (g**2 / (4 * m**2)) - k / m

    prefactor = g / (2 * m)

    plus = prefactor + np.sqrt(discriminant)
    minus = prefactor - np.sqrt(discriminant)

    Eigenvalues = namedtuple("Eigenvalues", ['plus', 'minus'])

    eigenvalues = Eigenvalues(plus, minus)

    return eigenvalues


def cValues(eigenvalues=eigenvalues(), timeStep = 0.001):
    """
    Calculates c values in the simulation of OTs according to Norrelykke et al.

    Computes: c = exp(-λ * ∆t)

    Args:
        eigenvalues (namedtuple): Eigenvalues of Langevin equation; (Fields: 'plus', 'minus')
        timeStep (float): Difference between two time points in [s]; delta time (Default: 0.001 s)

    Returns:
        cValues (namedtuple): C values in the simulation of an optical trap; (Fields: 'plus', 'minus')
    """
    cValues = namedtuple("cValues", ['plus', 'minus'])

    plus = np.exp(-eigenvalues.plus * timeStep)
    minus = np.exp(-eigenvalues.minus * timeStep)

    c = cValues(plus, minus)

    return c


def aValues(diffusionCoefficient=diffusion_coefficient(radius=1000, temperature=25, dynamicViscosity=1e-9),
            eigenvalues=eigenvalues(), cValues=cValues()):
    """
    Calculates A values in the simulation of OTs according to Norrelykke et al.

    Args:
        diffusionCoefficient (float): Diffusion coefficient in [nm^2 / s]
        eigenvalues (namedtuple): Eigenvalues of Langevin equation; (Fields: 'plus', 'minus')
        cValues (namedtuple): C values in the simulation of an optical trap; (Fields: 'plus', 'minus')

    Returns:
        aValues (namedtuple): A values in the simulation of an optical trap; (Fields: 'plus', 'minus')
    """

    l = eigenvalues
    c = cValues
    D = diffusionCoefficient

    factorA = (l.plus + l.minus) / (l.plus - l.minus)
    factorB = np.sqrt((1 - c.plus**2) * D / (2 * l.plus))
    factorC = np.sqrt((1 - c.minus**2) * D / (2 * l.minus))

    aValues = namedtuple("aValues", ['plus', 'minus'])

    plus = factorA * factorB
    minus = factorA * factorC

    A = aValues(plus, minus)

    return A


def alpha(eigenvalues=eigenvalues(), cValues=cValues()):
    """
    Calculates alpha in the simulation of OTs according to Norrelykke et al.

    Args:
        eigenvalues (namedtuple): Eigenvalues of Langevin equation; (Fields: 'plus', 'minus')
        cValues (namedtuple): C values in the simulation of an optical trap; (Fields: 'plus', 'minus')

    Returns:
        alpha (float): alpha in the simulation of an optical trap;
    """
    l = eigenvalues
    c = cValues

    factorA = 2 * np.sqrt(l.plus * l.minus) / (l.plus + l.minus)
    factorB = (1 - c.plus * c.minus) / np.sqrt((1 - c.plus**2) * (1 - c.minus**2))

    alphaValue = factorA * factorB

    return alphaValue


def exp_Matrix(eigenvalues=eigenvalues(), cValues=cValues()):
    """
    Calculates alpha in the simulation of OTs according to Norrelykke et al.

    Computes: exp(-M * ∆t)

    Args:
        eigenvalues (namedtuple): Eigenvalues of Langevin equation; (Fields: 'plus', 'minus')
        cValues (namedtuple): C values in the simulation of an optical trap; (Fields: 'plus', 'minus')

    Returns:
        M (np.matrix(2, 2)): exp(-M * ∆t) in the simulation of an optical trap;
    """
    l = eigenvalues
    c = cValues

    prefactor = 1 / (l.plus - l.minus)

    M = np.matrix(np.zeros((2, 2)))

    M[0, 0] = -l.minus * c.plus + l.plus * c.minus
    M[0, 1] = -c.plus + c.minus
    M[1, 0] = l.plus * l.minus * (c.plus - c.minus)
    M[1, 1] = l.plus * c.plus - l.minus * c.minus

    return prefactor * M


def step(eigenvalues=eigenvalues(), aValues=aValues(), alpha=alpha()):
    """
    Calculates single step in the simulation of OTs according to Norrelykke et al.

    Args:
        eigenvalues (namedtuple): Eigenvalues of Langevin equation; (Fields: 'plus', 'minus')
        aValues (namedtuple): A values in the simulation of an optical trap; (Fields: 'plus', 'minus')
        alpha (float): alpha in the simulation of an optical trap;

    Returns:
        step (np.array): Simulation step; step[0] - deltaX in [nm], step[1] - deltaV in [nm/s]
    """
    l = eigenvalues
    A = aValues

    v1 = np.array([-1, l.plus])
    v2 = np.array([1, -l.minus])

    partA = (A.plus * v1 + A.minus * v2) * np.sqrt(1 + alpha) * np.random.randn()
    partB = (A.plus * v1 - A.minus * v2) * np.sqrt(1 - alpha) * np.random.randn()

    # keep in mind that:
    # deltaX = partA[0] + partB[0]
    # deltaV = partA[1] + partB[1]
    step = partA + partB

    return step


def simulate_trap(dataPoints=1e3,
                  timeStep=0.001,
                  radius=1000,
                  viscosity=1e-9,
                  trapStiffness=0.1,
                  temperature=25):
    """
    Simulates the position time series of a sphere in an optical trap

    Args:
        dataPoints (int): Number of data points in final time series. (Default: 1e3)
        timeStep (float): Difference between two time points in [s]; delta time (Default: 0.001 s)
        radius (float): Radius of trapped sphere in [nm]. (Default: 1000)
        viscosity (float): Dynamic viscosity in [pN/nm^2 s] (Default: 1e-9 pN/nm^2 s)
        trapStiffness (float): Stiffness of simulated trap, k, in [pN/nm]. (Default: 0.1 pN/nm)
        temperature (float): Temperature in [°C]. (Default: 25 °C)

    Returns:
        state (pandas.DataFrame): Simulated state of an optical trap. (Columns: 't', 'x', 'v')

    """
    #TODO: add time index to pd.DataFrame; makes for easier plotting as time is inherent
    assert timeStep > 0
    assert radius > 0
    assert dataPoints > 0

    # boundary conditions
    drag = drag_sphere(radius=radius, dynamicViscosity=viscosity)
    mass = mass_sphere(radius=radius)
    D = diffusion_coefficient(radius=radius, temperature=temperature, dynamicViscosity=viscosity)

    l = eigenvalues(dragCoefficient=drag, massSphere=mass, trapStiffness=trapStiffness)

    c = cValues(eigenvalues=l, timeStep=timeStep)

    expM = exp_Matrix(eigenvalues=l, cValues=c).T

    A = aValues(diffusionCoefficient=D, eigenvalues=l, cValues=c)

    alphaValue = alpha(eigenvalues=l, cValues=c)

    state = np.zeros((dataPoints, 3))
    state[:, 0] = np.arange(0, timeStep * dataPoints, timeStep)

    for i in range(int(dataPoints) - 1):

        singleStep = step(eigenvalues=l, aValues=A, alpha=alphaValue)

        state[i + 1, 1:3] = state[i, 1:3] * expM + singleStep

    state = pd.DataFrame(state, columns=['t', 'x', 'v'])
    return state
