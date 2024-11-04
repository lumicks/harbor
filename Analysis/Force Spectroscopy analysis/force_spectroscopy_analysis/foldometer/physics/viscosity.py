#!/usr/bin/env python
#-*- coding: utf-8 -*-

__doc__ = """Calculate viscosity of water-glycerol mixtures according to \
[Cheng](http://www3.ntu.edu.sg/home/cnscheng/Publications/reprint/Glycerin-\
    water%20viscosity%20paper%20_sent%20to%20Ind%20Eng%20Chem%20Res.pdf)

The following variables are used:

Temperature T: [0, 100] in [C]

"""
import sys
import numpy as np
import pandas as pd

from math import exp
from itertools import chain


def dynamic_viscosity_of_mixture(waterVolume=1, glycerolVolume=0, temperature=25):
    """
    Power law equation for the dynamic viscosity of a water \
    to glycerol mixture according to:

    \mu = \mu_{w}^{\alpha} \mu_{g}^{1 - \alpha}

    """
    T = float(temperature)
    wV = float(waterVolume)
    gV = float(glycerolVolume)
    mu_w = water_dynamic_viscosity(temperature)
    mu_g = glycerol_dynamic_viscosity(temperature)
    Cm = calc_glycerol_fraction_by_mass(wV, gV, T)
    alpha = calc_alpha(Cm, temperature)

    mu = mu_w ** alpha * mu_g ** (1 - alpha)

    return mu


def calc_coef_a(temperature=25):
    """
    Coefficient *a* for temperature dependence of alpha

    Args:
        temperature (float): in Celsius in the range [0, 100]

    Returns:
        coefA (float): coefficient A
    """
    coefA = 0.705 - 0.0017 * temperature

    return coefA


def calc_coef_b(coefA, temperature=25):
    """
    Coefficient *b* for temperature dependence of alpha

    Args:
        temperature (float): in Celsius in the range [0, 100]
        coefA (float): Coefficient A

    Returns:
        coefB (float): coefficient B
    """
    coefB = (4.9 + 0.036 * temperature) * coefA ** (2.5)

    return coefB


def calc_alpha(glycerolMassFraction, temperature=25):
    """
    Calculates *alpha*, the exponent used in the power law of\
    the dynamic viscosity calculation

    Args:
        glycerolMassFraction (float): Cm in the range [0, 1]
        temperature (float): in Celsius in the range [0, 100]

    Returns:
        alpha (float): power law exponent
    """
    a = calc_coef_a(temperature)
    b = calc_coef_b(a, temperature)
    Cm = float(glycerolMassFraction)

    alpha = 1 - Cm + (a * b * Cm * (1 - Cm)) / (a * Cm + b * (1 - Cm))

    return alpha


def water_dynamic_viscosity(temperature=25):
    """
    Calculates *mu_w*, the dynamic viscosity of water using the\
    interpolation formula of Cheng

    Args:
        temperature (flaot): in Celsius in the range [0, 100]

    Returns:
        waterDynamicViscosity (float): Dynamic viscosity of water in [0.001 N s / m^2]
    """
    T = float(temperature)
    mu = 1.790 * exp(((-1230 - T) * T) / (36100 + 360 * T))
    waterDynamicViscosity = 0.001 * mu

    return waterDynamicViscosity


def glycerol_dynamic_viscosity(temperature=25):
    """
    Calculates *mu_g*, the dynamic viscosity of glycerol using the\
    interpolation formula of Cheng

    Args:
        temperature (float): in Celsius in the range [0, 100]

    Returns:
        glycerolDynamicViscosity (float): Dynamic viscosity of glycerol in [0.001 N s / m^2]
    """
    T = float(temperature)
    mu = 12100 * exp(((-1233 + T) * T) / (9900 + 70 * T))
    glycerolDynamicViscosity = 0.001 * mu

    return glycerolDynamicViscosity


def water_density(temperature):
    """
    Calculates the density of water from an interpolation by\
    Cheng (see viscosity docstring for reference.)

    Args:
        temperature (float): in Celsius in the range [0, 100]

    Returns:
        waterDensity (float): Density of water in kg/m^3
    """
    rho = 1000 * (1 - abs((temperature - 4) / (622.0)) ** (1.7))
    waterDensity = rho

    return waterDensity


def glycerol_density(temperature):
    """
    Calculates the density of glycerol from an interpolation by\
    Cheng (see viscosity docstring for reference.)

    Args:
        temperature (float): in Celsius in the range [0, 100]

    Returns:
        glycerolDensity (float): Density of Glycerol in kg/m^3
    """
    rho = 1277 - 0.654 * temperature
    glycerolDensity = rho

    return glycerolDensity


def calc_density_of_mixture(waterVolume, glycerolVolume, temperature=25):
    """
    Calculates the density of a glycerol-water mixture\
    from an interpolation by Cheng (see viscosity docstring for reference.)

    Args:
        waterVolume (float): volume of water in l
        glycerolVolume (float): volume of glycerol in l
        temperature (float): in Celsius in the range [0, 100]

    Returns
        rho (float): Description
    """
    T = float(temperature)
    rW = water_density(temperature)
    rG = glycerol_density(temperature)
    Cm = calc_glycerol_fraction_by_mass(waterVolume, glycerolVolume, T)

    rho = rG * Cm + rW * (1 - Cm)

    return rho


def calc_glycerol_fraction_by_volume(waterVolume, glycerolVolume):
    """
    Calculates the volume fraction of glycerol in a water - glycerol mixture

    args:
        waterVolume (float): volume of water in l
        glycerolVolume (float): volume of glycerol in l

    Returns:
        volumeFractionGlycerol (float): Fraction of glycerol by volume in [0, 1]
    """
    gV = float(glycerolVolume)
    gW = float(waterVolume)

    try:
        Cv = gV / (gW + gV)
    except ZeroDivisionError:
        Cv = 0.0

    volumeFractionGlycerol = Cv

    return volumeFractionGlycerol


def calc_water_fraction_by_volume(waterVolume, glycerolVolume):
    """
    Calculates the volume fraction of water in a water - glycerol mixture

    Args:
        waterVolume (float): volume of water in l
        glycerolVolume (float): volume of glycerol in l

    Returns:
        volumeFractionWater (float): Fraction of water by volume in [0, 1]
    """
    gV = float(glycerolVolume)
    wV = float(waterVolume)

    try:
        Cv = wV / (wV + gV)
    except ZeroDivisionError:
        Cv = 0.0

    volumeFractionWater = Cv

    return volumeFractionWater


def calc_glycerol_fraction_by_mass(waterVolume, glycerolVolume, temperature):
    """
    Calculates the mass fraction of glycerol in a water - glycerol mixture

    Args:
        waterVolume (float): volume of water in l
        glycerolVolume (float): volume of glycerol in l
        temperature (float): in Celsius in the range [0, 100]

    Returns:
        massFractionGlycerol (float): Fraction of glycerol by mass in [0, 1]
    """
    T = float(temperature)
    wM = calc_mass(waterVolume, water_density(T))
    gM = calc_mass(glycerolVolume, glycerol_density(T))

    try:
        Cm = gM / (wM + gM)
    except ZeroDivisionError:
        Cm = 0.0

    massFractionGlycerol = Cm

    return massFractionGlycerol


def calc_mass(volume, density):
    """
    Calculates the mass

    Args:
        volume (float): in m^3
         density (float): in kg/m^3

    Returns:
        mass (float): in kg
    """
    mass = volume * density
    return mass


def plot_dynamic_viscosity():
    """
    Plot the dynamic viscosity
    """
    # make a list of viscosities
    T = np.round(np.linspace(0, 100, 1001), 4)
    W = [n / 10.0 for n in range(0, 11)]
    G = [n / 10.0 for n in range(0, 11)]
    G = G.reverse()

    V = [(a, b, c, dynamic_viscosity_of_mixture(a, b, c)) for a, b, c
         in zip(list(chain(*([W] * 10))), list(chain(*([G] * 10))), T)]

    A = [n[0] for n in V]
    B = [n[1] for n in V]
    C = [n[2] for n in V]
    D = [n[3] for n in V]

    df = pd.DataFrame({'waterVolume': A, 'glycerolVolume': B,
                       'temperature': C, 'dynamicViscosity': D})

    # from ggplot import *

    # ggplot(data = df, aes())

    return V, df


def main():
    T = 20
    W = 4
    G = 15

    wString = "Water Volume [l]:\t{}".format(W)
    gString = "Glycerol Volume [l]:\t{}".format(G)
    tString = "Temperature [C]:\t{}".format(T)

    print("Parameter:\n{}\n{}\n{}\n".format(wString, gString, tString))

    print("Density of Water [kg/m^3]:")
    print("{}\n".format(water_density(T)))

    print("Density of Glycerol [kg/m^3]:")
    print("{}\n".format(glycerol_density(T)))

    print("Density of Mixture [kg/m^3]:")
    print("{}\n".format(calc_density_of_mixture(W, G, T)))

    print("Volume Fraction of Glycerol:")
    print("{}\n".format(calc_glycerol_fraction_by_volume(W, G)))

    print("Volume Fraction of Water:")
    print("{}\n".format(calc_water_fraction_by_volume(W, G)))

    print("Mass Fraction of Glycerol:")
    print("{}\n".format(calc_glycerol_fraction_by_mass(W, G, T)))

    print("Dynamic Viscosity of Water [Ns/m^2]:")
    print("{}\n".format(water_dynamic_viscosity(T)))

    print("Dynamic Viscosity of Glycerol [Ns/m^2]:")
    print("{}\n".format(glycerol_dynamic_viscosity(T)))

    print("Dynamic Viscosity of Mixture [Ns/m^2]:")
    print("{}\n".format(dynamic_viscosity_of_mixture(W, G, T)))

    print("Power law exponent:")
    Cm = calc_glycerol_fraction_by_mass(W, G, T)
    alpha = calc_alpha(Cm, T)
    print("{}\n".format(alpha))


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("\n\nCTRL-C detected, shutting down....")
        sys.exit(0)
