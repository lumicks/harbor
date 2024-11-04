# coding=utf-8

import numpy as np


def as_Kelvin(temperatureInCelsius=25):
    """
    Converts temperature from Celsius to Kelvin

    Args:
        temperatureInCelsius (float): Temperature in [°C]. (Default: 25 °C)

    Returns:
        temperatureInKelvin (float): Temperature in Kelvin [K]
    """
    assert temperatureInCelsius >= -273.15

    temperatureInKelvin = 273.15 + temperatureInCelsius

    return temperatureInKelvin


def as_Celsius(temperatureInKelvin=298):
    """
    Converts temperature from Kelvin to Celsius

    Args:
        temperatureInKelvin (float): Temperature in Kelvin [K]. (Default: 298 K)

    Returns:
        temperatureInCelsius (float): Temperature in [°C]. (Default: 25 °C)
    """
    assert temperatureInKelvin >= 0

    temperatureInCelsius = temperatureInKelvin - 273.15

    return temperatureInCelsius


def volume_sphere(radius: float=1000):
    """
    Volume of a sphere

    Args:
        radius (float): Radius in [nm]. (Default: 1000 nm)

    Returns:
        volume (float): Volume in [nm³]
    """
    assert radius > 0

    volume = 4 / 3 * np.pi * radius**3

    return volume


def mass_sphere(radius: float=1000, density: float=1e-21,
                verbose: bool=False, inferDensity: bool=True) -> float:
    """
    Calculates mass of a sphere

    Args:
        radius (float): Radius in [nm]. (Default: 1000 nm)
        density (float): Density, ϱ, in [g/nm³]. (Default: 1e-21 g/nm³)
        verbose (bool): Flag whether to show extra information. (Default: False)
        inferDensity (bool): If true inputs above 1e-5 are interpreted as g/cm³. (Default: True)

    Returns:
        mass (float):  Mass in [g]
    """
    assert radius > 0
    assert density > 0

    volume = volume_sphere(radius=radius)

    # check if density is provided in g/cm³.
    if inferDensity and density > 1e-5:
        print("Assuming that ϱ was provided in [g/cm³]. Converting to g/nm³")
        density /= 1e21
        print("ϱ is now: {} g/nm³".format(density))

    mass = volume * density

    if verbose:
        print("Input:")
        print("Radius of Sphere [nm]: {}".format(radius))
        print("Density of material [g/nm³]: {} \n".format(density))

        print("Output:")
        print("Mass of sphere [g]: {}\n".format(mass))

    return mass