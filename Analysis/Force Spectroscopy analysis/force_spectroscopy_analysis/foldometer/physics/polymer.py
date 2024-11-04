#!/usr/bin/env python
#-*- coding: utf-8 -*-

import numpy as np

from foldometer.physics.thermodynamics import thermal_energy
from foldometer.physics.utils import as_Kelvin
from collections import namedtuple


from statsmodels.formula.api import ols


CONTOUR_LENGTH = 960
PERSISTENCE_LENGTH = 45
STRETCH_MODULUS = 2100
ELASTIC_MODULUS = 800

def fit_wlc(fitData):
    """
    Fits the extensible worm-like chain model to data

    :param fitData:
    """

    WLC_FIT = namedtuple("WLC_FIT", ["L", "P", "S", "model", "fit"])

    kBT = thermal_energy()
    wlcModel = 'surfaceSeparation ~ 1 + np.sqrt(kBT/(force)) + force'

    try:
        wlcFit = ols(formula=wlcModel, data=fitData).fit()
    except:
        raise BaseException("Can't fit the data like this... ;-((")

    contourLength = wlcFit.params[0]
    persistenceLength = 1 / (-2 * (wlcFit.params[1] / wlcFit.params[0])) ** 2
    stretchModulus = 1 / (wlcFit.params[2] / wlcFit.params[0])

    fitResults = WLC_FIT(contourLength, persistenceLength, stretchModulus, wlcModel, wlcFit)

    return fitResults


def wlc(force, contourLength=CONTOUR_LENGTH, persistenceLength=PERSISTENCE_LENGTH, stretchModulus=STRETCH_MODULUS, T=21, A=0):
    """

    Args:
        force (array): array containing force data
        contourLength (float): contour length of the polymer
        persistenceLength (float): persistence length of the polymer.
        stretchModulus (float): stretch modulus of the polymer
        T (float): temperature in Celsius
        A (float): extension offset. Default 0

    Returns:
        extension (array): array containing the calculated extension according to the force
    """

    return contourLength * (1 - 0.5 * np.sqrt(thermal_energy(as_Kelvin(T)) / (force * persistenceLength))
                            + force / stretchModulus) + A



def extensible_WLC(extension, force, contourLength=CONTOUR_LENGTH, persistenceLength=PERSISTENCE_LENGTH,
                   elasticModulus=ELASTIC_MODULUS, T=21):
    #print((thermal_energy(as_Kelvin(T)) / persistenceLength))
    #print(np.mean(1 / (4 * (1 - extension / contourLength + force / elasticModulus) ** 2)))
    #print(np.mean(extension / contourLength))
    #print(np.mean(force / elasticModulus))

    return (thermal_energy(as_Kelvin(T)) / persistenceLength) * (
        (1 / (4 * (1 - extension / contourLength + force / elasticModulus) ** 2)
         - 1 / 4
         + extension / contourLength
         - force / elasticModulus))