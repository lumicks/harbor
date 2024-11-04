#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import bisect
import matplotlib.pyplot as plt
import math
from foldometer.physics.thermodynamics import thermal_energy
from foldometer.physics.hydrodynamics import drag_sphere as dg_sph
from foldometer.physics.utils import as_Kelvin
from foldometer.physics.viscosity import dynamic_viscosity_of_mixture
import pandas as pd
from scipy.signal import welch
from collections import OrderedDict as od
from scipy.optimize import curve_fit
from foldometer.ixo.binary import *
from foldometer.analysis.thermal_calibration import *
import seaborn as sns


def find_peak(rawPSD, dFreq, peakSurroundings=10):
    """
    Function to find the peak corresponding to the mechanical oscillation. It will look for the maximum value
    in an interval of width=peakSurroundings centered in dFreq

    Args:
        rawPSD (pandas.DataFrame): DataFrame with the frequency and the PSD of a trap in certain direction
        dFreq (float): the external driven frequency of the oscillation in [Hz].
        peakSurroundings (float): width of the interval (in number of data points) to look for the peak. Default: 10

    Returns:
        peak (tuple): tuple with the peak height in [V^2/Hz] and the peak width in [Hz]
    """
    peakWidth = rawPSD.ix[1, "f"] - rawPSD.ix[0, "f"]
    peakHeight = rawPSD[abs(rawPSD["f"]-dFreq) < peakWidth*peakSurroundings]["psd"].max()

    return peakHeight, peakWidth


def calculate_experimental_area(peakHeight, peakWidth):
    """
    Function to get the area under the peak corresponding to the mechanical oscillation

    Args:
        peakHeight (float): height of the peak in [V^2/Hz]
        peakWidth (float): width of the peak in Hz

    Returns:
        peakArea (float): experimental area under the oscillation peak [V^2]

    """
    peakArea = peakHeight * peakWidth
    return peakArea


def calculate_theoretical_area(cornerFreq, dFreq=32, amplitude=150e-9):
    """
    Function to calculate the theoretical power response of the mechanical oscillation, according to the paper

    Args:
        cornerFreq (float): corner frequency of the power spectrum in [Hz]
        dFreq (float): the external driven frequency of the oscillation in [Hz].
        amplitude (float): amplitude of the mechanical oscillation

    Returns:
        theoreticalArea (float): theoretical area under the peak [m^2]
    """
    theoreticalArea = 0.5 * amplitude**2 / (1 + (cornerFreq / dFreq)**2)
    return theoreticalArea


def calculate_beta_from_file(fileName, dFreq=32, amplitude=150e-9, peakSurroundings=20, columns=None):
    """
    Extract the beta factor as well as the rest of the thermal calibration parameters from a file

    Args:
        fileName (str): string with the path to the file to calibrate
        dFreq (float): the external driven frequency of the oscillation in [Hz].
        amplitude (float): amplitude of the mechanical oscillation
        peakSurroundings (float): width of the interval (in number of data points) to look for the peak. Default: 10
        columns (list): list of strings indicating which columns should be included in the calculation.
            If None, all four relevant columns will be taken (both traps, both directions)

    Returns:
        calibrationData (pandas.DataFrame): DataFrame with all calibration parameters,
        including the beta extracted from the peak
    """
    header, calibration, data, beadData = read_file(fileName)

    if columns is None:
        columns = ['PSD1VxDiff', 'PSD1VyDiff', 'PSD2VxDiff', 'PSD2VyDiff']

    calibrationParameters = od([("betaSinusoidal", []), ("D", []), ("fc", []), ("sigma", []), ("chiSqr", []),
                                ("beta", []), ("stiffness", []), ("errorBeta", []), ("errorStiff", [])])

    for column in columns:
        single_fit = calculate_single_beta(data[column], header["sampleFreq"], amplitude, dFreq, peakSurroundings)
        for key in list(calibrationParameters.keys()):
            calibrationParameters[key].append(single_fit[key])

    calibrationData = pd.DataFrame(calibrationParameters, index=[col[0:4] + col[5] for col in columns])
    return calibrationData


def calculate_single_beta(data, sFreq, dFreq=32, amplitude=150e-9, peakSurroundings=20):
    """
    Extract the beta factor as well as the rest of the thermal calibration parameters from a column of data

    Args:
        data (pandas.Series): series with the data from one direction of one trap
        sFreq (float): sampling frequency of the measurement in [Hz]
        dFreq (float): the external driven frequency of the oscillation in [Hz].
        amplitude (float): amplitude of the mechanical oscillation
        peakSurroundings (float): width of the interval (in number of data points) to look for the peak. Default: 10

    Returns:
        fit (pandas.DataFrame): DataFrame with all calibration parameters,
        including the beta extracted from the peak
    """

    psd = calculate_psd(data.astype(float), sFreq=sFreq, blockNumber=2)

    if dFreq < 150:
        eliminate_peak(psd, dFreq, peakSurroundings)
        fit = calibration_psd(psd, blockNumber=2)

    expArea = calculate_experimental_area(*find_peak(psd, dFreq, peakSurroundings))
    theoArea = calculate_theoretical_area(fit["fc"], amplitude, dFreq)
    fit["betaSinusoidal"] = theoArea / expArea

    return fit


def eliminate_peak(rawPSD, dFreq=30, peakSurroundings=20):
    """
    Function to eliminate the peak for the standard thermal calibration IN PLACE

    Args:
        rawPSD (pandas.DataFrame): DataFrame with the frequency and the PSD of a trap in certain direction
        dFreq (float): the external driven frequency of the oscillation in [Hz].
        peakSurroundings (float): width of the interval (in number of data points) to look for the peak. Default: 10
    """
    peakHeight, peakWidth = find_peak(rawPSD, dFreq, peakSurroundings)
    peakIndex = pd.Index(rawPSD["psd"]).get_loc(peakHeight)

    interpolationData = rawPSD.ix[peakIndex - 5 : peakIndex-2].\
        append(rawPSD.ix[peakIndex + 2 : peakIndex + 5], ignore_index=True)
    meanSurroundings = interpolationData.psd.mean()

    stdSurroundings = interpolationData.psd.std()
    print(meanSurroundings, stdSurroundings)

    rawPSD.ix[peakIndex-1:peakIndex+1, "psd"] = rawPSD.ix[peakIndex-1 : peakIndex+1, "psd"].\
        apply(lambda x: meanSurroundings + (np.random.random() - 0.5)*stdSurroundings)


