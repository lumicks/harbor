#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import math
from foldometer.ixo.binary import read_file
from foldometer.physics.thermodynamics import thermal_energy
from foldometer.physics.hydrodynamics import drag_sphere as dg_sph
from foldometer.physics.utils import as_Kelvin
import pandas as pd
from scipy.signal import welch
from collections import OrderedDict as od
from scipy.optimize import curve_fit
#from foldometer.ixo.binary import *
import seaborn as sns
from itertools import count



# <editor-fold desc="General Functions and Parameters">

VISCOSITY = 9e-10
TEMPERATURE = 21
RADIUS = 1050
RADIUS2 = 1050
LIMITS = (150, 15000)
BLOCKNUMBER = 64
TSLENGTH = 1000000
SAMPLEFREQUENCY = 100000
OVERLAP = 0
PLOT = False
MODE = "mle"
DEFAULTARGS = {"viscosity": VISCOSITY, "temperature": TEMPERATURE, "radius": RADIUS, "radius1": RADIUS, "radius2": RADIUS2,
               "limits": LIMITS, "blockNumber": BLOCKNUMBER, "TSLength": TSLENGTH, "sFreq": SAMPLEFREQUENCY,
               "overlap": OVERLAP, "columnLabel": None, "figure": None, "fmCalibration": None, "plot": PLOT, "mode": MODE}

def lorentzian(f, D, fc):
    """
    Lorentzian function

    Args:
        f (numpy.array): frequency in units of [Hz]
        D (float): diffusion constant in units of [V]
        fc (float): corner frequency in units of [Hz]

    Returns:
        lorentzian(numpy.array): Lorentzian function with the corresponding parameters
    """

    return D/(math.pi**2*(f**2 + fc**2))


def corrected_lorentzian(f, D, fc, radius=1050, kineticViscosity=1e12, dynamicViscosity=9e-10,
                         fluidDensity=0.9e-21, beadDensity=1.05e-21):
    """
    Corrected lorentzian function with antialiasing and hydrodynamic coupling

    Args:
        f (numpy.array): frequency in units of [Hz]
        D (float): diffusion constant in units of [V]
        fc (float): corner frequency in units of [Hz]
        radius (float): radius of the bead in units of [nm]
        kineticViscosity (float): kinetic viscosity of the fluid in units of [nm^2/s]
        dynamicViscosity (float): dynamic viscosity of the fluid in units of [pN s / nm]
        fluidDensity (float): density of the fluid in units of [g/nm^3]
        beadDensity (float): density of the bead in units of [g/nm^3]

    Returns:
        CorrectedLorentzianValues(numpy.array): Lorentzian function with the corresponding parameters
    """
    
    fv = kineticViscosity / (math.pi * radius**2)
    #fv = 330000
    effectiveMass = 4/3 * math.pi * radius**3 * beadDensity + 2 * math.pi * fluidDensity * radius **3 / 3
    fm = dg_sph(radius, dynamicViscosity) / (2 * math.pi * effectiveMass)
    #fm = 480000
    return D/(np.pi**2) * (1 + np.sqrt(f/fv)) / \
            ((fc - f**1.5/np.sqrt(fv) - f**2 / fm)**2 + (f + f**1.5 / np.sqrt(fv))**2)


def distance_calibration(D, radius=RADIUS, viscosity=VISCOSITY, temperature=TEMPERATURE):
    """Distance calibration factor (beta) in units of [V/nm]

    Args:
        D (float): diffusion constant in units of [V]
        radius (float): radius of the bead in units of [nm]
        viscosity (float): viscosity of the solution in units of [pN/nm^2s]
        temperature (float): temperature in units of [ºC]

    Returns:
        beta (float): distance calibration factor in units of [V/nm]
    """

    beta = 1/np.sqrt(thermal_energy(as_Kelvin(temperature))/(dg_sph(radius, viscosity)*D))

    return beta


def trap_stiffness(fc, radius=RADIUS, viscosity=VISCOSITY):
    """Trap stiffness in units of [pN/nm]

    Args:
        fc (float): corner frequency in units of [Hz]
        radius (float): radius of the bead in units of [nm]
        viscosity (float): viscosity of the solution in units of [pN/nm^2s]

    Returns:
        kappa (float): trap stiffness in units of [pN/nm]
    """

    kappa = 2*math.pi*fc*dg_sph(radius, viscosity)

    return kappa


def force_calibration(beta, kappa):
    """force calibration factor in units of [pN/V]

    Args:
        beta (float): distance calibration factor in units of [V/nm]
        kappa (float): trap stiffness in units of [pN/nm]

    Returns:
        alpha (float): force calibration factor in units of [pN/V]
    """

    alpha = kappa/beta

    return alpha


def residuals(psd, D, fc):
    """Performs a chi^2 test with the average of the residuals squared to test the best fitting limits

    Args:
        psd (pandas.DataFrame): two columns: frequency in [Hz] and power spectrum density in [V^2]
        D (float): diffusion constant in units of [V]
        fc (float): corner frequency in units of [Hz]
    Returns:
        res (float): the mean of the residuals squared
    """

    psdTheo = lorentzian(psd["f"], D, fc)
    res = np.mean([(exp-theo)**2 for exp, theo in zip(psd["psd"], psdTheo)])

    return res


def calculate_psd(data, blockNumber=BLOCKNUMBER, sFreq=SAMPLEFREQUENCY, overlap=OVERLAP):
    """Calculates the power spectral density of the data and sets the proper pandas format

    Args:
        data (pandas.DataFrame): data read from the file
        blockNumber (float): number of averaged blocks in which the time series is divided for calculating the psd
        sFreq (float): sample frequency in units of [Hz]
        overlap (float): number of points for overlapping blocks

    Returns:
        psd (pandas.DataFrame): power spectrum density of the data
    """

    #calculate the length of each block for the Welch algorithm
    blockLength = len(data.index)//blockNumber

    fRaw, psdRaw = welch(data.values.astype(float), nperseg=blockLength, window="hann", fs=sFreq, noverlap=overlap)

    # if there is a value for 0 Hz, delete it
    if fRaw[0] == 0:
        fRaw = np.delete(fRaw, 0)
        psdRaw = np.delete(psdRaw, 0)
    #set format to pandas
    psd = pd.DataFrame({"psd": psdRaw, "f": fRaw})

    return psd


def check_kwargs(kwargs):
    """Fill missing values in the arguments with the default values

    Args:
        viscosity (float): viscosity of the solution in units of [pN/nm^2s] (Default: 8.93e-10, pure water at 25 ºC)
        temperature (float): temperature in ºC
        radius (float): radius of the first bead in units of [nm]
        radius2 (float): radius of the second bead in units of [nm]
        limits (tuple): frequency limits for consideration in the fitting.
        blockNumber (float): number of averaged blocks in which the time series is divided for calculating the psd
        sFreq (float): sample frequency in units of [Hz]
        overlap (float): number of points for overlapping blocks
        mode (string): the method to use for the fitting:
            "lstsq": normal least squares fitting of a Lorentzian
            "wlstsq": weighted (theoretical weights: 1/PSD) least squares fitting of a Lorentzian
            "lstsq_corrected": hydrodynamic and alias corrections for the Lorentzian (Default)
            "mle": maximum likelihood estimator fitting of a Lorentzian (DEFAULT)
            "mean": the arithmetic mean of the parameters from LSTSQ and MLE fitting
    Returns:
        kwargs (dict): new keyword arguments with filled missing data with default values

    """
    for key in DEFAULTARGS:
        if key not in list(kwargs.keys()):
            kwargs[key] = DEFAULTARGS[key]

    return kwargs


def plot_calibration(psd, D, fc, limits=LIMITS, fmCalibration=None, columnLabel=None, figure=None, mode=MODE):
    """
    Function to plot the calibration

    Args:
        psd (pandas.DataFrame): two columns: frequency in [Hz] and power spectrum density in [V^2]
        psdLim (pandas.DataFrame): two columns within the set limits: frequency in [Hz] and power spectrum density in [
        V^2]
        D (float): diffusion constant in units of [V]
        fc (float): corner frequency in units of [Hz]
        columnLabel (str): name of the channel to be plotted
        figure (dict): dictionary with the axes of the plots (one for each PSD channel). None if only one channel.
        mode (str): the method to use for the fitting:
            "lstsq": normal least squares fitting of a Lorentzian
            "wlstsq": weighted (theoretical weights: 1/PSD) least squares fitting of a Lorentzian
            "lstsq_corrected": hydrodynamic and alias corrections for the Lorentzian (Default)
            "mle": maximum likelihood estimator fitting of a Lorentzian
            "mean": the arithmetic mean of the parameters from LSTSQ and MLE fitting

    """


    psdLim = psd[(psd["f"] > limits[0]) & (psd["f"] < limits[1])]

    if mode == "lstsq_corrected":
        fitData = [corrected_lorentzian(x, D, fc) for x in psdLim["f"]]
    else:
        fitData = [lorentzian(x, D, fc) for x in psdLim["f"]]

    if figure is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    else:
        ax = figure[columnLabel]

    ax.plot(psd["f"], psd["psd"], '.', alpha=0.4)
    ax.plot(psdLim["f"], fitData, color='orange', label='New Python fit')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('PSD (V^2)')
    ax.set_xlim(min(psd["f"]), max(psd["f"]))
    ax.set_title(columnLabel)

    if fmCalibration is not None:

        Dfm = thermal_energy(as_Kelvin(TEMPERATURE))/(dg_sph(RADIUS, VISCOSITY) * fmCalibration["beta"]**2)

        fitDataFoldometer = [lorentzian(x, D, fmCalibration["cornerFrequency"][columnLabel[:4]+columnLabel[5]]) for x in psdLim["f"]]
        ax.plot(psdLim["f"], fitDataFoldometer, color='red', label='Imported fit from machine')
    ax.legend()
    if figure is None:
        plt.show()
# </editor-fold>


# <editor-fold desc="Least squares calibration">
def lstsq_calibration(psd, n, p0, weights=None, mode="lstsq"):
    """Function to perform the least squares fitting of data

    Args:
        psd (pandas.DataFrame): two columns: frequency in [Hz] and power spectrum density in [V^2]
        n (float): number of averaged power spectra (total data points divided by the block length)
        p0 (list): the two initial guesses for the parameters of the fit (D and fc, respectively)
        weights (pandas.Series): array with the weight of each data point for the fit. Default: None
        mode (str): either "lstsq" or "lstsq_corrected"

    Returns:
        * **D** (float): fitted experimental diffusion coefficient in [V^2]
        * **fc** (float): fitted corner frequency in [Hz]
        * **errors** (float):  fitting errors for D anf fc int he same units as those
        * **chiSqr** (float):  Chi squared factor representing the goodness of the fit
    """

    if mode is "lstsq" or "lstsq_corrected":
        result, errors = curve_fit(lorentzian, psd["f"], psd["psd"], p0=p0)
    else:
        result, errors = curve_fit(corrected_lorentzian, psd["f"], psd["psd"], p0=p0)

    if result[1] < 0:
        result[1] = -result[1]

    # standard deviation errors
    errors = np.sqrt(np.diag(errors))

    D, fc = result

    if weights is not None:
        D = D*n/(n+1)

    chiSqr = residuals(psd, D, fc)

    return D, fc, errors, chiSqr
# </editor-fold>


# <editor-fold desc="Maximum Likelihood Estimator Functions">
def mle_factors(f, psd):
    """
    Calculation of the S coefficients related to the MLE, according to the paper of Norrelike

    Args:
        f (numpy.array): Frequency in [Hz]
        psd (numpy.array): Experimental PSD function in [V^2]

    Returns:
        s (list): matrix with the S coefficients
    """
    N = len(f)
    s01 = 1/N * np.sum(psd)
    s02 = 1/N * np.sum(np.power(psd, 2))
    s11 = 1/N * np.sum(np.multiply(np.power(f, 2), psd))
    s12 = 1/N * np.sum(np.multiply(np.power(f, 2), np.power(psd, 2)))
    s22 = 1/N * np.sum(np.multiply(np.power(f, 4), np.power(psd, 2)))
    s = [[0, s01, s02], [0, s11, s12], [0, s12, s22]]

    return s


def mle_ab(s, n):
    """
    Calculation of the pre-parameters a and b, according to the paper of Norrelike

    Args:
        s (list): matrix with the S coefficients
        n (float): number of averaged power spectra (total data points divided by the block length)

    Returns:
        * **a** (float): first pre-parameter for the calculation of D and fc
        * **b** (float): second pre-parameter for the calculation of D and fc
    """

    a = ((1+1/n)/(s[0][2]*s[2][2]-s[1][2]*s[1][2])) * (s[0][1]*s[2][2]-s[1][1]*s[1][2])
    b = ((1+1/n)/(s[0][2]*s[2][2]-s[1][2]*s[1][2])) * (s[1][1]*s[0][2]-s[0][1]*s[1][2])
    return a, b


def mle_parameters(a, b, n):
    """Calculate parameters from the factors of the MLE

    Args:
        a, b (float): pre-parameters for the calculation of D and fc
        n (float): number of averaged power spectra (total data points divided by the block length)

    Returns:
        D (float): diffusion constant in units of [V]
        fc (float): corner frequency in units of [Hz]

    """

    if a*b > 0:
        fc = math.sqrt(a/b)
    else:
        fc = 0
    D = (math.pi**2) / b

    return D, fc


def mle_errors(f, D, fc, a, b, n):
    """Function to get the standard deviation of the parameters according to the paper of Norrelyke

    Args:
        f (numpy.array): array of the frequencies in units of [Hz]
        D (float): diffusion constant in units of [V]
        fc (float): corner frequency in units of [Hz]
        a, b (float): pre-parameters for the calculation of D and fc
        n (float): number of averaged power spectra

    Returns:
        errorsMle (numpy.array): with sigma(D) and sigma(fc)
    """
    y = lorentzian(f, D, fc)
    s = mle_factors(f, y)
    sB = [[(n+1)/n*s[0][2], (n+1)/n*s[1][2]], [(n+1)/n*s[1][2], (n+1)/n*s[2][2]]]
    sError = 1/(len(f)*n)*(n+3)/n*np.linalg.inv(sB)

    sigmaFc = fc**2/4 * (sError[0][0]/a**2+sError[1][1]/b**2-2*sError[0][1]/(a*b))
    sigmaD = D**2*(sError[1][1]/b**2)
    errorsMle = [np.sqrt(sigmaD), np.sqrt(sigmaFc)]

    return errorsMle


def mle_calibration(psd, n):
    """Function to perform the Maximum Likelihood Estimator fitting of data

    Args:
        psd (pandas.DataFrame): two columns: frequency in [Hz] and power spectrum density in [V^2]
        n (float): number of averaged power spectra (total data points divided by the block length)


    Returns:
        * **D** (float): diffusion constant in units of [V]
        * **fc** (float): corner frequency in units of [Hz]
        * **errors** (list): standard deviation of D and fc
        * **chiSqr** (float): average of the squared residues (test of chi^2)
    """

    s = mle_factors(psd["f"], psd["psd"])
    a, b = mle_ab(s, n)
    D, fc = mle_parameters(a, b, n)

    errors = mle_errors(psd["f"], D, fc, a, b, n)
    #chiSqr = residuals(psd, D, fc)
    s0 = np.sum(np.divide(psd["psd"], lorentzian(psd["f"], D, fc)))
    chiSqr = math.erfc(abs(s0 - len(psd["psd"])) * np.sqrt(n/(len(psd["psd"]))))

    return D, fc, errors, chiSqr
# </editor-fold>


# <editor-fold desc="Maximum Likelihood Estimator Functions for Aliased PSD">
def mle_factors_aliased(f, psd, n):
    """
    Calculation of the S coefficients related to the MLE, according to the paper of Norrelike

    Args:
        f (numpy.array): Frequency in [Hz]
        psd (numpy.array): Experimental PSD function in [V^2]

    Returns:
        s (list): matrix with the S coefficients
    """
    K = len(f)

    f = np.arange(K)

    R01 = 1/K * np.sum(psd)
    R02 = 1/K * np.sum(np.power(psd, 2))
    R11 = 1/K * np.sum(np.multiply(np.cos(2*np.pi*f/(K*n)), psd))
    R12 = 1/K * np.sum(np.multiply(np.cos(2*np.pi*f/(K*n)), np.power(psd, 2)))
    R22 = 1/K * np.sum(np.multiply(np.power(np.cos(2*np.pi*f/(K*n)), 2), np.power(psd, 2)))
    R = [[0, R01, R02], [0, R11, R12], [0, R12, R22]]

    return R


def mle_ab_aliased(R, n):
    """
    Calculation of the pre-parameters a and b, according to the paper of Norrelike

    Args:
        s (list): matrix with the S coefficients
        n (float): number of averaged power spectra (total data points divided by the block length)

    Returns:
        * **a** (float): first pre-parameter for the calculation of D and fc
        * **b** (float): second pre-parameter for the calculation of D and fc
    """

    a = ((1+1/n)/(R[0][2]*R[2][2]-R[1][2]*R[1][2])) * (R[0][1]*R[2][2]-R[1][1]*R[1][2])
    b = ((1+1/n)/(R[0][2]*R[2][2]-R[1][2]*R[1][2])) * (R[1][1]*R[0][2]-R[0][1]*R[1][2])
    return a, b


def mle_parameters_aliased(a, b, fsample):
    """Calculate parameters from the factors of the MLE

    Args:
        a (float): first pre-parameters for the calculation of D and fc
        b (float): second pre-parameters for the calculation of D and fc
        fsample (float): sample frequency in [Hz]

    Returns:
        * **D** (float): diffusion constant in units of [V]
        * **fc** (float): corner frequency in units of [Hz]
    """

    u = np.arccosh(-a/b)
    fc = fsample * u / (2 * np.pi)

    D = fsample**2 * u / (a * np.tanh(u))

    return D, fc


def mle_errors_aliased(f, D, fc, a, b, n):
    """Function to get the standard deviation of the parameters according to the paper of Norrelyke

    Args:
        f (numpy.array): array of the frequencies in units of [Hz]
        D (float): diffusion constant in units of [V]
        fc (float): corner frequency in units of [Hz]
        a, b (float): pre-parameters for the calculation of D and fc
        n (float): number of averaged power spectra (total data points divided by the block length)

    Returns:
        errosMle (numpy.array): with sigma(D) and sigma(fc)
    """
    y = lorentzian(f, D, fc)
    s = mle_factors(f, y)
    sB = [[(n+1)/n*s[0][2], (n+1)/n*s[1][2]], [(n+1)/n*s[1][2], (n+1)/n*s[2][2]]]
    sError = 1/(len(f)*n)*(n+3)/n*np.linalg.inv(sB)

    sigmaFc = fc**2/4 * (sError[0][0]/a**2+sError[1][1]/b**2-2*sError[0][1]/(a*b))
    sigmaD = D**2*(sError[1][1]/b**2)
    errorsMle = [np.sqrt(sigmaD), np.sqrt(sigmaFc)]

    return errorsMle


def mle_calibration_aliased(psd, n):
    """Function to perform the Maximum Likelihood Estimator fitting of data

    Args:
        psd (pandas.DataFrame): two columns: frequency in [Hz] and power spectrum density in [V^2]
        n (float): number of averaged power spectra (total data points divided by the block length)


    Returns:
        D (float): diffusion constant in units of [V]
        fc (float): corner frequency in units of [Hz]
        errors (list): standard deviation of D and fc
        chiSqr (float): average of the squared residues (test of chi^2)
    """
    fsample = 100000

    s = mle_factors_aliased(psd["f"], psd["psd"], n)
    a, b = mle_ab_aliased(s, n)
    D, fc = mle_parameters_aliased(a, b, fsample)

    errors = mle_errors_aliased(psd["f"], D, fc, a, b, n)
    chiSqr = residuals(psd, D, fc)

    return D, fc, errors, chiSqr
# </editor-fold>


def calibration_fit(psd, **kwargs):
    """Function to fit a Power Spectrum Density distribution to a Lorentzian

    Args:
        psd (pandas.DataFrame): two columns: frequency in [Hz] and power spectrum density in [V^2]
        **kwargs: optional keyword arguments. Refer to the documentation of check_kwargs for possible arguments

    Returns:
        D (float): diffusion constant in units of V
        fc (float): corner frequency in units of Hz
        errors (list): standard deviation of the parameters
        chiSqr (float): average of the squared residues (test of chi^2)
    """

    kwargs = check_kwargs(kwargs)

    #eliminate the points out of the limit range
    psdLim = psd[(psd["f"] > kwargs["limits"][0]) & (psd["f"] < kwargs["limits"][1])]

    #calculate parameters with MLE to use them as initial guess for the rest
    D, fc, errors, chiSqr = mle_calibration(psdLim, kwargs["blockNumber"])

    if kwargs["mode"] is "lstsq" or "lstsq_corrected" or "wlstsq":
        if kwargs["mode"] is "wlstsq":
            weights = [1/x for x in lorentzian(psdLim["f"])]
        D, fc, errors, chiSqr = lstsq_calibration(psdLim, kwargs["blockNumber"], p0=[D, fc], mode=kwargs["mode"])
    elif kwargs["mode"] == "wlstsq":
        D, fc, errors, chiSqr = lstsq_calibration(psdLim, kwargs["blockNumber"], p0=[D, fc],
                                                  weights=[1/x for x in lorentzian(psdLim["f"], D, fc)])
    elif kwargs["mode"] == "mean":
        D2, fc2, errors2, chiSqr2 = lstsq_calibration(psdLim, kwargs["blockNumber"], p0=[D, fc])
        D = (D + D2)/2
        fc = (fc + fc2)/2
        errors = (errors + errors2)/2
        chiSqr = (chiSqr + chiSqr2)/2
    elif kwargs["mode"] == "aliased_mle":
        D, fc, errors, chiSqr = mle_calibration_aliased(psdLim, kwargs["blockNumber"])

    if kwargs["plot"]:
        plot_calibration(psd, D, fc, kwargs["limits"], kwargs["fmCalibration"], kwargs["columnLabel"], kwargs["figure"], kwargs["mode"])

    return D, fc, errors, chiSqr


def calibration_psd(psd, **kwargs):
    """Performs the fitting from a Power Spectrum Distribution series

    Args:
        data (pandas.DataFrame): single time series with data
        **kwargs: optional keyword arguments. Refer to the documentation of check_kwargs for possible arguments
    Returns:
        fitParameters (orderedDict): values of the fit: diffusion constant (D), corner frequency (fc), \
        errors of D and fc (sigma), distance calibration factor (beta), trap stiffness (kappa), \
        error of beta (eBeta), error of kappa (eKappa) and limits of the fit (limits)
    """

    kwargs = check_kwargs(kwargs)

    #calculate the number of averaged spectra (nBlocks)
    if kwargs["overlap"] != 0:
        kwargs["blockNumber"] *= kwargs["overlap"]

    D, fc, sigma, chiSqr = calibration_fit(psd, **kwargs)

    #create an ordered dictionary for the fitting parameters
    fitParameters = od()
    fitParameters["D"] = D
    fitParameters["fc"] = fc
    fitParameters["sigma"] = sigma
    fitParameters["chiSqr"] = chiSqr
    fitParameters["beta"] = distance_calibration(D, kwargs["radius"], kwargs["viscosity"], kwargs["temperature"])
    fitParameters["stiffness"] = trap_stiffness(fc, kwargs["radius"], kwargs["viscosity"])
    fitParameters["errorBeta"] = (sigma[0]/D)*fitParameters["beta"]
    fitParameters["errorStiff"] = (sigma[1]/fc)*fitParameters["stiffness"]
    fitParameters["beadDiameter"] = kwargs["radius"] * 2

    return fitParameters


def calibration_time_series(data, **kwargs):
    """Performs the fitting from a time series data

    Args:
        data (pandas.Series): single time series with data
        **kwargs: optional keyword arguments. Refer to the documentation of check_kwargs for possible arguments

    Returns:
        fitParameters (orderedDict): values of the fit: diffusion constant (D), corner frequency (fc), \
        errors of D and fc (sigma), distance calibration factor (beta), trap stiffness, \
        error of beta (eBeta), error of stiffness (eStiffness)
    """

    kwargs = check_kwargs(kwargs)

    #calculation of the power spectrum density function with the Welch algorithm
    psd = calculate_psd(data, kwargs["blockNumber"], kwargs["sFreq"], kwargs["overlap"])

    kwargs["columnLabel"] = data.name

    fitParameters = calibration_psd(psd, **kwargs)

    return fitParameters


def calibration_data(data, metadata=None, fmCalibration=None, **kwargs):
    """Performs the fitting from a time series data

    Args:
        data (pandas.Series): single time series with data
        metadata (pandas.Series): header values read from the file
        fmCalibration (pandas.DataFrame): parameters from the thermal calibration fit performed in Foldometer
        **kwargs: optional keyword arguments. Refer to the documentation of check_kwargs for possible arguments

    Returns:
        fitParameters (orderedDict): values of the fit: diffusion constant (D), corner frequency (fc), \
        errors of D and fc (sigma), distance calibration factor (beta), trap stiffness, \
        error of beta (eBeta), error of stiffness (eStiffness)

    Notes:
        If metadata and fmCalibration are not given, the default values will be used
    """
    kwargs = check_kwargs(kwargs)

    #Prepare relevant data like sample frequency and the PSd data
    if metadata is not None:
        kwargs["sFreq"] = metadata["sampleFreq"]
        if "radius" not in kwargs.keys():
            kwargs["radius1"] = metadata["beadRadius"]
        if "radius2" not in kwargs.keys():
            kwargs["radius2"] = metadata["beadRadius"]
        if "viscosity" not in kwargs.keys():
            kwargs["viscosity"] = metadata["viscosity"]
        if "temperature" not in kwargs:
            kwargs["temperature"] = metadata["temperature"]


    if fmCalibration is not None:
        kwargs["fmCalibration"] = fmCalibration

    #choose only the relevant columns
    psdData = data[['PSD1VxDiff', 'PSD1VyDiff', 'PSD2VxDiff', 'PSD2VyDiff']].astype(float)

    #define the ordered dictionary and the labels for the calibration parameters
    calibrationParameters = od()
    columnLabels = ["D", "fc", "sigma", "chiSqr", "beta", "stiffness", "errorBeta", "errorStiff", "beadDiameter"]
    for label in columnLabels:
        calibrationParameters[label] = []

    #Prepare the figure to plot the four calibrations in the same canvas
    if kwargs["plot"]:
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(20,12), sharex=True, sharey=True)
        kwargs["figure"] = {'PSD1VxDiff': ax1, 'PSD1VyDiff': ax2, 'PSD2VxDiff': ax3, 'PSD2VyDiff': ax4}

    #Perform the calibration for each of the four columns
    for label, column in psdData.iteritems():
        kwargs["columnLabel"] = label

        #Check which radius to use
        if "PSD1" in label:
            kwargs["radius"] = kwargs["radius1"]
        else:
            kwargs["radius"] = kwargs["radius2"]

        fitTimeSeries = calibration_time_series(data=column.astype(float), **kwargs)

        for key in calibrationParameters.keys():
            calibrationParameters[key].append(fitTimeSeries[key])

    if kwargs["plot"]:
        plt.show()

    #Create the pandas DataFrame
    calibrationData = pd.DataFrame(calibrationParameters, index=["PSD1x", "PSD1y", "PSD2x", "PSD2y"])

    return calibrationData


def calibration_file(file, limits=LIMITS, blockNumber=BLOCKNUMBER, overlap=OVERLAP, plot=PLOT, mode=MODE, **kwargs):
    """
    Function to calibrate a Foldometer data file

    Args:
        file (str): path of the file containing the data
        limits (list): frequency limits for consideration in the fitting. Default LIMITS (normally [150,15000] Hz)
        blockNumber (float): number of averaged blocks in which the time series is divided for calculating the psd
        overlap (float): number of points for overlapping blocks (Default: 0)
        plot (bool): if True, the fittings are plotted. (Default: False)
        mode (str): the method to use for the fitting:
            "lstsq": normal least squares fitting of a Lorentzian
            "wlstsq": weighted (theoretical weights: 1/PSD) least squares fitting of a Lorentzian
            "lstsq_corrected": hydrodynamic and alias corrections for the Lorentzian
            "mle": maximum likelihood estimator fitting of a Lorentzian (Default)
            "mean": the arithmetic mean of the parameters from LSTSQ and MLE fitting
        **kwargs: optional keyword arguments. Refer to the documentation of check_kwargs for possible arguments
    Returns:
        calibrationData (pandas.DataFrame): a pandas DataFrame with the fit parameters and errors for all channels

    """

    #Read the file
    metadata, fmCalibration, data, beadData = read_file(file, mode="calibration")

    calibrationData = calibration_data(data, metadata, fmCalibration, limits=limits, blockNumber=blockNumber,
                                       overlap=overlap, plot=plot, mode=mode, **kwargs)

    return calibrationData
