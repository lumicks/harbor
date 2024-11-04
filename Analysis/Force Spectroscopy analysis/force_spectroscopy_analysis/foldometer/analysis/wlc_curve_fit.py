#!/usr/bin/env python
# -*- coding: utf-8 -*-

#from foldometer.tools.plots import force_extension_curve
from foldometer.tools.misc import data_selection
from foldometer.physics.thermodynamics import thermal_energy
from foldometer.physics.utils import as_Kelvin
from foldometer.analysis.tweezers_parameters import PROTEIN_LENGTHS
from scipy import ndimage
from scipy.stats import norm
import lmfit
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from copy import deepcopy


FORCE_MARGIN = 0.01
CONTOUR_LENGTH_DNA = 906
PERSISTENCE_LENGTH_DNA = 40
STRETCH_MODULUS_DNA = 1000
CONTOUR_LENGTH_PROTEIN = 120
PERSISTENCE_LENGTH_PROTEIN = 0.5
TEMPERATURE = 25
PARAMETER_RANGE = 0.5
EXTENSION_OFFSET = 0.01
FORCE_OFFSET = 0
MIN_FORCE = 0.1
MAX_FORCE = 50
FIXED_PARAMETERS = ["contourLengthDNA", "persistenceLengthProtein", "temperature", "extensionOffset"]
DEFAULT_PROPERTIES = {"contourLengthDNA": CONTOUR_LENGTH_DNA, "persistenceLengthDNA": PERSISTENCE_LENGTH_DNA,
                      "stretchModulusDNA": STRETCH_MODULUS_DNA, "contourLengthProtein": CONTOUR_LENGTH_PROTEIN,
                      "persistenceLengthProtein": PERSISTENCE_LENGTH_PROTEIN, "temperature": TEMPERATURE,
                      "extensionOffset": EXTENSION_OFFSET, "parameterRange": PARAMETER_RANGE,
                      "fixedParameters": FIXED_PARAMETERS, "forceMargin": FORCE_MARGIN, "minForce": MIN_FORCE,
                      "maxForce": MAX_FORCE}


def check_kwargs(kwargs):
    """Fill missing values in the arguments with the default value

    Args:
        contourLengthDNA (float): contour length of the DNA handles in nm (Default: 906 nm)
        persistenceLengthDNA (float): persistence length of the DNA handles in nm (Default: 40 nm)
        stretchModulusDNA (float): stretch modulus of DNA handles in pN (Default: 1300 pN)
        contourLengthProtein (float): contour length of the protein in nm (Default: 120 nm)
        persistenceLengthProtein (float): persistence length of the protein in nm (Default: 0.5 nm)
        temperature (float): temperature in C
        extensionOffset (float): offset for the extension (Default: 0.01 for lmfit reasons)
        parameterRange (float): fraction of the default value for upper and lower limits in the fit (Default: 0.5)
        fixedParameters (list): list of labels for the parameters that are not to be fitted
        forceMargin (float): margin for subtracting the force in order to perform the fits
        minForce (float): minimum force to be considered in the WLC fit
        maxForce (float): maximum force to be considered in the WLC fit
    Returns:
        kwargs (dict): new keyword arguments with filled missing data with default values

    """
    for key in DEFAULT_PROPERTIES:
        if key not in list(kwargs.keys()):
            kwargs[key] = DEFAULT_PROPERTIES[key]

    return kwargs


def identify_individual_curves(data):
    """
    Function to identify and label each of the wlc regions in the data, including before and after an unfolding event

    Args:
        data (pandas.DataFrame): Measured data
    """

    assert ("unfolding" in data.columns), "In order to perform wlc fits, please analyse unfolding events first. This " \
                                          "can be done with the FoldoMeasurement class method analyse_data() or " \
                                          "manually with find_unfolding_events()"
    wlcLabels = ndimage.label([not value for value in data["unfolding"]])[0] - 1

    wlcLabels[wlcLabels == -1] = - data["pullingCycle"].max() - 1
    data["wlcRegion"] = data["pullingCycle"] + wlcLabels
    data.loc[data["region"] == "stationary", "wlcRegion"] = -1
    data.loc[data["wlcRegion"] < 0, "wlcRegion"] = -1


def calculate_force_offset(data, axis="X", forceChannel="force", margin=FORCE_MARGIN):
    """
    Calculates the force offset of each of the wlc curves in order to perform a better fit

    Args:
        data (pandas.DataFrame): Measured data
        axis (str): axis on which perform the analysis, either "X" or "Y"
        forceChannel (str): either "PSD1Force", "PSD2Force" or "force" (combined signal)
        margin (float): margin for subtracting the force in order to perform the fits

    Returns:
        forceOffsets (dict): dictionary with the offset for the force in [pN] for each wlc curve

    """
    forceOffsets = {}
    for wlcRegion in data["wlcRegion"].unique():
        forceOffsets[wlcRegion] = data.loc[data["wlcRegion"] == wlcRegion, forceChannel + axis].min() - margin
        print("forceOffset", data.loc[data["wlcRegion"] == wlcRegion, forceChannel + axis].min())
    return forceOffsets


def calculate_extension_offset(model, data, tmin, tmax, axis="X", forceChannel="force", distanceChannel="surfaceSep", **kwargs):
    """
    Calculates the global extension offset in order to perform the wlc fit, based on the curve corresponding to
    folded protein

    Args:
        model (lmfit.model.Model): model of worm-like-chain to be fitted
        data (pandas.DataFrame): Measured data
        axis (str): axis for which calculate the fit, either "X" or "Y"
        forceChannel (str): either "PSD1Force", "PSD2Force" or "force" (combined signal)
        distanceChannel (str): either "surfaceSep" (from PSDs) or "trackingSep" (from image tracking)
        **kwargs: optional keyword arguments. Refer to the documentation of wlc_curve_fit.check_kwargs for arguments
    Returns:
        extensionOffset (dict): dictionary with the offset for the extension in [nm] for each wlc curve

    """
    maskTime = (data["time"]>=tmin)*(data["time"]<=tmax)

    parameters = model.make_params()
    parameters["contourLengthDNA"].set(value=kwargs["contourLengthDNA"], vary=False)
    parameters["contourLengthProtein"].set(value=kwargs["contourLengthProtein"], vary=False)
    parameters["extensionOffset"].set(value=0, vary=True, min=-np.inf, max=np.inf)
    curve = data.loc[maskTime]

    initialFoldedFit = model.fit(curve[distanceChannel + axis], params=parameters, force=curve[forceChannel + axis],
                                 weights=curve[forceChannel + axis], nan_policy='omit')
    #print(curve["forceX"].min(), curve["forceX"].max())
    extensionOffset = initialFoldedFit.params["extensionOffset"].value

    parameters["contourLengthDNA"].set(value=kwargs["contourLengthDNA"], vary=True)
    parameters["contourLengthProtein"].set(value=kwargs["contourLengthProtein"], vary=True)
    parameters["extensionOffset"].set(value=kwargs["extensionOffset"], vary=False)

    return extensionOffset


def include_fit_curves(curve, fit, axis="X", forceChannel="force"):
    """
    Function to include the theoretical WLC in the dataset

    Args:
        curve (pandas.DataFrame) Subset of data of each identified WLC curve
        fit (lmfit.model.ModelResult): complete information of the fitting result for a particular curve
        axis (str): axis on which perform the analysis, either "X" or "Y"
        forceChannel (str): either "PSD1Force", "PSD2Force" or "force" (combined signal)

    Returns:
        curve (pandas.DataFrame): DataFrame containing the theoretical wlc force and extension from the fit
    """

    curve.loc[:, "wlcForce" + axis] = curve[forceChannel + axis]
    #curve.loc[:, "wlcExtension" + axis] = fit.best_fit
    curve.loc[:, "wlcExtension" + axis] = wlc_from_fit(curve[forceChannel + axis], fit)


    return curve.loc[:, ["wlcForce" + axis, "wlcExtension" + axis]]


def wlc_weights(curve, mean=30, sigma=10, plot=False):
    """

    Args:
        curve (pandas.Series):
        mean (float): where to center the normal distribution for weights
        sigma (float): width of the normal distribution
        plot (bool): if True, plot the weight distribution
    Returns:
        weights (numpy.array): array containing the weigth for each point of the curve
    """

    weights = norm(loc=mean, scale=sigma)
    if plot:
        plt.plot(curve, weights.pdf(curve))
        plt.show()
    return weights.pdf(curve)



def wlc_fit_single_curve(model, curve, axis="X", forceChannel="force", distanceChannel="surfaceSep",
                         construct="protein", upperContourLength=10, weightMode="norm", weightNormMean=30,
                         weightNormSigma=10):
    """
    Function to fit a single pulling or retracting curve to the WLC

    Args:
        model (lmfit.model.Model): model of worm-like-chain to be fitted
        curve (pandas.DataFrame): interval of the data corresponding to a single force-extension curve
        axis (str): axis for which calculate the fit, either "X" or "Y"
        forceChannel (str): either "PSD1Force", "PSD2Force" or "force" (combined signal)
        distanceChannel (str): either "surfaceSep" (from PSDs) or "trackingSep" (from image tracking)
        construct (str): either "DNA" or "protein"
        upperContourLength (float): margin given to the maximum possible contour length of protein
        weightMode (str): None will not apply weights, "high" will give higher weights to high forces and "norm" will
        use a normal distribution around weightNormMean with width weightNormSigma
        weightNormMean (float): the center of the weight distribution if weightMode is "norm"
        weightNormSigma (float): the width of the weight distribution is weightMode is "norm"
    Returns:
        wlcResult (lmfit.model.ModelResult): complete information of the fitting (check lmfit docs for more info)
    """

    parameters = model.make_params()
    if construct is "protein":
        proteinMaxLc = parameters["contourLengthProtein"].value + upperContourLength
        parameters["contourLengthProtein"].set(min=0, max=proteinMaxLc)
    if weightMode is "norm":
        weights = wlc_weights(curve.loc[:, forceChannel + axis], mean=weightNormMean, sigma=weightNormSigma)
    elif weightMode is "high":
        weights = curve.loc[:, forceChannel + axis]
    elif weightMode is None:
        weights = None
    wlcResult = model.fit(curve.loc[:, distanceChannel + axis], params=parameters,
                          force=curve.loc[:, forceChannel + axis], weights=weights)

    #
    return wlcResult


def wlc_fit_data(data, tmin, tmax, axis="X", forceChannel="force", distanceChannel="surfaceSep", joinFitCurves=True,
                 calculateForceOffset=False, calculateExtensionOffset=True, construct="protein",
                 protein="MBP", upperLc=10, weightMode="norm", weightNormMean=30, weightNormSigma=10, **kwargs):
    """
    Function to identify, classify and fit each pulling and retracting curve to the WLC model

    Args:
        data (pandas.DataFrame): Measured data
        axis (str): axis on which perform the analysis, either "X" or "Y"
        forceChannel (str): either "PSD1Force", "PSD2Force" or "force" (combined signal)
        distanceChannel (str): either "surfaceSep" (from PSDs) or "trackingSep" (from image tracking)
        calculateForceOffset (bool): if True, the force signal will be shifted
        calculateExtensionOffset (bool): if True, a fit is performed to calculate the offset in extension
        joinFitCurves (bool): offset in the extension in [nm]
        construct (str): either "protein" or "DNA", to generate a WLC model in accordance
        protein (str): name of the protein, to automatically set the proper contour length
        **kwargs: optional keyword arguments. Refer to the documentation of wlc_curve_fit.check_kwargs for arguments

    Returns:
        * **data** (pandas.DataFrame): new DataFrame in case of joinFitCurves=True, old one otherwise
        * **fits** (dict):  parameters from the wlc fits
    """

    if protein in PROTEIN_LENGTHS.keys() and "contourLengthProtein" not in kwargs:
        kwargs["contourLengthProtein"] = PROTEIN_LENGTHS[protein]
    else:
        print("Warning: the contour length of your protein is not registered in the database, falling back to MBP if "
              "not specified")
    kwargs = check_kwargs(kwargs)
    #kwargs["fixedParameters"] = list(set(kwargs["fixedParameters"].extend(FIXED_PARAMETERS)))
    kwargs["fixedParameters"].extend(FIXED_PARAMETERS)
    kwargs["fixedParameters"] = list(set(kwargs["fixedParameters"]))
    print(kwargs["fixedParameters"])


    #Create the lmfit model with the series WLC for protein or the single WLC model for DNA
    if construct is "protein":
        fittingParamsLabels = ["contourLengthDNA", "persistenceLengthDNA", "stretchModulusDNA",
                               "contourLengthProtein", "persistenceLengthProtein", "stretchModulusProtein",
                               "temperature", "extensionOffset"]
        wlcModel = lmfit.Model(wlc_series_accurate)
    elif construct is "DNA":
        fittingParamsLabels = ["contourLengthDNA", "persistenceLengthDNA", "stretchModulusDNA", "temperature",
                               "extensionOffset"]
        wlcModel = lmfit.Model(wlc)
        calculateExtensionOffset = False
    #Add fitting parameters to the model
    for key in kwargs.keys():
        if key in fittingParamsLabels:
            paramMin = kwargs[key] * (1 - kwargs["parameterRange"])
            paramMax = kwargs[key] * (1 + kwargs["parameterRange"])
            if key in kwargs["fixedParameters"]:
                wlcModel.set_param_hint(key, value=kwargs[key], vary=False, min=paramMin, max=paramMax)
            else:
                wlcModel.set_param_hint(key, value=kwargs[key], vary=True, min=paramMin, max=paramMax)

    #separate data in wlc single curves
    # identify_individual_curves(data)

    #calculate force and extension offsets and correct data
    if calculateForceOffset:
        forceOffsets = calculate_force_offset(data, axis, forceChannel, kwargs["forceMargin"])
        print("Force offsets: ", forceOffsets)
    else:
        forceOffsets = {}
        for wlcRegion in data["wlcRegion"].unique():
            forceOffsets[wlcRegion] = 0

    if calculateExtensionOffset:
        extensionOffset = calculate_extension_offset(wlcModel, data, tmin, tmax, axis, forceChannel, distanceChannel, **kwargs)
        print("Extension offset: ", extensionOffset)
    else:
        extensionOffset = 0

    data[forceChannel + axis] -= min(forceOffsets.values())
    # data.loc[:, distanceChannel + axis] -= extensionOffset

    fits = {}
    if joinFitCurves:
        data["wlcForce" + axis] = np.nan
        data["wlcExtension" + axis] = np.nan

    for wlcRegion in data["wlcRegion"].unique():
        if wlcRegion > -1:
            mask = (data["wlcRegion"] == wlcRegion) & (data["region"] != "stationary")
            maskedData = data[mask]
            forceMask = (maskedData[forceChannel + axis] > kwargs["minForce"]) & \
                        (maskedData[forceChannel + axis] < kwargs["maxForce"])
            #make the WLC fit only if the segment is longer than 4 data points
            if len(maskedData[forceMask]) > 4:
                fits[wlcRegion] = wlc_fit_single_curve(wlcModel, maskedData[forceMask], axis, forceChannel,
                                                       distanceChannel, construct, upperLc, weightMode, weightNormMean,
                                                       weightNormSigma)
            #include the fitted extension and force to the DataFrame
            if joinFitCurves and len(maskedData[forceMask]) > 4:
                data.loc[mask, ["wlcForce" + axis, "wlcExtension" + axis]] = \
                    include_fit_curves(data[mask], fits[wlcRegion], axis, forceChannel)

    return data, fits


def wlc(force, contourLengthDNA=CONTOUR_LENGTH_DNA, persistenceLengthDNA=PERSISTENCE_LENGTH_DNA,
        stretchModulusDNA=STRETCH_MODULUS_DNA, temperature=21, extensionOffset=0):
    """
    Function calculating a WLC model

    Args:
        force (numpy.array): array containing force data
        contourLength (float): contour length of the polymer
        persistenceLength (float): persistence length of the polymer.
        stretchModulus (float): stretch modulus of the polymer
        temperature (float): temperature in Celsius
        extensionOffset (float): offset in the distance between surfaces. Default 0

    Returns:
        extension (array): array containing the calculated extension according to the force
    """

    kbT = thermal_energy(as_Kelvin(temperature))

    return contourLengthDNA * (1 - 0.5 * np.sqrt(kbT / (force * persistenceLengthDNA)) + force / stretchModulusDNA) + \
           extensionOffset


def wlc_series(force, contourLengthDNA=CONTOUR_LENGTH_DNA, persistenceLengthDNA=PERSISTENCE_LENGTH_DNA,
               stretchModulusDNA=STRETCH_MODULUS_DNA, contourLengthProtein=CONTOUR_LENGTH_PROTEIN,
               persistenceLengthProtein=PERSISTENCE_LENGTH_PROTEIN, temperature=TEMPERATURE, extensionOffset=0):
    """
    Function calculating two WLC in series: one for the DNA handles (extensible WLC, according to Odijk 1995) and
    another for the protein, without any stretch modulus.

    Args:
        force (numpy.array): array containing force data
        contourLengthDNA (float): contour length of the DNA handles
        persistenceLengthDNA (float): persistence length of the DNA handles
        stretchModulusDNA (float): stretch modulus of the DNA handles
        contourLengthProtein (float): contour length of the protein
        persistenceLengthProtein (float): persistence length of the protein
        temperature (float): temperature in Celsius
        extensionOffset (float): offset in the distance between surfaces. Default 0

    Returns:
        extension (array): array containing the calculated extension from the force according to the WLC model
    """

    kbT = thermal_energy(as_Kelvin(temperature))

    return contourLengthDNA * (1 - 0.5 * np.sqrt(kbT / (force * persistenceLengthDNA)) + force / stretchModulusDNA) + \
           contourLengthProtein * (1 - 0.5 * np.sqrt(kbT / (force * persistenceLengthProtein))) + extensionOffset

def wlc_series_accurate(force, contourLengthDNA=CONTOUR_LENGTH_DNA, persistenceLengthDNA=PERSISTENCE_LENGTH_DNA,
                        stretchModulusDNA=STRETCH_MODULUS_DNA, contourLengthProtein=CONTOUR_LENGTH_PROTEIN,
                        persistenceLengthProtein=PERSISTENCE_LENGTH_PROTEIN, temperature=TEMPERATURE,
                        extensionOffset=EXTENSION_OFFSET):
    kbT = thermal_energy(as_Kelvin(temperature))
    A = (force * persistenceLengthDNA) / kbT
    B = np.exp((900 / A) ** (1/4))

    return contourLengthDNA * (4 / 3 * (1 - 1 / np.sqrt(A + 1)) - 10 * B / (np.sqrt(A) * (B -1) ** 2) +
                               A ** 1.62 / (3.55 + 3.8 * A ** 2.2) + force / stretchModulusDNA) + \
           contourLengthProtein * (1 - 0.5 * np.sqrt(kbT / (force * persistenceLengthProtein))) + extensionOffset



def wlc_from_fit(force, fit, accurate=True):
    args = (fit.params["contourLengthDNA"].value, fit.params["persistenceLengthDNA"].value,
              fit.params["stretchModulusDNA"].value, fit.params["contourLengthProtein"].value,
              fit.params["persistenceLengthProtein"].value, fit.params["temperature"].value, 0)
    if accurate:
        return wlc_series_accurate(force, *args)
    else:
        return wlc_series(force, *args)


def protein_contour_length(extension, force, contourLengthDNA=CONTOUR_LENGTH_DNA,
                           persistenceLengthDNA=PERSISTENCE_LENGTH_DNA, stretchModulusDNA=STRETCH_MODULUS_DNA,
                           contourLengthProtein=CONTOUR_LENGTH_PROTEIN,
                           persistenceLengthProtein=PERSISTENCE_LENGTH_PROTEIN, temperature=TEMPERATURE):
    """
    Function to calculate the protein contour length from the force and extension data
    Args:
        extension (numpy.array): array containing extemsopm data
        force (numpy.array): array containing force data
        contourLengthDNA (float): contour length of the DNA handles
        persistenceLengthDNA (float): persistence length of the DNA handles
        stretchModulusDNA (float): stretch modulus of the DNA handles
        persistenceLengthProtein (float): persistence length of the protein
        temperature (float): temperature in Celsius

    Returns:
        LcProtein (pandas.Series): array with the protein contour length
    """
    kbT = thermal_energy(as_Kelvin(temperature))

    return (extension - contourLengthDNA * (1 - 0.5 * np.sqrt(kbT / (force * persistenceLengthDNA)) +
                                            force / stretchModulusDNA)) / \
           (1 - 0.5 * np.sqrt(kbT / (force * persistenceLengthProtein)))


def protein_contour_length_accurate(extension, force, contourLengthDNA=CONTOUR_LENGTH_DNA,
                           persistenceLengthDNA=PERSISTENCE_LENGTH_DNA, stretchModulusDNA=STRETCH_MODULUS_DNA,
                           contourLengthProtein=CONTOUR_LENGTH_PROTEIN,
                           persistenceLengthProtein=PERSISTENCE_LENGTH_PROTEIN, temperature=TEMPERATURE):
    """
    Function to calculate the protein contour length from the force and extension data
    Args:
        extension (numpy.array): array containing extemsopm data
        force (numpy.array): array containing force data
        contourLengthDNA (float): contour length of the DNA handles
        persistenceLengthDNA (float): persistence length of the DNA handles
        stretchModulusDNA (float): stretch modulus of the DNA handles
        persistenceLengthProtein (float): persistence length of the protein
        temperature (float): temperature in Celsius

    Returns:
        LcProtein (pandas.Series): array with the protein contour length
    """

    kbT = thermal_energy(as_Kelvin(temperature))
    A = (force * persistenceLengthDNA) / kbT
    B = np.exp((900 / A) ** (1/4))
    return (extension - contourLengthDNA * (4 / 3 * (1 - 1 / np.sqrt(A + 1)) - 10 * B / (np.sqrt(A) * (B -1) ** 2) +
                               A ** 1.62 / (3.55 + 3.8 * A ** 2.2) + force / stretchModulusDNA)) / \
           (1 - 0.5 * np.sqrt(kbT / (force * persistenceLengthProtein)))