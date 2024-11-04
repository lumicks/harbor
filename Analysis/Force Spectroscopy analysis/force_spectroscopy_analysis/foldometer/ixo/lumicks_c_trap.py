#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from collections import OrderedDict as od
import pandas as pd
from nptdms import TdmsFile
import matplotlib.pyplot as plt
from foldometer.tools.maths import cross_correlation
from lmfit.models import GaussianModel, ExponentialModel
from skimage.filters import threshold_mean
import os

CHANNEL_LABEL_MAPPING = {"/'Sensor Data'/'Time (ms)'": "time",
                         "/'Sensor Data'/'Force Channel 0 (V)'": "trapSepX",
                         "/'Sensor Data'/'Force Channel 1 (V)'": "trapSepY",
                         "/'Sensor Data'/'Force Channel 2 (V)'": "2trapSepX",
                         "/'Sensor Data'/'Force Channel 3 (V)'": "2trapSepY",
                         "/'Power Spectrum Data'/'Fitted Power Spectrum Channel 0 (V^2/Hz)'": "PSD1ForceX",
                         "/'Power Spectrum Data'/'Fitted Power Spectrum Channel 1 (V^2/Hz)'": "PSD1ForceY",
                         "/'Power Spectrum Data'/'Fitted Power Spectrum Channel 2 (V^2/Hz)'": "PSD2ForceX",
                         "/'Power Spectrum Data'/'Fitted Power Spectrum Channel 3 (V^2/Hz)'": "PSD2ForceY"}

CALIBRATION_LABEL_MAPPING = {'Bead Diameter (um)': "beadDiameter", 'Viscosity (Pa*s)': 'viscosity',
                             'Temperature (C)': 'temperature', 'Lower Bound (Hz)': 'Lower Bound (Hz)',
                             'Upper Bound (Hz)': 'Upper Bound (Hz)', 'Exclusion Ranges (Hz)': 'Exclusion Ranges (Hz)',
                             'Corner Frequency (Hz)': "cornerFrequency",
                             'Force Response (pN/V)': "alpha", 'Distance Response (um/V)': "distanceResponse",
                             'Trap Stiffness (pN/m)': "stiffness", 'RMSE': "RMSE", 'Applied': "Applied"}

FLUORESCENCE_LABEL_MAPPING = {"/'Data'/'Actual position X (um)'": "positionX",
                              "/'Data'/'Actual position Y (um)'": "positionY",
                              "/'Data'/'Pixel ch 1'": "638nm",
                              "/'Data'/'Pixel ch 2'": "532nm",}


def lumicks_file(moleculeNumber, fileNumber=1, experiment="", folder="."):
    """
    Get the file name of a c-trap measurement from the exp and molecule number.
    Args:
        moleculeNumber (int): the number of the molecule as set in the program
        fileNumber (int): the number of the file as set in the program
        experiment (str): any string included in the file name
        folder (str): folder where to look for files

    Returns:
        lumicksFile (dict): dictionary containing the data file, the corresponding calibration file and fluo file (if)
    """

    lumicksFile = {"dataFile": None, "powerSpectrumFile": None, "fluoFile": None}
    for fileName in os.listdir(folder):
        if fileName.endswith(str(moleculeNumber).zfill(3) + "-" + str(fileNumber).zfill(3) + ".tdms") and experiment in fileName:
            lumicksFile["dataFile"] = os.path.join(folder, fileName)
    for spectrumFileName in os.listdir(folder):
        if spectrumFileName.endswith("Power Spectrum.tdms") and spectrumFileName[9:15] < os.path.split(lumicksFile["dataFile"])[1][9:15]:
            lumicksFile["powerSpectrumFile"] = os.path.join(folder, spectrumFileName)
    try:
        for fluoFileName in os.listdir(folder+"/fluo/"):
            if fluoFileName.endswith(str(moleculeNumber).zfill(3) + "-" + str(fileNumber).zfill(3) + ".tdms"):
                lumicksFile["fluoFile"] = os.path.join(folder, "fluo", fluoFileName)
    except:
        pass
    return lumicksFile

import streamlit as st
@st.cache()
def read_file_lumicks(dataFilePath, fluorescenceFilePath=None):
    """
    Read a data file from C-trap

    Args:
        dataFilePath (str): path of the file containing the force data
        fluorescenceFilePath (str): path of the file containing the fluorescence data

    Returns:
        data (pandas.DataFrame): relevant data for analysis
    """

    if fluorescenceFilePath is None:
        fluorescenceFilePath = str("./fluo" + dataFilePath)
        try:
            fluorescenceFile = TdmsFile(fluorescenceFilePath)
        except:
            Warning("No fluorescence file found")
    else:
        fluorescenceFile = TdmsFile(fluorescenceFilePath)

    data = read_data_file_lumicks(dataFilePath)

    return data


def extract_calibration_parameters(calibrationFilePath):
    """
    Get stiffness and diffusion coefficients from lumicks calibration file

    Args:
        calibrationFilePath (str): path for the calibration file

    Returns:
        calibrationFit (pandas.DataFrame): data of the calibrated values
    """

    calibrationFile = TdmsFile(calibrationFilePath)
    calibrationMap = {"PSD1x": "/'Power Spectrum Data'/'Fitted Power Spectrum Channel 0 (V^2/Hz)'",
                      "PSD1y": "/'Power Spectrum Data'/'Fitted Power Spectrum Channel 1 (V^2/Hz)'",
                      "PSD2x": "/'Power Spectrum Data'/'Fitted Power Spectrum Channel 2 (V^2/Hz)'",
                      "PSD2y": "/'Power Spectrum Data'/'Fitted Power Spectrum Channel 3 (V^2/Hz)'"}
    calibrationFit = pd.DataFrame({key: 4*[0] for key in CALIBRATION_LABEL_MAPPING.values()},
                                  index=["PSD1x", "PSD1y", "PSD2x", "PSD2y"])

    for key in calibrationMap:
        for parameter in CALIBRATION_LABEL_MAPPING:
            calibrationFit.loc[key, CALIBRATION_LABEL_MAPPING[parameter]] = \
                calibrationFile.objects[calibrationMap[key]].properties[parameter]

    calibrationFit.loc[:, "stiffness"] *= 1e-9

    return calibrationFit


def read_fluorescence_file_lumicks(fluorescenceFilePath):
    """
    Read and parse a C-trap fluorescence file, only kymographs

    Args:
        fluorescenceFilePath (str): path to the fluorescence file
    Returns:
        images (dict): dictionary with the scanning kymographs for both colors
        timeRes (float): time resolution of the kymograph
    """

    fluorescenceFile = TdmsFile(fluorescenceFilePath)
    metadata = fluorescenceFile.object().properties
    fluorescenceData = fluorescenceFile.as_dataframe().rename(index=str, columns=FLUORESCENCE_LABEL_MAPPING)
    diff = fluorescenceData["/'Data'/'Time (ms)'"].diff()
    timeRes = diff.loc[diff > 0].mean() / 1000
    pixelsPerLine = metadata["Pixels per line"]
    pixelsLength = len(fluorescenceData["638nm"]) // pixelsPerLine
    images = {"638nm": np.matrix.transpose(fluorescenceData["638nm"].values.reshape(pixelsLength, pixelsPerLine)),
              "532nm": np.matrix.transpose(fluorescenceData["532nm"].values.reshape(pixelsLength, pixelsPerLine))}

    return images, timeRes


def track_bead_fluorescence_LEGACY(image, startLine=0, correlationLength=1000, lineInterval=(66, -1), threshold=False):
    """
    Track the bead movement from the fluorescent data, using the cross-correlation of a reference line with the rest

    Args:
        image (numpy.ndarray): matrix containing the image information
        startLine (int): scanning line to take as a first reference
        referenceCount (int): number of additional reference lines to take for the average
        correlationLength (int): number of lines to cross-correlate to the reference
        halfLine (int): horizontal limit separating image in two to take only moving bead.
        threshold (bool): if True, perform a binary threshold to the image before the cross-correlation

    Returns:
        meanCorrelation (numpy.array): array with the averaged tracked bead movement based on cross-correlations
    """
    if threshold:
        thresh = threshold_mean(image)
        binary = image > thresh
        image = binary

    maxCorrelation = [0]


        #compare lines to reference and get the maximum position, which represents the pixel shift
    for line in np.arange(startLine, startLine + correlationLength):

        positiveCorr = cross_correlation(image[lineInterval[0]:lineInterval[1], line],
                                                       image[lineInterval[0]:lineInterval[1],line+1])
        negativeCorr = cross_correlation(image[lineInterval[0]:lineInterval[1], line + 1],
                                                       image[lineInterval[0]:lineInterval[1], line])
        corr = np.concatenate((np.flip(negativeCorr, 0), positiveCorr))


        x=np.arange(len(corr)) - (len(corr) - 1) / 2

        y=corr - min(corr)
        #plt.plot(x,y)
        gauss1  = GaussianModel(prefix='g1_')
        pars = gauss1.make_params()
        pars['g1_center'].set(0, min=-10, max=10)
        pars['g1_sigma'].set(1, min=0.1)
        pars['g1_amplitude'].set(100, min=2)
        mod = gauss1
        #init = mod.eval(pars, x=x)
        out = mod.fit(y, pars, x=x)
        #comps = out.eval_components(x=x)
#print(out.fit_report(min_correl=0.5))
        #plt.plot(x,y)
        #plt.plot(x, out.best_fit, 'r-')
        maxCorrelation.append(out.params["g1_center"].value + maxCorrelation[-1])
        #print(referenceLine + lineCount, line, out.params["g1_center"].value )
        #correlationPeaks = np.argmax(positiveCorr)
        #if correlationPeaks == 0:
        #    correlationPeaks = -np.argmax(negativeCorr)
        #correlation.append(correlationPeaks)
    #maxCorrelation.append(correlation)
    #print(x)
    #plt.show()
    #print(maxCorrelation)
    #meanCorrelation = np.mean(maxCorrelation, axis=0)

    return np.array(maxCorrelation)


def align_fluorescence_data_legacy(data, image, extensionTimeLength=20, fluoResolution=0.087, **kwargs):
    """
    Calculate the temporal shift of the fluorescence with respect to the PSD data by comparing the time position of the
    peak in the distance between surfaces.
    Args:
        data (pandas.DataFrame): whole dataframe with force extension data
        image (numpy.ndarray): matrix containing the fluorescence information
        extensionTimeLength (int): seconds to consider for the PSD data containing trap movement
        fluoResolution (float): estimated fluorescence time resolution acquisition
        **kwargs: keyword arguments for the function track_bead_fluorescence()

    Returns:
        fluorescenceTimeOffset (float): time shift in seconds between the PSD data and the fluorescence data
    """
    maxCorrelation = track_bead_fluorescence(image, **kwargs)

    peakFluo = np.argmax(maxCorrelation[np.count_nonzero(np.isnan(maxCorrelation)):]) + \
               np.count_nonzero(np.isnan(maxCorrelation))
    peakFluo *= fluoResolution

    mask = data["time"] < extensionTimeLength
    peakData = data.loc[np.argmax(data.loc[mask, "surfaceSepX"]), "time"]
    fluorescenceTimeOffset = peakData - peakFluo

    return fluorescenceTimeOffset

def read_data_file_lumicks(dataFilePath, compact=True):
    """

    Args:
        dataFilePath (str): path of the file containing the force data

    Returns:

    """
    CHANNEL_LABEL_MAPPING = {"Time (ms)": "time",
                             "Distance 1 (um)": "trapSepX",
                             "Distance 2 (um)": "trapSepY",
                             "Force Channel 0 (pN)": "PSD1ForceX",
                             "Force Channel 1 (pN)": "PSD1ForceY",
                             "Force Channel 2 (pN)": "PSD2ForceX",
                             "Force Channel 3 (pN)": "PSD2ForceY"}
    group_name = 'FD Data'
    with TdmsFile.open(dataFilePath) as tdms_file:
        channel = tdms_file[group_name]["Time (ms)"]
        n = len(channel)
        dictSlice = {}
        for channel_name in CHANNEL_LABEL_MAPPING:
            channel = tdms_file[group_name][channel_name]
            dictSlice[CHANNEL_LABEL_MAPPING[channel_name]] = channel[:]

        df = pd.DataFrame(dictSlice)
        df["time"] /= 1000
        df["trapSepX"] *= 1e3
        df["MirrorX"] = df.loc[:, "trapSepX"]
        df["forceX"] = (df["PSD2ForceX"] - df["PSD1ForceX"]) / 2
        data = pd.DataFrame(df)
    return data


import streamlit as st
@st.cache()
def process_lumicks_data(data, calibrationFitCTrap, calibrationFitPython=None):
    if calibrationFitPython is not None:
        print("Using Python calibration")
        correctionBeta = calibrationFitCTrap.loc[:, "beta"] / calibrationFitPython.loc[:, "beta"]

        data["PSD1xDisplacement"] = (data["PSD1ForceX"] / calibrationFitCTrap.loc["PSD1x", "stiffness"]) * \
                                    correctionBeta.loc["PSD1x"]
        data["PSD2xDisplacement"] = (data["PSD2ForceX"] / calibrationFitCTrap.loc["PSD2x", "stiffness"]) * \
                                     correctionBeta.loc["PSD2x"]
        data["PSD1ForceX"] = data["PSD1xDisplacement"] * calibrationFitPython.loc["PSD1x", "stiffness"]
        data["PSD2ForceX"] = data["PSD2xDisplacement"] * calibrationFitPython.loc["PSD2x", "stiffness"]
        data["forceX"] = (data["PSD2ForceX"] - data["PSD1ForceX"]) / 2
        calibrationFit = calibrationFitPython
    else:
        calibrationFit = calibrationFitCTrap
        data["PSD1xDisplacement"] = data["PSD1ForceX"] / calibrationFitCTrap.loc["PSD1x", "stiffness"]
        data["PSD2xDisplacement"] = data["PSD2ForceX"] / calibrationFitCTrap.loc["PSD2x", "stiffness"]


    data["surfaceSepX"] = data["trapSepX"] - calibrationFit.loc["PSD1x", "beadDiameter"] / 2 \
                          - calibrationFit.loc["PSD2x", "beadDiameter"] / 2 \
                          - data["PSD2ForceX"] / calibrationFit.loc["PSD2x", "stiffness"] \
                          + data["PSD1ForceX"] / calibrationFit.loc["PSD1x", "stiffness"]

    return data

def process_lumicks_data_old(data, calibrationFit):

    data["PSD1xDisplacement"] = data["PSD1ForceX"] / calibrationFit.loc["PSD1x", "stiffness"]
    data["PSD2xDisplacement"] = data["PSD2ForceX"] / calibrationFit.loc["PSD2x", "stiffness"]
    data["surfaceSepX"] = data["trapSepX"] - calibrationFit.loc["PSD1x", "beadDiameter"] / 2 \
                          - calibrationFit.loc["PSD2x", "beadDiameter"] / 2 \
                          - data["PSD2ForceX"] / calibrationFit.loc["PSD2x", "stiffness"] \
                          + data["PSD1ForceX"] / calibrationFit.loc["PSD1x", "stiffness"]

    return data