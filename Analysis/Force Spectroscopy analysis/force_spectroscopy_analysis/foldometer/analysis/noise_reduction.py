#!/usr/bin/env python
# -*- coding: utf-8 -*-


import pandas as pd
import numpy as np
import scipy.fftpack as fft
import os
import pkg_resources

P_CCD_DEFAULT_FILE_NAME = pkg_resources.resource_filename('foldometer', 'data/p_CCD.xlsx')
P_QPD_PIEZO_DEFAULT_FILE_NAME = pkg_resources.resource_filename('foldometer', 'data/p_QPD_piezo.xlsx')


def weight_from_ps(pCCD, pQPD, a=5):
    """
    Function to calculate the weights according the paper of Alireza and Vach:
    http://scitation.aip.org/content/aip/journal/rsi/82/11/10.1063/1.3658825

    Args:
        p1 (pandas.Series): first of the two noise spectra
        p2 (pandas.Series): second of the two noise spectra

    Returns:
        weight (numpy.array): the calculated weights
    """

    return np.asarray(1 / (np.exp(a) - 1) * ((np.exp(a) + 1) / (1 + np.exp(a * (pQPD - pCCD) / (pCCD + pQPD))) - 1))


def calculate_weights(data, pCCDFile=P_CCD_DEFAULT_FILE_NAME, pQPDFile=P_QPD_PIEZO_DEFAULT_FILE_NAME,
                      sampleFrequencyPiezo=50, sampleFrequencyCCD=50, a=5, join=True):
    """
    Calculate the weights for noise correction in the Fourier space using the noise power sepctra already recorded
    in the old setup. Read Alireza's and P. Vach paper for more information:
    http://scitation.aip.org/content/aip/journal/rsi/82/11/10.1063/1.3658825

    Args:
        data (pandas.DataFrame): frame containing all the data
        pCCDFile (str): path of the file containing the measured noise spectrum (p in the paper) of the CCD
        pQPDFile (str): path of the file containing the measured noise spectrum (p in the paper) of the QPD piezo
        sampleFrequencyPiezo (float): sample frequency of the data from the piezo
        sampleFrequencyCCD (float): sample frequency of the data from the camera
        a (float): factor of weighting
        join (bool): if True, there is no return and the weights are joined to the whole dataset. If False,
        nothing is joined and the function returns a pandas.Series with the weights

    """
    pCCD = pd.read_excel(pCCDFile, header=None)[0]
    pQPD = pd.read_excel(pQPDFile, header=None)[0]
    calibrationLength = len(pCCD)
    dataLength = len(data.forceX)

    weightsArray = weight_from_ps(pCCD, pQPD, a)
    #Mirror (by padding) the array to cover the symmetric Fourier transform of the data
    weightsArray = np.lib.pad(weightsArray, (0, len(weightsArray) - 1), 'reflect')
    #Create a pandas.Series of the array for convenience and rename it as "weights"
    weights = pd.Series(weightsArray)
    weights.name = "weights"

    if join:
        #Change the indexes to distribute the weights along the whole measured dataset
        weights.index = (np.arange(1, len(weightsArray) + 1) * (dataLength / 2) / calibrationLength).astype(int)
        data = data.join(weights)
        #Interpolate the weights to fill for missing values (if the measured dataset is bigger than the sample noise)
        data["weights"] = data["weights"].interpolate().fillna(method="bfill")

        return data

    else:
        return weights


def calculate_noise_fft(data):
    """
    Function to calculate the Fourier transform of the noise in the piezo and bead tracking data

    Args:
        data (pd.DataFrame): frame containing all the measured data
    """
    signalCCD = (data["xspt2"] - data["xspt1"]) * 0.0615
    signalQPD = data["surfaceSepX"] / 1000

    noiseCCD = signalCCD - signalCCD.mean()
    noiseQPD = signalQPD - signalQPD.mean()

    data["fourierNoiseCCD"] = fft.fft(noiseCCD)
    data["fourierNoiseQPD"] = fft.fft(noiseQPD)


def apply_weights_to_noise(data):
    """
    Function to apply the weights to the noise in the Fourier space

    Args:
        data (pandas.DataFrame): frame containing all the measured data

    """

    data["fourierNoiseCCD"] *= (1 - data["weights"])
    data["fourierNoiseQPD"] *= data["weights"]



def subtract_noise_from_data(data):
    """
    Function to correct the noise signal from the extension extracted from the CCD tracking

    Args:
        data (pandas.DataFrame): frame containing all the measured data
    Returns:


    """

    data["fourierNoiseCombined"] = data["fourierNoiseCCD"] + data["fourierNoiseQPD"]
    data.rename(columns={"surfaceSepX": "surfaceSepXQPD"}, inplace=True)
    data["surfaceSepX"] = np.real(fft.ifft(data["fourierNoiseCombined"]) + np.mean((data.xspt2 - data.xspt1) * 0.0615))
    data.drop(["fourierNoiseCCD", "fourierNoiseQPD", "fourierNoiseCombined"], axis=1, inplace=True)


def correct_signal_noise(data, **kwargs):
    """
    Function to reduce the noise in the signal of the extension, using the data from the piezo mirror and the camera
    tracking. It is mostly thought to be used for old setup data. The code is based on the paper by Alireza and Vach;
    http://scitation.aip.org/content/aip/journal/rsi/82/11/10.1063/1.3658825
    Args:
        data (pandas.DataFrame): frame containing all the measured data
        **kwargs: optional arguments for calculate_weights function. Check that function for more information
    Returns:
        data (pandas.DataFrame): frame containing all the measured data and the additional corrected extension
    """
    data = calculate_weights(data, **kwargs)
    calculate_noise_fft(data)
    apply_weights_to_noise(data)
    subtract_noise_from_data(data)
    data.index = data["time"]
    data.surfaceSepX *=1000
    return data