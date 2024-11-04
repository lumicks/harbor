#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from foldometer.tools.misc import data_selection
import numpy as np


def track_gaussian_feature(fluoData, startLine, endLine=None, lineEdges=None, fluoTimeRes=0.08596, shiftCenter=0,
                           averageLines=0, invert=False):
    """
    Track a Gaussian-resembling feature in the fluorescence using a Gaussian fit. The feature can be the edge of a
    bead, a bright spot in the middle or a whole bead (by using invert)/

    Args:
        fluoData (numpy.array): data containing the fluorescence fluoData
        startLine (int): line to start tracking
        endLine (int): line to end the tracking
        lineEdges (tuple): limits for the fit. If None, select a region from the profile
        fluoTimeRes (float): time resolution of the kymograph lines
        shiftCenter (int): displacement of the center if the Gaussian-like feature is not symmetric
        averageLines (int): number of scanning lines to add to the standard for smoother tracking
        invert (bool): to use for tracking dark features, instead of bright

    Returns:
        trackedFluorescence (numpy.array): array with the location of the tracked feature (in pixels)

    """
    if not lineEdges:
        if averageLines is 0:
            profile = pd.DataFrame({"x": np.arange(len(fluoData[:, startLine])), "y": fluoData[:, startLine]})
        else:
            profile = pd.DataFrame({"x": np.arange(len(fluoData[:, startLine])),
                                    "y": np.mean(fluoData[:, startLine - averageLines : startLine + averageLines + 1],
                                                 axis=1)})
            print(len(fluoData[0, startLine - averageLines : startLine + averageLines + 1]))
        if invert:
            profile["y"] = np.max(profile["y"]) - profile["y"]
        selection = data_selection(profile, columnX="x", columnY="y")
        lineEdges = (selection["x"].min(), selection["x"].max())

    if not endLine:
        endLine = len(fluoData[0])
    guessSigma = 1
    guessAmplitude = 1000

    def gauss(x, A, mean, sigma):
        """
        Simple Gaussian function
        """
        return A*np.exp(-(x-mean)**2/(2*sigma**2))

    guessCenter = (lineEdges[0] + lineEdges[1]) // 2 + shiftCenter
    halfWidth = (lineEdges[1] - lineEdges[0]) // 2
    halfWidthLeft = halfWidth + shiftCenter
    halfWidthRight = halfWidth - shiftCenter
    trackingFits = []

    for line in np.arange(startLine + averageLines, endLine - averageLines):

        x = np.arange(guessCenter - halfWidthLeft, guessCenter + halfWidthRight)
        if averageLines is 0:
            signal = fluoData[guessCenter - halfWidthLeft:guessCenter + halfWidthRight, line]
        else:
            signal = np.mean(fluoData[:, line - averageLines: line + averageLines], axis=1)[guessCenter -
                                                                                             halfWidthLeft:guessCenter + halfWidthRight]
        #plt.plot(signal)
        try:
            if invert:
                y = max(signal) - signal
            else:
                y = signal - min(signal)

            fittedParams, fittedCov = curve_fit(gauss, x, y, p0=[guessAmplitude, guessCenter, guessSigma])
            trackingFits.append(fittedParams[1])
            guessCenter = int(trackingFits[-1])
            guessSigma = fittedParams[2]
            guessAmplitude = fittedParams[0]
        except:
            trackingFits.append(np.nan)
            guessCenter = int(np.mean(pd.Series(trackingFits).dropna()))

    trackedFluorescence = pd.Series(trackingFits, index=(np.arange(len(trackingFits)) + startLine) * fluoTimeRes)
    #plt.show()
    return trackedFluorescence


def align_fluorescence_data(data, fluoData, startLine, endLine=None, fluoTimeRes=0.085356, maximumOffset=4,
                            defaultOffset=2, channel="trapSepX", **kwargs):
    """
    Automatically align the fluorescence and force data by tracking a bead, resampling high-res data and minimizing
    residues between signals at different time shifts
    Args:
        data (pandas.DataFrame): high-res data
        fluoData (numpy.array): array containing the fluorescent image of a given color
        startLine (int): line to start the tracking for alignment
        endLine (int): line to end the tracking for alignment
        fluoTimeRes (float): time resolution of the fluorescence data
        maximumOffset (float): the maximum shift to look for (in seconds)
        **kwargs: keyword arguments for the tracking function

    Returns:
        fluoOffset (float): time offset in seconds between fluorescence and force signal
    """
    trackedFluorescence = track_gaussian_feature(fluoData, startLine, endLine, fluoTimeRes=fluoTimeRes, **kwargs)

    lowData = data[["time", "trapSepX", "surfaceSepX"]].copy()
    lowData.index = lowData["time"]
    fluo = pd.DataFrame({"fluo": trackedFluorescence.values * 80}, index=trackedFluorescence.index)
    fluo["time"] = fluo.index
    print(channel)
    lowData = lowData.merge(fluo, on="time", how="outer", sort=True)
    lowData.iloc[:, 0:-1] = lowData.iloc[:, 0:-1].interpolate()

    newData = lowData.loc[lowData["fluo"].dropna().index, :]
    newData.index = np.arange(len(newData))
    shifts = []
    residues = []
    average = (newData[channel] - newData["fluo"]).mean()
    for shift in range(-int((maximumOffset + defaultOffset) // fluoTimeRes), int((maximumOffset) // fluoTimeRes)):
        shifts.append(shift)

        difference = (newData[channel].shift(shift) - newData["fluo"]-average).dropna()**2
        residues.append(difference.sum() / difference.count())

    def gauss(x, A, mean, sigma, B):
        """
        Simple Gaussian function
        """
        return A*np.exp(-(x-mean)**2/(2*sigma**2)) + B
    #plt.plot(shifts, max(residues) - residues, '.')
    fitParams, fitCov = curve_fit(gauss, shifts, max(residues) - residues)
    #plt.plot(shifts, gauss(shifts, fitParams[0], fitParams[1], fitParams[2], fitParams[3]))
    fluoOffsetRough = -shifts[residues.index(min(residues))] * fluoTimeRes
    fluoOffset = -fitParams[1] * fluoTimeRes
    plt.plot(newData["time"] + fluoOffset, newData["fluo"] + average, '.')
    plt.plot(newData["time"], newData[channel])
    #plt.show()
    print(fluoOffsetRough, fluoOffset)
    return fluoOffset



def plot_contour_fluorescence(data, fluoData, color, signalRegion, topChannel="proteinLc",
                              backgroundRegion=None, plotSignalRegion=False, fluoTimeRes=0.08596, fluoOffset=0,
                              savgolWindow=13, savgolOrder=1, timeUnits="s", normedFluo=False, vmaxKymo=1000,
                              cmap="inferno", proteinLength=120, lw=1, fontsize=15):
    """
    Plot the PSD data (either force or protein contour length) together with the fluorescence data

    Args:
        data (pandas.DataFrame): data from the PSD (force, bead displacement...)
        fluoData (dict): data containing the fluorescence images
        color (str): either "532nm" or "638nm"
        signalRegion (tuple): top and bottom pixels for calculating the signal region
        topChannel (str): which channel to plot on top, normally contour length ("proteinLc") or force ("force")
        backgroundRegion (tuple): region to calculate background. Should be same length as signalRegion. Default None
        plotSignalRegion (bool): if True, the moving signal region is plotted along with the kymograph
        fluoTimeRes (float): the time resolution of fluorescence data
        fluoOffset (float): seconds of delay between the fluorescenece data and the PSD data
        savgolWindow (int): window for the Savgol filter for the signal
        savgolOrder (int): order for the Savgol filter for the signal
        timeUnits (str): either "s" (seconds) or "lines" (scanning lines)
        normedFluo (bool): if True, normalize the fluorescence signal to the background
        vmaxKymo (int): maximum value for the kymograph display
        cmap (str): colormap for the kymograph
        proteinLength (float): length of the protein
        lw (float): linewidth for plots
        fontsize (float): fontsize for plots

    """

    sav_fil = lambda x: savgol_filter(x, savgolWindow, savgolOrder)
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(6,8))
    if timeUnits is "s":
        unit = 1
        ax3.set_xlabel("Time (s)")
    elif timeUnits is "lines":
        unit = fluoTimeRes
        ax3.set_xlabel("Scanning line")
    newData = data.copy()
    newData.index = pd.to_timedelta(newData.index, unit='s')
    resDataHigh = newData.resample(rule='5ms').mean()
    resData = newData.resample(rule='50ms').mean()

    if topChannel is "proteinLc":
        ax1.plot(resDataHigh["time"] / unit, resDataHigh["proteinLc"], alpha=0.2, zorder=-1, color="gray", lw=lw)##acdb99")
        ax1.plot(resData["time"] / unit, resData["proteinLc"], color="#bf311a", lw=lw)
        ax1.set_ylim(-50, proteinLength + 50)
        ax1.set_ylabel("Protein Lc (nm)", fontsize=fontsize)
    elif topChannel is "force":
        ax1.plot(resDataHigh["time"] / unit, resDataHigh["forceX"], alpha=0.2, zorder=-1, color="gray", lw=lw)##acdb99")
        ax1.plot(resData["time"] / unit, resData["forceX"], color="#bf311a", lw=lw)
        ax1.set_ylabel("Force (pN)", fontsize=fontsize)


    extensionAverage = int(resDataHigh["surfaceSepX"].mean() / 80)

    resDataHigh["beadMovement"] = (resDataHigh["surfaceSepX"].interpolate().astype(int) / 80 - extensionAverage) //2

    time = np.arange(fluoData[color].shape[1]) * fluoTimeRes + fluoOffset

    signal = np.sum(fluoData[color][signalRegion[0]: signalRegion[1], :], axis=0)
    laserData = resDataHigh.copy()
    laserData.index = laserData["time"]
    fluo = pd.DataFrame({"signal": signal}, index=time)
    fluo["time"] = fluo.index
    laserData = laserData.merge(fluo, on="time", how="outer", sort=True)

    laserData.iloc[:, 0:-1] = laserData.iloc[:, 0:-1].interpolate()

    last = laserData.loc[laserData["signal"].dropna().index, :].dropna()

    last.index = np.arange(len(last))
    #ax2.plot(time / unit, last["beadMovement"], '.')
    movingSignal = []
    #print(last.dropna())
    for scanLine in np.arange(len(time)):
        #print(last.loc[scanLine, "beadMovement"])
        try:
            movingSignal.append(np.sum(fluoData[color][signalRegion[0]+int(last.loc[scanLine, "beadMovement"]):
            signalRegion[1]+int(last.loc[scanLine, "beadMovement"]), scanLine], axis=0))
        except:
            movingSignal.append(np.nan)

    if not backgroundRegion:
        backgroundRegion = (5, 5 + signalRegion[1] - signalRegion[0])
    background = np.sum(fluoData[color][backgroundRegion[0]: backgroundRegion[1], :], axis=0)
    if normedFluo:
        signal /= np.mean(background)
        background /= np.mean(background)

    #print(movingSignal)
    ax2.plot(time / unit, movingSignal, color="#005b96", alpha=0.3, lw=lw)
    ax2.plot(time / unit, background, zorder=-1, color="gray", alpha=0.3, lw=lw)
    ax2.plot(time / unit, sav_fil(movingSignal), color="#005b96", lw=lw)
    ax2.plot(time / unit, sav_fil(background), color="gray", lw=lw)

    #ax2.set_ylabel("I (a.u.)", fontsize=fontsize)
    #ax3.set_xlabel("Time (" + timeUnits + ")", fontsize=fontsize)

    for ax in [ax1, ax2, ax3]:
        ax.tick_params(direction="in", length=8, width=lw, pad=10, top='on', bottom='on')
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(lw)
    kymoExtent = np.array([fluoOffset / unit, (fluoTimeRes*fluoData[color].shape[1]+fluoOffset) / unit,
                                               fluoData[color].shape[0], 0])
    ax3.imshow(fluoData[color], extent=kymoExtent, aspect='auto', vmax=vmaxKymo,  cmap=cmap)

    if plotSignalRegion:
        ax3.plot(last["time"]/unit, signalRegion[0] + last["beadMovement"], color="yellow", lw=0.5)
        ax3.plot(last["time"]/unit, signalRegion[1] + last["beadMovement"], color="yellow", lw=0.5)

    ax3.set_ylabel("Pixels", fontsize=fontsize)
    ax3.set_xlim(time[0]/unit, time[-1]/unit)
    ax3.set_ylim(fluoData[color].shape[0], 0)
    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.95, wspace=None, hspace=0.1)
    plt.show()