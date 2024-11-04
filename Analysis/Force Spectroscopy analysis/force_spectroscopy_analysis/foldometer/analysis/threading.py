#!/usr/bin/env python
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.ndimage import label
from scipy.signal import savgol_filter, welch

from foldometer.physics.thermodynamics import thermal_energy
from foldometer.physics.utils import as_Kelvin
from foldometer.tools.misc import data_selection, column_indexer
import numpy as np
from matplotlib.widgets import SpanSelector
from scipy.stats import linregress


def calculate_rolling_slopes(data, proteinLength, dt=None, window=600, columnY="surfaceSepX", selectData=False,
                             ylim=None, center=True, residuals=False):
        if ylim is None and columnY is "proteinLc":
            ylim = (-50, proteinLength + 50)
        if selectData:
            data = data_selection(data, columnX="time", columnY=columnY, ylim=ylim)
            data = data[["time", columnY]]
            data.index = data["time"]

        data = data[["time", columnY]]
        data.index = data["time"]
        if not dt:
            dt = data["time"].diff().mean()

        def summing_window(data):
            return (data * np.arange(1, len(data)+1)).sum()
        #print(self.metadata["sampleFreq"])
        denominator = dt / 12 * (2 * window + window ** 3)
        slopes = (data[columnY].rolling(window, center=center).apply(summing_window) - window /
                  2 * data[columnY].rolling(window, center=center).mean() * (window + 1)) / denominator

        data["transSpeed"] = np.nan
        data.loc[slopes.index, "transSpeed"] = slopes
        intercepts = data[columnY].rolling(window, center=center).mean() - \
                     data["transSpeed"] * data["time"].rolling(window, center=center).mean()
        residuals = data[columnY] - (data["transSpeed"] * data["time"] + intercepts)
        #residualsSTD =
        return data["transSpeed"]


def pairwise_calculation(data, proteinLength, selectData=False, channel="proteinLc", mode="nm", order=5, window=21,
                         bins=None, plot=True):

    if channel is "surfaceSepX":
        mode = "nm"

    if selectData:
        if channel is "proteinLc":
            ylim = (-50, 400)
        else:
            ylim = None
        data = data_selection(data, columnX="time", columnY=channel, ylim=ylim)
        data = data.loc[:, ["time", channel]]
        data.index = data["time"]
    sav_fil = lambda x: savgol_filter(x, window, order)
    pairs = []
    if plot:
        fig = plt.figure(figsize=(8,6))
        ax1 = plt.subplot(121)
        ax2 = plt.subplot(222)
        ax3 = plt.subplot(224)

    if channel is "proteinLc":
        if mode is "aa":
            stepsData = pd.Series(sav_fil((proteinLength - data[channel]) / 0.35))
            if plot:
                ax1.plot(data["time"], (proteinLength - data[channel]) / 0.35, color="gray", alpha=0.2)
                ax1.plot(data["time"], sav_fil((proteinLength - data[channel]) / 0.35), color="#bf311a")
                print(min(data[channel])/0.35, max(data[channel])/0.35)
                major_ticks = np.arange(int(min(stepsData)), int(max(stepsData)), 5)
                ax1.set_yticks(major_ticks)
        else:
            stepsData = pd.Series(sav_fil((proteinLength - data[channel])))
            if plot:
                ax1.plot(data["time"], (proteinLength - data[channel]), color="gray", alpha=0.2)
                ax1.plot(data["time"], sav_fil((proteinLength - data[channel])), color="#bf311a")

    else:
        stepsData = pd.Series(sav_fil(data[channel]))
        if plot:
            ax1.plot(data["time"], data[channel], color="gray", alpha=0.2)
            ax1.plot(data["time"], sav_fil(data[channel]), color="#bf311a")
            major_ticks = np.arange(data[channel].min(), data[channel].max(), 10)

            ax1.set_yticks(major_ticks)

    for shift in np.arange(len(data[channel]))[1:]:
        pairs.extend(stepsData.diff(shift).dropna().values)
    if not bins:
        bins = len(pairs) // 500
    hist, binEdges = np.histogram(pairs, bins=bins, density=True)

    blockLength = len(hist) // 2
    fRaw, psdRaw = welch(hist, nperseg=blockLength, window="hann", fs=1/np.diff(binEdges).mean(), noverlap=8)
    print("Time interval: ", data["time"].min(), data["time"].max())
    if plot:
        ax1.set_xlabel("Time (s)")
        ax2.hist(pairs, bins=bins, normed=True)
        ax2.tick_params(labelleft="off")
        ax3.plot(fRaw, psdRaw, '-')
        ax3.set_ylim(0, max(psdRaw))
        ax3.set_xlim(0,0.1)
        ax3.tick_params(labelleft="off")
        ax1.set_ylabel("Translocated length (" + mode + ")")
        ax2.set_xlabel("Pairwise length (" + mode + ")")
        ax3.set_xlabel("Spatial frequency (" + mode + "$^{-1}$)")
        ax1.grid(b=True, axis="y")
        if mode =="aa":
            ax3.set_xlim(0,0.2)
        elif mode == "nm":
            ax3.set_xlim(0,2)
        plt.subplots_adjust(hspace=0.5, wspace=0.1, top=0.95, right=0.95, left=0.15, bottom=0.15)
        plt.show()
        #return (hist, (fRaw, psdRaw))
    return fRaw, psdRaw


def manual_translocation_speed_fit(data, smooth=True, columnX="time", columnY="forceX", window=101, order=1):
    """
    Opens a plot and allows to select a portion of the data by dragging the mouse. Pressing any key will confirm the
    selection and close the window

    Args:
        data (pandas.DataFrame): DataFrame with the data source to select from
        columnX (str): column to display in the X axis
        columnY (str): column to display in the Y axis

    Returns:
        subData (pandas.DataFrame): selected subset of original data
    """

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)

    x = data[columnX]
    y = data[columnY]
    ax.plot(x, y, '-')
    if smooth:
        window =101
        order = 1
        sav_fil = lambda x: savgol_filter(x, window, order)
        ax.plot(sav_fil(data[columnX]), sav_fil(data[columnY]))
    if columnY is "proteinLc":
        ax.set_ylim(-250, 500)
    fig.suptitle("Press any key to confirm selection and close this window", fontsize=14, fontweight='bold')
    ax.set_title("Select data by dragging the mouse")
    labels = ("slopes", "intercepts", "r_values", "p_values", "std_err")
    fits = {key: [] for key in labels}
    fits["force"] = []
    fits["startTime"] = []
    fits["endTime"] = []
    if "proteinLc" in data.columns:
        fits["startLc"] = []
        fits["endLc"] = []
    initCoords = []
    endCoords = []
    def onselect(xmin, xmax):
        """
        Sub-function to assign the minimum and maximum of the range and plot the selection"

        Args:
            xmin (int): minimum index of selection
            xmax (int): maximum index of selection
        """
        global indmin
        global indmax
        global txt
        global span

        #txt = ax.text(0,0, ' ', transform=ax.transAxes)
        try:
            txt.remove()
            span.remove()
            fig.canvas.draw()
        except:
            pass

        indmin, indmax = np.searchsorted(x, (xmin, xmax))
        indmax = min(len(x) - 1, indmax)
        span = ax.axvspan(x.values[int(indmin)], x.values[int(indmax)], facecolor="red", alpha=0.5)
        fit = linregress(x.values[int(indmin): int(indmax)],y.values[int(indmin): int(indmax)])
        for label, fitOut in zip(labels, fit):
            fits[label].append(fitOut)
        fits["force"].append(np.mean(data["forceX"].values[int(indmin): int(indmax)]))
        fits["startTime"].append(data["time"].values[int(indmin)])
        fits["endTime"].append(data["time"].values[int(indmax)])
        if "proteinLc" in data.columns:
            fits["startLc"].append(data["proteinLc"].values[int(indmin)])
            fits["endLc"].append(data["proteinLc"].values[int(indmax)])
        initCoords.append(int(indmin))
        endCoords.append(int(indmax))
        txt = ax.text(0,0.9, fit[0], transform=ax.transAxes, color="red")
        fig.canvas.draw()


    # set useblit True on gtkagg for enhanced performance
    span = SpanSelector(ax, onselect, 'horizontal', useblit=True, rectprops=dict(alpha=0.5, facecolor='red'))
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized()
    plt.show()

    while True:
        wait = plt.waitforbuttonpress(timeout=20000)
        if wait:
            plt.close('all')
            break
        else:
            pass

    fitsDF = pd.DataFrame(fits)
    fitsDF["start"] = initCoords
    fitsDF["end"] = endCoords

    fitsDF["slopes"] *= -1
    print(fitsDF[["slopes", "force"]])
    idCol = column_indexer(data)
    mask = (data[columnX] >= data.iloc[indmin, idCol[columnX]]) & \
           (data[columnX] <= data.iloc[indmax, idCol[columnX]])

    return fitsDF