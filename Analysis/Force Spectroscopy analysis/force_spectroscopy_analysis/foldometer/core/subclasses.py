#!/usr/bin/env python
# -*- coding: utf-8 -*-

from foldometer.core.main import Folding
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.ndimage import label, binary_closing
from scipy.signal import savgol_filter
from foldometer.physics.thermodynamics import thermal_energy
from foldometer.physics.utils import as_Kelvin
from foldometer.analysis.wlc_curve_fit import protein_contour_length
from foldometer.tools.misc import data_selection
from foldometer.analysis.threading import calculate_rolling_slopes
from matplotlib.widgets import SpanSelector
from matplotlib.widgets import Cursor
from matplotlib.widgets import Button
import numpy as np
import os
from scipy.stats import linregress
from scipy.signal import find_peaks_cwt
import gc

class Thread(Folding):
    """
    Class for analyzing threading experiments by AAA+ protein like ClpB, Hsp104 or ClpG
    """
    def __init__(self, filePath, protein="unknown", condition={}, setup="DT", fluorescenceFilePath=None):
        """
        Main class of the Foldometer package, containing all necessary attributes and automated methods to analyse the
        pulling experiments data for any of the setups at AMOLF (and easily expandable to more setups).

        Args:
            filePath (str): Path to the file containing the data
            protein (str): Protein used for the experiment. In the newer Foldometer version this is in metadata
            condition (dict): Dictionary containing different conditions like chaperones, buffer, temperature...
            setup (str): which setup was used for the data. Either "ST", "DT" or "LT"
        """
        Folding.__init__(self, filePath, protein, condition, setup, fluorescenceFilePath)

    def plot_contour_length(self, timeRes=50, threaded=True, normalized=False, lengthOffset=0, ax=None, lw=1,
                            show=True):
        if not ax:
            fig, ax = plt.subplots(1,1)
        plotData = self.data.copy()
        plotData.index = pd.to_timedelta(plotData.index, unit='s')
        plotData["proteinLc"] -= lengthOffset
        if threaded:
            plotData["proteinLc"] = self.proteinLength - plotData["proteinLc"]
        ax.set_ylim(-50, self.proteinLength + 50)
        if normalized:
            plotData["proteinLc"] /= self.proteinLength
            ax.set_ylim(-0.1, 1.1)

        resLcProtein = plotData["proteinLc"].resample(rule=str(timeRes)+'ms').mean()
        ax.plot(resLcProtein.index.total_seconds(), resLcProtein, color="#bf311a", lw=lw)
        ax.plot(plotData.index.total_seconds(), plotData["proteinLc"], alpha=0.2, zorder=-1, color="gray", lw=lw)

        ax.set_ylabel("$L_t$ (nm)")
        ax.set_xlabel("Time (s)")
        if show:
            plt.subplots_adjust(bottom=0.15, left=0.15, top=0.95, right=0.95)
            plt.show()
        del plotData
        gc.collect()


    def translocation_speed(self, window=600, columnY="surfaceSepX", selectData=False, ylim=None, center=True):
        """
        Calculate the translocation speed of the threading events
        Args:
            window (int): number of data points to calculate the slope from
            columnY (str): for translocation speed, either surfaceSepX or proteinLc
            selectData (bool): if True, select a subset of the data
            ylim (tuple): limits for plotting and selecting data. If None, automatically adapt to protein total length
            center (bool): if True, the window is set in the center of the step
        """
        dt = 1 / self.metadata["sampleFreq"]
        self.data["transSpeed"] = np.nan
        self.data["transSpeed"] = calculate_rolling_slopes(self.data, self.proteinLength, dt, window=window,
                                                           columnY=columnY, selectData=selectData, ylim=ylim,
                                                           center=center)

    def identify_threading_events(self, window=101, order=1, slopeThreshold=10, closingLength=100,
                                  closingLengthSmall=10):
        """
        Automatically identify individual translocation runs looking at the time derivative of the contour length

        Args:
            window (int): number of points to use for the Savitzky-Golay filter
            order (int): order of the polynomial to use for the Savitzky-Golay filter
            slopeThreshold (float): from which slope to consider as a running event
            closingLength (int): the length used for the dilation in order to close artifact gaps
            closingLengthSmall (int): the length used to remove artifact gaps from first dilation

        Returns:

        """
        savgol_fil = lambda x: savgol_filter(x, window, order, deriv=1)
        plt.plot(self.data["time"], self.data["proteinLc"], color="gray", alpha=0.2)
        speeds = savgol_fil(self.data["proteinLc"])*(-self.metadata["sampleFreq"])
        self.data["threading"] = label(binary_closing(
                ~binary_closing(speeds < slopeThreshold, structure=closingLength*[1]),
                structure=closingLengthSmall*[1]))[0]

        for lab in sorted(self.data["threading"].unique())[1:]:
            mask = self.data["threading"] == lab
            plt.plot(self.data.loc[mask, "time"], self.data.loc[mask, "proteinLc"])
        plt.show()

    def fit_translocation_runs(self, window=125, order=1, speedThreshold=120, closingLength=100, fitQualityLimit=0.75,
        plot=False, channel="proteinLc"):
        """
        Identify subruns and fit each of them using a linear regression to calculate the average run speed

        Args:
            window (int): number of points to use for the Savitzky-Golay filter
            order (int): order of the polynomial to use for the Savitzky-Golay filter
            speedThreshold (float): from which slope to consider as a running event
            closingLength (int): the length used for the dilation in order to close artifact gaps
            fitQualityLimit (float: minimum p value required to consider the fit
            plot (bool): whether or not to plot the fits
            channel (str): either "proteinLc" or "trapSepX"

        Returns:
            fits (pandas.DataFrame): linear fits of the speed of individual runs
        """
        variables = {"slope": [], "errorSlope": [], "intercept": [], "startTime": [], "endTime": [], "startForce": [],
                     "endForce": [], "rValue": [], "event": []}

        if plot:
            plt.plot(self.data["time"], self.data[channel], color="gray", alpha=0.3)
        for threadingEvent in sorted(self.data["threading"].unique())[1:]:
            eventMask = (self.data["threading"] == threadingEvent)
            _data = self.data.loc[eventMask, ["time", channel, "forceX"]].copy()
            if np.count_nonzero(eventMask) > window:
                if plot:
                    pass

                savgol_fil2 = lambda x: savgol_filter(x, window, order, deriv=1, mode="interp")
                bins = pd.cut(savgol_fil2(_data[channel])*(-self.metadata["sampleFreq"]),
                              [-np.inf, 20, speedThreshold, np.inf], labels=[0, 1, 2])
                slopesMask = []
                for _bin in bins.unique():
                    closing = binary_closing(bins == _bin, structure=[1]*closingLength)
                    bins[closing] = _bin
                _data["region"] = np.zeros_like(bins)
                for _bin in [1.0, 2.0]:
                    newLabels = label(bins==_bin)[0]
                    newLabels[bins==_bin] = newLabels[bins==_bin] + max(_data["region"])
                    _data["region"] += newLabels
                    try:
                        newEvent = max(variables["event"]) + 1
                    except:
                        newEvent = 0
                for region in sorted(_data["region"].unique())[:]:
                    regionMask = _data["region"] == region
                    if plot:
                        #pass
                        plt.plot(_data.loc[regionMask, "time"], _data.loc[regionMask, channel])
                    x = _data.loc[regionMask, "time"]
                    fit = linregress(x, _data.loc[regionMask, channel])

                    if -fit[0] > 0 and -fit[0] < 300 and -fit[2]> fitQualityLimit:
                        variables["slope"].append(-fit[0])
                        variables["errorSlope"].append(fit[4])
                        variables["intercept"].append(-fit[1])
                        variables["startTime"].append(_data.loc[regionMask, "time"].min())
                        variables["endTime"].append(_data.loc[regionMask, "time"].max())
                        variables["startForce"].append(_data.loc[regionMask, "forceX"].min())
                        variables["endForce"].append(_data.loc[regionMask, "forceX"].max())
                        variables["rValue"].append(fit[2])
                        variables["event"].append(newEvent)
                        if plot:
                            plt.plot(x, x*fit[0] + fit[1], color="red", lw=3)
            else:
                try:
                    newEvent = max(variables["event"]) + 1
                except:
                    newEvent = 0
                if plot:
                    plt.plot(_data["time"], _data[channel])
                x = _data["time"]
                fit = linregress(x, _data[channel])

                if -fit[0] > 0 and -fit[0] < 300 and -fit[2]> fitQualityLimit:

                    variables["slope"].append(-fit[0])
                    variables["errorSlope"].append(fit[4])
                    variables["intercept"].append(-fit[1])
                    variables["startTime"].append(_data["time"].min())
                    variables["endTime"].append(_data["time"].max())
                    variables["startForce"].append(_data["forceX"].min())
                    variables["endForce"].append(_data["forceX"].max())
                    variables["rValue"].append(fit[2])
                    variables["event"].append(newEvent)
                    if plot:
                        plt.plot(x, x*fit[0] + fit[1], color="red", lw=3)

        fits = pd.DataFrame(variables).sort_values("startTime").reset_index(drop=True)
        fits["subEvent"] = 0
        for event in fits["event"].unique():
            fits.loc[fits["event"] == event, "subEvent"] = np.arange(np.count_nonzero(fits["event"] == event))
        fits["startLength"] = -(fits["startTime"] * fits["slope"] + fits["intercept"])
        fits["endLength"] = -(fits["endTime"] * fits["slope"] + fits["intercept"])
        fits["fileName"] = self.filePath
        return fits

    def identify_threading_events_legacy(self, window=5, slopeThreshold=10):
        """
        Look for breaking events to identify individual runs. Similar code as to identify unfolding events

        Args:
            window (int): number of points to calculate the difference
            slopeThreshold (int): threshold to consider data as a breaking event
        """
        if "transSpeed" not in self.data.columns:
            self.translocation_speed()
        _dataEvents = self.data.copy()
        _dataEvents.index = pd.to_timedelta(_dataEvents.index, unit='s')
        _dataEvents = _dataEvents.resample(rule='20ms').mean()

        def mean_diff(dataSet, n):
            """
            Subfunction to calculate the rolling mean of the difference

            Args:
                dataSet (numpy.array): array with the data subset (window) to calculate the moving mean of the difference
                n (integer): number of datapoints to calculate the diff from
            Returns:
                mean of the difference of the subset
            """
            return np.mean(np.diff(dataSet, n))

        meanChange = _dataEvents["proteinLc"].rolling(window, center=True).apply(func=mean_diff, args=(1,))
        meanChange[meanChange < slopeThreshold] = 0
        _dataEvents["runBreak"] = False
        _dataEvents.loc[meanChange > 0, "runBreak"] = True
        _dataEvents["runBreak"] = label(_dataEvents["runBreak"])[0]
        _dataEvents.loc[_dataEvents["runBreak"] == 0, "runBreak"] = np.nan
        if "runBreak" in self.data.columns:
            self.data.drop("runBreak", axis=1, inplace=True)
        self.data = self.data.merge(_dataEvents.loc[_dataEvents["runBreak"] != np.nan, ["time","runBreak"]], on="time",
                             how="outer", sort=True)
        firstEvent = self.data.loc[ : , "transSpeed"].first_valid_index()
        self.data.loc[0: firstEvent, "runBreak"] = -1
        self.data["runBreak"] = self.data["runBreak"].fillna(method="bfill").fillna(method="ffill")
        self.data.dropna(subset=["trapSepX"], inplace=True)
        self.data["threading"] = self.data["runBreak"].values
        self.data.index = self.data["time"]
        del _dataEvents
        gc.collect()

    def export_threading_events(self):
        """
        Save a file containing the event information: time, translocated length, translcoation speed, event number,
        force and file name. Allows for several exported files per data file.
        """
        if self.setup is "DT":
            fileName = self.filePath[:-4]
            fileNumber = self.filePath[-7:-4]
        elif self.setup is "CT":
            fileName = self.filePath[:-5]
            fileNumber = self.filePath[-12:-5]
        counter = 1
        for file in os.listdir(os.path.split(fileName)[0]):
            if fileNumber in file and "threading_events" in file:
                counter += 1
            filePath = fileName + "_threading_events_" + str(counter).zfill(2) + ".txt"
        columns = ["time", "proteinLc", "threading", "file", "forceX", "transSpeed"]
        try:
                self.data["file"] = fileNumber
                self.data.loc[:,columns].to_csv(filePath)
                self.data.drop("file", axis=1, inplace=True)
        except:
            print("Calculate contour length and identify threading runs first!")

class LifeTimeMeasurement(Folding):
    def __init__(self, filePath, protein="unknown", condition={}, setup="DT", fluorescenceFilePath=None):
        """
        Class initialization
        Args:
            filePath (str): path to the file containing the data
            protein (str): protein used for the experiment. In the newer Foldometer version this is in metadata
            condition (dict): dictionary containing different conditions like chaperones, buffer, temperature...
        """
        Folding.__init__(self, filePath, protein, condition, setup, fluorescenceFilePath)

    def select_noise_region(self):
        """

        Returns:

        """
        noiseData = data_selection(self.data)
        noisePolyFitForce = np.polyfit(noiseData["surfaceSepX"], noiseData["forceX"], deg=5)
        polynomialForce = np.polyval(noisePolyFitForce, self.data["surfaceSepX"])
        self.data.loc[:, "forceX"] -= polynomialForce

    def identify_binding_events(self, timeRes=20, minimumRegionLength=5, minForce=1, maxForce=50,
                                timeLimit=60, plot=True):


        self.data.index = pd.to_timedelta(self.data.index, unit='s')
        stationaryMask = self.data["region"] == "stationary"
        resData = self.data.loc[stationaryMask].resample(rule=str(timeRes) + 'ms').mean()
        forceMask = (resData["forceX"] > minForce) & (resData["forceX"] < maxForce)
        labeled, num = label(forceMask)
        resData["region"] = labeled
        print(num)
        if num is 0:
            num=1
        if plot:
            plt.plot(self.data["time"], self.data["forceX"], color="gray", alpha=0.5)
        resData[resData["region"] ==0] = np.nan
        lifes = []
        forces = []

        with sns.hls_palette(num, l=.3, s=.8):
            for region in list(resData["region"].unique()):
                regionMask = (resData["region"] == region)
                #print(len(resData.loc[regionMask, "pullingCycle"]))
                if len(resData.loc[regionMask, "pullingCycle"]) < minimumRegionLength:
                    #make nan the statioanry regions that are too sort (for fill methods)
                    resData.loc[regionMask, "region"] = np.nan
                    #fill first forward half of the region and the other half backwards
                else:
                    lifes.append(len(resData.loc[regionMask, "region"])*timeRes/1000)
                    forces.append(resData.loc[regionMask, "forceX"].mean())
                    if plot:
                        plt.plot(resData.loc[regionMask, "time"], resData.loc[regionMask, "forceX"], )
        resData.dropna(inplace=True)
        lifeTimes = np.asarray(lifes)
        self.lifeTimes = pd.DataFrame({"lifeTime": lifeTimes, "force": forces, "file": self.filePath})
        self.lifeTimes["break"] = True
        self.lifeTimes.loc[self.lifeTimes["lifeTime"] > timeLimit, "break"] = False

        if plot:
            #plt.plot(self.data["time"], self.data["forceX"])
            #plt.plot(resData["time"], resData["forceX"], '.')
            plt.show()

    def append_lifetimes(self, file):
        try:
            with open('lifetimes.csv', 'a') as f:
                self.lifeTimes.to_csv(f, header=False)
        except:
            print("error!")


class DisaggMeasurement(Folding):
    def __init__(self, filePath, protein, condition={}):
        """
        Class initialization
        Args:
            filePath (str): path to the file containing the data
            protein (str): protein used for the experiment. In the newer Foldometer version this is in metadata
            condition (dict): dictionary containing different conditions like chaperones, buffer, temperature...
        """
        Folding.__init__(self, filePath, protein, condition)

    def chaperoneTime(self):
        """
        The SpanSelector is a mouse widget to select a xmin/xmax range and plot the
        detail view of the selected region in the lower axes
        """

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, axisbg='#FFFFCC')

        x = self.data["time"]
        y = self.data["forceX"]
        ax.plot(x, y, '-')

        #def close(event):
        #    plt.close('all')

        #axClose = plt.axes([0.4, 0.92, 0.2, 0.05])
        #closeButton = Button(axClose, 'Ok, close!')
        #closeButton.on_clicked(close)

        def on_click(event):
            # get the x and y coords, flip y from top to bottom
            x, y = event.x, event.y
            if event.button == 1:
                if event.inaxes is not None:
                    print('data coords %f %f' % (event.xdata, event.ydata))
                    self.chaperoneIndex = event.xdata

        plt.connect('button_press_event', on_click)

        # set useblit True on gtkagg for enhanced performance
        cursor = Cursor(ax, useblit=True, color='red', horizOn=False, linewidth=2)

        plt.show()
        self.data["chaperone"] = False
        self.data.loc[self.chaperoneIndex:, "chaperone"] = True
        print(self.chaperoneIndex)



