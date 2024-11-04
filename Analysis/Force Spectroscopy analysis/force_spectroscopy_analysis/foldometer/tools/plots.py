#!/usr/bin/env python
# -*- coding: utf-8 -*-

import seaborn as sns
import matplotlib.pyplot as plt
import time
import os
import pandas as pd
import numpy as np
import string
from foldometer.tools.maths import cross_correlation
from foldometer.tools.misc import format_pulling_cycles
#from foldometer.physics.polymer import wlc, wlc_series
from foldometer.analysis.wlc_curve_fit import wlc, wlc_series, wlc_from_fit, wlc_series_accurate
from foldometer.tools.misc import column_indexer, resample_data
from scipy.signal import savgol_filter
from scipy.signal import medfilt
import pkg_resources

sns.set_style("whitegrid")
import scipy.spatial as spatial
import copy

sMBP_1300_DIG_BIOTIN_RULER = pkg_resources.resource_filename('foldometer', 'data/sMBP_ruler_1300_dig_biotin.txt')


def event_info_display(meanSeparation, unfoldingData):
    """
    What to display in the message when hovering the mouse over an event

    Args:
        meanSeparation (float): the nearest middle point of an event to the mouse position
        unfoldingData (pandas.DataFrame); dataframe with the unfolding events information

    Returns:
        string with what to show in the pop up box of the event (str)
    """

    eventMask = unfoldingData["meanSeparation"] == meanSeparation
    eventID = unfoldingData.loc[eventMask, "force"].index[0]
    force = abs(unfoldingData.loc[eventMask, "force"].values[0])
    extensionChange = abs(unfoldingData.loc[eventMask, "extensionChange"].values[0])

    return 'event ID:{ID}\n$F_u$: {fu:0.2f} pN\n$\Delta L$: {deltaL:0.2f} nm'.format(ID=eventID, fu=force,
                                                                                     deltaL=extensionChange)


class FollowDotCursor(object):
    """Display the event characteristics of the nearest event to the mouse.
    http://stackoverflow.com/a/4674445/190597 (Joe Kington)
    http://stackoverflow.com/a/13306887/190597 (unutbu)
    http://stackoverflow.com/a/15454427/190597 (unutbu)
    """
    def __init__(self, ax, x, y, data, unfoldingData, axis='x', pulls=None, tolerance=5, formatter=event_info_display,
                 offsets=(-20, 20)):
        """
        Class initiation

        Args:
            ax (matplotlib.pyplot.axes): axes where to plot the event info
            x (np.array): set of x points where events happen (mean separation of the event, (end-start)/2)
            y (np.array): set of y points where events happen (mean force of the event, (end-start)/2)
            data (pandas.DataFrame): DataFrame with the force and separation info
            unfoldingData (pandas.DataFrame): DataFrame with the unfolding events information
            pulls (tuple): tuple with the initial and the last pulling cycles to be plotted
            tolerance (int): how far from the set points the mouse can hover to display certain point
            formatter (function): the function that returns the text of the event information
            offsets (tuple): offsets of the displaying box with respect to the point of the event

        """

        try:
            x = np.asarray(x, dtype='float')
        except (TypeError, ValueError):
            x = np.asarray(mdates.date2num(x), dtype='float')
        y = np.asarray(y, dtype='float')
        mask = ~(np.isnan(x) | np.isnan(y))
        x = x[mask]
        y = y[mask]
        self.times = x
        self._points = np.column_stack((x, y))
        self.offsets = offsets
        y = y[np.abs(y-y.mean()) <= 3*y.std()]
        self.scale = x.ptp()
        self.scale = y.ptp() / self.scale if self.scale else 1
        self.tree = spatial.cKDTree(self.scaled(self._points))
        self.formatter = formatter
        self.tolerance = tolerance
        self.ax = ax
        self.fig = ax.figure
        self.ax.xaxis.set_label_position('top')
        self.unfoldingData = unfoldingData
        self.eventPlots = []
        self.filteredSeparations = []
        margin = 15
        iCol = {columnName: columnIndex for columnIndex, columnName in enumerate(data.columns)}
        newMargin = copy.deepcopy(margin)
        for startTime, endTime in zip(unfoldingData.startTime, unfoldingData.endTime):
            startIndex = data.index.get_loc(startTime)
            endIndex = data.index.get_loc(endTime)
            if data.loc[endTime, "pullingCycle"] in pulls:
                if startIndex < newMargin:
                    margin = data.index.get_loc(startTime) - 1
                else:
                    margin = newMargin
                self.filteredSeparations.append(unfoldingData.loc[unfoldingData["startTime"] == startTime,
                                                                  "meanSeparation"].values[0])
                slicer = slice(startIndex - margin, endIndex + margin)
                self.eventPlots.append(ax.plot(data.iloc[slicer, iCol["smoothSeparation"]],
                     data.iloc[slicer, iCol["smoothForce"]], color="#B9F408", linewidth = 8, alpha=0, zorder=-1,
                                               label="_"))

        self.annotation = self.setup_annotation()
        plt.connect('motion_notify_event', self)

    def scaled(self, points):
        points = np.asarray(points)
        return points * (self.scale, 1)

    def __call__(self, event):
        ax = self.ax
        unfoldingData = self.unfoldingData
        # event.inaxes is always the current axis. If you use twinx, ax could be
        # a different axis.
        if event.inaxes == ax:
            x, y = event.xdata, event.ydata
        elif event.inaxes is None:
            return
        else:
            inv = ax.transData.inverted()
            x, y = inv.transform([(event.x, event.y)]).ravel()
        annotation = self.annotation
        x, y = self.snap(x, y)
        annotation.xy = x, y
        annotation.set_text(self.formatter(x, unfoldingData))
        annotation.set_color("white")
        for plotIndex, t in enumerate(self.filteredSeparations):

            if t == x:

                self.eventPlots[plotIndex][0].set_alpha(1)
                self.eventPlots[plotIndex][0].set_zorder(50)
            else:
                self.eventPlots[plotIndex][0].set_alpha(0)
                self.eventPlots[plotIndex][0].set_zorder(-1)

        bbox = ax.viewLim
        event.canvas.draw()

    def setup_annotation(self):
        """Draw and hide the annotation box."""
        annotation = self.ax.annotate(
            '', xy=(0, 0), ha = 'right',
            xytext = self.offsets, textcoords = 'offset points', va = 'bottom',
            bbox = dict(
                boxstyle='round,pad=0.5', fc='#1a3245', alpha=0.75),
            arrowprops = dict(
                arrowstyle='->', connectionstyle='arc3,rad=0'), zorder=1000)
        return annotation

    def snap(self, x, y):
        """Return the value in self.tree closest to x, y."""
        dist, idx = self.tree.query(self.scaled((x, y)), k=1, p=1)
        try:
            return self._points[idx]
        except IndexError:
            # IndexError: index out of bounds
            return self._points[0]


def force_extension_curve(data, unfoldingData=None, wlcData=None, axis="x", forceChannel="force",
                          distanceChannel="surfaceSep", timeRes=50, pulls=None, wlcParameters=None,
                          palettes=("PuBu_d", "YlOrRd_d", "Greens_d"), regions=("pulling", "retracting", "stationary"),
                          alpha=0.7, legend=True, dataPoints=True, interactive=False, eventMargin=15, xLimits=None,
                          yLimits=None, displayUnfoldingEvents=False, rulers=None, show=True, figure=None, axes=111,
                          fontSize=20, lineWidth=5, **kwargs):
    """
    Function to plot force extension curves in an automated way, distinguishing between different pulling,
    retracting and stationary cycles and also with possibility of interactive display of identified unfolding events

    Args:
        data (pandas.DataFrame): DataFrame containing all the processed data, including the regions (check assign_regions() for more details)
        unfoldingData (pandas.DataFrame): DataFrame containing the unfolding events data
        wlcData (dictionary): dictionary containing the fits of the WLC model for each wlc region
        axis (str): direction to plot the curves, either "x" or "y"
        forceChannel (str): either "PSD1Force", "PSD2Force" or "force" (combined signal)
        distanceChannel (str): either "surfaceSep" (from PSDs) or "trackingSep" (from image tracking)
        timeRes (int): time window in milliseconds to resample (and therefore smooth) the data
        pulls (tuple): tuple with the initial and the last pulling cycles to be plotted
        palettes (tuple): tuple with the names of the palettes used for the regions
        regions (tuple): tuple with regions to be plotted
        alpha (float): for setting the transparency of the main line data
        legend (bool): if True, a legend will be drawn
        dataPoints (bool): if True, the scattered data is plotted in the background
        interactive (bool): whether to show interactive information of the unfolding events. Default: False
        eventMargin (int): number of data points to consider before and after the event region for correct display
        xLimits (tuple): set the x limits of the plot
        yLimits (tuple): set the y limits of the plot
        displayUnfoldingEvents (bool): if True, the events are shown in highlighted black
        rulers (bool): if True, show custom rulers (in development still)
        show (bool): if True, the plot will be shown
        figure (matplotlib.figure): if needed to save in an already created figure. Default None creates new figure
        axes (int): in which axes to plot
        fontSize (float): the fontsize of the labels
        lineWidth (float): linewidth of the main data plotting
    """

    axis = axis.upper()

    if "region" not in data.columns.values:
        raise ValueError("Regions have not been assigned, please first run assign_regions() function")

    if axis not in ["x", "X", "y", "Y"]:
        raise ValueError("Axis should be 'x', 'y'")

    if figure is None:
        figure = plt.figure(figsize=(12,11))


    axForceExtension = figure.add_subplot(axes)

    axForceExtension.set_ylabel("Force (pN)", fontsize=fontSize)
    axForceExtension.set_xlabel("Extension (nm)", fontsize=fontSize)
    #set the data and the smooth data according to the axis
    separation = distanceChannel + axis
    force = forceChannel + axis

    #set the number of pull cycles to be plotted
    pulls = format_pulling_cycles(pulls, data.iloc[0, column_indexer(data)["pullingCycle"]],
                                  len(data.pullingCycle.unique()))
    dataPulls = data[data["pullingCycle"].isin(pulls)]

    #Prepare data for filtering with resampling instead of savgol_filter
    dataPullsResample = dataPulls.copy(deep=True)
    dataPullsSmooth = resample_data(dataPullsResample, timeRes)

    #Plot the three different regions with different palettes:
    for region, palette in zip(regions, palettes):
        pullingCycles = dataPulls.pullingCycle[dataPulls.region == region].unique()
        for regionN, color in zip(pullingCycles, sns.color_palette(palette, len(pullingCycles))):
            condition = (dataPulls.pullingCycle == regionN) & (dataPulls.region == region)
            conditionSmooth = (dataPullsSmooth.pullingCycle == regionN) & (dataPullsSmooth.region == region)
            if dataPoints:
                axForceExtension.plot(dataPulls.loc[condition, separation], dataPulls.loc[condition, force],
                                      '.', alpha=0.02, color=color, label="_")
            axForceExtension.plot(dataPullsSmooth.loc[conditionSmooth, separation],
                                  dataPullsSmooth.loc[conditionSmooth, force],
                                  color=color, alpha=alpha, label=region + str(int(regionN)), lw=lineWidth)

    events = False

    if wlcData is not None:

        force = np.arange(0.8, 45, 0.2)
        for wlcRegion in dataPulls["wlcRegion"].unique():
            if wlcRegion > -1 and wlcRegion in wlcData.keys():
                axForceExtension.plot(wlc_from_fit(force, wlcData[wlcRegion]), force, "--", zorder=0, color="gray")

    if unfoldingData is not None and displayUnfoldingEvents:
        meanSeparations = []
        meanForces = []
        iCol = {columnName: columnIndex for columnIndex, columnName in enumerate(data.columns)}
        bufferedMargin = copy.deepcopy(eventMargin)
        unfoldingDataPulls = unfoldingData[unfoldingData["pullingCycle"].isin(pulls)]
        for startTime, endTime in zip(unfoldingDataPulls.startTime, unfoldingDataPulls.endTime):
            startIndex = data.index.get_loc(startTime)
            endIndex = data.index.get_loc(endTime)

            meanSeparations.append((data.iloc[startIndex, iCol[separation]] + data.iloc[endIndex, iCol[separation]])/ 2)
            meanForces.append((data.iloc[startIndex, iCol[force]] + data.iloc[endIndex, iCol[force]]) / 2)

            if startIndex < bufferedMargin:
                margin = data.index.get_loc(startTime) - 1
            else:
                margin = bufferedMargin

            slicer = slice(startIndex - margin, endIndex + margin)
            axForceExtension.plot(data.iloc[slicer, iCol["smoothSeparation"]],
                                  data.iloc[slicer, iCol["smoothForce"]], color="black", linewidth=8, label="_")

        for pullingCycle in pulls:
            if True in data.loc[data["pullingCycle"] == pullingCycle, "unfolding"].values:
                events = True

        if interactive and events:
            unfoldingDataPulls.loc[:, "meanSeparation"] = meanSeparations
            unfoldingDataPulls.loc[:, "meanForce"] = meanForces
            FollowDotCursor(axForceExtension, meanSeparations, meanForces, data, unfoldingDataPulls,
                            axis=axis, pulls=pulls, tolerance=20)
    if rulers is not None and wlcParameters is not None:
        print("ye")
        forceArray = np.arange(0.5, 60, 0.2)
        for length in rulers:
            axForceExtension.plot(wlc_series_accurate(forceArray,
                                                      wlcParameters["contourLengthDNA"],
                                                      wlcParameters["persistenceLengthDNA"],
                                                      wlcParameters["stretchModulusDNA"],
                                                      length, wlcParameters["persistenceLengthProtein"],
                                                      ), forceArray, color="gray")
        #ruler = pd.read_csv(sMBP_1300_DIG_BIOTIN_RULER, sep="\t", decimal=",", dtype=float)
        #ruler = ruler.ix[ruler.iloc[:, 3] < 40, :]
        #ruler.iloc[:,[0, 1, 2]] *= 1e3

        #axForceExtension.plot(ruler.iloc[:,2],ruler.iloc[:,3])
        #axForceExtension.plot(ruler.iloc[:,1],ruler.iloc[:,3])
        #axForceExtension.plot(ruler.iloc[:,0],ruler.iloc[:,3])
    if legend:
        legend = plt.legend(loc="upper left")
        legend.get_frame().set_facecolor('white')
    if xLimits is not None:
        axForceExtension.set_xlim(xLimits)
    if yLimits is not None:
        axForceExtension.set_ylim(yLimits)
    for tick in axForceExtension.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontSize)
    for tick in axForceExtension.yaxis.get_major_ticks():
        tick.label.set_fontsize(fontSize)
    if show:
        #mng = plt.get_current_fig_manager()
        #mng.window.showMaximized()
        plt.show()

    if interactive and events:
        unfoldingDataPulls.drop(["meanSeparation", "meanForce"], axis=1, inplace=True)


def grid_plotting(data, events=None, wlcFits=False, colWrap=4, xlim=(600,1000), extensionOffset=0, forceOffset=0):
    """
    Function to plot in a grid of subplots the different force extension curves

    Args:
        data (pandas.DataFrame): data containing all the information (also the unfolding)
        events (pandas.DataFrame): data containing the information of the unfolding events
        wlcFits (pandas.DataFrame): whether or not to plot the wlc fits to the data
        colWrap (int): maximum number of subplots in columns
        xlim (tuple): range for the x axis in the plots
    """
    if "region" not in data.columns.values:
        raise ValueError("Regions have not been assigned, please first run assign_regions() function")
    if "unfolding" not in data.columns.values:
        raise ValueError("Unfolding events have not been classified, please first run unfolding_data()")

    condition = (data["region"] != "stationary") & (data["region"] != "")
    #data.loc[data["unfolding"] == True, "region"] = "unfolding"
    #in case there are less than the colWrap number of regions, wrap to the number of regions
    if len(data.loc[condition, "pullingCycle"].unique()) < colWrap:
        colWrap = len(data.loc[condition, "pullingCycle"].unique())
    offsetedData = copy.copy(data)
    offsetedData["surfaceSepX"] += extensionOffset
    offsetedData["forceX"] += forceOffset

    g = sns.FacetGrid(offsetedData[condition], hue="region", col="pullingCycle", col_wrap=colWrap, palette="Set2")

    g.map_dataframe(force_extension_grid)
    g.set(xlim=xlim)
    #g.add_legend(loc=0)
    #plt.legend(loc='upper left')
    sns.plt.show()


def force_extension_grid(data, **kwargs):
    """
    Function to very simply plot force extension curves. Intended for the more complex grid plotting with seaborn
    Args:
        data (pandas.DataFrame): data to be plotted

    """
    axis = "x"
    filterParam=(51, 1)

    sns.set_style("white")
    if axis in ["x", "X", 0]:
        separation = data.surfaceSepX
        force = data.forceX
    elif axis in ["y", "Y", 1]:
        separation = data.surfaceSepY
        force = data.forceY
    if data["region"] is "unfolding":
        lw=20
    else:
        lw=2
    if "wlcForce" + str.upper(axis) in data.columns:
        plt.plot(data["wlcExtension" + str.upper(axis)] , data["wlcForce" + str.upper(axis)], '.', ms=5,
                 color="black")
    plt.plot(separation, force, '.', alpha=0.7, lw=lw, **kwargs)
    if data.iloc[0, 15] != "unfolding" and len(data)>filterParam[0]*2:
        plt.plot(savgol_filter(separation, *filterParam), savgol_filter(force, *filterParam), lw=2, **kwargs)
    else:
        plt.plot(separation, force, alpha=1, lw=5, color="red")
    #print(data.eventID)
    ruler = pd.read_csv(sMBP_1300_DIG_BIOTIN_RULER, sep="\t", decimal=",", dtype=float)
    ruler = ruler.ix[ruler.iloc[:, 3] < 40, :]
    ruler.iloc[:,[0, 1, 2]] *= 1e3

    #plt.plot(ruler.iloc[:,2],ruler.iloc[:,3])
    #plt.plot(ruler.iloc[:,1],ruler.iloc[:,3])
    #plt.plot(ruler.iloc[:,0],ruler.iloc[:,3])

    for unfoldingEvent in list(data.eventID.unique()[1:]):
        #print(unfoldingEvent)
        condition = data["eventID"] == unfoldingEvent
        plt.plot(separation.loc[condition], force.loc[condition], lw=5, color="black")


def force_time(data, smooth="both", regions=True):

    colors = [(0.3333333333333333, 0.6588235294117647, 0.40784313725490196),
              (0.7686274509803922, 0.3058823529411765, 0.3215686274509804),
              (0.5058823529411764, 0.4470588235294118, 0.6980392156862745)]
    for region, color in zip(data.region.unique(), colors):

        for pullingCycle in data.pullingCycle.unique():
            condition = (data.region == region) & (data.pullingCycle == pullingCycle)
            plt.plot(data.time[condition], data.forceX[condition], '.', alpha=0.3, color=color)
            plt.plot(data.time[condition], medfilt(data.forceX[condition], 51), color=color)

    plt.show()


def trap_separation_time(data, axis="x", palettes=["PuBu_d", "YlOrRd_d", "Greens_d"], alpha=0.7, legend=False):

    #Plot the three different regions with different palettes:
    for region, palette in zip(sorted(data.region.unique()), palettes):
        pullingCycles = data.pullingCycle[data.region == region].unique()

        for regionN, color in zip(pullingCycles, sns.color_palette(palette, len(pullingCycles))):
            condition = (data.pullingCycle == regionN) & (data.region == region)
            plt.plot(data.loc[condition, "time"], data.loc[condition, "trapSep" + axis.upper()]-2100,
                         color=color, alpha=alpha, label=region + str(int(regionN)), lw=4)
            plt.plot(data.loc[condition, "time"], data.loc[condition, "surfaceSep" + axis.upper()],
                         color=color, alpha=alpha, lw=4)

    if legend:
        legend = plt.legend(loc="upper left")#, bbox_to_anchor=(1,1))
        legend.get_frame().set_facecolor('white')


    plt.show()


def plot_fluorescence_signal(fluorescenceData, data=None, color="532nm", signalRegion=(70,76), backgroundRegion=(0,20),
                             timeResolution=0.087, timeOffset=0):

    if data is not None:
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
    else:
        fig, (ax2, ax3) = plt.subplots(2, 1, sharex=True)

    image = fluorescenceData[color]
    time = np.arange(image.shape[1]) * timeResolution + timeOffset
    signal = np.mean(image[signalRegion[0]: signalRegion[1], :], axis=0)
    background = np.mean(image[backgroundRegion[0]: backgroundRegion[1], :], axis=0)
    ax2.plot(time, signal)
    ax2.plot(time, background)
    ax3.imshow(image, extent=[timeOffset, timeResolution*image.shape[1], image.shape[0], 0])
    if data is not None:
        ax1.plot(data["time"], data["LcProtein"])
        ax1.set_ylim(0,500)
    plt.show()


def plot_xcorr(x, y, sFreq=100000, maxlag=None, lagStep=1, show=True):
    if maxlag is None:
        maxlag = len(x) // 4
    plt.plot(np.arange(0, maxlag + 1, lagStep)/sFreq*1000, cross_correlation(x, y, maxlag, lagStep))
    plt.plot(np.arange(0, maxlag + 1, lagStep)/sFreq*1000, cross_correlation(x, y, maxlag, lagStep))
    plt.ylabel("Cross-correlation")
    plt.xlabel("Lag time (ms)")
    if show:
        plt.show()


def save_figure(metadata, extension=".pdf"):
    """
    Function to save the figure and the metadata properly. The name of the figure is the date and the first keyword.
    The figures are saved in the specified folder inside the function (hardcoded)

    Args:
        metadata (dict): dictionary with the metadata of the figure. At least, it should include:
            filePath (list): path(s) of the data file(s) used for the plot ||
            variables (list): variables plotted ||
            keywords (list): a list with the keywords representing the plot ||
            varParameters (dict): comments or parameters used for each variable

    Example:
        >>> plt.plot(np.random.rand(10), np.random.rand(10))
        >>> plt.xlabel("Random variable X")
        >>> plt.ylabel("Random variable Y")
        >>> metadata = {"filePath": myPath, "calibrationFilePath": myCalibrationPath, "variables": ["randomX", "randomY"], "keywords": ["randomX", "randomY", "usageExample"]}
        >>> save_figure(metadata)

    """

    figuresFolder = "D:/data/figures/" + time.strftime("%Y_%m") + "/"
    #Check folder and if not create
    if not os.path.exists(figuresFolder):
        os.makedirs(figuresFolder)


    metadata["date"] = time.strftime("%Y/%m/%d")
    figurePath = figuresFolder + time.strftime("%Y_%m_%d") + "_" + metadata["keywords"][0] + extension

    #check if there is a figure with the same name and if so iterate over the alphabet to rename the new figure
    if os.path.exists(figurePath):
        appendices = iter(string.ascii_lowercase)
        currentAppendix = appendices.__next__()
        while os.path.exists(figurePath[:-4] + "_" + currentAppendix + extension):
            currentAppendix = appendices.__next__()
        figurePath = figurePath[:-4] + "_" + currentAppendix + extension
    metadata["figurePath"] = figurePath

    #Update the record file with all the plots
    recordsFile = "D:/data/figures/record.csv"
    records = pd.read_csv(recordsFile, sep="\t", index_col=0)
    #Append an empty row
    records = records.append(pd.DataFrame({key: None for key in metadata.keys()},
                                          dtype=object, index=[0]), ignore_index=True)
    #Complete the row with the new data
    for key in list(metadata.keys()):
        records[key].iloc[-1] = metadata[key]
    records.to_csv(recordsFile, sep="\t")

    #Save the figure to the figures folder
    #plt.tight_layout()
    plt.savefig(figurePath)


#def plot_thermal_calibration(data, fit, *args, **kwargs):
#    if "psd" not in data.columnsvalues.tolist():


