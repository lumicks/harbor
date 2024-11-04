#!/usr/bin/env python
# -*- coding: utf-8 -*-

from foldometer.analysis.region_classification import *
from foldometer.tools.plots import force_extension_curve
from scipy import ndimage
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor
import numpy as np
import pandas as pd
from copy import deepcopy


def select_event_threshold(data, forceSTD, forceMeanChange, axis="X", forceChannel="force"):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, axisbg='#FFFFCC')

    x, y = 4*(np.random.rand(2, 100)-.5)
    ax.plot(data["time"], data[forceChannel + axis] / (data[forceChannel + axis].max() / forceSTD.max()), alpha=0.5)
    #ax.set_xlim(-2, 2)
    #ax.set_ylim(-2, 2)
    ax.plot(data["time"], forceSTD)
    #ax.plot(data["time"], forceMeanChange)

    # set useblit = True on gtkagg for enhanced performance
    cursor = Cursor(ax, useblit=True, color='red', linewidth=2, vertOn=False )

    plt.show()


def identify_unfolding_events(data, axis="x", forceChannel="force", window=15, STDthreshold=0.8, forceThreshold=5):
    """
    Function to identify potential unfolding events using the rolling standard deviation and rolling mean of the
    difference

    Args:
        data (pandas.DataFrame): pandas DataFrame containing the data
        axis (str): axis to perform the identification, either "x" or "y"
        forceChannel (str): either "PSD1Force", "PSD2Force" or "force" (combined signal)
        window (int): number of points to compute the rolling std and mean of the difference
        STDthreshold (float): threshold in the standard deviation for choosing unfolding events
        forceThreshold (float): threshold in the force for choosing unfolding events (too much noise at low forces)

    Returns:
        data (pandas.DataFrame): the same data including the column with the potential unfolding events
    """

    axis = axis.upper()

    data["unfolding"] = False

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

    force = data[forceChannel + axis]
    #Calculate the rolling standard deviation
    forceSTD = force.rolling(window, center=True).std()
    #print(window)
    #Calculate the rolling mean of the difference
    forceMeanChange = force.rolling(window, center=True).apply(func=mean_diff, args=(1,))
    #select_event_threshold(data, forceSTD, forceMeanChange, axis)
    #forceMeanChange2 = pd.rolling_apply(force, window, func=mean_diff, args=(2,), center=True)

    #plt.plot(forceMeanChange)
    #plt.plot(forceSTD)
    #plt.axhline(STDthreshold)
    #plt.show()
    forceSTD[forceSTD < STDthreshold] = 0
    forceMeanChange[abs(forceMeanChange) < STDthreshold / 3] = 0
    #forceMeanChange2[abs(forceMeanChange2) < 0.3] = 0


    forceSTD = forceSTD.fillna(0)
    forceMeanChange = forceMeanChange.fillna(0)
    #Create a mapping mask to filter the data and assign unfolding events as True
    mask = forceSTD + abs(forceMeanChange)# + abs(forceMeanChange2)
    mask[(mask < STDthreshold/2) | (data[forceChannel + axis] < forceThreshold) | (data[forceChannel + axis] > 60)] = 0
    mask = (mask != 0)
    # mask = (((forceSTD - forceMeanChange) > threshold) & (data["force" + axis] > 7))
    data.loc[mask, "unfolding"] = True

    return data


def calculate_force_change(data, axis="x", forceChannel="force", distanceChannel="surfaceSep", window=15):
    """
    Function to calculate the force change in unfolding events with the mean before and after the identified event

    Args:
        data (pandas.DataFrame): pandas DataFrame containing the data
        axis (str): axis to perform the identification, either "x" or "y"
        forceChannel (str): either "PSD1Force", "PSD2Force" or "force" (combined signal)
        distanceChannel (str): either "surfaceSep" (from PSDs) or "trackingSep" (from image tracking)
        window (int): number of points to be considered in the force calculation before and after the event

    Returns:
        unfoldingData (pandas.DataFrame): pandas DataFrame with information about the unfolding events
    """
    axis = axis.upper()

    #check if the regions have been assigned
    if "unfolding" not in data.columns.values.tolist():
        raise ValueError("The unfolding events have not yet been identified. See function identify_unfolding_events")

    #Label the different isolated events using scipy.ndimage
    data["eventID"], eventsNumber = ndimage.label(data["unfolding"])

    #Start the counting in 0
    data["eventID"] -= 1
    #Show how many events were identified
    print(eventsNumber, "events identified")

    def averaged_values(column, startT, endT, window=5):
        start = column.index.get_loc(startT)
        end = column.index.get_loc(endT)
        averagedBefore = column.iloc[start-window: start - 3].mean()
        averagedAfter = column.iloc[end + 3: end + window].mean()
        diffAverage = averagedAfter - averagedBefore
        return averagedBefore, averagedAfter, diffAverage

    startForce = []
    forceChange = []

    pullingCycle = []
    #Take the first and last times point of each unfolding event, discarding the first point because it is the
    # unclassified regions
    times = {"startTimes": data.groupby("eventID").time.first()[1:], "endTimes": data.groupby("eventID").time.last()[1:]}
    newWindow = deepcopy(window)
    import streamlit as st
    for startTime, endTime in zip(times["startTimes"], times["endTimes"]):
        if (startTime > data["time"].min())*(startTime < data["time"].max())*(endTime > data["time"].min())*(endTime < data["time"].max()):
            if data.index.get_loc(startTime) < newWindow:
                window = data.index.get_loc(startTime) - 1
            else:
                window = newWindow
            forceBefore, forceAfter, forceDifference = averaged_values(data[forceChannel+axis], startTime, endTime, window)
            startForce.append(forceBefore)
            forceChange.append(forceDifference)
            pullingCycle.append(data.loc[startTime, "pullingCycle"])

    unfoldingData = pd.DataFrame({"startTime": times["startTimes"], "endTime": times["endTimes"],
                                    "force": startForce, "forceChange": forceChange, "pullingCycle": pullingCycle})

    return unfoldingData


import streamlit as st
@st.cache()
def find_unfolding_events(data, axis="x", forceChannel="force", distanceChannel="surfaceSep", unfoldingThreshold=1,
                          forceThreshold=5, rollingWindow=5, unfoldingWindow=15, plot=False, **kwargs):
    """
    Function to calculate the force change in unfolding events with the mean before and after the identified event

    Args:
        data (pandas.DataFrame): pandas DataFrame containing the data
        axis (str): axis to perform the identification, either "x" or "y"
        forceChannel (str): either "PSD1Force", "PSD2Force" or "force" (combined signal)
        distanceChannel (str): either "surfaceSep" (from PSDs) or "trackingSep" (from image tracking)
        order (int): the order wanted for the event finder algorithm
        unfoldingThreshold (float): minimum value for choosing unfolding events
        forceThreshold (float): threshold in the force for choosing unfolding events (too much noise at low forces)
        rollingWindow (int): number of points to be considered in the standard deviation and diff calculation
        unfoldingWindow (int): number of points to be considered in the force calculation before and after the event
        plot (bool): whether or not to plot the unfolding events on the force-extension curve. Default is True
        **kwargs (): arguments for the plot_unfolding_events function (palettes, alpha and legend)

    Returns:
        unfoldingForces (pandas.DataFrame): pandas DataFrame with information about the unfolding events force,
        extension change and time
    """

    axis = axis.upper()
    #print(rollingWindow)
    data = identify_unfolding_events(data, axis, forceChannel, rollingWindow, unfoldingThreshold, forceThreshold)

    #check if there are any unfolding events
    if True in data["unfolding"].unique():
        unfoldingData = calculate_force_change(data, axis, forceChannel, distanceChannel, unfoldingWindow)

        if plot:
            force_extension_curve(data, unfoldingData, None, axis, forceChannel, distanceChannel, **kwargs)
        #re-arrange the dataframe
        unfoldingData = unfoldingData[["startTime", "endTime", "force", "forceChange", "pullingCycle"]]
        return unfoldingData

    else:
        if plot:
            figure = plt.figure()
            ax = figure.add_subplot(111)
            force_extension_curve(data, axis=axis, show=False, figure=figure)
            ax.text(0.5, 0.75, "No unfolding events detected", horizontalalignment='center', verticalalignment='top',
                    transform=ax.transAxes, fontsize=15,
                    bbox={'facecolor':'#FFE3B9', 'alpha':0.9, 'boxstyle':'round,pad=0.5'})
            plt.show()
        print("No unfolding events detected, check plot for potential missing events")
