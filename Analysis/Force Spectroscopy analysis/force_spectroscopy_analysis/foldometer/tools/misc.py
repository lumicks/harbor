#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector
import warnings

warnings.filterwarnings('ignore')


def column_indexer(data):
    """
    Creates a dictionary with column labels as keys and column integer index as value, in order to use iloc method
    with labels

    Args:
        data (pandas.DataFrame): dataframe with the columns to be indexed
    Returns:
        idCol (dict): dictionary with labels mapped to indices
    """
    idCol = {label: index for index, label in enumerate(data.columns)}
    return idCol


def format_pulling_cycles(pulls, initialPull, totalPulls):
    """
    Function to unify the format of input pulls by the user

    Args:
        pulls: pulls to be processed, input from the user in the form of int, tuple or list
        totalPulls (int): total number of pulls in the data
    Returns:
        pulls (list): reformatted pulls in the form of a list
    """
    if pulls is None:
        pulls = list(np.arange(initialPull, initialPull + totalPulls))
    elif isinstance(pulls, int):
        pulls = [pulls]
    elif isinstance(pulls, tuple):
        pulls = list(np.arange(pulls[0], pulls[1]))
    if not isinstance(pulls, (int, tuple, list)):
        raise ValueError("Incorrect input for pulls, use an integer, a list or a tuple")

    return pulls


def data_selection(data, columnX="time", columnY="forceX", ylim=None):
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
    if ylim is not None:
        ax.set_ylim(ylim[0], ylim[1])
    x = data[columnX]
    y = data[columnY]
    ax.plot(x, y, '-')

    fig.suptitle("Press any key to confirm selection and close this window", fontsize=14, fontweight='bold')
    ax.set_title("Select data by dragging the mouse")

    def onselect(xmin, xmax):
        """
        Sub-function to assign the minimum and maximum of the range and plot the selection"

        Args:
            xmin (int): minimum index of selection
            xmax (int): maximum index of selection
        """
        global indmin
        global indmax
        global span

        try:
            span.remove()
            fig.canvas.draw()
        except:
            pass

        indmin, indmax = np.searchsorted(x, (xmin, xmax))
        indmax = min(len(x) - 1, indmax)
        ax.axvspan(x.values[int(indmin)], x.values[int(indmax)], facecolor="red", alpha=0.5)
        fig.canvas.draw()

    # set useblit True on gtkagg for enhanced performance
    span = SpanSelector(ax, onselect, 'horizontal', useblit=True, rectprops=dict(alpha=0.5, facecolor='red'))
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized()
    plt.show()

    while True:
        wait = plt.waitforbuttonpress(timeout=60)
        if wait:
            plt.close('all')
            break
        else:
            pass

    idCol = column_indexer(data)
    mask = (data[columnX] >= data.iloc[indmin, idCol[columnX]]) & \
           (data[columnX] <= data.iloc[indmax, idCol[columnX]])
    subData = data.loc[mask, :]

    return subData


def data_deletion(data, columnX="time", columnY="forceX"):
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

    subset = data_selection(data, columnX, columnY)

    data = data.drop(subset.index)
    return data


def resample_data(data, time, unit='s'):
    """
    Function to resample and decimate the data to certain time resolution using the mean

    Args:
        data (pandas.DataFrame): data containing
        time (int): time resolution in milliseconds
    """
    if not isinstance(data.index, pd.core.indexes.timedeltas.TimedeltaIndex):
        data.index = pd.to_timedelta(data.index, unit=unit)

    if "region" in data.columns:
        regionDict = {"pulling": 0, "retracting": 1, "stationary": 2}
        data["region"] = data["region"].map(regionDict)
        dataResampled = data.resample(rule=str(time)+'ms').mean()
        dataResampled["region"] = dataResampled["region"].map({v: k for k, v in regionDict.items()})
    else:
        dataResampled = data.resample(rule=str(time)+'ms').mean()

    return dataResampled