#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.signal import savgol_filter


def assign_rough_regions(data, axis="x", **kwargs):
    """
    Function to assign and include ROUGHLY the regions (stationary, pulling or retracting) of the traps, according
    to criteria related to the first and the second derivative of the trap separation with time

    Args:
        data (pandas.DataFrame): pandas DataFrame containing the data
        threshold (list): list in which the first element indicates the maximum value accepted for the first derivative
        and the second the maximum value for the second derivative
        axis (str): axis of the movement of the trap, either "x" or "y"


    Returns:
        data (pandas.DataFrame): the same data including the column with the regions
    """
    pd.set_option('mode.chained_assignment',None)
    axis = axis.upper()
    data["region"] = ""

    mirrorDiff = pd.Series(savgol_filter(data["trapSep" + axis], 71, 1), index=data.index).diff().fillna(
        method="backfill")

    sigma = np.max([mirrorDiff.max(), mirrorDiff.min()]) / 20

    #set the criteria for assigning region using the first and second derivative of the trapSepX
    data.loc[mirrorDiff <= sigma, "region"] = "stationary"
    data.loc[mirrorDiff > sigma, "region"] = "pulling"
    data.loc[mirrorDiff < -sigma, "region"] = "retracting"
    #Column conversion to category
    data.region.astype("category")

    return data


def label_regions(data, minRegionLength, verbose=True):
    """
    Function to improve the region assignment and to add region number

    Args:
        data (pandas.DataFrame): pandas DataFrame containing the data
        minRegionLength (int): minimum number of data points to consider a region as valid (otherwise it is removed)
        verbose (bool): if True, region information will be printed
    Returns:
        data (pandas.DataFrame): the same data including the column with the regions
    """

    #check if the regions have been assigned
    if "region" not in data.columns.values.tolist():
        raise ValueError("The regions have not been assigned yet. Please see the function assign_rough_regions")
    #define the column pullingCycle
    data["pullingCycle"] = 0.0
    #label the isolated regions
    for region in list(data["region"].unique()):
        data["pullingCycle"] += ndimage.label(data["region"] == region)[0]

    #remove and rearrange the false stationary regions (transition from pulling to retracting or viceversa)
    data = remove_false_stationary(data, minRegionLength)
    #Start the region counter at 0, not at 1 (labelling starts at 1)
    data["pullingCycle"] -= 1
    if verbose:
        for region in list(data["region"].unique()):
            print(ndimage.label(data["region"] == region)[1], region, "regions identified")

    return data


def remove_false_stationary(data, minRegionLength=10):
    """
    Function to remove and rearrange the false stationary regions (transition from pulling to retracting or viceversa)

    Args:
        data (pandas.DataFrame): pandas DataFrame containing the data
        minRegionLength (int): minimum number of data points to consider a region as valid (otherwise it is removed)
    Returns:
        data (pandas.DataFrame): the same data with fixed stationary regions
    """

    for pullingCycle in list(data.loc[data["region"] == "stationary", "pullingCycle"].unique()):
        regionMask = (data["region"] == "stationary") & (data["pullingCycle"] == pullingCycle)
        if len(data.loc[regionMask, "pullingCycle"]) < minRegionLength:
            #make nan the statioanry regions that are too sort (for fill methods)
            data.loc[regionMask, ["region", "pullingCycle"]] = np.nan
            #fill first forward half of the region and the other half backwards
            data = data.fillna(method="ffill", limit=len(data.loc[regionMask])//2+1)
            data = data.fillna(method="bfill", limit=len(data.loc[regionMask])//2+1)
    #relabel the stationary regions
    data.loc[data["region"] == "stationary", "pullingCycle"] = 0
    data["pullingCycle"] += ndimage.label(data["region"] == "stationary")[0]

    return data


import streamlit as st
@st.cache()
def assign_regions(data, axis="X", minRegionLength=10, verbose=True, **kwargs):
    """
    Function to assign and include the regions (stationary, pulling or retracting) of the traps, based on mirror
    movement

    Args:
        data (pandas.DataFrame): pandas DataFrame containing the data
        axis (str): axis of the movement of the trap, either "x" or "y"
        minRegionLength (int): minimum number of data points to consider a region as valid (otherwise it is removed)
        verbose (bool): if True, region information will be printed

    Returns:
        data (pandas.DataFrame): the same data including the column with the regions
    """

    if axis not in ["x", "X", "y", "Y"] and axis is not None:
        raise ValueError("Axis should be 'x', 'y' or None")

    regionData = assign_rough_regions(data, axis, **kwargs)
    data = label_regions(regionData, minRegionLength, verbose)

    return data

