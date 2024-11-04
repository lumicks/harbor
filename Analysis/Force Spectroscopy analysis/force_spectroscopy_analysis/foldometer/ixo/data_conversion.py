#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Example:
    >>> import foldometer as fm
    >>> fileName = "fileName.dat"
    >>> calName = "calibrationFileName.dat"
    >>> data = fm.analyse_file(fileName, calName, plotEvents=True)
"""

from foldometer.ixo.binary import read_file
from foldometer.analysis.thermal_calibration import calibration_file
from foldometer.analysis.tweezers_parameters import MIRRORVOLTDISTANCEFACTOR, CCD_PIXEL_NM, CCD_FREQUENCY, \
    get_mirror_values
from foldometer.analysis.region_classification import assign_regions
from foldometer.analysis.event_classification import find_unfolding_events
from foldometer.tools.misc import data_selection
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from copy import deepcopy
import warnings
import pkg_resources

noiseFile = pkg_resources.resource_filename('foldometer', 'data/noise.dat')


def calculate_psd_offset(psdData):
    """
    Function to get the offset of the psd position signal in V

    Args:
        psdData (pandas.Series): single column with PSD position information

    Returns:
        offset (pandas.Series): offset of the PSD signals to be removed from data

    """
    offset = psdData.iloc[0:10].mean()
    return offset


def single_volts_to_distance(rawData, beta, offset, normData=None, normalized=True):
    """
    Function to convert a single column of raw PSD data (in V) into usable data (in nm)

    Args:
        rawData (pandas.Series): data containing the PSD signal (displacement only) in [V]
        beta (float): distance calibration factor in units of [V/nm]
        offset (float): offset from a single channel of the PSD signal

    Returns:
        data (pandas.Series): converted data from [V] to [nm]

    """
    data = (rawData - offset) / beta
    if normalized is False:
        data = (rawData - offset) / (beta * normData)
    return data


def psd_volts_to_distance(rawData, betaFactors, offsets=None, normData=None, normalized=True):
    """
    Function to convert raw PSD data (in V) into usable data (in nm)

    Args:
        rawData (pandas.DataFrame): data containing the PSD signal (displacement only) in [V]
        betaFactors (dict): dictionary with the four distance calibration factor in units of [V/nm]
        offsets (dict): offsets from both directions of both PSD signals

    Returns:
        displacement (pandas.DataFrame): displacement of the bead from the center of the traps in units of [nm]

    """

    if not isinstance(rawData, pd.core.frame.DataFrame):
        raise TypeError("For single column conversion use the function single_volts_to_distance")
    for label in rawData.columns.values.tolist():
        assert ("PSD1VxDiff" in label or "PSD1VyDiff" in label or "PSD2VxDiff" in label or
                "PSD2VyDiff" in label), "At least one of the columns does not contain PSD information"
    if isinstance(betaFactors, (int, float)):
        raise TypeError("For single column conversion use the function single_volts_to_distance")

    if offsets is None:
        offsets = {}
        for key in ["PSD1x", "PSD1y", "PSD2x", "PSD2y"]:
            offsets[key] = calculate_psd_offset(rawData[key[:4]+"V"+key[-1:]+"Diff"])

    displacement = deepcopy(rawData)

    for label in displacement.columns.values.tolist():
        labelParameter = label[0:4]+label[5:6]
        displacement[label] = single_volts_to_distance(displacement[label], betaFactors[labelParameter],
                                                       offsets[labelParameter], normData=normData[label[0:6]+"Sum"],
                                                       normalized=normalized)
        displacement.rename(columns={label: labelParameter+"D"}, inplace=True)

    return displacement


def single_calculate_force(psdDisplacementData, stiffness):
    """
    Function to calculate force (in [pN]) exerted on the bead from the displacement

    Args:
        psdDisplacementData (pandas.Series): single column of the data with the displacement (in [nm])
        of the bead from the center of the trap

        stiffness (float): stiffness of the trap as calculated from the thermal calibration in units of [pN/nm]

    Returns:
        force (pandas.Series): column with the calculated force in units of [pN]

    """

    force = psdDisplacementData * stiffness

    return force


def psd_calculate_force(psdDisplacementData, stiffness):
    """
    Function to calculate force exerted on the beads from the displacements

    Args:
        psdDisplacementData (pandas.DataFrame): sdata with the displacement (in [nm])
        of the beads from the center of the traps
        stiffness (dict): stiffness of the traps as calculated from the thermal calibration in units of [pN/nm]

    Returns:
        force (pandas.DataFrame): new DataFrame with added information of the force

    """

    if not isinstance(psdDisplacementData, pd.core.frame.DataFrame):
        raise TypeError("For single column conversion use the function single_volts_to_distance")
    for label in psdDisplacementData.columns.values.tolist():
        assert ("PSD1xD" in label or "PSD1yD" in label or "PSD2xD" in label or
                "PSD2yD" in label), "At least one of the columns does not contain PSD information"
    if isinstance(stiffness, (int, float)):
        raise TypeError("For single column conversion use the function single_volts_to_distance")

    force = deepcopy(psdDisplacementData)

    for label in force.columns.values.tolist():

        labelParameter = label[0:5]
        force[label] = single_calculate_force(force[label], stiffness[labelParameter])
        force.rename(columns={label: label[0:4]+"Force" + str.upper(label[4])}, inplace=True)

    averagedForceX = (force["PSD1ForceX"] - force["PSD2ForceX"]) / 2
    averagedForceY = (force["PSD1ForceY"] - force["PSD2ForceY"]) / 2
    averagedForceX.name = "forceX"
    averagedForceY.name = "forceY"

    force = pd.concat([force[["PSD1ForceX", "PSD2ForceX", "PSD1ForceY", "PSD2ForceY"]],
                       averagedForceX,
                       averagedForceY], axis=1)

    return force


def trap_separation(mirrorPosition, header=None):
    """
    Function to calculate the separation between the center of the traps

    Args:
        mirrorPosition (pandas.DataFrame): absolute position (in x and y) of the steerable mirror in [V]
    Returns:
        trapSeparation (pandas.DataFrame): x and y separation between the traps in [nm]

    """

    # The minus sign in the trap separation depends on the direction of the movement of the trap
    STATIONARYMIRRORX, STATIONARYMIRRORY = get_mirror_values(header)
    trapSeparationX = (-STATIONARYMIRRORX + mirrorPosition["MirrorX"]) * MIRRORVOLTDISTANCEFACTOR
    trapSeparationY = (-STATIONARYMIRRORY + mirrorPosition["MirrorY"]) * MIRRORVOLTDISTANCEFACTOR
    trapSeparationX.name = "trapSepX"
    trapSeparationY.name = "trapSepY"
    trapSeparation = pd.concat([trapSeparationX, trapSeparationY], axis=1)

    return trapSeparation


def surface_separation(displacement, trapSeparation, radii=(1050, 1050)):
    """
    Function to calculate the separation between the surface of the beads

    Args:
        displacement (pandas.DataFrame): displacement of the beads from the center of the traps [nm]
        trapSeparation (pandas.DataFrame): x and y separation between the traps in [nm]
    Returns:
        surfaceSeparation (pandas.DataFrame): x and y separation between the surfaces of the beads in [nm]

    """

    surfaceSeparationX = trapSeparation.trapSepX - (radii[0] + radii[1]) + displacement["PSD2xD"] - displacement[
        "PSD1xD"]
    surfaceSeparationY = trapSeparation.trapSepY - 2110 + displacement["PSD2yD"] - displacement["PSD1yD"]
    surfaceSeparationX.name = "surfaceSepX"
    surfaceSeparationY.name = "surfaceSepY"
    surfaceSeparation = pd.concat([surfaceSeparationX, surfaceSeparationY], axis=1)
    return surfaceSeparation


def remove_parasitic_noise(rawData, data, header):
    try:
        noiseFileName = 'data/' + str(header["fileName"][-34:-27], 'utf-8') + '_noise.dat'

        noiseFile = pkg_resources.resource_filename('foldometer', noiseFileName)
        noiseMetadata, noiseFit, noiseRawData, noiseBead = read_file(noiseFile)
    except:
        noiseFile = pkg_resources.resource_filename('foldometer', 'data/noise.dat')
        noiseMetadata, noiseFit, noiseRawData, noiseBead = read_file(noiseFile)

    print(noiseFile)

    offset = noiseFit["offset"]

    noiseData = process_data(noiseRawData, offset, noiseFit, noiseBead, header=header,
                        radii=[header["beadRadius1"], header["beadRadius2"]], normalized=True, noiseRemoval=False)
    noisePolyFitPSD1 = np.polyfit(noiseRawData["MirrorX"], noiseData["PSD1ForceX"], deg=10)
    polynomialPSD1 = np.polyval(noisePolyFitPSD1, rawData["MirrorX"])
    noisePolyFitPSD2 = np.polyfit(noiseRawData["MirrorX"], noiseData["PSD2ForceX"], deg=10)
    polynomialPSD2 = np.polyval(noisePolyFitPSD2, rawData["MirrorX"])

    data.loc[:, "PSD1ForceX"] -= polynomialPSD1
    data.loc[:, "PSD2ForceX"] -= polynomialPSD2


    data.loc[:, "PSD1ForceX"] -= data["PSD1ForceX"].min() - 0.5
    data.loc[:, "PSD2ForceX"] -= data["PSD2ForceX"].max() + 0.5
    data.loc[:, "forceX"] = (data.loc[:, "PSD1ForceX"] - data.loc[:, "PSD2ForceX"]) / 2
    data.loc[:, "forceX"] -= data["forceX"].min() - 0.5

import streamlit as st
@st.cache()
def process_data(rawData, offset, calibrationParameters, beadData=None, header=None, radii=(1050, 1050),
                 normalized=True, beadTracking=False, noiseRemoval=False):
    """
    Function to convert a set of raw data (already read from a file) into useful distance and force information

    Args:
        rawData (pandas.DataFrame): data recorded from Foldometer
        offset (dict): offset from both directions of both PSD signals
        calibrationParameters (pandas.DataFrame): parameters from the thermal calibration
        radii (tuple): radius of each bead
    Returns:
        data (pandas.DataFrame): processed data including separation and force

    """

    psdColumns = ["PSD1VxDiff", "PSD1VyDiff", "PSD2VxDiff", "PSD2VyDiff"]
    normColumns = ["PSD1VxSum", "PSD1VySum", "PSD2VxSum", "PSD2VySum"]

    #offsets = 10
    displacement = psd_volts_to_distance(rawData[psdColumns], calibrationParameters.beta,
                                         offsets=offset, normData=rawData[normColumns], normalized=normalized)
    force = psd_calculate_force(displacement, calibrationParameters.stiffness)
    trapSeparation = trap_separation(rawData[["MirrorX","MirrorY"]], header)
    surfaceSeparation = surface_separation(displacement, trapSeparation, radii=radii)
    #create instance of the DataFrame and join everything together
    data = displacement
    data = data.join([force, trapSeparation, surfaceSeparation])
    data.insert(0, "time", rawData["time"])

    if noiseRemoval:
        remove_parasitic_noise(rawData, data, header)


    #data["forceX"] -= data["forceX"].min() - 0.5

    #Bead Tracking
    if len(beadData.index) is 0 and beadTracking:
        warnings.warn("WARNING: there is no bead tracking data in this file, proceeding without merging")
        beadTracking = False

    if beadTracking:
        beadDataMerged = data_beadData_merging(beadData, data)
        beadDataMerged.index = beadDataMerged["time"]
        return beadDataMerged
    else:
        return data


def process_file(fileName, calibrationFileName=None, beadTracking=False, fitParameters=False,
                 noiseRemoval=False, mode=None, normalized=True, **kwargs):
    """
    Function to convert the data of a file into useful distance and force information

    Args:
        fileName (string): file with data recorded from Foldometer
        calibrationFileName (string): calibration file with data recorded from Foldometer. If None, takes the fit stored
         in the metadata
        beadTracking (bool): whether or not to use and merge the beadTracking data. Default is False
        fitParameters (bool): if True, function returns a tuple with first the processed data
         and second the thermal calibration parameters. Default is False
         mode (str): if None, all data is read. If "compact[AXIS]" only the Diff columns to the corresponding axis
        **kwargs: all optional arguments accepted by the function calibrate_file
    Returns:
        data (pandas.DataFrame): processed data including separation and force

    """

    header, calibrationParameters, rawData, beadData = read_file(fileName, mode=mode)
    offset = calibrationParameters["offset"]

    if calibrationFileName is not None:
        calibrationParameters = calibration_file(calibrationFileName, **kwargs)

    data = process_data(rawData, offset, calibrationParameters, header=header,
                        radii=[header["beadRadius"], header["beadRadius"]], normalized=normalized,
                        noiseRemoval=noiseRemoval)

    if len(beadData.index) is 0 and beadTracking:
        warnings.warn("WARNING: there is no bead tracking data in this file, proceeding without merging")
        beadTracking = False

    if beadTracking:
        beadDataMerged = data_beadData_merging(beadData, data)

    if fitParameters:
        if beadTracking:
            return data, beadDataMerged, calibrationParameters
        else:
            return data, calibrationParameters
    else:
        if beadTracking:
            return beadDataMerged
        else:
            return data


def data_beadData_merging(beadData, data, factor=0.0074, columns=None):
    """
    Function to merge the PSD and the bead tracking data by interpolating PSD values.

    Args:
        beadData (pandas.DataFrame): data from the bead tracking
        data (pandas.DataFrame): data from the PSDs
        columns (list): list of the columns from the PSD to be merged. If it is None, all columns are merged (default)

    Returns:
        mergedData (pandas.DataFrame): merged DataFrame with as many entries as the bead tracking data
    """

    sampleFreq = 1/data.time.diff().mean()

    if columns is None:
        if sampleFreq < CCD_FREQUENCY:
            columns = list(beadData.columns)
            beadData = beadData[columns]
        else:
            columns = list(data.columns)
            data = data[columns]
    #ToDo: Add the possibility of choosing the columns

    mergedData = data.merge(beadData, how="outer").sort_values(by="time")
    mergedData.index = np.arange(0, len(mergedData))
    mergedData[columns] = mergedData[columns].interpolate(method="index")
    mergedData.dropna(inplace=True)
    mergedData["trackingSepX"] = (mergedData["Bead2X"] - mergedData["Bead1X"]) / CCD_PIXEL_NM - 2100

    return mergedData


def analyse_file(fileName, calName=None, unfoldingEvents=True, beadTracking=False, fitParameters=False,
                 select_data=True, mode=None, unfoldingThreshold=1, unfoldingWindow=15, plotEvents=True,
                 normalized=True, **kwargs):

    if fitParameters and beadTracking:
        data, beadData, fit = process_file(fileName, calName, beadTracking=beadTracking, fitParameters=fitParameters,
                                    mode=mode, **kwargs)
    elif fitParameters is False and beadTracking:

        data, beadData = process_file(fileName, calName, beadTracking=beadTracking, fitParameters=fitParameters,
                                    mode=mode, **kwargs)
    elif fitParameters and beadTracking is False:
        data, fit = process_file(fileName, calName, beadTracking=beadTracking, fitParameters=fitParameters,
                                    mode=mode, **kwargs)

    else:
        data = process_file(fileName, calName, beadTracking=beadTracking, fitParameteres=fitParameters,
                            normalized=normalized, **kwargs)
    if select_data and beadTracking is False:
        data = data_selection(data)
    elif select_data and beadTracking:
        beadData = data_selection(beadData)
    if beadTracking:
        regions = assign_regions(beadData)
    else:
        regions = assign_regions(data)
    if unfoldingEvents:
        events = find_unfolding_events(regions, unfoldingThreshold=unfoldingThreshold, unfoldingWindow=unfoldingWindow,
                            plot=plotEvents)
        return regions, events
    else:
        return regions