#!/usr/bin/env python
# -*- coding: utf-8 -*-

from struct import unpack
import numpy as np
from collections import OrderedDict as od
import pandas as pd
import array


CCD_FREQUENCY = 90


def read_measurement_block(dataFile, measurementVariables, variableMask, sampleFrequency):
    """Function to read a measurement data block, always 1 second of data (the number of data points depend on the
    Sample Frequency). For a description of the data contained in this block, refer to the non-updated Manual
    (\\SUPPORTSRV\Projects\BioPhysics\Foldometer\Manuals & Datasheets\FoldometerDataFileFormat.pdf)

    Args:
        dataFile (file): the opened file from which the data is read
        measurementVariables (list): list containing all the variables to be extracted from the file.
            It also contains the variable "time", not read but computed just adding a time step from the very first
            start time of the first block (and only the first one)
        variableMask (list): list containing only the indices of variables in measurementVariables to be appended to
            the dictionary (and DataFrame)
        sampleFrequency (float): data sampling frequency in order to compute the timestamp

    """
    # Define an extend function for increased speed and performance (extend is like appending but for arrays)
    appendFunctions = [variable.extend for variable in measurementVariables]

    # Read the timestamp at which the block recording started
    startTimeFile = unpack('d', dataFile.read(8))[0]
    #Number of samples in the block
    nSamples = unpack('i', dataFile.read(4))[0]
    #Due to a finite precision in the computer time generation, only consider the first timestamp, and calculate the
    #rest using the sample frequency
    if len(measurementVariables[0]) == 0:
        startTimeReal = startTimeFile
    else:
        startTimeReal = measurementVariables[0][-1] + 1 / sampleFrequency

    #Read the whole block, much more efficient than struct (around 10x faster)
    block = array.array('i')
    block.fromfile(dataFile, nSamples*12)
    blockArray = np.asarray(block).reshape((nSamples, 12))
    #Append time
    appendFunctions[0](np.arange(0, nSamples)/sampleFrequency + startTimeReal)
    #Append the rest of variables depending on the mode chosen, subtract 1 to the index because time is not in the file
    for variableIndex in variableMask:
        appendFunctions[variableIndex](tuple(blockArray[:, variableIndex-1]))


def read_bead_block(dataFile, beadVariables):
    """Function to read a bead data block, containing the information for a single frame
    For a description of the data contained in this block, refer to the Manual
    (\\SUPPORTSRV\Projects\BioPhysics\Foldometer\Manuals & Datasheets\FoldometerDataFileFormat.pdf)

    Args:
        dataFile (file): the opened file from which the data is read
        beadVariables (list): list containing all the variables regarding the bead
            to be extracted from the file. It also contains the variable "time", not read but computed

    """

    startTimeFile = unpack('d', dataFile.read(8))
    beadVariables[1].append(startTimeFile[0])
    beadVariables[0].append(unpack('Q', dataFile.read(8))[0])  # Frame number

    for variable in beadVariables[2:]:
        variable.append(unpack('d', dataFile.read(8))[0])


def read_header(fileObject, version="1.9"):
    """
    Function to read the header of a binary file from Foldometer software.
    For a description of the data contained in this kind of file, refer to the Manual
    (\\SUPPORTSRV\Projects\BioPhysics\Foldometer\Manuals & Datasheets\FoldometerDataFileFormat.pdf)

    Args:
        fileObject(file): opened object from wehere toe xtract the header

    Returns:
        header (dict): dictionary with the metadata of the file
    """

    fileVersion = unpack('i', fileObject.read(4))[0]
    strLen = unpack('B', fileObject.read(1))
    fileName = unpack(str(strLen[0]) + 's', fileObject.read(strLen[0]))[0]
    strLen = unpack('B', fileObject.read(1))
    sampleFreqString = unpack(str(strLen[0]) + 's', fileObject.read(strLen[0]))[0]
    sampleFreq = unpack('d', fileObject.read(8))[0]
    strLen = unpack('B', fileObject.read(1))
    userName = unpack(str(strLen[0]) + 's', fileObject.read(strLen[0]))[0]
    strLen = unpack('B', fileObject.read(1))
    machineName = unpack(str(strLen[0]) + 's', fileObject.read(strLen[0]))[0]
    date = unpack('Q', fileObject.read(8))[0]

    header = {"fileVersion": fileVersion, "fileName": fileName, "sampleFreqString": sampleFreqString,
              "sampleFreq": sampleFreq, "userName": userName, "machineName": machineName, "date": date}
    calibrationParametersLabels = ["beadRadius", "temperature", "viscosity", "diffCoeffTheo", "dragCoeffTheo"]

    if "TestingDataFile" not in str(fileName):
        if fileVersion > 26:
            #date is only there after certain version
            strLen = unpack('B', fileObject.read(1))
            date = unpack(str(strLen[0]) + 's', fileObject.read(strLen[0]))[0]
            header["date"] = date
        strLen = unpack('B', fileObject.read(1))
        protein = unpack(str(strLen[0]) + 's', fileObject.read(strLen[0]))[0]
        strLen = unpack('B', fileObject.read(1))
        #strLen=(308,)
        #print(protein)
        #print(strLen)
        condition = unpack(str(strLen[0]) + 's', fileObject.read(strLen[0]))[0]

        header["protein"] = protein
        header["condition"] = condition
        #print(condition)
        calibrationParametersLabels = ["temperature", "viscosity", "DrivingFrequency", "DrivingAmplitudeX","DrivingAmplitudeY"]

    for parameter in calibrationParametersLabels:
        header[parameter] = unpack('d', fileObject.read(8))[0]

    if "TestingDataFile" in str(fileName):
        header["beadRadius1"] = 1e9 * header["beadRadius"]
        header["beadRadius2"] = 1e9 * header["beadRadius"]
    header["viscosity"] *= 1e-6
    header["temperature"] -= 273

    return header


def read_calibration_fit_values(fileObject, header):
    """
    Function to read the thermal calibration fit values performed live of a binary file from Foldometer software.
    For a description of the data contained in this kind of file, refer to the Manual
    (\\SUPPORTSRV\Projects\BioPhysics\Foldometer\Manuals & Datasheets\FoldometerDataFileFormat.pdf)

    Args:
        fileObject(file): opened object from wehere toe xtract the header

    Returns:
        calibrationFit (pandas.DataFrame): DataFrame with the calibration fit values performed by Foldometer live

    """
    if "TestingDataFile" not in str(header["fileName"]):
        calibrationFitLabels = ["alphaSinusoidal", "beta", "calibrationQuality", "cornerFrequency",
                                "diffCoeffExp", "dragCoeffExp", "forceFactor", "peakHeight", "trapStiffnessSinusoidal",
                                "stiffness", "beadRadius", "dragCoeff", "diffCoeff"]
    else:
        calibrationFitLabels = ["alphaSinusoidal", "beta", "calibrationQuality", "cornerFrequency",
                                "diffCoeffExp", "dragCoeffExp", "drivingAmplitude", "drivingFrequency",
                                "forceFactor", "peakHeight", "trapStiffnessSinusoidal", "stiffness"]
    #1/beta is called in the foldometer alphaStationary
    calibrationFitDict = {key: [] for key in calibrationFitLabels}
    for i in range(4):
        for key in calibrationFitLabels:
            calibrationFitDict[key].append(unpack('d', fileObject.read(8))[0])

    calibrationFit = pd.DataFrame(calibrationFitDict, index=["PSD1x", "PSD1y", "PSD2x", "PSD2y"])
    calibrationFit.ix[:, "stiffness"] = calibrationFit.ix[:, "stiffness"] / 1000
    calibrationFit.ix[:, "beta"] = 0.001 / calibrationFit.ix[:, "beta"]
    if "TestingDataFile" not in str(header["fileName"]):
        calibrationFit.loc[:, "beadRadius"] *= 1e9
        header["beadRadius1"] = calibrationFit.loc["PSD1x", "beadRadius"]
        header["beadRadius2"] = calibrationFit.loc["PSD2x", "beadRadius"]


    offsets = pd.Series([unpack('d', fileObject.read(8))[0] for i in range(4)], index=["PSD1x", "PSD1y", "PSD2x",
                                                                                       "PSD2y"])
    offsets.name = "offset"

    calibrationFit = calibrationFit.join(offsets)

    return calibrationFit


def read_file(fileName, mode=None):
    """
    Function to read a binary file from Foldometer software.
    For a description of the data contained in this kind of file, refer to the Manual
    (\\SUPPORTSRV\Projects\BioPhysics\Foldometer\Manuals & Datasheets\FoldometerDataFileFormat.pdf)

    Args:
        fileName(str): name of the file to be opened (including the path)
        mode (str): None will read the whole file, "calibration" only the 4 Diff columns and "compact[AXIS]" the
        2 Diff columns corresponding to that axis. Default is None

    :returns: * **header** (dict) -- dictionary with the metadata of the file
              * **fitParameters** (pandas.DataFrame) --  parameters from the fit performed in Foldometer
              * **data** (pandas.DataFrame) --  information from the PSDs
              * **beadTrack** (pandas.DataFrame) --  spatial coordinates of both beads

    """

    with open(fileName, 'rb') as f:

        header = read_header(f)
        calibrationFit = read_calibration_fit_values(f, header)


        #Measurement Data Block Variables Definition
        #=======================================================
        measurementVariablesLabels = ["time", "index", "PSD1VxDiff", "PSD1VxSum", "PSD1VyDiff", "PSD1VySum",
                                      "PSD2VxDiff", "PSD2VxSum", "PSD2VyDiff", "PSD2VySum", "MirrorX",
                                      "MirrorY", "Status"]

        if mode == "calibration":
            finalVariablesLabels = ["time", "PSD1VxDiff", "PSD1VyDiff", "PSD2VxDiff", "PSD2VyDiff"]
        elif mode == "compactX":
            finalVariablesLabels = ["time", "PSD1VxDiff", "PSD2VxDiff", "MirrorX"]
        else:
            finalVariablesLabels = measurementVariablesLabels
        varDict = od([(variable, []) for variable in measurementVariablesLabels])
        measurementVariables = list(varDict.values())  #for performance improvement in the iteration
        #Create a mask with the columns to append, in order to save time
        variableMask = [value for value, key in enumerate(measurementVariablesLabels) if key in
                        finalVariablesLabels[1:]]

        #=======================================================

        #Bead Data Block Variables Definition
        #=======================================================
        beadVariablesLabels = ["frame", "time", "Bead1X", "Bead1Y", "Bead2X", "Bead2Y"]
        beadDict = od([(variable, []) for variable in beadVariablesLabels])
        beadVariables = list(beadDict.values())
        #=======================================================
        #print(header)
        #print(calibrationFit)

        #Read the file, create the header, data and bead tracking dataframes
        #===================================================================


        checkByte = f.read(4)
        while checkByte != b'':
            #print(checkByte)
            blockType = unpack('i', checkByte)[0]
            if blockType == 17:
                read_measurement_block(f, measurementVariables, variableMask, header["sampleFreq"])

            else:
                read_bead_block(f, beadVariables)

            checkByte = f.read(4)

    finalDict = od([(key, value) for key, value in varDict.items() if key in finalVariablesLabels])

    data = pd.DataFrame(finalDict, index=varDict["time"])
    beadTrack = pd.DataFrame(beadDict, index=beadDict["time"])
    #===================================================================

    return header, calibrationFit, data, beadTrack


def binary_to_csv(fileName, newFileName):
    """
    Function to convert a binary file to csv

    Args:
        fileName(str): name of the file to be opened (including the path)
        newFileName (str): name of the file where to save the data (including the path)

    Returns:
        data (pandas.DataFrame):  information from the PSDs

    """
    read_file(fileName)[2].to_csv(newFileName)