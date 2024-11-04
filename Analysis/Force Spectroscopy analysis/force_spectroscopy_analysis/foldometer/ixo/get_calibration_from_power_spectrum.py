from foldometer.ixo.lumicks_c_trap import extract_calibration_parameters
import os
from scipy import constants
import pickle
import streamlit as st
import pandas as pd

def get_calibration_from_power_spectrum(folderData, fileName, extension=".csv"):
    if extension==".csv":
        calibrationFilePath = folderData + "/" + fileName + " Power Spectrum" + extension
        thermalCalibration = pd.read_csv(calibrationFilePath, index_col=0)
    elif extension==".tdms":
        calibrationFilePath = folderData + "/" + fileName + " Power Spectrum" + extension
        thermalCalibration = extract_calibration_parameters(calibrationFilePath)
        pd.set_option("display.max_rows", None, "display.max_columns", None)
        print(calibrationFilePath)
        print(thermalCalibration)
        thermalCalibration["beadDiameter"] *= 1000
        thermalCalibration["beta"] = 0.001 / thermalCalibration["distanceResponse"]
        thermalCalibration["diffusionTheory"] = ((constants.Boltzmann * (thermalCalibration["temperature"] + 273.15)) \
                                                / (3 * constants.pi * thermalCalibration["viscosity"] * 1000 * thermalCalibration["beadDiameter"] * 10 ** -9))
        thermalCalibration["diffusionExp"] = thermalCalibration["diffusionTheory"] / (
                (thermalCalibration["distanceResponse"] * 10 ** -6) ** 2)
    else:
        print("calibrationPath:" + folderData + "/" + fileName + " Power Spectrum" + extension)
        print("thermalCalibration not defined because extension doesn't match the handled case: extension='.csv' extension='.tdms'")
        thermalCalibration = None
    return thermalCalibration