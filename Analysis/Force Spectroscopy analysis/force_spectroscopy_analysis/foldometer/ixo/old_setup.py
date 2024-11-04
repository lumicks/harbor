#!/usr/bin/env python
# -*- coding: utf-8 -*-


import pandas as pd
import numpy as np
from scipy import stats
from foldometer.tools.misc import data_selection

def extract_calibration_parameters(commentsFileName, fileName):
    """

    Args:
        commentsFileName (str): path the file containing the comments of the measurements
        fileName (str): path to the file containing the data

    Returns:

    """
    pd.options.mode.chained_assignment = None
    comments = pd.read_csv(commentsFileName, header=None, sep="\t")
    comments.columns = ["date", "id", "comment"]
    comments.loc[:, "comment"] = comments.comment.astype(str)
    commentsCalibration = comments[comments["comment"].str.contains("\*\*k ")]
    commentsCalibration.loc[:, "date"] = commentsCalibration["date"].astype(float)
    values = commentsCalibration["comment"].str.split().values
    for value in values:
        value.pop(0)
    commentsCalibration.loc[:, "comment"] = values
    commentsCalibration["minimum"] = abs(commentsCalibration.loc[:, "date"] - float(fileName[:-4]))
    calibrationValues = commentsCalibration.loc[commentsCalibration["minimum"] == min(commentsCalibration[
                                                                                             "minimum"]),
                                     "comment"].values[-1]
    print(calibrationValues)
    calibrationParameters = {"kx": float(calibrationValues[0]), "ky": float(calibrationValues[1]),
                             "betax": float(calibrationValues[2]), "betay": float(calibrationValues[3])}
    print(calibrationParameters)
    return calibrationParameters


def compute_beta_factor(data):
    slope, intercept, r_value, p_value, std_err = stats.linregress(data["xspt2"] * 0.0615, data["vx"])
    return slope, intercept

def read_file_old_setup(fileName, commentsFileName="comment.txt"):
    """
    Function to read a text file from the old setup

    Args:
        fileName(str): name of the file to be opened (including the path)
        commentsFileName (str): file containing metadata such as trap stiffnes and beta factor
    Return:
        data (pandas.DataFrame): DataFrame containing the whole measured set of data

    """
    data = pd.read_csv(fileName, sep="\t")
    fit = extract_calibration_parameters(commentsFileName, fileName)
    data["xpz"] = - data["xpz"]
    data = data.rename(columns={"xpz": "MirrorX", "t": "time"})

    fit["betax"], intercept = compute_beta_factor(data)
    #fit["betax"] = -2
    print(fit["betax"])
    data["surfaceSepX"] = 1000 * (data.MirrorX + (data.vx) / fit["betax"])

    data["forceX"] = - (data.vx) * fit["kx"] / fit["betax"]
    #data = data.rename(columns={"xpz": "MirrorX", "t": "time"})
    data["trapSepX"] = data["MirrorX"]
    data.index = data.time

    data = data_selection(data)

    return data



