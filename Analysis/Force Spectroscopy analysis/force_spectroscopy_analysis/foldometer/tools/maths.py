#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from scipy.signal import savgol_filter


def column_derivative(columnY, columnX, smooth=True, returnAll=False, **kwargs):
    """
    Function to calculate the derivative of a certain columnY of the data with respect to a columnX

    Args:
        columnY (pandas.Series): column of the independent variable column
        columnX (pandas.Series): column dependent variable column
        smooth (bool): if True, a Savitzky-Golay filter is applied before the computation
        returnAll (bool): if True, a DataFrame is returned with the columnX, columnY and its derivative
        kwargs (dict): keyword arguments for the savgol_filter (most important "window_length" and "polyorder")
            (http://docs.scipy.org/doc/scipy-dev/reference/generated/scipy.signal.savgol_filter.html#scipy.signal.savgol_filter)
    Returns:
        derivative (pandas.DataFrame): series with the derivative of the columnY with respect to the columnX
    """

    if not isinstance(columnY, pd.core.series.Series) or not isinstance(columnX, pd.core.series.Series):
        raise TypeError("The function needs a label for the column, "
                        "please convert the input to a pandas.Series and label it")

    # check if arguments for the Svitzky-Golay filter have been passed
    if not kwargs:
        kwargs = {"window_length": 551, "polyorder": 2}

    if smooth:
        derivativeColumn = pd.Series(savgol_filter(columnY, **kwargs)).diff().fillna(method="backfill") / \
                           pd.Series(savgol_filter(columnX, 101, 3)).diff().fillna(method="backfill")
        derivativeColumn = derivativeColumn.fillna(derivativeColumn.iloc[1])

    else:
        derivativeColumn = columnY.diff().fillna(method="backfill") / columnX.diff().fillna(method="backfill")
        derivativeColumn = derivativeColumn.fillna(derivativeColumn.iloc[1])

    if returnAll:
        derivative = pd.DataFrame({columnX.name: columnX, columnY.name: columnY,
                                   (columnY.name + "Derivative" + columnX.name): derivativeColumn})
    else:
        derivative = derivativeColumn
        derivative.name = columnY.name + "Derivative" + (columnX.name[0].upper() + columnX.name[1:])

    return derivative


def include_derivative(data, columnY, columnX="time", smooth=True, force=False, **kwargs):
    """
    Function to calculate and append the derivative of a certain columnY of the data with respect to a columnX

    Args:
        data (pandas.DataFrame): DataFrame containing the foldometer data
        columnY (str): label of the independent variable column
        columnX (str): label dependent variable column. Default is "time"
        smooth (bool): if True, a Savitzky-Golay filter is applied before the computation
        kwargs (dict): keyword arguments for the savgol_filter (most important "window_length" and "polyorder")
            (http://docs.scipy.org/doc/scipy-dev/reference/generated/scipy.signal.savgol_filter.html#scipy.signal.savgol_filter)

    Returns:
        data (pandas.DataFrame): the same data including the column with the calculated derivative

    Example:
        >>> import pandas as pd
        >>> import numpy as np
        >>> distance = np.random.random(10)
        >>> time = np.arange(0, 10)
        >>> data = pd.DataFrame({"distance" : distance, "time": time})
        >>> data = include_derivative(data, columnY="distance", columnX="time")

    """

    # check if the column has already been added
    nameDerivativeColumn = columnY + "Derivative" + (columnX[0].upper() + columnX[1:])
    if nameDerivativeColumn in data.columns.values.tolist() and force is False:
        raise ValueError("The derivative is already calculated and included in the DataFrame. If you want to overwrite,"
                         " set force=True")

    data[nameDerivativeColumn] = column_derivative(data[columnY], data[columnX], smooth=smooth, **kwargs).values

    return data


def cross_correlation(x, y, maxLag=None, lagStep=1):
    """
    Compute the correlation coefficients of two one dimensional arrays with a given maximum lag.\
    The input vectors are converted to numpy arrays.
    http://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient
    Or see: Chatfield, 2004, The analysis of time series

    Args:
        x (array_like): one dimensional input array
        y (array_like): one dimensional input array
        maxLag (int): maximum lag (default is len(x)//4)
        lagStep (int): number of data points between two consecutive lags

    Returns:
        res (numpy.array): A one dimensional vector of length maxLag + 1 (to account for zero lag).
    """

    if maxLag is None:
        maxLag = len(x) // 4

    # convert input to numpy array for speed
    x = np.array(x)
    y = np.array(y)

    # prepare output array
    xcorr = []
    # check length of y
    if len(y) < maxLag + 1:
        raise ValueError('y array is too short.')

    # pre-computing mean and denominator
    mx = np.mean(x)
    my = np.mean(y)
    denominator = np.sqrt(np.sum((x - mx) ** 2) * np.sum((y - my) ** 2))

    # compute correlation coefficient for each lag
    for i in np.arange(0, maxLag + 1, lagStep):
        dx = x[:len(x) - i] - mx
        dy = y[i:] - my
        nom = np.sum(dx * dy)
        xcorr.append(nom / denominator)
    return xcorr


def goodness_fit(y, residuals):
    """
    Calculation of the R^2 factor, or goodness of the fit,
    according to: https://en.wikipedia.org/wiki/Coefficient_of_determination
    Args:
        y:
        residuals:

    Returns:
        R2 (float): R^2 factor, or goodness of the fit
    """

    ssErr = np.sum(np.power(residuals, 2))
    print(ssErr)
    ssTot = np.sum(np.power(y - np.mean(y), 2))
    print(ssTot)
    R2 = 1 - (ssErr / ssTot)
    return R2
