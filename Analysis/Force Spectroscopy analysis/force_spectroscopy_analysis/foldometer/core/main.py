#!/usr/bin/env python
# -*- coding: utf-8 -*-

from foldometer.ixo.binary import read_file, read_header, read_calibration_fit_values
from foldometer.analysis.thermal_calibration import calibration_file, calibration_data
from foldometer.analysis.region_classification import assign_regions
from foldometer.analysis.event_classification import find_unfolding_events
from foldometer.analysis.noise_reduction import correct_signal_noise
from foldometer.analysis.fluorescence import plot_contour_fluorescence, align_fluorescence_data
from foldometer.ixo.data_conversion import process_data, data_beadData_merging
from foldometer.ixo.old_setup import read_file_old_setup
from foldometer.ixo.lumicks_c_trap import *
from foldometer.tools.misc import data_selection, column_indexer, data_deletion
from foldometer.analysis.tweezers_parameters import MIRRORVOLTDISTANCEFACTOR, get_mirror_values, CCD_PIXEL_NM, PROTEIN_LENGTHS
from foldometer.analysis.wlc_curve_fit import wlc_fit_data, protein_contour_length, protein_contour_length_accurate
from foldometer.tools.plots import force_extension_curve
import gc

from nptdms import TdmsFile
from nptdms import TdmsFile
from scipy import ndimage
from scipy.signal import savgol_filter
from copy import copy, deepcopy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pkg_resources
import simplejson as json
import os
from pprint import pprint

def lazyprop(fn):
    attr_name = '_lazy_' + fn.__name__
    @property
    def _lazyprop(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, fn(self))
        return getattr(self, attr_name)
    return _lazyprop


class Folding(object):
    """
    Main class of the Foldometer package, containing all necessary attributes and automated methods to analyse the
    pulling experiments data for any of the setups at AMOLF (and easily expandable to more setups).

    Args:
        filePath (str): Path to the file containing the data
        protein (str): Protein used for the experiment. In the newer Foldometer version this is in metadata
        condition (dict): Dictionary containing different conditions like chaperones, buffer, temperature...
        setup (str): which setup was used for the data. Either "ST", "DT" or "LT"

    Attributes:
        axis (str): Axis to perform operations. Default when creating the class is X axis
        metadata (dict): Dictionary containing information about bead radius, viscosity and other common parameters
        foldometerCalFit (pandas.DataFrame): Fitting parameters from the Foldometer software routine
        allRawData (pandas.DataFrame): DataFrame containing all raw data generated in the setup
        beadTrack (pandas.DataFrame): DataFrame containing all bead tracking data (empty if tracking was disable)
        thermalCalibration (pandas.DataFrame): DataFrame with the thermal calibration parameters used by the class
    """

    def __init__(self, filePath, protein="unknown", condition={}, setup="DT", fluorescenceFilePath=None):
        self.protein = protein
        self.filePath = filePath
        self.condition = condition
        self.setup = setup
        self.axis = "X"
        self.paramsWLC = {}
        try:
            self.proteinLength = PROTEIN_LENGTHS[self.protein]
        except:
            self.proteinLength = None

        #========New setup data format========
        if self.setup is "DT":
            self.metadata, self.foldometerCalFit, self.allRawData, self.beadTrack = read_file(self.filePath)
            self.rawData = deepcopy(self.allRawData)
            self.thermalCalibration = deepcopy(self.foldometerCalFit)
            self.offset = self.thermalCalibration["offset"]
            if "TestingDataFile" not in filePath:
                #self.thermalCalibration = deepcopy(self.foldometerCalFit)0
                self.protein = self.metadata["protein"].decode('utf-8')
                try:
                    self.condition = json.loads(self.metadata["condition"].decode('utf-8'))
                except:
                    self.condition = self.metadata["condition"].decode('utf-8')
            if self.metadata["fileVersion"] > 40 and len(self.beadTrack["Bead1X"]) > 0:
                self.beadTrack.loc[:, "Bead1X":"Bead2Y"] *= CCD_PIXEL_NM * 1e9

        #========Old setup data format========
        elif self.setup is "ST":
            self.allRawData = read_file_old_setup(self.filePath)
            self.rawData = deepcopy(self.allRawData)
            self.allData = deepcopy(self.allRawData)
            self.data = deepcopy(self.allData)

        #========Lumicks C-Trap data format========
        elif self.setup is "CT":

            self.allRawData = read_file_lumicks(self.filePath)
            self.rawData = deepcopy(self.allRawData)
            self.allData = deepcopy(self.allRawData)
            self.data = deepcopy(self.allData)
            self.metadata = {}
            for spectrumFileName in os.listdir(os.path.split(self.filePath)[0]):
                if spectrumFileName.endswith("Power Spectrum.tdms") and spectrumFileName[9:15] <= os.path.split(
                        self.filePath)[1][9:15]:
                    calibrationFilePath = os.path.join(os.path.split(self.filePath)[0], spectrumFileName)
            self.pythonCalFit = None
            self.foldometerCalFit = extract_calibration_parameters(calibrationFilePath)
            self.foldometerCalFit["beadDiameter"] *= 1000
            self.foldometerCalFit["beta"] = 0.001 / self.foldometerCalFit["distanceResponse"]
            self.thermalCalibration = self.foldometerCalFit
            self.metadata["beadRadius1"] = self.thermalCalibration.loc["PSD1x", "beadDiameter"] / 2
            self.metadata["beadRadius2"] = self.thermalCalibration.loc["PSD2x", "beadDiameter"] / 2
            if fluorescenceFilePath is not None:
                self.fluoData, self.fluoTimeRes = read_fluorescence_file_lumicks(fluorescenceFilePath)

        self.timeInterval = [self.rawData["time"].min(), self.rawData["time"].max()]
        gc.collect()


    def __str__(self):
        return "Foldometer measurement corresponding to the file: %s" % (self.filePath)

    def __repr__(self):
        return "Foldometer measurement %s" % (self.filePath)

    def set_axis(self, axis):
        """
        Set the axis to perform operations. Default when creating the class is X axis

        Args:
            axis (str): axis to perform operations, either "X" or "Y"
        """
        if axis not in ["x", "X", "y", "Y"] and axis is not None:
            raise ValueError("Axis should be 'x'or 'y'")
        self.axis = str.upper(axis)

    # <editor-fold desc="Thermal calibration methods">
    def new_thermal_calibration(self, calibrationPath, **kwargs):
        """
        Perform a thermal calibration using Python routines instead of the Foldometer values. Requires a calibration
        file.

        Args:
            calibrationPath (str): path to the calibration file
            kwargs: keyword arguments for the thermal calibration function (see documentation of
            foldometer.analysis.calibration.calibration_file)
        """
        self.calibrationPath = calibrationPath
        if self.setup is "DT":
            self.pythonCalFit = calibration_file(calibrationPath, **kwargs)
        elif self.setup is "ST":
            raise AttributeError("ERROR: still not implemented for the old setup")
        elif self.setup is "CT":
            calibrationFile = TdmsFile(calibrationPath)
            _calibrationData = pd.DataFrame({"time": calibrationFile.objects["/'Sensor Data'/'Time (ms)'"].data,
                               "PSD1VxDiff": calibrationFile.objects["/'Sensor Data'/'Force Channel 0 (V)'"].data,
                               "PSD1VyDiff": calibrationFile.objects["/'Sensor Data'/'Force Channel 1 (V)'"].data,
                               "PSD2VxDiff": calibrationFile.objects["/'Sensor Data'/'Force Channel 2 (V)'"].data,
                               "PSD2VyDiff": calibrationFile.objects["/'Sensor Data'/'Force Channel 3 (V)'"].data,})
            self.pythonCalFit = calibration_data(_calibrationData, fmCalibration=self.foldometerCalFit, sFreq=50000,
                                                 mode="lstsq_corrected", limits=(10, 15000), **kwargs)


        self.thermalCalibration = self.pythonCalFit

    def choose_thermal_calibration(self, method):
        """
        Choose between Foldometer or python calibrations

        Args:
            method (str): either "python" or "foldometer". If a python calibration was not performed yet,
            an error occurs.
        """
        if method == "foldometer":
            self.thermalCalibration = self.foldometerCalFit
        elif method == "python":
            try:
                self.thermalCalibration = self.pythonCalFit
            except AttributeError:
                print("ERROR: please first perform a thermal calibration with the method 'new_thermal_calibration()'")

    # </editor-fold>

    def remove_drift(self):
        """
        Removes the drift introduced by the machine in long term
        """
        if self.setup is "DT":
            STATIONARYMIRRORX, STATIONARYMIRRORY = get_mirror_values(self.metadata)
            pd.set_option('mode.chained_assignment', None)

            trapSeparationX = (-STATIONARYMIRRORX + self.rawData["MirrorX"]) * MIRRORVOLTDISTANCEFACTOR

            self.rawData["trapSepX"] = trapSeparationX
            self.rawData = assign_regions(self.rawData, verbose=False)
            self.averageData = data_selection(self.rawData, columnX="MirrorX", columnY="PSD1VxDiff")

            grouped1 = self.averageData.groupby(["region", "pullingCycle"])["PSD1VxDiff"]
            grouped2 = self.averageData.groupby(["region", "pullingCycle"])["PSD2VxDiff"]

            plt.plot(self.rawData["trapSepX"], (self.rawData["PSD1VxDiff"] - self.rawData["PSD2VxDiff"])/2, color="blue")
            for region in ["pulling", "retracting"]:
                for pullingCycle in self.averageData.loc[self.averageData["region"]==region, "pullingCycle"].unique():
                    mask = (self.rawData["region"] == region) & (self.rawData["pullingCycle"] == pullingCycle)
                    firstPoint = self.averageData.loc[self.averageData["region"]==region, "pullingCycle"].unique()[0]
                    diff1 = grouped1.get_group((region, pullingCycle)).mean() - grouped1.get_group((region,
                                                                                                    firstPoint)).mean()
                    self.rawData.loc[mask, "PSD1VxDiff"] -= diff1
                    diff2 = grouped2.get_group((region, pullingCycle)).mean() - grouped2.get_group((region,
                                                                                                    firstPoint)).mean()
                    self.rawData.loc[mask, "PSD2VxDiff"] -= diff2

            #print(self.rawData.loc[self.averageData.index, "PSD1VxDiff"].groupby("region").mean())
            self.rawData.loc[self.rawData["region"]=="retracting", "PSD1VxDiff"] -= \
                self.rawData.loc[self.averageData.index, :].groupby("region").get_group("retracting").PSD1VxDiff.mean()-\
                    self.rawData.loc[self.averageData.index, :].groupby("region").get_group("pulling").PSD1VxDiff.mean()
            self.rawData.loc[self.rawData["region"]=="retracting", "PSD2VxDiff"] -= \
                self.rawData.loc[self.averageData.index, :].groupby("region").get_group("retracting").PSD2VxDiff.mean()-\
                    self.rawData.loc[self.averageData.index, :].groupby("region").get_group("pulling").PSD2VxDiff.mean()

            plt.plot(self.rawData["trapSepX"], (self.rawData["PSD1VxDiff"] - self.rawData["PSD2VxDiff"])/2, color="red")
            plt.show()

        elif self.setup is "ST":
            self.data = correct_signal_noise(self.data)
            self.data["forceX"] -= self.data["forceX"].min() - 0.2

    def select_data(self, raw=True, subset=False, columnX="time", columnY="PSD1VxDiff"):
        """
        Method to select a portion of either the raw or the processed data

        Args:
            raw (bool): if True, the selection is of the raw data directly read from file. If False, select from
            processed data
            columnX (str): name of the column to display in the X axis
            columnY (str): name of the column to display in the Y axis
        """
        if raw:
            if subset:
                self.rawData = data_selection(self.rawData, columnX=columnX, columnY=columnY)
            else:
                self.rawData = data_selection(self.allRawData, columnX=columnX, columnY=columnY)
        else:
            if columnY not in self.allData.columns:
                columnY = "forceX"
            try:
                if subset:
                    self.data = data_selection(self.data, columnX=columnX, columnY=columnY)
                else:
                    self.data = data_selection(self.allData, columnX=columnX, columnY=columnY)
                    self.unfoldingEvents = find_unfolding_events(self.data, self.axis, plot=False)

                self.timeInterval = [self.data["time"].min(), self.data["time"].max()]
            except AttributeError:
                print("ERROR: please first process your raw data with 'process_data()' or 'analyse_data()' methods")
        gc.collect()

    def delete_data(self, raw=True, columnX="time", columnY="PSD1VxDiff"):
        """
        Method to select and remove a portion of either the raw or the processed data

        Args:
            raw (bool): if True, the selection is of the raw data directly read from file. If False, select from
            processed data
            columnX (str): name of the column to display in the X axis
            columnY (str): name of the column to display in the Y axis
        """

        if raw:
            self.rawData = data_deletion(self.rawData, columnX=columnX, columnY=columnY)
        else:
            if columnY not in self.allData.columns:
                columnY = "forceX"
            try:
                self.allData = data_deletion(self.allData, columnX=columnX, columnY=columnY)
                self.data = deepcopy(self.allData)
                self.unfoldingEvents = find_unfolding_events(self.data, self.axis, plot=False)
            except AttributeError:
                print("ERROR: please first process your raw data with 'process_data()' or 'analyse_data()' methods")

    # <editor-fold desc="Data analysis methods">
    def mergeBeadTrack(self):
        """
        Merge and synchronize the PSD and the bead signals. Warning: all the high temporal resolution of the PSDs is
        downgraded to the limiting camera 90 frames per second resolution.
        """
        self.rawData = data_beadData_merging(self.beadTrack, self.rawData)

    def process_data(self, radii=None, beadTracking=False, noiseRemoval=False, calibrationPath=None,):
        """
        Method to convert raw data to forces and distances usign the thermal calibration parameters

        Args:
            radii (tuple): radius of both beads (PSD1 first)
            beadTracking (bool): if True, resample data to include the bead tracking
            noiseRemoval (bool): if True, use file to remove proximity noise
            calibrationPath (str): the path to the thermal calibration file
        """
        self.forceOffset = 0
        self.extensionOffset = 0
        if radii is None:
            radii = (self.metadata["beadRadius1"], self.metadata["beadRadius2"])

        if self.setup == "DT":
            self.allData = process_data(self.rawData, self.offset, self.thermalCalibration, self.beadTrack, self.metadata,
                                        beadTracking=beadTracking, radii=radii, noiseRemoval=noiseRemoval)

        elif self.setup == "CT":
            _rawDataCopy = self.rawData.copy()
            self.allData = process_lumicks_data(_rawDataCopy, self.foldometerCalFit, self.pythonCalFit)
            self.metadata["sampleFreq"] = 1 / self.allData["time"].diff().mean()
            print(self.metadata)
            del _rawDataCopy


        self.data = deepcopy((self.allData))
        gc.collect()


    def assign_regions(self):
        """
        Method to find and assign pulling, retracting and stationary regions
        """
        try:
            self.data = assign_regions(self.data)
        except AttributeError:
                print("ERROR: please first process your raw data with 'process_data()' or 'analyse_data()' methods")

    def find_unfolding_events(self, **kwargs):
        """
        Method to find unfolding events
        """
        try:
            self.unfoldingEvents = find_unfolding_events(self.data, self.axis, **kwargs)
        except AttributeError:
                print("ERROR: please first find regions of the data with 'assign_regions()' method")

    def analyse_data(self, beadTracking=False, radii=None, noiseRemoval=False, calibrationPath=None,
                     windowFactor=20, **kwargs):
        """
        Method to directly assign regions and find unfolding events

        Args:
            beadTracking (bool): if True, the bead data is merged with the PSD data
            radii (tuple): the radius of each of the beads. Default the values in the metadata of the file
            noiseRemoval (bool): if True, take a base noise dataset and subtract it from the data
            **kwargs: Keyword arguments for the find_unfolding_events() force_extension_curve() plotting functions
        """
        self.forceOffset = 0
        self.extensionOffset = 0
        if radii is None:
            radii = (self.metadata["beadRadius1"], self.metadata["beadRadius2"])

        if self.setup == "DT":
            self.allData = process_data(self.rawData, self.offset, self.thermalCalibration, self.beadTrack, self.metadata,
                                        beadTracking=beadTracking, radii=radii, noiseRemoval=noiseRemoval)

        elif self.setup == "CT":
            _rawDataCopy = self.rawData.copy()
            self.allData = process_lumicks_data(_rawDataCopy, self.foldometerCalFit, self.pythonCalFit)
            self.metadata["sampleFreq"] = 1 / self.allData["time"].diff().mean()
            print(self.metadata)
            del _rawDataCopy

        self.allData = assign_regions(self.allData)

        #Use a 50 ms window to identify events
        window = int(self.metadata["sampleFreq"] // windowFactor)
        #print(window)
        # self.unfoldingEvents = find_unfolding_events(self.data, self.axis, rollingWindow = window, **kwargs)

        self.data = copy(self.allData)
        if self.unfoldingEvents is not None:
            self.refoldingRate = len(self.unfoldingEvents["pullingCycle"].unique()) / \
                                   len(self.data.loc[self.data["region"] == "pulling", "pullingCycle"].unique())
            self.pullingCycles = len(self.unfoldingEvents["pullingCycle"].unique())
        gc.collect()

    def remove_force_offset(self, forceOffset=None, forceChannel="force", axis="X"):
        """
        Shift the data in the force channel by substracting forceOffset

        Args:
            forceOffset (float): amount to be substracted to the force daata
            forceChannel (str): either "force", "PSD1Force" or "PSD2Force"
            axis (str): either "x" or "y"
        """
        forceChannel += axis
        if forceOffset is None:
            self.data[forceChannel] -= self.forceOffset
        else:
            self.data[forceChannel] -= forceOffset
        self.forceOffset += forceOffset

    def to_timedelta(self):
        """
        Converts the index to timedelta
        """
        self.data.index = pd.to_timedelta(self.data.index, unit='s')

    def to_time(self):
        """
        Converts the index to timedelta
        """
        self.data.index = self.data["time"]

    def fit_wlc(self, tmin, tmax, recalculateExtensions=False, **kwargs):
        """
        Method to fit the different pulling curves to the WLC using a series of two WLC models by Odjik (1995)

        Args:
            recalculateExtensions (bool): add the contour length changes to unfoldingEvents
            kwargs: keyword arguments for the fitting function, see fm.analysis.wlc_curve_fit.wlc_fit_data for info
        """

        def recalculate_extension_change(pullingCycles):
            """
            Based on the fitted WLC regions, estimate the change in contour length of each unfolding event
            Args:
                pullingCycles (list): list with the pulling cycles to consider
            """
            extensions = []
            extensionErrors = []
            for pullingCycle in pullingCycles.unique():
                wlcRegions = [region for region in self.data.loc[self.data["pullingCycle"]==pullingCycle,
                                                              "wlcRegion"].unique() if region != -1]
                for wlcRegion in np.arange(len(wlcRegions) - 1):
                    if wlcRegions[wlcRegion] in self.wlcData.keys() and wlcRegions[wlcRegion + 1] in \
                            self.wlcData.keys():
                        param1 = self.wlcData[wlcRegions[wlcRegion + 1]].params["contourLengthProtein"]
                        param0 = self.wlcData[wlcRegions[wlcRegion]].params["contourLengthProtein"]
                        extensions.append(abs(param1.value - param0.value))
                        extensionErrors.append(np.sqrt((param1.stderr)**2 + (param0.stderr)**2))
                    else:
                        extensions.append(np.nan)
                        extensionErrors.append(np.nan)


            return extensions, extensionErrors

        temporal = self.data["surfaceSepX"].mean()
        self.data, self.wlcData = wlc_fit_data(self.data, tmin, tmax, protein=self.protein, **kwargs)
        self.extensionOffset += temporal - self.data["surfaceSepX"].mean() #accumulate extension offsets


        if recalculateExtensions:
            self.unfoldingEvents["extensionChange"], self.unfoldingEvents["extensionChangeError"] = \
                recalculate_extension_change(self.unfoldingEvents["pullingCycle"])
        Lc = []
        LcErrors = []
        for key in self.wlcData:
            Lc.append(self.wlcData[key].params["contourLengthProtein"].value)
            LcErrors.append(self.wlcData[key].params["contourLengthProtein"].stderr)
        self.contourLengths = pd.DataFrame(data = {"wlcRegion": list(self.wlcData.keys()), "contourLength" : Lc,
                                                   "contourLengthError": LcErrors})
        pullingCycles = pd.DataFrame(self.data.groupby("wlcRegion")["pullingCycle"].max())
        pullingCycles["wlcRegion"] = pullingCycles.index
        pullingCycles.index = pullingCycles.index.values
        self.contourLengths = self.contourLengths.merge(pullingCycles, how="inner", on="wlcRegion")

        self.contourLengths[["maximumForce", "region"]] = self.data.groupby(["wlcRegion"])["forceX", "region"].max()[0:]

        labelsWLC = ["persistenceLengthProtein", "persistenceLengthDNA", "contourLengthDNA", "stretchModulusDNA"]
        self.paramsWLC = {}
        for paramLabel in labelsWLC:
            self.paramsWLC[paramLabel] = np.mean([self.wlcData[fit].params[paramLabel].value
                                                  for fit in self.wlcData])
            kwargs[paramLabel] = self.paramsWLC[paramLabel]
        gc.collect()

    def calculate_protein_contour_length(self, proteinLength=None, accurate=True, **kwargs):
        """
        Function to calculate the protein contour length from the force and extension data

        Args:
            proteinLength (float): length of the protein in nanometers
            **kwargs (dict): WLC parameters, if not specified, will try to take averages from WLC fits.
                extension (numpy.array): array containing extemsopm data
                force (numpy.array): array containing force data
                contourLengthDNA (float): contour length of the DNA handles
                persistenceLengthDNA (float): persistence length of the DNA handles
                stretchModulusDNA (float): stretch modulus of the DNA handles
                persistenceLengthProtein (float): persistence length of the protein
                temperature (float): temperature in Celsius

        Returns:
            LcProtein (pandas.Series): array with the protein contour length
        """

        #Check for WLC fitted parameters if they are not specified in kwargs and
        # assign the average of the fits to calculate the contour length
        if not kwargs:
            try:
                labelsWLC = ["persistenceLengthProtein", "persistenceLengthDNA", "contourLengthDNA", "stretchModulusDNA"]
                for paramLabel in labelsWLC:
                    kwargs[paramLabel] = self.paramsWLC[paramLabel]
            except:
                print("Warning: no WLC fit was performed and no WLC parameters passed, falling back to default")
        if proteinLength is not None:
            self.proteinLength = proteinLength
        self.data["proteinLc"] = np.nan
        extensionMask = self.data["surfaceSepX"] > self.proteinLength
        if accurate:
            contour_function = protein_contour_length_accurate
        else:
            contour_function = protein_contour_length
        self.data.loc[extensionMask, "proteinLc"] = contour_function(self.data.loc[extensionMask, "surfaceSepX"],
                                                                           self.data.loc[extensionMask, "forceX"],
                                                                           **kwargs)

    # </editor-fold>


    def force_extension_curve(self, rulers=None, **kwargs):
        """
        Method to plot the force extension curve

        Args:
            **kwargs: Keyword arguments for the main plotting function
        """
        if rulers is not None:
            try:
                force_extension_curve(self.data, self.unfoldingEvents, rulers=rulers, wlcParameters=self.paramsWLC,
                **kwargs)
            except AttributeError:
                print("ERROR: please first perform the WLC fit to save WLC parameters")

        else:
            try:
                force_extension_curve(self.data, self.unfoldingEvents, self.wlcData, **kwargs)

            except:

                try:
                    force_extension_curve(self.data, self.unfoldingEvents, **kwargs)
                except:
                    try:
                        force_extension_curve(self.data, **kwargs)
                    except AttributeError:
                        print("ERROR: please first process your raw data with 'process_data()' or 'analyse_data()' methods")
        gc.collect()

    def align_fluorescence(self, color="532nm", startLine=0, endLine=None, maximumOffset=3, channel="trapSepX", **kwargs):
        """
        Calculate the temporal shift of the fluorescence with respect to the PSD data by comparing the time position of the
        peak in the distance between surfaces.

        Args:
            color (str): either "532nm" or "638nm"
            **kwargs: check fm.ixo.lumicks_c_trap for keyword arguments
        """
        #try:
        self.fluoOffset = align_fluorescence_data(self.data, self.fluoData[color], startLine, endLine,
                                                  self.fluoTimeRes, maximumOffset, channel=channel, **kwargs)

        #except:
        #    print("Error!")
        gc.collect()

    def plot_fluorescence_contour_length(self, color, signalRegion, topChannel="proteinLc", backgroundRegion=None,
                                         plotSignalRegion=False, **kwargs):
        #try:
        plot_contour_fluorescence(self.data, self.fluoData, color, signalRegion, topChannel,
                                  backgroundRegion, plotSignalRegion, self.fluoTimeRes, self.fluoOffset,
                                  proteinLength=self.proteinLength, **kwargs)
        #except:
        #    print("Error: some missing data, calculate contour length and fluorescence offset")
        gc.collect()

    @property
    def unfoldingForce(self):
        print("Unfolding force: ", '{:0.0f}'.format(self.unfoldingEvents.force.mean()), "+/-",
              '{:0.0f}'.format(self.unfoldingEvents.force.std()), "pN")
        return (self.unfoldingEvents.force.mean(), self.unfoldingEvents.force.std())

    @property
    def summary(self):
        summaryDict = {}
        summaryDictLabels = ["filePath", "protein", "condition", "stiffness", "dataRange", "unfoldingEvents",
                             "extensionOffset", "pullingCycles", "refoldingRate", "contourLengths"]
        for label in summaryDictLabels:
            if label in self.__dir__():
                try:
                    summaryDict[label] = json.dumps(self.__getattribute__(label))
                except:
                    summaryDict[label] = self.__getattribute__(label).to_json()


        #summaryDict = {"filePath": self.filePath,
        #               "protein": self.protein,
        #               "condition": self.condition,
        #               "stiffness": self.thermalCalibration["stiffness"].to_json(),
        #               "dataRange": [self.data.index[0], self.data.index[-1]],}
        #try:
        #    summaryDict["events"] = self.unfoldingEvents.to_json()
        #except:
        #    pass
        return summaryDict

    def export_data_parameters(self, path):
        """
        Function to export in a text file some analysis parameters, including the extension offset, the force offset,
        the protein, which time interval was chosen for analysis and the average of the WLC fit parameters
        """

        self.dataParameters = {}
        dataParamsLabels = ["protein", "timeInterval", "forceOffset", "extensionOffset", "paramsWLC", "fluoOffset"]
        for label in dataParamsLabels:
            if label in self.__dir__():
                self.dataParameters[label] = json.dumps(self.__getattribute__(label))

        print(self.dataParameters)

        jsonParameters = json.dumps(self.dataParameters)
        with open(path, 'w') as f:
            f.write(jsonParameters)

    def load_data_parameters(self, parametersFile=None, cropData=False):
        """
        Function to load analysis parameters calculated previously.
        Args:
            parametersFile (str): path to the file containing the analysis parameters. If None, look for same name as data
            cropData (bool): if True, crop the data according to the time interval loaded
        """

        if parametersFile is None:
            if self.setup == "DT":
                parametersFile = str(self.filePath[: -4] + "_parameters.txt")
            elif self.setup == "CT":
                parametersFile = str(self.filePath[: -5] + "_parameters.txt")
        dataParams = json.load(open(parametersFile))
        print("Loaded parameters:")
        pprint(dataParams, indent=4)
        self.protein = json.loads(dataParams["protein"])
        self.loadedTimeInterval = json.loads(dataParams["timeInterval"])
        try:
            self.fluoOffset = json.loads(dataParams["fluoOffset"])
        except:
            pass
        if cropData:
            self.data = self.data.loc[self.loadedTimeInterval[0]: self.loadedTimeInterval[1]]
        self.forceOffset = json.loads(dataParams["forceOffset"])
        self.extensionOffset = json.loads(dataParams["extensionOffset"])
        try:
            self.data["forceX"] -= self.forceOffset
            self.data["surfaceSepX"] -= self.extensionOffset
        except:
            print("Warning: the force and extension offsets were loaded but not applied. First process your data. ")

        self.paramsWLC = json.loads(dataParams["paramsWLC"])

    def to_json(self):
        #print(self.summary)
        jsonSummary = json.dumps(self.summary)
        summaryFile = str(self.filePath[2: 12] + "_" + self.filePath[-7:-4] + "_" + self.protein + "_" + "summary" +
                          ".txt")
        with open(str(summaryFile), 'w') as f:
            f.write(jsonSummary)


