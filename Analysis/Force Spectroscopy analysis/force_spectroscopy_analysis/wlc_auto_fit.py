import foldometer as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib
import streamlit as st
import os
from scipy.signal import savgol_filter
from foldometer.analysis.wlc_curve_fit import wlc_series_accurate
st.set_option('deprecation.showPyplotGlobalUse', False)

folderProject = "D:/projects/test/raw_data"
folderFigure = "D:/projects/test/figure"
folderDataSave = "D:/projects/test/fitted_data"
if not os.path.exists(folderFigure):
    os.mkdir(folderFigure)
if not os.path.exists(folderDataSave):
    os.mkdir(folderDataSave)

names = [name[:-5] for name in os.listdir(folderProject) if (not "Power Spectrum" in name) and (not ".tdms_index" in name)]
name = st.sidebar.selectbox("Select data", names)

rulersDefault = [0,100]
rulersStr = st.sidebar.text_input("rulers:", str(rulersDefault))
rulers = [float(r) for r in rulersStr.replace("]","").replace("[","").replace(" ","").split(",")]

contourLengthProtein = st.sidebar.number_input("contourLengthProtein",None,None, 100.)
contourLengthDNA = st.sidebar.number_input("contourLengthProtein",None,None, 3000.)

filePath = folderProject+"/"+name+".tdms"
calibrationFilePath = folderProject+"/"+name+" Power Spectrum.tdms"

#This function will load the data and read the calibration file
measurement = fm.Folding(filePath, protein="your_protein", setup="CT")

#Select the time interval you want to analyse (remove double tethers and tether breaking)
st.header("Select the time interval you want to analyse (remove double tethers and tether breaking)")
tmin = st.sidebar.number_input("tmin = ",None,None, measurement.data["time"].min())
tmax = st.sidebar.number_input("tmax = ",None,None, measurement.data["time"].max())
maskTime = (measurement.data["time"]>=tmin)*(measurement.data["time"]<=tmax)


measurement.process_data(calibrationPath=calibrationFilePath)
measurement.data = measurement.data.loc[maskTime]

plt.plot(measurement.allData["time"], measurement.allData["forceX"], "b")
plt.plot(measurement.data["time"], measurement.data["forceX"], "r")
st.pyplot()
#
# st.text("time surfaceSepX")
# plt.plot(measurement.data["time"], measurement.data["surfaceSepX"])
# st.pyplot()

measurement.assign_regions()


# #You can change the force offset with (it will SUBTRACT the number):


measurement.data = measurement.data.loc[maskTime]
measurement.find_unfolding_events()
measurement.analyse_data()
measurement.data = measurement.data.loc[maskTime]

tminFit = st.sidebar.number_input("tminFit",None,None,measurement.data["time"].min())
tmaxFit = st.sidebar.number_input("tmaxFit",None,None,measurement.data["time"].max())
maskTimeFit = (measurement.data["time"]>=tminFit)*(measurement.data["time"]<=tmaxFit)
maskForce = measurement.data["forceX"]>=2
measurement.data["wlcRegion"] = maskTimeFit
plt.plot(measurement.data["surfaceSepX"], measurement.data["forceX"], "lightgrey")
plt.plot(measurement.data["surfaceSepX"].loc[maskTimeFit], measurement.data["forceX"].loc[maskTimeFit], "r")
plt.plot(savgol_filter(measurement.data["surfaceSepX"].loc[maskTimeFit*maskForce],51,1), measurement.data["forceX"].loc[maskTimeFit*maskForce], "k", linewidth=0.2)
plt.xlabel("Extension (nm)")
plt.ylabel("Force (pN)")
st.pyplot()

fOffset = st.sidebar.number_input("-force offfset:",None,None, 0.)
measurement.data["forceX"]+=fOffset
maskForce20 = (measurement.data["forceX"]>=19)*(measurement.data["forceX"]<=21)
ext20 = wlc_series_accurate(20.,contourLengthDNA,35,800, 0., 0.75)
measurement.data["surfaceSepX"] += ext20-measurement.data["surfaceSepX"].loc[maskForce20].mean()
measurement.extensionOffset = 0
measurement.fit_wlc(tminFit, tmaxFit, contourLengthProtein=contourLengthProtein, parameterRange=0.5, persistenceLengthDNA=35, persistenceLengthProtein=0.75, stretchModulusDNA=800, contourLengthDNA=contourLengthDNA, minForce=2)
measurement.extensionOffset = 0

maskForce20 = (measurement.data["forceX"]>=19)*(measurement.data["forceX"]<=21)
ext20 = wlc_series_accurate(20.,float(measurement.wlcData[True].params["contourLengthDNA"].value),
                    float(measurement.wlcData[True].params["persistenceLengthDNA"].value),
                    float(measurement.wlcData[True].params["stretchModulusDNA"].value),
                    0., float(measurement.wlcData[True].params["persistenceLengthProtein"].value))
extOffsetEstimate = ext20-measurement.data["surfaceSepX"].loc[maskForce20].mean()
extOffsetCorrection = st.sidebar.number_input("extension offset: ",-100000.,100000., 0.)
measurement.data["surfaceSepX"] += extOffsetEstimate + extOffsetCorrection

measurement.extensionOffset = 0
measurement.fit_wlc(tminFit, tmaxFit, contourLengthProtein=contourLengthProtein, parameterRange=0.5, persistenceLengthDNA=35, persistenceLengthProtein=0.75, stretchModulusDNA=800, contourLengthDNA=contourLengthDNA, minForce=2)
measurement.extensionOffset = 0


# #You can check the fitted parameters with:
# for fit in measurement.wlcData:
#     print(measurement.wlcData[fit].params["contourLengthProtein"])
#     print(measurement.wlcData[fit].params["persistenceLengthDNA"])
#     print(measurement.wlcData[fit].params["stretchModulusDNA"])
#
# #or the averaged (which will be exported) with:
# measurement.paramsWLC
#

st.sidebar.header("WLC parameters")
contourLengthDNA = st.sidebar.number_input("contourLengthDNA"+ str(measurement.wlcData[True].params["contourLengthDNA"].value),None,None, float(measurement.wlcData[True].params["contourLengthDNA"].value))
persistenceLengthDNA = st.sidebar.number_input("persistenceLengthDNA"+ str(measurement.wlcData[True].params["persistenceLengthDNA"].value),None,None, float(measurement.wlcData[True].params["persistenceLengthDNA"].value))
stretchModulusDNA = st.sidebar.number_input("stretchModulusDNA"+ str(measurement.wlcData[True].params["stretchModulusDNA"].value),None,None, float(measurement.wlcData[True].params["stretchModulusDNA"].value))
persistenceLengthProtein = st.sidebar.number_input("persistenceLengthProtein"+ str(measurement.wlcData[True].params["persistenceLengthProtein"].value),None,None, float(measurement.wlcData[True].params["persistenceLengthProtein"].value))
measurement.wlcData[True].params["contourLengthDNA"].value = contourLengthDNA
measurement.wlcData[True].params["persistenceLengthDNA"].value = persistenceLengthDNA
measurement.wlcData[True].params["stretchModulusDNA"].value = stretchModulusDNA
measurement.wlcData[True].params["persistenceLengthProtein"].value = persistenceLengthProtein
measurement.paramsWLC["contourLengthDNA"] = contourLengthDNA
measurement.paramsWLC["persistenceLengthDNA"] = persistenceLengthDNA
measurement.paramsWLC["stretchModulusDNA"] = stretchModulusDNA
measurement.paramsWLC["persistenceLengthProtein"] = persistenceLengthProtein

boolSave = st.sidebar.button("save figures and data")

# #For plotting, there are many arguments to choose what to plot, check documentation
measurement.force_extension_curve(pulls=(0,10), rulers=rulers)
plt.plot(savgol_filter(measurement.data["surfaceSepX"].loc[maskTime],51,1), savgol_filter(measurement.data["forceX"].loc[maskTime],51,1), "k", linewidth=0.2)
if boolSave:
    plt.savefig(folderFigure+"/"+name+"_ExtForce"+".png")
    plt.savefig(folderFigure+"/"+name+"_ExtForce"+".pdf")
st.pyplot()




# #calculate the instantaneous contour length of the protein
measurement.calculate_protein_contour_length(proteinLength=120)

for length in rulers:
    plt.plot([length, length], [0,60], "k")
plt.plot(measurement.data["proteinLc"].loc[maskTime], measurement.data["forceX"].loc[maskTime], "r")
plt.plot(savgol_filter(measurement.data["proteinLc"].loc[maskTime].loc[maskForce],51,1), savgol_filter(measurement.data["forceX"].loc[maskTime].loc[maskForce],51,1), "k", linewidth=0.2)
plt.xlabel("Extension (nm)")
plt.ylabel("Force (pN)")
plt.xlim((rulers[0]-20, rulers[-1]+20))
if boolSave:
    plt.savefig(folderFigure+"/"+name+"_LcForce"+".png")
    plt.savefig(folderFigure+"/"+name+"_LcForce"+".pdf")
st.pyplot()

if boolSave:
    measurement.data.to_csv(folderDataSave+"/"+name+".csv", index=None)
    st.text(measurement.paramsWLC)
    df = pd.DataFrame(measurement.paramsWLC, index=[0])
    df.to_csv(folderDataSave+"/"+name+"_PowerSpectrum"+".csv")
