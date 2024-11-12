"""
BSD 2-Clause License

Copyright (c) 2024, Yasaman Madraki
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

# This code could be used for the analysis of the cell membrane force measurement obtained using three steps. The goal is to export the force corresponding to the last stage of the experiment. 
# This codes reads the force measurement file (.h5) which includes three stage of the experiments (approaching+pausing+detaching)
# This code identifies the last step of the experiment (detaching) and exports the data correlating to this step as an excel file

import numpy as np
import lumicks.pylake as lk
import matplotlib.pyplot as plt
import scipy.optimize

file = lk.File("20211029-155038 Marker Cell_N3_3_trapmove_P5_L15_V3.h5")
#NanoStage_x = file["Trap stage position"]["X"]
TrapPosition_x=file["Trap position"]["1X"]
print(dir(file.force1x.seconds))


f1_x_data=file.force1x.downsampled_by(1000).data
f1_y_data=file.force1y.downsampled_by(1000).data

f1_data=np.sqrt(f1_x_data**2+f1_y_data**2)
f1x_timestamps = file.force1x.downsampled_by(1000).timestamps
f1x_time = file.force1x.downsampled_by(1000).seconds
f1_NanoStage_x=TrapPosition_x.downsampled_by(1000).data
f1_time_nanostage=TrapPosition_x.downsampled_by(1000).seconds
#-------finding the right range of the position of the stage at steady state------
BinRegulaingData=f1_NanoStage_x[len(f1_NanoStage_x)-200:len(f1_NanoStage_x)-10]
Nanostage_position_mean=np.mean(BinRegulaingData)
Nanostage_position_min=np.amin(BinRegulaingData)
Nanostage_position_max=np.amax(BinRegulaingData)
#--------------------------------------------------------------------------------------------------------------------------
#-----------------Defining one singel bin: min and max that is obtained from analyzing last 200 data points of position stage
bins = np.array([Nanostage_position_min, Nanostage_position_max])
inds = np.digitize(f1_NanoStage_x, bins) # This will create 0 & 1 array. 0: data outside of the bin; 1: position inside the bin
print(inds)
target_indices=np.argwhere(inds==1) # just finds the index of data points which were shown as "1"; index of data points that are inside the bin
Data=f1_NanoStage_x[target_indices] # final array of nanostage position 
print(target_indices[0])
time=f1x_time[target_indices]-f1x_time[target_indices[0]]
force=f1_data[target_indices]

plt.plot(time, force)
plt.show()

np.savetxt('09_29_S2_Cell_N3_3_force.csv', force, delimiter=' ', fmt='%s')
np.savetxt('09_29_S2_Cell_N3_3_time.csv', time, delimiter=' ', fmt='%s')


