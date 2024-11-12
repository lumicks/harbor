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

# This code includes automation of the tether formation from a cell by a trapped bead. 
# In this code we are using a single bead
# 


import bluelake
from bluelake import trap1, trap2, trap12z, nanostage, timeline, power
force_1x=timeline["Force HF"]["Force 1x"]
force_1y=timeline["Force HF"]["Force 1y"]
force_1x_LF=timeline["Force LF"]["Force 1x"]
force_1y_LF=timeline["Force LF"]["Force 1y"]

# information needed for approaching step
approaching_speed=1     #um/s
approaching_step=0.02   #um
max_force=80            #pN

#information needed for detaching step
detaching_speed=5       #um/s
detaching_length=15     #length of the tether (um)
pause=5                 #sec, this is the pause time to let the bead gets attached to the cell
counter=1               #this is the counter of tethers formed from the same cell
Exp_name="Cell_N1_"+str(counter)+"_P"+str(pause)+"_L"+str(detaching_length)+"_V"+str(detaching_speed)


# creating "readme.txt" file to record any important information about the tether formation
file_readme = open('Z:/Data/Yasaman/20210702_SUM_Bead1p87_Fasudil40uM/readme.txt' , 'w')
file_readme.write("Experiment name: "+Exp_name+'\n')
file_readme.write("Movie in: x\n")
file_readme.write("Trapping Laser Power: "+str(power.trapping_laser)+'\n')
file_readme.write("Overall Power: "+str(power.overall_trapping_power)+'\n')
file_readme.write("Trap 1 split: "+str(power.trap1_split)+'\n')
file_readme.write("approaching speed: "+str(approaching_speed)+'\n')
file_readme.write("detaching speed: "+str(detaching_speed)+'\n')
file_readme.write("The length of tether: "+str(detaching_length)+'\n')
cal=bluelake.get_force_calibration()
stiffness_1x=cal["Force 1x"]["kappa (pN/nm)"]
stiffness_1y=cal["Force 1y"]["kappa (pN/nm)"]
file_readme.write("Trap stiffness, 1x: "+str(stiffness_1x)+'\n')
file_readme.write("Trap stiffness, 1y: "+str(stiffness_1y)+'\n')

timeline.mark_begin(Exp_name)   # marker is started

# Step 1: Approaching
# This step might be tricky as the "max_force" is not always constant for different cells and conditions. In my experience you can move the nanostage manually to the position which bead is in contact with the cell and then run the script without this "Approaching" step.
#NOTE: This step could be removed and performed manually
while (abs(trap1.current_force)<max_force):
    nanostage.move_by(dx=-approaching_step, speed=approaching_speed) # This is the command if bead is on the left and cell on the right in the BFV
    #nanostage.move_by(dx=approaching_step, speed=approaching_speed) # This is the command if bead is on the right and cell on the left in the BFV

# Step 2: Pausing    
bluelake.pause(pause)           # pause for a few seconds allowing the bead to attach to the cell
    
# Step 3: Detaching
nanostage.move_by(dx=detaching_length, speed=detaching_speed)       # This is the command if bead is on the left and cell on the right in the BFV
#nanostage.move_by(dx=-detaching_length, speed=detaching_speed)     # This is the command if bead is on the right and cell on the left in the BFV


bluelake.pause(60)              # we record the tether force for 60 seconds


timeline.mark_end()             # marker is ended

file_readme.close()
