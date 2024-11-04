# BSD 2-Clause License

# Copyright (c) [2020], [Olivia Yang]
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
   # list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright notice,
   # this list of conditions and the following disclaimer in the documentation
   # and/or other materials provided with the distribution.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#======

# imaging automation for RPA
# automated shutoff upon DNA breaks or everything photobleaches

# 1) move trap 1 into 'measure point' position and hang out in the 'protein'
# 2) move back into the buffer channel
# 3) start image
# 4) record starting force and photon count
# 5) if force goes below cutoff, wait n time and stop scan
#    if photon count goes below cutoff, wait n time and stop scan, 
#     then return to the protein channel

from bluelake import pause, timeline, confocal, Trap, stage, fluidics, excitation_lasers
import numpy as np

greencutoff = 11 #photon count cutoff

#%% FUNCTIONS ====================
trap2 = Trap("2","XY")

def autoshutoffkymo(): #both force and photon counts
    greenphoton = timeline["Photon count"]["Green"]
    confocal.start_scan("RPA Kymo")
    
    f1=trap2.current_force
    print("initial force: " + str(round(f1,2)))

    n=0; photonbin=20
    while 1:
        pause(1/78125)
        n+=1
        photonbin += greenphoton.latest_value
        if n%(3000) == 0:
            f2 = trap2.current_force
            print("periodic check #"+str(round(n/3000)))
            print("  photon count: " + str(photonbin))
            print("  force: " + str(round(f2,2)))
            if (photonbin < greencutoff) or (f2 < 0.5*f1):
                pause(7)
                break
            else: 
                photonbin=0
    confocal.abort_scan()
    excitation_lasers.green = 0
    print("!!!!!!!!! IT'S DONE !!!!!!!!!")

def setpressure(pres):
    pt=0.1
    curpres=fluidics.pressure
    if curpres<pres:
        while fluidics.pressure < pres:
            fluidics.increase_pressure()
            pause(pt) #important!
    else:
        while fluidics.pressure > pres:
            fluidics.decrease_pressure()
            pause(pt) #important!

def gogetprotein(hangouttime,flowon,initforce):
    print("getting protein for " + str(hangouttime) + " sec")
    trap = Trap("1", "XY")
    trap.move_to(waypoint="measure point",speed=4)
    stage.move_to("J1")
    stage.move_to("Ch1")
    setpressure(0.12)
    if flowon:
        fluidics.open(1,2,3,4,6)
    pause(hangouttime)
    
    fluidics.close(1,2,3,4,6)
    stage.move_to("J1")
    stage.move_to("buffer")


def stretchitout(stretchtime):
    print("pre stretch time: "+ str(stretchtime) + " sec")
    f1=trap2.current_force
    print(" initial force: " + str(round(f1,2)))
    trap = Trap("1", "XY")
    trap.move_to(waypoint="more streatch",speed=4)
    pause(stretchtime)
    trap.move_to(waypoint="measure point",speed=4)
    pause(2)
    f2=trap2.current_force
    print(" final force: " + str(round(f2,2)))
    if f2<(0.5*f1) and f1 > 1:
        print("it broke")
        return 0
    return 1

def parrot():
    parrot = """            .cc;.  ...  .;c.                     
         .,,cc:cc:lxxxl:ccc:;,.                   
        .lo;...lKKklllookl..cO;                   
      .cl;.,:'.okl;..''.;,..';:.                  
     .:o;;dkd,.ll..,cc::,..,'.;:,.                
     co..lKKKkokl.':lloo;''ol..;dl.               
   .,c;.,xKKKKKKo.':llll;.'oOxl,.cl,.             
   cNo..lKKKKKKKo'';llll;;okKKKl..oNc             
   cNo..lKKKKKKKko;':c:,'lKKKKKo'.oNc             
   cNo..lKKKKKKKKKl.....'dKKKKKxc,l0:             
   .c:'.lKKKKKKKKKk;....lKKKKKKo'.oNc             
     ,:.'oxOKKKKKKKOxxxxOKKKKKKxc,;ol:.           
     ;c..'':oookKKKKKKKKKKKKKKKKKk:.'clc.         
   ,xl'.,oxo;'';oxOKKKKKKKKKKKKKKKOxxl:::;,.      
  .dOc..lKKKkoooookKKKKKKKKKKKKKKKKKKKxl,;ol.     
  cx,';okKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKl..;lc.   
  co..:dddddddddddddddddddddddddddddddddl::',::.  
  co...........................................  """
    print(parrot)



#%% MAIN ========================================
excitation_lasers.green = 0
goodstretch = stretchitout(25)
#goodstretch=1
if goodstretch:
    f1=trap2.current_force
    gogetprotein(10,1,f1) #sec to hang out, 1=flow on, initial force
    pause(1)
    f2 = trap2.current_force
    if f2 > 0.5*f1 or f1 < 2:
        excitation_lasers.green=40
        pause(2)
        autoshutoffkymo()
        parrot()
    else:
        print("something is wrong!!")