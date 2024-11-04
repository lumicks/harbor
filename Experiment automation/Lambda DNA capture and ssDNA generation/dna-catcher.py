
# Catch beads, catch lambda DNA, generate ssDNA
# for bluelake 1.7

# Assume configuration:
# CH1 = 4.2um bead
# CH2 = lambda DNA
# CH3 = buffer
# waypoint 'catch DNA' ~7.75um distance between traps
# waypoint 'stretched' ~17.75um distance between traps

#from bluelake import Trap, stage, fluidics, pause, reset_force, timeline, traps
from bluelake import mirror1, microstage, fluidics,pause,reset_force, timeline, shutters
import numpy as np

trap = mirror1
trap2 = timeline["Force LF"]["Trap 2"]

#user controls =============================
echannel=0 #change this to 1 if you want to end in CH4
dnadwell=1

 #how long to dwell in dna channel (sec)
makessDNA=0 #if just want to catch lambda, change to 1

#functions =================================
def pingpongforce3(): #false positive success threshold
    trap.move_by(dx=5,speed=6)
    f1=trap2.latest_value
    trap.move_by(dx=5.5,speed=5)
    f2=trap2.latest_value
    trap.move_by(dx=-10.5,speed=5)
    print(". stretch: "+str(round(f1,2))+" to "+str(round(f2,2))+"pN")
    if f2>(2*f1) and f2>20:
        return 1
    else:
        return 0

def gohome():
    trap.move_to(waypoint="catch dna",speed=7)

pt=0.1 #pause time for pressure changes
def setpressure(pres):
    curpres=fluidics.pressure
    if curpres<pres:
        while fluidics.pressure < pres:
            fluidics.increase_pressure()
            pause(pt) #important!
    else:
        while fluidics.pressure > pres:
            fluidics.decrease_pressure()
            pause(pt) #important!

def movetoch4():
    fluidics.stop_flow()
    microstage.move_to("J1")
    microstage.move_to("Ch1")
    setpressure(0.1)
    fluidics.open(1,2,3,4,6)
    reset_force()

def beadtest():
    bead1 = timeline['Tracking Match Score']['Bead 1']
    bead2 = timeline['Tracking Match Score']['Bead 2']
    bead_scores = [bead1, bead2]
    if all(bead.latest_value > 0.85 for bead in bead_scores):
        print("you already have beads")
        return 1
    else:
        print("gotta get beads")
        return 0

def catch_beads(min_score = 30) : #from Zsombor
    microstage.move_to("beads")
    shutters.clear(1,2)

    bead1 = timeline['Tracking Match Score']['Bead 1']
    bead2 = timeline['Tracking Match Score']['Bead 2']
    bead_scores = [bead1, bead2]

    while any(bead.latest_value < min_score for bead in bead_scores):
        if any(0 < bead.latest_value < min_score for bead in bead_scores):
            shutters.clear(1,2)  # bad beads
            print(' bad beads')
            pause(1)
    print("  beads should be caught")

def false_positive(): #1 if it was ok, 0 if false positive
    fluidics.close(1,2,3,4,5,6)
    gohome()
    pause(0.8)
    reset_force()
    trap.move_by(dx=5,speed=6)
    f1=trap2.latest_value
    trap.move_by(dx=5.5,speed=6)
    f2=trap2.latest_value
    trap.move_by(dx=-10.5,speed=6)
    #print("  -FP"+str(round(f1,2))+" to "+str(round(f2,2)))
    if f2>(4*f1):
        return 1
    else:
        print(".   just kidding, need to catch dna again")
        return 0

def themainloop(): #will return 1 if caught, 0 after 10 attempts
    #check for beads first
    bt=beadtest()
    fluidics.close(1,2,3,4,5,6)
    if bt==0:
        setpressure(0.3)
        fluidics.open(1,2,3,6)
        microstage.move_to("beads")
        catch_beads(86)
    else:
        #turn on fluidics
        fluidics.open(1, 2, 3, 6)

    # set up the pressure
    microstage.move_to("buffer")
    setpressure(0.22)

    #set traps at initial position and reset force
    gohome()
    reset_force()
    print("Initialized position for force checking")
    caught=pingpongforce3()
    attempts=1

    #moving into dna channel and starting pingpong
    while caught!=1:
        print(". attempted DNA catching "+str(attempts)+" times")
        microstage.move_to("DNA")
        pause(dnadwell)
        microstage.move_to("buffer")
        caught=pingpongforce3()
        if attempts==5:
            return 0

        else:
            attempts+=1

    print(". hopefully we finished")
    return 1

def thelooploop(): #return 1 if done
    chk=0
    while chk==0:
        chk=themainloop()
        if chk==1:
            chk=false_positive()
        else:
            shutters.clear(1,2)

    #hopefully caught something, move to buffer channel and reset all to check
    microstage.move_to("buffer")
    gohome()
    fluidics.stop_flow()
    reset_force()
    
    #0 for no ssDNA, 1 for good, 2 for failed to pull
    if makessDNA==0:
        print("Time to try for ssDNA")
        s_attempts=0
        s=0
        while s != 1:
            s = ssDNAstretch()
            if s == 0:
                setpressure(0.3)
                fluidics.open(3,6)
                shutters.clear(1,2)
                fluidics.close(3,6)
                return s
            if s==2:
                s_attempts+=1
                continue
    else:
        trap.move_to(waypoint="mid stretch",speed=5)
        pause(1)
        if trap2.latest_value < 10:
            print(". probably broke")
            setpressure(0.2)
            fluidics.open(1,2,3,6)
            shutters.clear(1,2)
            return 0
        trap.move_to(waypoint="dsdna stretch",speed=5)
        return 1

def ssDNAstretch(): #0 if nothing/broke, 1 if good, 2 if unsuccessful pull
    success=0
    trap.move_to(waypoint="dsdna stretch",speed=6)
    trap.move_by(dx=12.5,speed=1)
    pause(3)
    highforce=trap2.latest_value
    trap.move_by(dx=-12.5,speed=1)
    lowforce = trap2.latest_value
    print(". stretch: "+str(round(lowforce,2))+" to "+str(round(highforce,2))+"pN")
    if highforce>50:
        print("   extension ok")
        success = 1
        if lowforce>20:
            print("   looks like its still dsDNA")
            success = 2
    if success==0:
        print(".   ssDNA broke ;-;")
        fluidics.open(1,2,3,6)
    if success==1:
        print(".   ssDNA done !!!!")
    return success


def bigbird():
    bb = """        )_)_)__)_)_)_)__)_.
         )_))_))_ _))_ ))) )
          )_))  )~)  )_))~)  )
           )  ) )  )_) )  ) ) )  )
           )  )))  ) )  ) )~_)   ))
           | _(~   ((~ C(C( (  _(~(
          (~  (    C(   ((~ ((C(~ C(
         _(  _(   (~   (  _=~  _C(~
         (   (   (        _C=(~~
         (_    __   ~~~==(X
         _(_  _@@_  _     ~(C_
       C(~ (_ ~((~   (_      ~(C_
      (~  _C(_       _(         (C_
      (C_C(C_(     _=~            (C_
        ~(  C(                      (C
           _(                        ~(_
           (                    _      (C
           (                    ~C      ~(_
           (                     (_       (C
           (_                     ~C       ~(
            (_                      (_      ~(
    """
    print(bb)




#THIS IS THE START OF THE MAIN =========================
finished = 0
while finished == 0:
    finished = thelooploop()

#if want to move to ch4 to ready for experiment start
if echannel==1:
    movetoch4()
    trap.move_to(waypoint="stretched",speed=3)
    print('ready to image!')

bigbird()
