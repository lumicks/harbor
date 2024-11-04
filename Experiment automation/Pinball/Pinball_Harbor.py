"""
BSD 2-Clause License

Copyright (c) 2022, LUMICKS B.V.
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


### Description

Play pinball on the C-Trap using two optical traps with beads caught.
Developed by: Michael Bugiel

Trap 1, the 'ball', bounces of the walls of its ROM. The user moves Trap 2,
the 'racket', vertically along the right wall of the ROM
and has to prevent Trap 1 from 'escaping' through the right wall.

The game ends when the player misses the ball.

Initialize the game as follows:
1) Trap 2 beads.
2) Choose trapping powers that allow you to move the traps without loosing the beads
3) Place Trap 2 just to the right of the ROM of Trap 1. The beads should not touch.
4) Lock the motion of Trap 2 in x
5) Enable video tracking by drawing the bead templates. Choose a small template, a tight box around the bead is best.
6) Draw a minimal ROI in brightfield and corresponding tracking region, where both beads can be always tracked.
Frame rate should be >35 fps.
7) Set threshold for bead tracking to 75%.
8) Enter bead sizes in Bluelake.

Set parameters of the script
1) Define the rectangle in which Trap 1 is 'trapped' by setting x_max and y_max; move Trap 1 to the max x
value and round down to the decimal for x_max. Similar for y_max.
2) x_min and y_min can best be set to 0.1.
3) Set 'do_init = True', if you run the script for the first time. When started, the script will move Trap 1 along
the right 'wall' and determine the minimal distance between the Traps.
This distance is then used to determine whether the racket hit the ball. This minimal distance is printed and can
be put in as 'minimal_distance_manual' for the next run. When using minimal_distance_manual, set do_init = False.

Start the script.
If do_init=True, the minimal distance between the Traps will be determined.
Next, Trap 1 moves to the start position on the right wall.
After a count down of 5 seconds, Trap 1 will start moving.
Bounce Trap 1 and try to get the next high score!


Full list of parameters that can be tuned

## Input parameters:
x_min, x_max, y_min, y_max:     Floats. Minimal/maximal positions in um of the region in which Trap 1 moves,
                                must be above/below the physical ROM of Trap 1, i.e. don't use
                                0 / the maximal possible position, respectively.
beta_min, beta_max:             Floats. Allowed range of the start angle in deg. Required:
                                0 < beta_min < beta_max < 90
speed_calib:                    Float. Speed in um/s with which Trap 1 moves to the starting position
                                at the beginning and returns to the starting position at the end.
                                Required: smaller than the maximal possible speed of mirror.
speed_reflect:                  Float. Speed in um/s with which Trap 1 moves between the reflection points
                                on the edge of the ROM.
                                Required: smaller than the maximal possible speed of mirror.
save_video:                     Bool. Save a video of the game or not.
path_video:                     String. Path were to save it. The user must ensure that the path is valid.
tilted_blow:                    Bool. Blow the 'ball' with the 'racket' under a tilted angle, which is
                                calculated according to their positions (i.e. beads 1 and 2).
                                If False, the ball is just reflected at the right wall as long as the blow was successful.
do_init:                        Bool. Make an initialization measurement or not.
minimal_distance_manual:        Float. Manual value for the minimal distance between beads 1 and 2 in um.
                                To be used when no initialization is done before the game.
speed_init:                     Float. Speed in um/s of Trap 1 for initialization measurement and movement
                                to start position. Required: smaller than the maximal possible speed of mirror.
t_check:                        Float. Time in s to wait after motion during initialization to make TMS check.
points_per_blow:                Float. Points per successful blow.
var_distance:                   Float. Allowed variation in distance for a successful blow.
"""

import bluelake as bl
import numpy as np
import random as rnd
import math
## Parameters

# Region of motion (ROM), in um
#limits = bl.mirror1.motion_limits      # Read out from Bluelake, only for newer BL versions
#x_max = math.floor(limits.upper.x * 100)/100.0
#y_max = math.floor(limits.upper.y * 100)/100.0
# Manual input
x_max = 40.7
y_max = 28.8
x_min = 0.1
y_min = 0.1

# Range of start angle, in deg
beta_min = 5
beta_max = 70

# Speeds for reflections per level, in um/s
speed_reflect = 25      # Speed for reflections

# Saving video
save_video = False      # Save video? True or False
path_video = "D:/"      # Path, make sure it is valid!

# Tilt blown ball according to mutual position of ball and racket
tilted_blow = True

# Points per successful blow.
points_per_blow = 100

# Make initialization measurement: Checks if beads are trapped, tracked, and well positioned,
# measures their minimal distance
do_init = True
minimal_distance_manual = 4.547          # If no init measurement, use manual value for minimal distance, in um
# Speed for initialization
speed_init = 25         # Speed for going to start position and calibration, in um/s
# Check times for TMS
t_check = 0.2           # in s

# Start position
x0 = x_max              # Start always at "right" side
# Rest (y, angle) is random
y0 = rnd.uniform(y_min, y_max)
beta0 = rnd.uniform(beta_min, beta_max) * np.sign(rnd.uniform(-1, +1)) *np.pi/180   # in rad

trap_numbers = np.array(["1", "2"])
axes = np.array(["x", "y"])

# Some channels that are used
distance = bl.timeline["Distance"]["Distance 1"]
match_score = {}
for trap_number in trap_numbers:
    match_score[trap_number] = bl.timeline["Tracking Match Score"][f"Bead {trap_number}"]

bead_position = {}
for trap_number in trap_numbers:
    bead_position[trap_number] = {}
    for axis in axes:
        bead_position[trap_number][axis] = bl.timeline["Bead position"][f"Bead {trap_number} {axis.upper()}"]

bead_diameter_1 = bl.timeline["Bead diameter"]["Template 1"].latest_value
bead_diameter_2 = bl.timeline["Bead diameter"]["Template 2"].latest_value
var_distance = 0.25*min(bead_diameter_1, bead_diameter_2)  # Allowed variation of distance for successful blow

## Define functions

def deg_to_rad(alpha):
    """Converts an angle alpha in degrees into units of radians"""
    return np.pi/180*alpha

def TMS_check():
    """Check TMS if tracking is enabled and beads are trapped for both traps."""
    for trap_number in trap_numbers:
        if match_score[trap_number].latest_value==0:
            raise ValueError(f"Trap {trap_number}: Tracking not enabled and/or no bead tracked.")    

def countdown(t_up, t_low, t_step=1):
    """Count down, display numbers on screen.
    Input:
        t_up, t_low,    Floats. Start and stop time in s,
                        Required: t_up > t_low
        t_step          Float. Interval time in s.
    """
    times = np.flip(np.arange(t_low, t_up+t_step, t_step))
    for time in times:
        print(time)
        bl.pause(t_step)

def count_digits(float):
    """Counts entered digits after decimal. Be careful how you use that!"""
    try:
        n = len(str(float).split(".")[1])
    except:
        n = 0
    return n
    
def find_position(x, y, x_min, x_max, y_min, y_max):
    """Find position on rectangle, as left, right, top, bottom. 
    Test whether position is on edge of rectangle.
    Input:
        x, y            Floats. Recent position on rectangle.
        x_min, x_max,   Floats. Minimal/maximal position in x/y,
        y_min, y_max    defining the rectangle.
    Output:
        position        String. Describing the recent position on the 
                        edge of the rectangle. Possible values:
                        'left', 'right', 'top', 'bottom'
    """
    # digits after decimal
    n_digits = np.max([count_digits(x_min),
                       count_digits(x_max),
                       count_digits(y_min),
                       count_digits(y_max),
                       ])
    position = "init"
    if round(x, n_digits) not in(x_min, x_max) and round(y,n_digits) not in(y_min, y_max):
        raise ValueError("Point (x,y) not on edge of rectangle!")
    if round(x, n_digits) == x_min:
        position = "left"
    else:
        if round(x, n_digits) == x_max:
            position = "right"
    if round(y,n_digits) == y_min:
        position = "bottom"
    else:
        if round(y,n_digits) == y_max:
            position = "top"
    if position == "init":
        raise ValueError("No correct position determined!")
    return(position)

def find_direction(position, beta):
    """Find direction of outgoing reflected beam on edge of rectangle.
    Input:
        position:       String. Recent position on edge of rectangle.
                        Possible value: 'left', 'right', 'top', 'bottom'
        beta:           Float. Outgoing reflection angle in rad.
    Output:
        new_direction:  Array of 2 strings. Describing the direction
                        of the outgoing beam. Possible values:
                        'left' or 'right' for 1st string,
                        'up' or 'down' for 2nd string.
    """
    new_direction = np.array(["inits", "inits"])
    if position=="right" or (position=="top" and beta>0) or  (position=="bottom" and beta<0):
        new_direction[0] = "left" 

    if position=="left" or (position=="top" and beta<0) or  (position=="bottom" and beta>0):
        new_direction[0] = "right" 

    if position=="top" or (position=="left" and beta<0) or (position=="right" and beta>0):
        new_direction[1] = "down"
        
    if position=="bottom" or (position=="left" and beta>0) or (position=="right" and beta<0):
        new_direction[1] = "up"
    
    if "inits" in new_direction:
        raise ValueError("No correct direction found!")
    return new_direction

## A note on directions and positions:
# In the C-Trap, the coordinate system, used in BF and for the mirror1 ROM is turned upside down.
# Thus, what is "up" (in directions) or "top" (in positions) in the opposite in Bluelake.

def find_next_position(x, y, beta, 
                       x_min, x_max,
                       y_min, y_max):
    """Find next connection point of outgoing reflected beam with rectangle, also calculates next outgoing angle.
    Input:
        x, y            Floats. Recent position on rectangle.
        beta            Float. Outgoing reflection angle in rad.
        x_min, x_max,   Floats. Minimal/maximal position in x/y,
        y_min, y_max    defining the rectangle.
    Output:
        x_next, y_next  Floats. Next position on rectangle.
        beta_next       Float. Next outgoing reflection angle in rad.
    """
    # Inits
    lin_func_paras = np.zeros(2)
    beta_next = np.nan
    # Current position on reactangle as "left", "right", "top", "bottom"
    position = find_position(x, y, x_min, x_max, y_min, y_max)
    
    # Linear function to get new position
    if position=="top" or position=="bottom":
        lin_func_paras[0] = 1/np.tan(beta)
        lin_func_paras[1] = -1*x/np.tan(beta) + y
    else:
        lin_func_paras[0] = np.tan(beta)
        lin_func_paras[1] = -1*x*np.tan(beta) + y
        
    # Next y to test
    if position=="top" or (position=="left" and beta<0) or (position=="right" and beta>0):
        y_test = y_min
    else:
        y_test = y_max
        
    # Next x to test, from linear extrapolation
    x_test = (y_test - lin_func_paras[1]) / lin_func_paras[0]
    
    # Find new direction
    new_direction = find_direction(position, beta)
    
    # Next x
    if new_direction[0] == "left":
        x_next = max(x_min, x_test)
    else:
        x_next = min(x_max, x_test)
        
    # Next y
    y_next = lin_func_paras[0]*x_next + lin_func_paras[1]

    # Next position x,y
    position_next = find_position(x_next, y_next, x_min, x_max, y_min, y_max)

    # Next beta
    # to opposite side
    if position in ("top","bottom") and position_next in ("top","bottom"):
        beta_next = -1*beta
    if position in ("left","right") and position_next in ("left","right"):
        beta_next = -1*beta
    # around corner
    if position in ("top","bottom") and position_next in ("left","right"):
        beta_next = -1*np.sign(beta)*(np.pi/2 - np.abs(beta))
    if position in ("left","right") and position_next in ("top","bottom"):
        beta_next = -1*np.sign(beta)*(np.pi/2 - np.abs(beta))
    # If none of these cases occured, something went wrong!
    if beta_next == np.nan:
        raise ValueError("No correct next angle calculated.")

    return x_next, y_next, beta_next
       
def check_for_good_blow_simple(minimal_distance):
    """Checks if the racket hits the ball by simply checking the bead center distance:
        If the distance at the moment when the ball hits the right side of the
        ROM was below a certain threshold, the blow was successful, i.e. good.
    Input:
        minimal_distance    Float. Threshold for bead center distance for the blow to be successful.
    Output:
        Bool
    """
    return distance.latest_value < minimal_distance

def print_in_box(text):
    """Prints text in a nice box. Don't use '\n' in string 'text'!
    Input:
        text:       A string. Must not contain '\n'.
    """
    if '\n' in text:
        raise ValueError("Do not use linebreaks (backslash+n) in input text.")
    length = len(text)
    print( "   " + "~"*(length+3))
    print(f"  | {text} |")
    print( "   " + "~"*(length+3))    

    
### Start
## Some inits
n_digits = np.max([count_digits(x_min),
                   count_digits(x_max),
                   count_digits(y_min),
                   count_digits(y_max),
                   ])
x = x0
y = y0
beta = beta0
score = 0
nn_blow = 0
good_blow = True

print("WELCOME TO PINBALL ON THE C-TRAP!\n")

if do_init:
    ## Initialization and checks for game
    print("Initialization starts.")
    # Check TMS if tracking is enabled and beads are trapped
    TMS_check()
        
    # Calibration run of trap 1 to get minimal position distance between both traps.
    # Also check if they are placed correctly.
    print("Going to right upper corner for initialization.")
    bl.mirror1.move_to(x=x_max, y=y_min, speed=speed_init)     # Move trap 1 to right upper corner
    bl.pause(t_check)
    # Another TMS check if beads were not lost
    TMS_check()
    
    # Check if trap 2 is right of trap 1
    if bead_position["2"]["x"].latest_value - bead_position["1"]["x"].latest_value < sum((bead_diameter_1, bead_diameter_2))/2:
        raise ValueError("Bead 2 is not right enough from bead 1. Please locate bead 2 right of ROM of mirror 1 by at least the sum of both radii and start again.")
    # Check if trap 2 is too much right of trap 1
    if bead_position["2"]["x"].latest_value - bead_position["1"]["x"].latest_value > 4*sum((bead_diameter_1, bead_diameter_2)):
        raise ValueError("Bead 2 is too right from bead 1. Please locate bead 2 right of ROM of mirror 1 by maximum the sum of both diameters and start again.")
    # Check if trap 2 is too much above ROM
    if bead_position["2"]["y"].latest_value - bead_position["1"]["y"].latest_value < 0:
        raise ValueError("Bead 2 is too much above ROM. Please locate bead 2 in the middle of the ROM of mirror 1 in y and start again.")
    
    # Go down and measure
    print("Going to right lower corner for minimal distance check.")  
    t1 = bl.timeline.current_time  # Starting point of measurement, ensure a high enough frame rate
    bl.mirror1.move_to(x=x_max, y=y_max, speed=speed_init)      # Move trap 1 to right lower corner
    t2 = bl.timeline.current_time  # Ending point of measurement
    bl.pause(t_check)
    # Check if tracking is not undermined by proximity of beads
    for trap_number in trap_numbers:
        match_score_init = match_score[trap_number][t1:t2].data
        if 0 in match_score_init:
            raise ValueError("Track matching score was 0 during motion down. Possible reason: Bead(s) lost or beads too close for simultaneous video tracking. Move bead 2 more to the right and/or use smaller template.")
    # Another TMS check if beads were not lost
    TMS_check()
    # Check if trap 2 is too much below ROM
    if bead_position["2"]["y"].latest_value - bead_position["1"]["y"].latest_value > 0:
        raise ValueError("Bead 2 is too much below ROM. Please locate bead 2 in the middle of the ROM of mirror 1 in y.")
    # Get minimal distance between both beads while trap 1 moved down
    distance_init = distance[t1:t2].data
    # This distance will be used to check if the back-reflected bead 1 ("the ball") is close enough to the player-controlled bead 2 ("the racket") to be blown back.
    minimal_distance = np.min(distance_init)  
    print(f"Minimal distance between beads is {minimal_distance:.3f} um.\nInitialization ended.\n\n")
else:
    minimal_distance = minimal_distance_manual      # If no init measurement, use manual value for minimal distance

## Move to random start position
print("Going to random start position.\n")
bl.mirror1.move_to(x=x0, y=y0, speed=speed_init)

## Start the actual game
print("Let's start!\n")
countdown(t_up=5, t_low=1)
print("BEGIN!\n")

# Basic game: 
# No levels, game just continues until missed blow

# Begin
try:
    if save_video:
        bl.brightfield.start_recording(path_video)
    while good_blow:
        # Calculate next position
        x_next, y_next, beta_next = find_next_position(x, y, beta,
                                                       x_min, x_max,
                                                       y_min, y_max)
        # Go there    
        bl.mirror1.move_to(x = x_next, y = y_next, speed = speed_reflect)
        # Init for next reflection
        x = x_next
        y = y_next
        beta = beta_next
        # Check if this next position is "right"
        if round(x_next,n_digits) == x_max:
            # bl.pause(0.1)       # Wait a bit to measure tracked distances and positions
            # Check for a good blow
            good_blow = check_for_good_blow_simple(minimal_distance + var_distance)
            if good_blow:
                score = score + points_per_blow
                print(f"+ {points_per_blow}, total score = {score}")
                nn_blow += 1
            else:
                # If player missed the ball, the game stops
                print("You missed the ball :-(")    
            if tilted_blow:
                # Alternative: angle of stroken ball gets titled according to position of racket
                dy_12 = bead_position["1"]["y"].latest_value - bead_position["2"]["y"].latest_value  # Distance between both beads along y
                alpha = np.arctan(dy_12/minimal_distance)  # Tilt angle in rad
                # Calculate tilted reflection angle in rad 
                beta = beta_next - 2*alpha
                # Make sure beta is not outside of allowed range
                beta = np.sign(beta) * min((np.abs(beta), deg_to_rad(beta_max)))       
finally:
    if save_video:
        bl.brightfield.stop_recording()       
    # Print score
    print("\n")
    print_in_box(f"Game over! You made {nn_blow} blow(s). Your score is {score}.")
    