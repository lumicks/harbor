"""
BSD 2-Clause License

Copyright (c) 2021, LUMICKS B.V.
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
---------------------------
System requirements: Bluelake 1.7
This script creates 1 or more DNA knots. The traps start in one of the two following configurations:
            1  3   or 1  2
            4  2      4  3
When knotting, x and y positions change roughly as follows:
                     2
         1   3  -> 1   3 -> 1   3
         4   2     4        4   2

Before starting the script the following steps are required:

1) Make sure the initialization parameters match your experiment. This script assumes
    a) DNA length of ~11 um.
    b) Trap 1 has a different piezo stage. That is why speeds and distances are different between trap 1
    and the other traps.
2) Find the z-focus using one of the following 2 approaches and store it in the variables focus12 and focus34
    a) Use the auto fluorescence of the beads prior to the experiment. Step T1+2 and T3+4 until
    they show the sharpest bead image.
    b) Once T1+2 plane is found, you can also use the BF along with the template matching score and
    visual inspection to step T3+4 until the score is highest, i.e. until the beads
    in Trap 3 and Trap 4 resemble the beads in Trap 1 and Trap 2
3) Check the range of motion of the telescopes (z-motion) above and below the focus and whether the script stays withing
that range. You may have to adjust the knotting steps if that is not the case.
4) Capture two single tethers of DNA, depending on the force detection on the system, this will be on the trap pairs
1-3/4-2, or 1-2/4-3
5) Place beads roughly at the position from where you would like to start knotting
"""

from bluelake import mirror1, mirror2, mirror3, mirror4, stage, fluidics, pause, timeline, reset_force, telescope12, telescope34
import numpy as np
from itertools import combinations

number_of_knots = 1
bead_size = 4.42
dna_hold = 12.5 + bead_size   # Spacing between two beads holding DNA in um for starting configuration
trap_space = 9 + bead_size    # Vertical spacing in um for starting configuration

focus12 = 3.9749  # z focus position for Traps 1 and 2
focus34 = 4.0752   # z focus position for Traps 3 and 4

bead_config = "1-2/4-3"  # Choose the starting configuration for knotting, you can choose from "1-3/4-2" and "1-2/4-3"

force = timeline["Force LF"]["Trap 2"]
match_score1 = timeline["Tracking Match Score"]["Bead 1"]
match_score2 = timeline["Tracking Match Score"]["Bead 2"]
match_score3 = timeline["Tracking Match Score"]["Bead 3"]
match_score4 = timeline["Tracking Match Score"]["Bead 4"]


def throw_if_beads_lost(match_threshold):
    """Raise an exception if we lose the beads.

    Parameters
    ----------
    match_threshold : float
        raise exception if match with template is below match_threshold
    """
    if match_score1.latest_value < match_threshold or match_score2.latest_value < match_threshold\
            or match_score3.latest_value < match_threshold or match_score4.latest_value < match_threshold:
        raise RuntimeError("Lost beads")


def throw_if_templates_overlap(min_distance=bead_size):
    """
    Raise exception if two templates overlap, implying that one of the beads has disappeared or looks bad.

    Parameters
    ----------
    min_distance : float
        the minimal distance between templates if they are occupying different beads

    """
    x_positions = np.array([timeline["Bead position"]["Bead 1 X"].latest_value,
                                timeline["Bead position"]["Bead 2 X"].latest_value,
                                timeline["Bead position"]["Bead 3 X"].latest_value,
                                timeline["Bead position"]["Bead 4 X"].latest_value])

    y_positions = np.array([timeline["Bead position"]["Bead 1 Y"].latest_value,
                                timeline["Bead position"]["Bead 2 Y"].latest_value,
                                timeline["Bead position"]["Bead 3 Y"].latest_value,
                                timeline["Bead position"]["Bead 4 Y"].latest_value])

    comb = combinations([0, 1, 2, 3], 2)
    distances = np.zeros(6)

    for n, it in enumerate(list(comb)):
        template_distance = np.sqrt((x_positions[it[0]]-x_positions[it[1]])**2+(y_positions[it[0]]-y_positions[it[1]])**2)
        distances[n] = template_distance
    if any(distances < min_distance):
        raise RuntimeError("Two or more templates overlap. Make sure all beads look good, "
                           "eg by improving the focus or catching new beads ")


def find_bead(trap,  moving_speed, trap_num, moving_distance, match_threshold=80):
    """
    Find which trap corresponds to which bead
    """
    x_positions_pre = np.array([timeline["Bead position"]["Bead 1 X"].latest_value,
                                timeline["Bead position"]["Bead 2 X"].latest_value,
                                timeline["Bead position"]["Bead 3 X"].latest_value,
                                timeline["Bead position"]["Bead 4 X"].latest_value])
    trap.move_by(dx=moving_distance, speed=moving_speed)

    pause(2)
    throw_if_beads_lost(match_threshold)
    x_positions_post = np.array([timeline["Bead position"]["Bead 1 X"].latest_value,
                                 timeline["Bead position"]["Bead 2 X"].latest_value,
                                 timeline["Bead position"]["Bead 3 X"].latest_value,
                                 timeline["Bead position"]["Bead 4 X"].latest_value])
    diff_x_position = np.abs(x_positions_pre - x_positions_post)
    trap.move_by(dx=-moving_distance, speed=moving_speed)
    pause(2)
    idx = np.argmax(diff_x_position) + 1
    print(f"Trap {trap_num} has bead number {idx}")

    return idx


def position_trap(bead_number):
    """
    Parameters
    ----------
    bead_number : int
        number of the bead of which the position is computed

    Returns
    -------
    pos : tuple
        position of the trap, given the bead number, [x,y]
    """
    pos = [timeline["Bead position"][f"Bead {bead_number} X"].latest_value,
           timeline["Bead position"][f"Bead {bead_number} Y"].latest_value]
    return pos


def calibrate(trap, bead_number, dx=0.1, dy=0.1):
    """
    Calculates calibration factor for trap movement. Assumes that the number of the bead in the trap is known.
    If the trap is moved r amount, and the tracking system records r' movement through tracking the bead, then there is
    an C = r' / r factor between intended movement (r) and real movement (r'). To make the tracking system record r as
    real movement, we need to move the trap by r / C amount.
    """
    pos_init_x = timeline["Bead position"][f"Bead {bead_number} X"].latest_value
    trap.move_by(dx=dx, speed=1)
    pause(1)
    pos_new_x = timeline["Bead position"][f"Bead {bead_number} X"].latest_value

    pos_init_y = timeline["Bead position"][f"Bead {bead_number} Y"].latest_value
    trap.move_by(dy=dy, speed=1)
    pause(1)
    pos_new_y = timeline["Bead position"][f"Bead {bead_number} Y"].latest_value

    conv_factor_x = (pos_new_x - pos_init_x) / dx
    conv_factor_y = - (pos_new_y - pos_init_y) / dy
    trap.move_by(dx=-dx, dy=-dy, speed=1)
    print(f"Trap {bead_number} conv factors: {conv_factor_x,conv_factor_y}")

    return conv_factor_x, conv_factor_y


def move_to_start_position(trap, final_pos, start_pos, conversion_factor, moving_speed, bead_num,
                           match_threshold=80, max_attempts=2, max_distance=5, min_movement=0.2):
    """
    Move traps into position for knotting
    """
    attempts = 0
    distance = 10000
    current_pos = start_pos
    while distance > max_distance:
        if attempts > max_attempts:
            raise Exception(f'Required more than {max_attempts} attempts to reach start position. Try placing the '
                            f'traps closer to the starting position before starting the script.')
        throw_if_beads_lost(match_threshold)
        throw_if_templates_overlap()
        trap.move_by(dx=(final_pos[0] - current_pos[0]) / conversion_factor[0],
                     dy=-(final_pos[1] - current_pos[1]) / conversion_factor[1], speed=moving_speed)
        throw_if_beads_lost(match_threshold)
        current_pos = position_trap(bead_num)
        distance_old = distance
        distance = np.sqrt((final_pos[0] - current_pos[0])**2 + (final_pos[1] - current_pos[1])**2)
        if distance_old < distance:
            raise Exception(f'Bead moved away from destination, check if the setup moved, or if the template'
                            f'positioning changed.')
        if (np.abs(distance - distance_old) < min_movement) & (distance > max_distance):
            raise Exception(f'Bead did not move, check if the template positioning changed')
        attempts += 1


print("Setting focus")
telescope12.move_to(z=focus12, speed=0)  # this will move T1+2 to the corresponding best widefield height
telescope34.move_to(z=focus34, speed=0)  # this will move T3+4 to the corresponding best widefield height
pause(1)
print("Setting focus DONE")


print("Calibration")
bead_per_trap = [
    find_bead(mirror1, 15, 1, 2),
    find_bead(mirror2, 1, 2, 0.2),
    find_bead(mirror3, 1, 3, 0.2),
    find_bead(mirror4, 1, 4, 0.2)
]

# Use the current coordinates of Trap 1 as the starting position:
start_x_trap1 = position_trap(bead_per_trap[0])[0]
start_y_trap1 = position_trap(bead_per_trap[0])[1]

if bead_config == "1-3/4-2":
    trap_pos_refs = [
        # Choose configuration:
        # 1   3
        # 4   2
        [start_x_trap1, start_y_trap1],
        [start_x_trap1 + dna_hold, start_y_trap1 + trap_space],
        [start_x_trap1 + dna_hold, start_y_trap1],
        [start_x_trap1, start_y_trap1 + trap_space]]
    mirrorM = mirror2
    c_index = 1
elif bead_config == "1-2/4-3":
    trap_pos_refs = [
        # Choose configuration:
        # 1   2
        # 4   3
        [start_x_trap1, start_y_trap1],
        [start_x_trap1 + dna_hold, start_y_trap1],
        [start_x_trap1 + dna_hold, start_y_trap1 + trap_space],
        [start_x_trap1, start_y_trap1 + trap_space]]
    mirrorM = mirror3
    c_index = 2
else:
    raise Exception("Wrong bead config, please choose '1-3/4-2' or '1-2/4-3'.")

for i in np.array([1, 2, 3, 4]):
    print(f'trap {i} starting position is {trap_pos_refs[i-1]}')

conversion_factors = [
    (1, 1),
    calibrate(mirror2, bead_per_trap[1], dx=-1.5, dy=1.5),
    calibrate(mirror3, bead_per_trap[2], dx=1.5, dy=-1.5),
    calibrate(mirror4, bead_per_trap[3], dx=-1.5, dy=-1.5)
]

print("Calibration DONE")

# Find all trap positions:
trap_positions = [
    position_trap(bead_per_trap[0]),
    position_trap(bead_per_trap[1]),
    position_trap(bead_per_trap[2]),
    position_trap(bead_per_trap[3])
]

print('Move traps into starting position')

traps = [mirror1, mirror2, mirror3, mirror4]
speeds = [10, 1, 1, 1]
for i in np.arange(4):
    move_to_start_position(traps[i], trap_pos_refs[i], trap_positions[i], conversion_factors[i], speeds[i],
                           bead_per_trap[i])
    pause(1)

print('Traps are in starting position, start knotting')

mirror4.move_by(dx=0, dy=4 / conversion_factors[3][1], speed=2)    # Moves Trap 4 slightly up to allow Trap 2
# to move across 1-3 without applying too much force on the 4-2 tether

mirrorM.move_by(dx=-7 / conversion_factors[c_index][0], dy=0, speed=2)
mirrorM.move_by(dx=0, dy=6 / conversion_factors[c_index][1], speed=2)

for _ in np.arange(number_of_knots):
    telescope34.move_by(dz=-3, speed=0)                            # moves 3+4
    for _ in np.arange(8):  # Take large movements in small steps, such that Bluelake can be stopped during motion
        mirrorM.move_by(dx=0, dy=2 / conversion_factors[c_index][1], speed=2)   # Moves Trap M up
    telescope34.move_by(dz=6, speed=0)                            # moves up Trap 3 and 4
    for _ in np.arange(8):
        mirrorM.move_by(dx=0, dy=-2 / conversion_factors[c_index][1], speed=2)  # Moves Trap M down and below tether 1-3
    telescope34.move_to(z=focus34, speed=0)                              # moves Trap 3 and 4 back to the focus position
    telescope12.move_to(z=focus12, speed=0)                              # moves Trap 1 and 2 back to the focus position


mirror4.move_by(dx=0, dy=-1.5 / conversion_factors[3][1], speed=2)   # Moves Trap 4 slightly down to have more space
# between beads in Trap 1 and 2

pause(1)
print('Knotting script is done!')
