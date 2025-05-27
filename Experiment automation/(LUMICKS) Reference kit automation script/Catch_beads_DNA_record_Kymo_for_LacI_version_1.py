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
--------------------------------
This script Iteratively catches beads, fishes for DNA and makes a kymograph at 10 pN.
The script is optimized for DNA tethers with a length of 37.8 kbp (12.85 micron).
Before running the script
    1) The traps have to be aligned horizontally.
    2) Calibrate the force
    3) Set the scan options to a continuous kymograph and configure the kymograph,
    ie make sure the correct location is being scanned when making a kymograph
    4) Make sure the traps are not too close together for catching beads

Workflow of the script
    1) Catch beads
    2) Catch a DNA tether
    3) Check if it is only a single tether, if not, try to break tethers until there is 1 left.
    4) Move using the microstage to 'J1' at 60% speed
    5) Move using the microstage to protein channel at 60% speed
    6) Set the force to 10 pN
    7) Start a kymograph.
    8) Stop the kymograph.
    9) Restart step 1 to 8
"""
from datetime import datetime

import bluelake as bl
import time
import os
import numpy as np


name_bead_channel = "beads"  # This has to match the correct channel name in the UI
name_dna_channel = "DNA"  # This has to match the correct channel name in the UI
name_protein_channel = "Ch4"
name_junction = "J1"
name_buffer_channel = "buffer"
time_to_wait_for_flow = 5.0
dna_length = 0.34*37.8  # Length of DNA in micrometers

distance = bl.timeline["Distance"]["Distance 1"]
force = bl.timeline["Force LF"]["Trap 2"]
match_score1 = bl.timeline["Tracking Match Score"]["Bead 1"]
match_score2 = bl.timeline["Tracking Match Score"]["Bead 2"]


def throw_if_beads_lost(match_threshold):
    """Raise an exception if we lose the beads."""
    if match_score1.latest_value < match_threshold or match_score2.latest_value < match_threshold:
        raise RuntimeError("Lost beads")


def set_pressure(target):
    """Increase pressure until we are above a certain target level"""
    while bl.fluidics.pressure < target:
        bl.fluidics.increase_pressure()
        bl.pause(1.0)


def start_flow(pressure):
    print("Starting Flow.")
    set_pressure(pressure)
    bl.fluidics.open(1, 2, 3, 6)


def stop_flow():
    print("Stopping flow.")
    bl.fluidics.close(1, 2, 3, 6)


def catch_beads(match_threshold, pressure=0.25):
    """Starts the flow and attempts to catch beads. Toggles the shutters when match score is too low."""
    start_flow(pressure)

    print(f"Waiting {time_to_wait_for_flow} sec for flow to begin.")
    bl.pause(time_to_wait_for_flow)

    print("Moving to bead channel.")
    bl.microstage.move_to(name_bead_channel)
    bl.pause(1.0)
    bl.mirror1.clear()
    bl.mirror2.clear()

    print("Trapping beads.")
    start_time = bl.time.time()
    n = 0
    while match_score1.latest_value < match_threshold or match_score2.latest_value < match_threshold:
        """Drop beads that do not fulfill the template"""
        if 0 < match_score1.latest_value < match_threshold:
            bl.mirror1.clear()
        if 0 < match_score2.latest_value < match_threshold:
            bl.mirror2.clear()
        if n > 50:
            stop_flow()
            raise Exception("Need more than 750 seconds to catch beads, stop script")
        """If it's taking too long, maybe something is stuck in the trap. Clear both traps."""
        if bl.time.time() - start_time > 15:
            bl.mirror1.clear()
            bl.mirror2.clear()
            start_time = bl.time.time()
            n += 1
        bl.pause(1.0)

    print("Got beads!")


def goto_distance(target, match_threshold, speed=1, tolerance=0.2):
    """Move trap 1 until it reaches the `target` distance from trap 2.

    Note: This throws an error if the beads are lost (since we will not have a reliable distance then either)"""
    dx = target - distance.latest_value
    throw_if_beads_lost(match_threshold)

    while abs(dx) > tolerance:  # um
        bl.mirror1.move_by(dx=0.1 if dx > 0 else -0.1, speed=speed)
        dx = target - distance.latest_value
        throw_if_beads_lost(match_threshold)


def goto_force(target, match_threshold, speed=1, tolerance=1, tether_lost_threshold=5):
    """Move trap 1 until it reaches the `target` force on trap 2.

    Note: This throws an error if the beads are lost"""
    df = target - force.latest_value
    throw_if_beads_lost(match_threshold)

    while abs(df) > tolerance:  # um
        bl.mirror1.move_by(dx=0.05 if df > 0 else -0.05, speed=speed)
        df = target - get_force(dt=0.2)
        throw_if_beads_lost(match_threshold)
    check_tether_breakage(tether_lost_threshold)


def get_force(dt=0.5):
    t0 = bl.timeline.current_time
    bl.pause(dt)
    t1 = bl.timeline.current_time
    f = np.mean(force[t0:t1].data)

    return f


def catch_dna(min_distance, max_distance, match_threshold, force_threshold, fishing_speed, fishing_attempts, dt):
    """Moves to the DNA channel and starts oscillating the trap until a prescribed force threshold is reached."""
    print("Moving to DNA channel")
    bl.microstage.move_to(name_dna_channel)
    bl.pause(2.0)
    print("reset force")
    bl.reset_force()
    f_current = get_force(dt=dt)

    attempts = 0
    max_retries = fishing_attempts*fishing_attempts
    while f_current < force_threshold:
        if attempts > max_retries:
            raise RuntimeError(f"Max DNA fishing attempts {max_retries} reached.")

        if attempts % fishing_attempts == 0:  # If the number of attempts is an integer times max_retries, do this:
            print("Moving to DNA channel")
            bl.microstage.move_to(name_dna_channel)
            print(f"Moving beads a distance {dna_length} um apart")
            goto_distance(dna_length, match_threshold, speed=fishing_speed)
            bl.pause(3.0)
            print("Moving to dna channel")
            bl.microstage.move_to(name_dna_channel)
            bl.pause(1.0)
        print(f"Fishing for DNA: attempt {attempts}/{max_retries}")
        goto_distance(min_distance, match_threshold, speed=fishing_speed)
        bl.pause(0.5)
        goto_distance(max_distance, match_threshold, speed=fishing_speed)
        f_current = get_force(dt=dt)

        attempts += 1


def make_kymograph(force_kymo, match_threshold, name, path):
    goto_force(force_kymo, match_threshold)
    try:
        goto_force(force_kymo, match_threshold)
        bl.timeline.mark_begin(name)
        bl.confocal.start_scan()
        print("start kymograph")
        bl.pause(10.0)
        # for df in [10]:
        print(f"go to a force of {force_kymo} pN")
        goto_force(force_kymo, match_threshold, speed=1)
        bl.pause(30)       ## needs to be 300 s in final version
    finally:
        print("Stop kymograph")
        bl.confocal.abort_scan()
        bl.timeline.mark_end(export=True if path else False, filepath=f"{path}/{name}.h5")
        # exit()
    bl.pause(1)


def validate_dir(path):
    """Make sure the user entered a valid path"""
    path = path if path[-1] == '/' else path + '/'
    if os.path.isdir(path):
        return path
    else:
        raise RuntimeError("Invalid path specified. Did you pass a directory name?")


def check_tether_breakage(test_force):
    if get_force() < test_force:
        raise RuntimeError("Lost tether")


def check_multiple_tethers(max_distance, multiple_tether_force, force_threshold, match_threshold, max_retries,speed):
    goto_distance(0.5 * dna_length, match_threshold, speed=speed)
    bl.reset_force()
    goto_distance(0.95 * dna_length, match_threshold, speed=speed)
    if get_force() > multiple_tether_force:
        print("Caught multiple DNA tethers, let's try to break tethers until there is only one left")
        step_size = 0.1 * dna_length
        one_tether_remains = False
        n = 0
        while not one_tether_remains:
            if n == max_retries:
                raise RuntimeError(f"Required more than {max_retries} attempts to break tethers")
            """Go to a large distance, trying to brake tethers:"""
            goto_distance(dna_length + n * step_size, match_threshold, speed=speed)
            bl.pause(0.5)
            """Go to smaller distance, to check if there are still multiple tethers"""
            goto_distance(0.95 * dna_length, match_threshold, speed=speed)
            bl.pause(0.5)
            if get_force() < multiple_tether_force:
                """Check whether there is still a tether left, as it can happen that multiple tethers break at once"""
                goto_distance(max_distance, match_threshold, speed=speed)
                if get_force() > force_threshold:
                    print("Looks like only one tether is left, continue with the rest of the protocol")
                    one_tether_remains = True
                else:
                    raise RuntimeError("Lost all tethers")
            n += 1


def fd_workflow(experiment_name, path, match_threshold, dna_fishing_speed, min_distance_fishing, max_distance_fishing,
                force_threshold, tether_lost_threshold, max_kymos, force_kymo):
    path = validate_dir(path)

    kymo_count = -1
    while kymo_count < max_kymos:
        try:
            kymo_count += 1
            catch_beads(match_threshold=match_threshold)
            catch_dna(min_distance=min_distance_fishing, max_distance=max_distance_fishing,
                      fishing_speed=dna_fishing_speed, match_threshold=match_threshold, force_threshold=force_threshold,
                      fishing_attempts=4, dt=0.5)
            print("Stop flow")
            stop_flow()
            bl.pause(3.0)
            check_multiple_tethers(max_distance=max_distance_fishing, multiple_tether_force=10, force_threshold=25,
                                   match_threshold=60, max_retries=7, speed=dna_fishing_speed)
            print("Go to distance of 0.75 x tether length")
            goto_distance(0.75 * dna_length, match_threshold, speed=dna_fishing_speed)
            bl.pause(1.0)
            print("reset force")
            bl.reset_force()
            print("Go to distance of 1 x tether length")
            goto_distance(dna_length, match_threshold, speed=dna_fishing_speed)
            check_tether_breakage(tether_lost_threshold)
            print(f"go to a force of {force_kymo} pN")
            goto_force(force_kymo, match_threshold)
            print("Go to channel junction")
            bl.microstage.move_to(name_junction, speed=60)
            bl.pause(2.0)
            print("Go to protein channel")
            bl.microstage.move_to(name_protein_channel, speed=60)
            bl.pause(2.0)
            print(f"Recording kymograph (kymo: {kymo_count}).")
            current_time = datetime.now().strftime("%H%M%S")
            name = f"{current_time}-{experiment_name}_tether={kymo_count}"
            check_tether_breakage(tether_lost_threshold)
            make_kymograph(force_kymo, match_threshold, name, path)
            print("Script finished!")

        except RuntimeError as e:
            print(e)
            if kymo_count == max_kymos:
                print(f"Reached max number of {max_kymos} attempts, stop script")
            else:
                print("Restarting protocol.")


fd_workflow(experiment_name="kymo",  # Name of the experiment (goes into the filename)
            path="D:/",  # Location to output the experiment.
            match_threshold=60,  # Minimal template match threshold
            dna_fishing_speed=4,  # DNA fishing speed [-]
            min_distance_fishing=7.5,  # Minimal distance when fishing for DNA [um]
            max_distance_fishing=dna_length,  # Maximal distance when fishing for DNA [um]
            force_threshold=20,  # Force threshold when fishing for DNA [pN]
            tether_lost_threshold=5.0,  # When force falls below this threshold, tether is assumed lost [pN]
            max_kymos=15,  # Maximal number of attempted kymographs [#]
            force_kymo=10)  # Force at which the kymograph is recorded, which should be 10 pN to find the reference values provided by LUMICKS
