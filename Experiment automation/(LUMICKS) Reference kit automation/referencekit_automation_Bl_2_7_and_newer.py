"""
BSD 2-Clause License

Copyright (c) 2025, LUMICKS B.V.
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
This script has been written for Bluelake >= 2.7

This script Iteratively catches beads, fishes for DNA and makes a kymograph at 10 pN.
The script is optimized for DNA tethers with a length of 37.8 kbp (12.85 micron).
Before running the script
    1) The traps have to be aligned horizontally.
    2) Calibrate the force
    3) Set the scan options to a continuous kymograph and configure the kymograph,
    ie make sure the correct location is being scanned when making a kymograph
    4) Check that the variables at the top and at the bottom of the script have the right values (eg, the channel names and the pathname where you want to store the data)

Workflow of the script
    1) Catch beads
    2) Catch a DNA tether
    3) Check if it is only a single tether, if not, try to break tethers until there is 1 left
    4) Move the microstage to 'J1' 
    5) Move the microstage to the protein channel 
    6) Set the force to 10 pN
    7) Start a kymograph (the variable for setting the duration of the kymograph is 'kymo_duration')
    8) Stop the kymograph and export the data 
    9) Restart step 1 to 8 n times, where n is the number of kymographs you want to make (the variable for setting the number of kymographs is 'max_kymos')
"""
import os
import time
from datetime import datetime
from pathlib import Path

import bluelake as bl
import numpy as np

name_bead_channel = "beads"  # This has to match the correct channel name in the UI
name_dna_channel = "DNA"  # This has to match the correct channel name in the UI
name_protein_channel = "Ch4"
name_junction = "J1"
name_buffer_channel = "buffer"
time_to_wait_for_flow = 5.0 # [s]
pressure_flow = 0.25  # Pressure when catching beads and fishing for DNA [bar]
dna_length_um = 0.34*37.8  # Length of DNA in micrometers [um]
kymo_duration = 300  # [s]
pathname = "D:/"  # Folder where the exported kymographs will be stored

distance = bl.timeline["Distance"]["Distance 1"]
force = bl.timeline["Force LF"]["Trap 2"]
match_score1 = bl.timeline["Tracking Match Score"]["Bead 1"]
match_score2 = bl.timeline["Tracking Match Score"]["Bead 2"]


def throw_if_beads_lost(match_threshold):
    """Raise an exception if we lose the beads."""
    if match_score1.latest_value < match_threshold or match_score2.latest_value < match_threshold:
        raise RuntimeError("Lost beads")

def start_flow(pressure):
    print("Starting Flow.")
    bl.fluidics.set_pressure(pressure)
    bl.fluidics.open(1, 2, 3, 6)


def close_valves():
    print("Stopping flow.")
    bl.fluidics.close(1, 2, 3, 6)


def catch_beads(match_threshold, pressure=pressure_flow):
    """Starts the flow and attempts to catch beads. Toggles the shutters when match score is too low."""
    start_flow(pressure)

    print(f"Waiting {time_to_wait_for_flow} sec for flow to begin.")
    bl.pause(time_to_wait_for_flow)

    print("Moving to bead channel.")
    bl.microstage.move_to(name_bead_channel)
    bl.pause(1.0)
    bl.shutters.clear(1,2)

    print("Trapping beads.")
    start_time = time.time()
    n = 0
    while match_score1.latest_value < match_threshold or match_score2.latest_value < match_threshold:
        """Drop beads that do not fulfill the template"""
        if 0 < match_score1.latest_value < match_threshold:
            bl.shutters.clear(1)
        if 0 < match_score2.latest_value < match_threshold:
            bl.shutters.clear(2)
        if n > 50:
            close_valves()
            raise RuntimeError("Needed more than 750 seconds to catch beads, stopping script")
        # If it's taking too long, maybe something is stuck in the trap. Clear both traps.

        if time.time() - start_time > 15:
            bl.shutters.clear(1,2)
            start_time = time.time()
            n += 1
        bl.pause(1.0)

    print("Got beads!")


def goto_distance(target, match_threshold, speed=1.0, tolerance=0.2):
    """Move trap 1 until it reaches the `target` distance from trap 2.

    Raises
    ------
    RuntimeError
        If the beads are lost (since we will not have a reliable distance then either)"""
    dx = target - distance.latest_value
    throw_if_beads_lost(match_threshold)

    while abs(dx) > tolerance:  # um
        bl.mirror1.move_by(dx=0.1 if dx > 0 else -0.1, speed=speed)
        dx = target - distance.latest_value
        throw_if_beads_lost(match_threshold)


def goto_force(target, match_threshold, speed=1, tolerance=1, tether_lost_threshold=5):
    """Move trap 1 until it reaches the `target` force on trap 2.

    Parameters
    ----------
    target : float
        The target force to reach [pN]
    match_threshold: float
        The threshold for the match score to check whether the beads are still caught [%]
    speed : float
        The speed at which to move the trap [um/s]
    tolerance : float
        The tolerance for the force measurement [pN]
    tether_lost_threshold : float
        The force threshold below which the tether is assumed lost [pN]

    Raises
    ------
    RuntimeError
        If the beads are lost (since we will not have a reliable distance then either)"""

    df = target - force.latest_value
    throw_if_beads_lost(match_threshold)

    while abs(df) > tolerance:  # pN
        bl.mirror1.move_by(dx=0.05 if df > 0 else -0.05, speed=speed)
        df = target - get_force(dt=0.2)
        throw_if_beads_lost(match_threshold)
    check_tether_breakage(tether_lost_threshold)


def get_force(dt=0.5):
    t0 = bl.timeline.now
    bl.pause(dt)
    t1 = bl.timeline.now
    return np.mean(force[t0:t1].data)


def catch_dna(min_distance, max_distance, match_threshold, force_threshold, fishing_speed, fishing_attempts, dt):
    """Moves to the DNA channel and starts oscillating the trap until a prescribed force threshold is reached.

    Parameters
    ---------
    min_distance : float
        The minimum distance to oscillate the trap to when fishing for DNA (in microns)
    max_distance : float
        The maximum distance to oscillate the trap to (in microns)
    match_threshold : float
        The threshold for the match score to check whether the beads are still caught [%]
    force_threshold : float
        The force threshold to check whether we have caught DNA [pN]
    fishing_speed : float
        The speed at which to fish for DNA [um/s]
    fishing_attempts : float
        The number of attempts to fish for DNA
    dt : float
        The time interval over which the force is averaged [s]
    """
    print("Moving to DNA channel")
    bl.microstage.move_to(name_dna_channel)
    bl.pause(2.0)
    print("reset force")
    bl.reset_force()
    bl.pause(1)
    f_current = get_force(dt=dt)
    attempts = 0
    while f_current < force_threshold:
        if attempts > fishing_attempts:
            raise RuntimeError(f"Max DNA fishing attempts {fishing_attempts} reached.")

        print(f"Fishing for DNA: attempt {attempts}/{fishing_attempts}")
        goto_distance(min_distance, match_threshold, speed=fishing_speed)
        bl.pause(0.5)
        goto_distance(max_distance, match_threshold, speed=fishing_speed)
        f_current = get_force(dt=dt)

        attempts += 1


def make_kymograph(force_kymo, match_threshold, name, path):
    try:
        goto_force(force_kymo, match_threshold)

        t0 = bl.timeline.now
        bl.confocal.start_scan()
        print("start kymograph")
        bl.pause(10.0)
        print(f"go to a force of {force_kymo} pN")
        goto_force(force_kymo, match_threshold, speed=1)
        bl.pause(kymo_duration)
    finally:
        print("Stop kymograph")
        bl.confocal.abort_scan()
        bl.pause(2)
        t1 = bl.timeline.now
        marker = bl.Marker(name, start = t0, stop = t1)
        marker.export(Path(path))
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


def check_multiple_tethers(max_distance, multiple_tether_force, force_threshold, match_threshold, max_retries, speed, dna_length):
    """Check whether there are multiple tethers by going to a large distance and checking the force.

    Parameters
    ---------
    max_distance : float 
        The maximum distance to go to check for multiple tethers [um]
    multiple_tether_force : float
        The force threshold above which we assume there are multiple tethers [pN]
    force_threshold : float 
        The force threshold above which we assume there is still a tether [pN]
    match_threshold : float 
        The threshold for the match score to check whether the beads are still caught [%]
    max_retries : integer
        The maximum number of attempts to break tethers
    speed : float 
        The speed at which to move the trap [um/s]
    """
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
            """Go to a large distance, trying to break tethers:"""
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
                force_threshold, tether_lost_threshold, max_kymos, force_kymo, dna_length):

    """
    Main workflow of the script. Catches beads, fishes for DNA and makes a kymograph at force_kymo pN.

    Parameters
    ----------
    experiment_name : string
        Name of the experiment (goes into the filename)
    path : string
        Where to store the data
    match_threshold : float
        The threshold for the match score to check whether the beads are still caught [%]
    dna_fishing_speed : float
        The speed at which the Trap should move when fishing for DNA [um/s]
    min_distance_fishing : float
        The minimum distance to oscillate the trap to when fishing for DNA (in microns)
    max_distance_fishing : float
        The maximum distance to oscillate the trap to (in microns)
    force_threshold : float
        The force threshold to check whether we have caught DNA [pN]
    tether_lost_threshold : float
        The force threshold to check whether the tether has been lost [pN]
    max_kymos : integer
        The maximum number to record
    force_kymo : float
        The force at which to record a kymograph [pN]
    dna_length : float
        The length of the DNA tether [micron]
    """
    path = validate_dir(path)

    kymo_count = -1

    while kymo_count < max_kymos:
        try:
            kymo_count += 1
            catch_beads(match_threshold=match_threshold)
            catch_dna(min_distance=min_distance_fishing, max_distance=max_distance_fishing,
                      fishing_speed=dna_fishing_speed, match_threshold=match_threshold, force_threshold=force_threshold,
                      fishing_attempts=10, dt=0.5)
            close_valves()
            bl.pause(3.0)
            bl.microstage.move_to(name_buffer_channel)
            check_multiple_tethers(max_distance=max_distance_fishing, multiple_tether_force=10, force_threshold=force_threshold,
                                   match_threshold=match_threshold, max_retries=7, speed=dna_fishing_speed, dna_length=dna_length)
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
            bl.microstage.move_to(name_junction, speed=1)  # speed in mm/s
            bl.pause(2.0)
            print("Go to protein channel")
            bl.microstage.move_to(name_protein_channel, speed=1)  # speed in mm/s
            bl.pause(2.0)
            print(f"Recording kymograph (kymo: {kymo_count}).")
            current_time = datetime.now().strftime("%H%M%S")
            name = f"{current_time}-{experiment_name}_tether={kymo_count}"
            check_tether_breakage(tether_lost_threshold)
            make_kymograph(force_kymo, match_threshold, name, path)

        except RuntimeError as e:
            if kymo_count == max_kymos:
                print(f"Reached max number of {max_kymos} attempts, stop script")
            else:
                print("Restarting protocol.")

        finally:
            print("Stop the flow and vent")
            bl.fluidics.stop_flow()

    print("Script finished!")


fd_workflow(experiment_name = "kymo",  # Name of the experiment (goes into the filename)
            path = pathname,  # Location to output the experiment.
            match_threshold = 70,  # Minimal template match threshold
            dna_fishing_speed = 4,  # DNA fishing speed [-]
            min_distance_fishing = 6,  # Minimal distance when fishing for DNA [um]
            max_distance_fishing = dna_length_um,  # Maximal distance when fishing for DNA [um]
            force_threshold = 20,  # Force threshold when fishing for DNA [pN]
            tether_lost_threshold = 5.0,  # When force falls below this threshold, tether is assumed lost [pN]
            max_kymos = 15,  # Maximal number of attempted kymographs [#]
            force_kymo = 10,  # Force at which the kymograph is recorded, which should be 10 pN to find the reference values provided by LUMICKS
            dna_length = dna_length_um)
