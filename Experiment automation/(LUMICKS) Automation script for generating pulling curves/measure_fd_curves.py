"""
Copyright 2020, LUMICKS B.V.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:
 
1. Redistributions of source code must retain the above copyright notice, this 
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, 
this list of conditions and the following disclaimer in the documentation 
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

"""Iteratively catch beads, fish for DNA and make F,d curves"""
from bluelake import trap1, trap2, stage, fluidics, pause, power, timeline,reset_force, confocal
import time
import os
import numpy as np

distance = timeline["Distance"]["Distance 1"]
force = timeline["Force LF"]["Trap 2"]
match_score1 = timeline["Tracking Match Score"]["Bead 1"]
match_score2 = timeline["Tracking Match Score"]["Bead 2"]
name_bead_channel = "beads"  # This has to match the correct channel name in the UI
name_dna_channel = "DNA"  # This has to match the correct channel name in the UI
time_to_wait_for_flow = 5.0


def throw_if_beads_lost(match_threshold):
    """Raise an exception if we lose the beads."""
    if match_score1.latest_value < match_threshold or match_score2.latest_value < match_threshold:
        raise RuntimeError("Lost beads")


def set_pressure(target):
    """Increase pressure until we are above a certain target level"""
    while fluidics.pressure < target:
        fluidics.increase_pressure()
        pause(1.0)


def start_flow(pressure):
    print("Starting Flow.")
    set_pressure(pressure)
    fluidics.open(1, 2, 3, 6)


def stop_flow():
    print("Stopping flow.")
    fluidics.close(1, 2, 3, 6)


def catch_beads(match_threshold, pressure=0.2):
    """Starts the flow and attempts to catch beads. Toggles the shutters when match score is too low."""
    start_flow(pressure)

    print(f"Waiting {time_to_wait_for_flow} sec for flow to begin.")
    pause(time_to_wait_for_flow)

    print("Moving to bead channel.")
    stage.move_to(name_bead_channel)

    pause(1.0)
    trap1.clear()
    trap2.clear()

    print("Trapping beads.")
    start_time = time.time()
    while match_score1.latest_value < match_threshold or match_score2.latest_value < match_threshold:
        """Drop beads that do not fulfill the template"""
        if 0 < match_score1.latest_value < match_threshold:
            trap1.clear()
        if 0 < match_score2.latest_value < match_threshold:
            trap2.clear()

        """If it's taking too long, maybe something is stuck in the trap. Clear both traps."""
        if time.time() - start_time > 15:
            trap1.clear()
            trap2.clear()
            start_time = time.time()

        pause(1.0)

    print("Got beads!")


def goto_distance(target, match_threshold, speed=1, tolerance=0.2):
    """Move trap 1 until it reaches the `target` distance from trap 2.

    Note: This throws an error if the beads are lost (since we will not have a reliable distance then either)"""
    dx = target - distance.latest_value
    throw_if_beads_lost(match_threshold)

    while abs(dx) > tolerance:  # um
        if dx > 0:
            trap1.move_by(dx=+0.1, speed=speed)
        else:
            trap1.move_by(dx=-0.1, speed=speed)

        dx = target - distance.latest_value
        throw_if_beads_lost(match_threshold)


def get_force():
    t0 = timeline.current_time
    pause(0.1)
    t1 = timeline.current_time
    f = np.mean(force[t0:t1].data)

    return f


def catch_dna(min_distance, max_distance, match_threshold, force_threshold, fishing_speed, max_retries=15):
    """Moves to the DNA channel and starts oscillating the trap until a prescribed force threshold is reached."""
    goto_distance(min_distance, match_threshold, speed=0)
    pause(1.0)

    print("Moving to DNA channel")
    stage.move_to(name_dna_channel)

    reset_force()
    pause(1.0)

    attempts = 0
    while get_force() < force_threshold:
        goto_distance(min_distance, match_threshold, speed=fishing_speed)
        reset_force()
        pause(1.0)
        goto_distance(max_distance, match_threshold, speed=fishing_speed)
        print(f"Fishing for DNA: attempt {attempts}/{max_retries}")

        if attempts > max_retries:
            raise RuntimeError(f"Max retries {max_retries} reached.")

        attempts += 1

    print("Moving to buffer.")
    stage.move_to("buffer")
    stop_flow()
    pause(5.0)


def make_fd_curve(min_distance, max_distance, force_threshold, match_threshold, fd_speed, name, path):
    """Measure an F, d curve."""
    goto_distance(min_distance, match_threshold, speed=5)
    reset_force()
    pause(1)

    try:
        timeline.mark_begin(name)
        trap1.move_by(dx=max_distance - min_distance, speed=fd_speed)
        pause(1.0)
    finally:
        timeline.mark_end(export=True if path else False, filepath=f"{path}/{name}.h5")

    if get_force() < force_threshold:
        raise RuntimeError(f"Lost tether.")

    pause(1)


def validate_dir(path):
    """Make sure the user entered a valid path"""
    path = path if path[-1] == '/' else path + '/'
    if os.path.isdir(path):
        return path
    else:
        raise RuntimeError("Invalid path specified. Did you pass a directory name?")


def fd_workflow(experiment_name, path, match_threshold, min_distance_fishing, max_distance_fishing, dna_fishing_speed,
                force_threshold, min_distance_fd, max_distance_fd, tether_lost_threshold, fd_speed, replicates,
                max_tethers):

    path = validate_dir(path)

    tether_count = -1
    while tether_count < max_tethers:
        try:
            tether_count += 1

            catch_beads(match_threshold=match_threshold)
            catch_dna(min_distance=min_distance_fishing, max_distance=max_distance_fishing,
                      fishing_speed=dna_fishing_speed, match_threshold=match_threshold, force_threshold=force_threshold)

            for replicate in np.arange(replicates):
                print(f"Recording F,d curve (tether: {tether_count}, replicate: {replicate}).")
                name = f"{experiment_name}_tether={tether_count:04}_replicate={replicate:04}"
                make_fd_curve(min_distance=min_distance_fd, max_distance=max_distance_fd,
                              force_threshold=tether_lost_threshold, match_threshold=match_threshold,
                              fd_speed=fd_speed, name=name, path=path)

        except RuntimeError as e:
            print(e)
            print("Restarting protocol.")


# This command runs the entire workflow.
# Note that the traps have to be aligned horizontally before starting the script.

fd_workflow(experiment_name="exp",              # Name of the experiment (goes into the filename)
            path="C:/Test",                     # Location to output the experiment. Note the slash
            match_threshold=80,                 # Minimal template match threshold
            dna_fishing_speed=10,               # DNA fishing speed [-]
            min_distance_fishing=10.0,          # Minimal distance when fishing for DNA [um]
            max_distance_fishing=16.0,          # Maximal distance when fishing for DNA [um]
            force_threshold=15.0,               # Force threshold when fishing for DNA [pN]
            min_distance_fd=12.5,               # Minimal distance of F,d curve [um]
            max_distance_fd=25.0,               # Approximate maximal distance of F,d curve [um]
            tether_lost_threshold=5.0,          # When force falls below this threshold, tether is assumed lost [pN]
            fd_speed=0.5,                       # Speed at which to pull [-]
            replicates=5,                       # Maximal number of replicates per tether (assuming it isn't lost) [#]
            max_tethers=1000)                   # Maximal number of attempted tethers [#]
