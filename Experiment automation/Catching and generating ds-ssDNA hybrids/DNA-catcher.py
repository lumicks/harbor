"""

Copyright 2022, LUMICKS B.V.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that 
the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the 
following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the 
following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, 
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""


import bluelake as bl
import lumicks.pylake as lk
import time
import numpy as np
from dataclasses import dataclass


"""
To catch beads and DNA, then overstretch to generate ssDNA
ends in buffer channel, stretched to flow_stretch * Lc and flow off
"""

# ===== user variables
total_nt = 17853
ssDNA_nt = 5005  # to be unwound

# ===== maybe variables
## tether catching x * Lc
flow_stretch = 0.4

## waypoint names
main_channel = ["beads", "DNA", "buffer"]

## flow rates for each step
bead_catching_pressure = 0.3
dna_catching_pressure = 0.2


# ===== functions and setup for functions
@dataclass
class DNA_model:
    """This dataclass contains a set of parameters used for DNA model calculations

    Returns
    -------
    dataclass : DNA_model
        Contains a set of parameters used for DNA model calculations
    """

    def __init__(self, ds_DNA, ss_DNA):
        if ds_DNA < 0:
            raise ValueError("final ssDNA will not exceed total nucleotides")

        self.ds_DNA = ds_DNA
        self.ss_DNA = ss_DNA

    @property
    def nt_DNA(self) -> int:
        return self.ds_DNA + self.ss_DNA

    model_ewlc = lk.odijk("ds")  # eWLC model for dsDNA
    model_fjc = lk.freely_jointed_chain("ss")  # FJC model for ssDNA
    model_hybrid = model_ewlc + model_fjc

    # FIXED - don't use defaults
    Lp_ds: float = 42.0  # persistent length (nm)
    Lp_ss: float = 0.84  # persistence length ssDNA (nm)
    bp_ss: float = 0.56  # basepair length ssDNA (nm)
    bp_ds: float = 0.34  # basepair length dsDNA (nm)

    # constant - use defaults
    kT: float = model_hybrid.defaults["kT"].value  # Boltzmann*temperature, default is 4.11 pN nm
    S: float = model_hybrid.defaults["ds/St"].value  # stretching modulus, default is 1500

    @property
    def Lc_ds(self) -> float:
        return self.ds_DNA * self.bp_ds  # contour length double stranded (nm)

    @property
    def Lc_ss(self) -> float:
        return self.ss_DNA * self.bp_ss  # contour length single stranded (nm)

    @property
    def Lc_DNA(self) -> float:
        return (self.Lc_ds + self.Lc_ss) * 1e-3  # contour length of total DNA (um)!!!!

    def DNA_length_force(self, force_DNA) -> float:
        dna_length_force = (
            self.model_hybrid(
                force_DNA,
                {
                    "ds/Lp": self.Lp_ds,
                    "ds/Lc": self.Lc_ds,
                    "ss/Lp": self.Lp_ss,
                    "ss/Lc": self.Lc_ss,
                    "kT": self.kT,
                    "ss/St": self.S,
                    "ds/St": self.S,
                },
            )
        ) * 1e-3
        return dna_length_force  # contour length of DNA at a certain force (um)

    def DNA_force_length(self, length_DNA) -> float:  # length (um?) in force out
        dna_force_length = self.model_hybrid.invert()(
            [length_DNA * 1e3],
            {
                "ds/Lp": self.Lp_ds,
                "ds/Lc": self.Lc_ds,
                "ss/Lp": self.Lp_ss,
                "ss/Lc": self.Lc_ss,
                "kT": self.kT,
                "ss/St": self.S,
                "ds/St": self.S,
            },
        )[0]
        return dna_force_length


match_score1: float = bl.timeline["Tracking Match Score"]["Bead 1"]
match_score2: float = bl.timeline["Tracking Match Score"]["Bead 2"]
distance: float = bl.timeline["Distance"]["Distance 1"]
force: float = bl.timeline["Force LF"]["Trap 2"]
dna_fishing_speed: float = 10  # DNA fishing speed [μm/s]
dna_stretching_speed: float = 500  # DNA stretching speed [nm/s]
match_threshold = 60
default_model = DNA_model(total_nt, 0)
if not ssDNA_nt == 0:  # ssDNA model parameters from workflow parameters
    unwound_model = DNA_model(total_nt - ssDNA_nt, ssDNA_nt)
global dna_contour_length
dna_contour_length = default_model.Lc_ds * 1e-3


def set_pressure(target):
    if target <= 0:
        bl.fluidics.start_venting()
    else:
        while bl.fluidics.pressure <= target:
            old_press = bl.fluidics.pressure
            bl.fluidics.increase_pressure()
            new_press = bl.fluidics.pressure
            if abs(new_press - old_press) < 0.01:
                raise RuntimeError("!! Target pressure cannot be reached!")
        while bl.fluidics.pressure > target:
            bl.fluidics.decrease_pressure()


def set_flow(pressure, channels):
    print(f"Changing pressure to {pressure:.2f}")
    set_pressure(target=pressure)
    open_channels = [i + 1 for i, x in enumerate(channels) if x]
    if open_channels:
        print("Opening channels " + str(open_channels))
        bl.fluidics.open(*open_channels)
    closed_channels = [i + 1 for i, x in enumerate(channels) if not x]
    if closed_channels:
        print("Closing channels " + str(closed_channels))
        bl.fluidics.close(*closed_channels)


def catch_beads():
    countbeaddrop = 0
    bead_drop_time = 30
    bead_pressure_cycle_time = 15

    print("Trapping beads.")
    set_flow(pressure=bead_catching_pressure, channels=[True, True, True, False, False, True])
    next_bead_drop = time.time() + bead_drop_time
    next_pressure_cycle = time.time() + bead_pressure_cycle_time
    while (
        match_score1.latest_value < match_threshold or match_score2.latest_value < match_threshold
    ):
        """Drop beads that do not fulfill the template"""
        if 0 < match_score1.latest_value < match_threshold:
            bl.shutters.clear(1, delay_ms=100)
        if 0 < match_score2.latest_value < match_threshold:
            bl.shutters.clear(2, delay_ms=100)

        set_pressure(bead_catching_pressure)

        if time.time() > next_bead_drop:
            countbeaddrop += 1
            bl.shutters.clear(1, 2, delay_ms=100)
            print(f"{countbeaddrop*bead_drop_time}s since last beads")
            next_bead_drop += bead_drop_time

        if time.time() > next_pressure_cycle:
            print(f"pressure cycling to flush more beads")
            set_pressure(0.7)
            set_pressure(0.0)
            bl.pause(3)
            next_pressure_cycle += bead_pressure_cycle_time

        if countbeaddrop > 20:
            print("It's been too long with no beads; stopping protocol.")
            exit()

        bl.pause(1.0)

    print("Got beads!")
    return 1


def get_distance():
    t0 = bl.timeline.current_time
    bl.pause(0.1)
    t1 = bl.timeline.current_time
    d = np.mean(distance[t0:t1].data)
    return d


def throw_if_beads_lost():
    if match_score1.latest_value < match_threshold or match_score2.latest_value < match_threshold:
        raise RuntimeError("Lost beads.")


def throw_if_tether_lost():
    tether_lost_threshold = 10
    if get_force() < tether_lost_threshold and get_distance() > 1.1 * dna_contour_length:
        raise RuntimeError("Lost tether.")


def goto_distance(target_distance):
    dx = target_distance - distance.latest_value
    throw_if_beads_lost()

    while abs(dx) > 0.2:
        if dx > 0:
            bl.mirror1.move_by(dx=+0.1, speed=dna_fishing_speed)
        else:
            bl.mirror1.move_by(dx=-0.1, speed=dna_fishing_speed)

        dx = target_distance - distance.latest_value
        throw_if_beads_lost()


def get_force():
    t0 = bl.timeline.current_time
    bl.pause(0.1)
    t1 = bl.timeline.current_time
    f = np.mean(force[t0:t1].data)

    return f


def catch_dna():
    max_retries = 15
    force_threshold = 15

    set_flow(pressure=dna_catching_pressure, channels=[True, True, True, False, False, True])
    goto_distance(target_distance=dna_contour_length * flow_stretch)
    bl.pause(1.0)
    bl.reset_force()
    bl.pause(0.5)

    attempts = 0
    while get_force() < force_threshold:
        goto_distance(target_distance=dna_contour_length * flow_stretch)
        bl.pause(1.0)
        bl.reset_force()
        bl.pause(0.5)
        goto_distance(target_distance=dna_contour_length)
        print(f"Fishing for DNA: attempt {attempts}/{max_retries}")

        if attempts > max_retries:
            raise RuntimeError(f"Max retries {max_retries} reached.")

        attempts += 1


def set_force_feedback(target_force, speed):
    # config force feedback
    bl.force_feedback.enabled = False
    bl.force_feedback.set_device("1")
    bl.force_feedback.set_detector("Trap 2")
    bl.force_feedback.set_frequency(30)
    # check these?
    bl.force_feedback.set_angle(0)
    bl.force_feedback.set_lock_motion_angle()
    bl.force_feedback.set_pid_settings(
        kp=15,
        ki=100,
        kd=0,
        max_step=int(speed / 30),
        reversed=False,
    )
    bl.force_feedback.set_target(target_force)


def goto_force(target_force):
    speed = dna_stretching_speed
    set_force_feedback(target_force, speed)
    bl.force_feedback.enabled = True
    print("moving to force:" + str(target_force) + " pN with speed of:" + str(speed) + " nm/s.")
    while abs(get_force() - target_force) > 0.5:
        throw_if_tether_lost()
        throw_if_beads_lost()
    bl.force_feedback.enabled = False


def check_dna():
    force_check = 15
    distance_threshold = 2.5  # difference extension measured - model

    # relax tether to zero force
    goto_distance(target_distance=dna_contour_length * 0.4)
    bl.pause(1.0)
    bl.reset_force()
    bl.pause(1.0)

    # check for single tether
    goto_distance(target_distance=dna_contour_length * 1)
    if get_force() > 80:  # Check for multiple tethers; TODO is 300 ok? more/less?
        print("multiple tethers found, trying to break them.")
        set_force_feedback(300, dna_stretching_speed)
        bl.force_feedback.enabled = True
        print(f"moving to force: {300} pN with speed of: {dna_stretching_speed:.0f} nm/s.")
        while abs(get_force() - 300) > 0.5 and get_force() > 75:
            throw_if_tether_lost()
            throw_if_beads_lost()
        bl.force_feedback.enabled = False
        print("Possible single tether, relaxing.")

    # check if extension is between default and unwound models (if exists) at target force
    goto_force(force_check)
    default_length = default_model.DNA_length_force(force_check)
    if not ssDNA_nt == 0:
        unwound_length = unwound_model.DNA_length_force(force_check)
    else:
        unwound_length = default_length

    current_dist = get_distance()
    if (
        current_dist >= default_length - distance_threshold
        and current_dist <= unwound_length + distance_threshold
    ):
        print("tether passed check")
    else:
        dx = default_length - current_dist
        dx2 = unwound_length - current_dist
        raise RuntimeError(
            f"Measured DNA tether length difference with model is too large (>{distance_threshold}um)"
            + f"{dx:.2f} from default, {dx2:.2f} from unwound"
        )

    # go back to relaxed distance
    goto_distance(target_distance=dna_contour_length * 0.4)


def overstretching():
    tether_lost_threshold = 10

    strained = unwound_model.Lc_DNA * 1.25
    relax = default_model.Lc_DNA * 0.4
    probe = default_model.Lc_DNA * 1.1

    # at strain distance its completely unwound
    strain_force = unwound_model.DNA_force_length(strained)
    # after unwind, this is the target force
    probe_force = unwound_model.DNA_force_length(probe)

    print(f".expected probe force = {round(probe_force,1)}")
    print(f".expected strain force = {round(strain_force,1)}")

    # internal variables and things actually start from here
    acceptable_probe_force_diff = 15
    max_attempts = 5  # can prob go more than 5, but how patient are you
    overstretchtime = 2
    attempts = 0
    success = 0
    while success == 0:
        # 0 means keep going (still trying because return (to relaxed) was still dsDNA-like)
        # 1 means too many attempts, 2 means broken?
        if attempts > max_attempts:
            # should stop the whole overstretching attempting
            success = 4
            print(f"Tried enough times ({max_attempts})")
        attempts += 1

        goto_distance(target_distance=probe)
        unwprobeforce = get_force()
        print(f"..unwprobeforce = {round(unwprobeforce,1)}")

        goto_distance(target_distance=strained)
        bl.pause(overstretchtime)
        measstrainedforce = get_force()
        print(f"..measstrainedforce = {round(measstrainedforce,1)}")

        goto_distance(target_distance=probe)
        rewprobeforce = get_force()
        print(f"..rewprobeforce = {round(rewprobeforce,1)}")

        goto_distance(target_distance=relax)

        # these are the checks
        if measstrainedforce < tether_lost_threshold:
            # tether broke
            success = 2

        else:
            # tether intact
            # check if the force at 'probe' distance is at the right value
            if abs(rewprobeforce - probe_force) < acceptable_probe_force_diff:
                # seems right
                success = 1
            else:
                # check if hysterisis?
                if measstrainedforce > strain_force:
                    print("! measured strain exceeds predicted, maybe multiple tether")
                    print(f"  {round(measstrainedforce/strain_force,1)}x more than expected")
                elif abs(rewprobeforce - unwprobeforce) < acceptable_probe_force_diff:
                    print("looks about the same, try again")
                else:
                    print("something else went wrong; stopping for this tether")
                    success = 3

    return success, unwound_model.Lc_DNA


def catch_beads_n_tether(dna_contour_length):  # the main script
    tetherpass = False
    while not tetherpass:
        try:
            # Step 0
            print(">>> Resetting Flow.")
            bl.fluidics.stop_flow()

            # Step 1 start in bead channel
            print(">> Moving to bead channel.")
            bl.microstage.move_to(waypoint=main_channel[0], speed=1)

            # Step 2 turn on flow
            print(">> Starting Flow.")
            set_flow(
                pressure=bead_catching_pressure, channels=[True, True, True, False, False, True]
            )
            bl.pause(1)

            # Step 3 discard old things (might be junk, or previously used)
            print(">> Discarding trapped things.")
            bl.shutters.clear(1, 2, delay_ms=100)

            # Step 4 catch new beads
            print(">> Catching beads.")
            catch_beads()

            # Step 5 move to dna after catching
            print(">> Moving to DNA channel.")
            bl.microstage.move_to(waypoint=main_channel[1], speed=1)

            # Step 6 check for tether
            print(">> Checking for DNA tether.")
            catch_dna()
            goto_distance(target_distance=dna_contour_length * flow_stretch)

            # Step 7 move out after tethered
            print(">> Moving to Buffer channel.")
            bl.microstage.move_to(waypoint=main_channel[2], speed=1)

            # Step 8 stop flow for calibration
            print(">> Stopping Flow for calibration")
            bl.fluidics.stop_flow()
            bl.pause(3)

            # Step 10 zero force
            bl.reset_force()
            bl.pause(1)

            # Step 11 - check dna single tether
            check_dna()

            # Step 12 - overstretching
            print(">> making ssDNA")
            if not ssDNA_nt == 0:
                ssDNA, dna_contour_length = overstretching()
                if not ssDNA == 1:
                    raise RuntimeError("ssDNA generation failed")
            tetherpass = True

        except RuntimeError as e:
            bl.force_feedback.enabled = False
            print("Runtime error: " + e.args[0])
            print(">>> Restarting protocol.")
    print("ssDNA generation succesful, please continue your experiment!")


# Code execution
try:
    catch_beads_n_tether(dna_contour_length)
finally:
    bl.fluidics.stop_flow()
