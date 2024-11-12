# coding=utf-8
"""
BSD 2-Clause License

Copyright (c) 2021, SB FMH-lab (CNB-CSIC)
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

--------------------------------------------------------------------------------

PROGRAMMABLE FORCE SEQUENCE
Automation script for Bluelake 1.7.1 (C-Trap, Lumicks)
Sara de Bragança, sbraganca@cnb.csic.es
2021-09-21

This automation script allows the user to program a sequence of forces to be applied to a tethered molecule.
Using this script to apply forces increases the reproduciblity of the experiments. Firstly, because the user,
instead of defining a particular distance between the traps, defines a target force. Then, the machine finds
the experimental distance that correlates to the target force, reducing the error related to the variablity
between different molecules. In addition, the user can choose to repeat the exact same force sequence several
times to the same molecule or even for two different experimental conditions (by choosing to repeat in two
different flow channels).

Workflow:
1 - Saves the initial position of Trap1.
2 - Validates the parameters given by the user.
3 - Finds the absolute positions for Trap1 corresponding to each target force given by the user.
4 - Applies the sequence of movements and pauses as predefined by the user.
5 - Repeats the sequence, if it applies.
6 - If it applies: Changes to the second channel.
                   Pauses for incubation.
                   Proceeds to apply the force sequence and repetitions.
"""

from bluelake import mirror1, mirror2, pause, timeline, microstage
import numpy as np


def apply_steps_workflow(name_id, forces, times, speeds, repetitions, repeat_in_another_channel, channel_info):
    """
    Run the force sequence designed by the user.

    :param name_id: Identification label for the experiment. This will be the exported markers filename.
    :param forces: Target forces values for each step. [pN]
    :param times: Time to pause at each step. [s]
    :param speeds: Velocity of the movement from one target force value to the next one in the sequence. [um/s]
    :param repetitions: Number of times the sequence should be repeated.
    :param repeat_in_another_channel: Boolean. Repeat execution in another channel or not.
    :param channel_info: Contains channel_name and incubation_time.
                         The channel_name must be a waypoint in Bluelake.
                         incubation_time in seconds.
    """
    print('Starting the execution.')
    safe_position = mirror1.position.x
    print('Safe position saved.')
    inputs = rearrange_inputs(name_id, forces, times, speeds, repetitions)
    validate_inputs(inputs)
    inputs = find_positions(inputs, safe_position)
    apply_steps(inputs, f'{name_id}')
    if repeat_in_another_channel:
        move_to_channel(channel_info)
        apply_steps(inputs, f'{name_id}_{channel_info[0]}')
    print('Finished!')


def rearrange_inputs(name_id, forces, times, speeds, repetitions):
    """
    Rearrange inputs given by the user to facilitate code readability.
    :return: dictionary with the inputs from the user
    """
    inputs = {'name_id': name_id, 'reps': repetitions}
    quick_speed = 10  # add speed for the first movement, always quick 10 um/s
    speeds = [quick_speed] + speeds
    seq = np.array([forces, times, speeds])
    seq = seq.transpose()
    inputs['seq'] = seq
    return inputs


def validate_inputs(inputs):
    """
    Find if any values are missing or are out of the expected ranges.
    :param inputs: dictionary with the inputs from the user
    """
    print('Validating inputs.')
    inputs = inputs['seq']
    for step in inputs:
        force, time, speed = step
        validate_force(force)
        validate_time(time)
        validate_speed(speed)


def validate_force(value):
    assert 0 <= value < 300, 'Execution aborted. The accepted force range goes from 0 pN to 300 pN.'


def validate_time(value):
    assert 0 <= value < 60, 'Execution aborted. The accepted time range goes from 0 s to 60 s.'


def validate_speed(value):
    assert 0 <= value < 20, 'Execution aborted. The accepted speed range goes from 0 um/s to 20 um/s.'


def find_positions(inputs, safe_position):
    """
    Find the absolute positions of the Trap1 for a sequence of forces measured in the Trap2.
    :param inputs: dictionary with the inputs from the user
    :param safe_position: position to return to if anything goes wrong [um]
    :return: dictionary with the inputs from the user and the positions for the target forces
    """
    print('Finding the absolute position of Trap1 for each target force.')
    print('Press Stop button if the molecule breaks.')
    seq = []
    steps = inputs['seq']
    for step in steps:
        target, time, speed = step
        position, force = find_target_position(target, safe_position)
        seq.append([target, time, speed, position.x])
    inputs['seq'] = seq
    return inputs


def find_target_position(target, safe_position, increments=0.1, acceptable_error=0.5, max_movement=20):
    """
    Move Trap1 by small increments while evaluating the force measured in the Trap2, until a target force is reached.
    Get the absolute position of Trap1.
    :param target: force value in pN
    :param safe_position: position to return to if anything goes wrong
    :param increments: x-axis position increments [um]
    :param acceptable_error: accepted margin of error for the force value
    :param max_movement: maximum relative displacement allowed
    :return: position [um] and real force [pN]
    """
    force = get_force()
    slope = 1 if force < target else -1
    dx = slope * increments  # [um]
    df = slope * round(target - force, 1)  # [pN]
    verify_beads()
    n = 0
    n_max = int(max_movement / increments)
    while df >= acceptable_error:
        mirror1.move_by(dx=dx)
        verify_beads()
        force = get_force()
        df = slope * round(target - force, 1)
        verify_relative_displacement(n, n_max, safe_position)
        n += 1
    position = mirror1.position
    print(f'F = {round(force, 2)} pN -> x = {round(position.x, 2)} um')
    return position, force


def get_force():
    """
    Get the last force value measured in Trap2. Low frequency value and in the x-direction.
    :return: force value [pN]
    """
    force = timeline['Force LF']['Force 2x'].latest_value
    return force


def verify_beads(acceptable_score=70):
    """
    Verify that there are beads in the traps.
    :param acceptable_score: threshold used to evaluate if there is a good bead in the trap [%]
    """
    score_bead1 = timeline['Tracking Match Score']['Bead 1'].latest_value
    score_bead2 = timeline['Tracking Match Score']['Bead 2'].latest_value
    assert score_bead1 > acceptable_score and score_bead2 > acceptable_score, 'Execution aborted. It seams that one of the beads is gone.'


def verify_relative_displacement(n, n_max, safe_position):
    """
    Verify how much has Trap1 moved from previous iterations.
    :param n: number of iterations
    :param n_max: maximum number of iterations allowed
    :param safe_position: position to return to if anything goes wrong
    """
    if n > n_max:
        mirror1.move_to(x=safe_position, speed=5)
        assert False, 'Something went wrong. Maybe the molecule is broken.'


def apply_steps(inputs, label):
    """
    Apply the sequence of movements in Trap1.
    :param inputs: dictionary with the input from the user and the positions for the target forces
    :param label: string with label to include in the markers name
    """
    print('Starting the sequence.')
    reps = inputs['reps']
    steps = inputs['seq']
    print('Do not press the Stop button now... Otherwise you may stop the code without properly ending a marker.')
    for n in range(0, reps):
        verify_beads()
        print(f'Repetition {n + 1}...', end='')
        timeline.mark_begin(f'{label}_rep{n + 1}')
        for step in steps:
            _, time, speed, position = step
            mirror1.move_to(x=position, speed=speed)
            pause(time)
        timeline.mark_end(
            export=True)  # remove export=True for Bluelake v1.6.1 as it is not implemented and raises an error.
        print('Done.')


def move_to_channel(channel_info):
    """
    Move to another flow channel in the chamber.
    :param channel_info: Contains channel_name and incubation_time.
                         The channel_name must be a waypoint in Bluelake.
                         incubation_time in seconds.
    """
    channel_name, incubation_time = channel_info
    print(f'Moving to channel {channel_name}...', end='')
    microstage.move_to(channel_name)
    print('Done.')
    print(f'Pausing for {incubation_time} seconds...', end='')
    pause(incubation_time)
    print('Done.')
