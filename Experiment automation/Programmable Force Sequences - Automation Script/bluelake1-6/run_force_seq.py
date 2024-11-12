# coding=utf-8
"""
PROGRAMMABLE FORCE SEQUENCE
Automation script for Bluelake 1.6 (C-Trap, Lumicks)
Sara de Bragança, sbraganca@cnb.csic.es
2021-09-03

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

''' INPUT FROM THE USER '''
# Identification label for your experiment. Will appear in the name of the markers.
name_id = 'lambda_10nM-protein'

# Define the parameters for the force sequence:
forces = [0, 55, 0]  # [pN] Target forces (F)
times = [0, 0, 0]  # [s] Time to pause at each F
speeds = [1, 1]  # [um/s] The speed to go from F1 to F2, F2 to F3, ...
repetitions = 5  # Number of times to repeat the sequence

# Repeat the exact same sequence (and repetitions) in another channel:
repeat_in_another_channel = True  # True to repeat in another channel. Otherwise False.
channel_name = 'Protein'  # Exact name of the waypoint (Bluelake user interface)
incubation_time = 60  # [s] Time to wait in the new channel before continuing execution.

"""
Use your creativity to design the best sequence for your experiment.
Some examples:

RAMP-UP (for example: to study the resistance of a protein to force)
forces = [0, 50]
times = [5, 5]
speeds = [1] or even slower

RAMP-DOWN (for example: to study condensation)
forces = [50, 0]
times = [5, 0]
speeds = [1] or even slower

FORCE CYCLES
forces = [10, 50]
times = [10, 10]
speeds = [10]
repetitions = 20

GRADIENT (for example: to evaluate the working force ranges)
forces = [10, 20, 30, 40, 50]
times = [5, 5, 5, 5, 5]
speeds = [2, 2, 2, 2]

"""

''' RUN THE WORKFLOW '''
import mod_force_seq as mod
mod.apply_steps_workflow(name_id, forces, times, speeds, repetitions,
                         repeat_in_another_channel, channel_info=(channel_name, incubation_time))
