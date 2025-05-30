# BSD 2-Clause License

# Copyright (c) [2023], [Olivia Yang]
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

# ======
import bluelake as bl
from datetime import datetime

# user inputs =================================
cycles = 1  # number of cycles.
high_t = 20  # minutes of high pressure
low_t = 1  # minutes of low pressure
high_p = 2.0  # high pressure in bar
low_p = 0.4  # low pressure in bar


# functions ===================================
def setpressure(pres):
    if pres == 0:
        bl.fluidics.start_venting()
    else:
        curpres = bl.fluidics.pressure
        if curpres < pres:
            while bl.fluidics.pressure < pres:
                bl.fluidics.increase_pressure()
        else:
            while bl.fluidics.pressure > pres:
                bl.fluidics.decrease_pressure()


def timeprogression(time, numsteps):
    stepsize = round(time / numsteps, 2)
    for i in range(1, numsteps + 1):
        bl.pause(stepsize)
        print(
            f" {str(round(i * stepsize * 100 / time))}% of time has elapsed ({str(round(i * stepsize / 60, 1))} min)"
        )


# MAIN ------------------------------------------
print("> started at " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
print(high_t, " mins at", high_p, " high pressure")
print(low_t, " mins at", low_p, " low pressure")
print(f"{cycles} cycles.")
try:
    for cycle in range(cycles):

        print(f"Starting cycle {cycle+1} out of {cycles}.")
        print("starting high pressure")
        setpressure(high_p)
        bl.fluidics.open(1, 2, 3, 4, 5, 6)
        timeprogression(60 * high_t, int(high_t))

        print("high done, starting low pressure")
        setpressure(low_p)
        timeprogression(60 * low_t, int(low_t))
except:
    print("Script aborted unexpectedly.")
finally:
    bl.fluidics.stop_flow()
    print("Done and vented")
    print("> finished at " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
