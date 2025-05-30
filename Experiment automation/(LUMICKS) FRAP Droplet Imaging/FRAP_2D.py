"""
BSD 2-Clause License

Copyright (c) 2024, Rachel Leicher
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
"""

from bluelake import timeline, excitation_lasers, pause, confocal

# just to make sure the marker can start fresh
confocal.abort_scan()
timeline.mark_end()
pause(0.5)


# setting up all the parameters for the FRAP experiment
photobleaching_power = 50 # percent
photobleacing_time = 2 # seconds

frap_imaging_power = 10 # percent
frap_imaging_time = 30 # seconds


# Starting the experiment, first with photobleaching at high laser power
timeline.mark_begin('')
pause(0.5)

# Take 2D scan

print('Taking initial scan at'.format(frap_imaging_power) )

excitation_lasers.red = frap_imaging_power
confocal.start_scan("FRAP_2D_Scan")
pause(0.5)

while confocal.is_scanning:
    pause(0.1)

confocal.abort_scan()

# choose point scan area, then start point scan
excitation_lasers.red = photobleaching_power
print (excitation_lasers.red)

confocal.start_scan("FRAP_Point_Scan")

print('Photobleaching at high laser power ({:.0f}%) for {:.0f} seconds'.format(photobleaching_power, photobleacing_time) )

pause (photobleacing_time)

confocal.abort_scan()

# Here you turn the laser power down and look for the fluorescence recovery

print('Turning laser power down to {:.0f}% for {:.0f} seconds'.format(frap_imaging_power, frap_imaging_time) )

excitation_lasers.red = frap_imaging_power
print (excitation_lasers.red)

confocal.start_scan("FRAP_2D_Scan_Continuous")

pause(frap_imaging_time)

print(excitation_lasers.red)

confocal.abort_scan()

# automatically exports the scan
# timeline.mark_end(export=True())

timeline.mark_end()

# make sure the scan is off
excitation_lasers.red = 0
print(excitation_lasers.red)

print('Done')