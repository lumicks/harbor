"""
BSD 2-Clause License

Copyright (c) 2024, Yasaman Madraki
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

# This code creates a single multipage scan by doing z-Scanning. 
# please set your scanning window the way you like and set the stage position at the lowest level you want to scan. 
# please set the scanning mode on "continuous" in the scan setting window
# please enter the scanning time, step size and the number of steps you would like to have in this script
# when you run the code, scanning starts and after scanning each plane the nanostage moves upward and does another scan
# NOTICE: please be mindful of the total height you are scanning. Nanostage should not move too high to hit the condenser.
import bluelake
from bluelake import trap1, trap2, trap12z, nanostage, timeline, confocal


step_size=0.1           # um
step_number=20          # number of planes you want to scan
scan_time=4             # set your scanning window and enter the scantime from scan setting window here
confocal.start_scan()   # make sure the scanning mode is set on the "continuous" mode in the scan setting window
bluelake.pause(scan_time)
for i in range(1,step_number,1):
    nanostage.move_by(dz=step_size,speed=0)
    bluelake.pause(scan_time)
confocal.abort_scan()
