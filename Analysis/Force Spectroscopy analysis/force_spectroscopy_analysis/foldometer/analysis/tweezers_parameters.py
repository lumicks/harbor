#!/usr/bin/env python
# -*- coding: utf-8 -*-

from datetime import date
import operator
from collections import OrderedDict as od

#Default values, corresponding to two 1.3kb DNA handles and single MBP
PROTEIN_LENGTHS = {"MBP": 120, "4MBP": 480, "4MBP_casein":560, "dmMBP": 120, "LUCI": 200, "LUCIFERASE": 200, "GR": 90,
                   "P53": 80, "2MBP": 240, "2MBP_pB1": 240, "casein":70, "rubisco":138, "Rubisco":138, "YPet":105}

CCD_FREQUENCY = 90
MIRRORVOLTDISTANCEFACTOR = 0.175#0.1725#0.165753#0.123697 #0.10604#0.1176# or 0.00018155565
# MIRRORVOLTDISTANCEFACTOR=1/(CCD_PIXEL_NM * VOLTS_PIXELS_SLOPE)
CCD_PIXEL_NM = 0.00804
VOLTS_PIXELS_SLOPE = 710 #D(V) = VOLTS_PIXELS_SLOPE * D(pixel) - 255650

STATMIRRORXUNORDDICT = {date(2014, 10, 1): -31283,
                        date(2016, 11, 21): -32059,
                        date(2017, 1, 17): -27583,
                        date(2017, 8, 25): -24199,
                        }
STATMIRRORYUNORDDICT = {date(2014, 10, 1): -13551,
                        date(2016, 11, 21): -10417,
                        date(2017, 1, 17): -5195,
                        date(2017, 8, 25): 1055,
                        }

STATMIRRORXDICT = od(sorted(STATMIRRORXUNORDDICT.items(), key=operator.itemgetter(0)))
STATMIRRORYDICT = od(sorted(STATMIRRORYUNORDDICT.items(), key=operator.itemgetter(0)))


def get_mirror_values(header):
    """
    Function to find the correct stationary mirror coordinates depending on date of measurement

    Args:
        header: the metadata from the file, containing the date of the file

    Returns:
        stationaryMirrorCoords (tuple): a tuple with the x and the y coordinates from the correct date
    """
    fileDateStr = str(header["fileName"][19:29], 'utf-8')
    fileDate = date(int(fileDateStr[0:4]), int(fileDateStr[5:7]), int(fileDateStr[8:10]))
    for index, key in enumerate(STATMIRRORXDICT.keys()):
        if fileDate >= key:
             stationaryMirrorX = STATMIRRORXDICT[key]
             stationaryMirrorY = STATMIRRORYDICT[key]
    #print(stationaryMirrorX, stationaryMirrorY)
    return (stationaryMirrorX, stationaryMirrorY)

