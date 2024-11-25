<h1>Correction of Tweezer Calibration Factors (CTCF)</h1><br>
<b>Author</b>: Dieter Kamp<br>
<b>Research subjects</b>: Other applications<br>
<br>
<p>This script is available on the StiglerLab <a href=https://github.com/StiglerLab/CTCF>GitHub page</a>.</p>
<br>
This code provides functions to correct miscalibration of optical tweezer measurements using the force noise signal.  
A detailed description of the theory behind this code has been published in Freitag et al. J. Chem. Phys. 155, 175101 (2021). https://doi.org/10.1063/5.0063690


## Use
The main class provided by CTCF is Trace. It contains all data required for correction, which is performed using the correct() method.
```python
from CTCF import Trace
example_trace = Trace()

##################################
### Add data to example_trace  ###
##################################

example_trace.correct()
```

## Adding data
Data is added to the Trace object by assigning it to the following attributes
+ force: Differential force signal
+ force_mob: Force signal of mobile trap
+ force_fix: Force signal of fixed trap
+ ext_orig: Extension
+ dist: Distance
+ stdev: Noise of differential force signal
+ stdev_mob: Noise of mobile force signal
+ stdev_fix: Noise of fixed force signal
+ bead_diameter: Bead diameter, default can be set at the top of trace.py (DIAM_BEAD)
+ bead_diameter1: Bead diameter in mobile trap: defaults to DIAM_BEAD
+ bead_diameter2: Bead diameter in fixed trap: defaults to DIAM_BEAD
+ k_mobile: Trap stiffness of mobile trap
+ k_fixed: Trap stiffness of fixed trap
+ bead: Bead selection, defaults to 0 (both beads), other options 1 (mobile) and 2 (fixed)
+ filters: string of filters
+ parameters: dictionary of filter parameters


By default the following filters are implemented in the psd_filter.py file.

| Filter | Filter command | Parameter |
| ------ | ------ | ------ |
| Bessel filter (Order 8) | bessel8 | Cutoff frequency|
| Boxcar filter | boxcar | averaging window size|
| Butterworth filter | butterworth | Cutoff frequency|
| Subsampling | subsample, sample, ss | Downsampling factor|
| Response of quadrant photodiode | qpd | - |
| Response of National Instruments NI 447x filter| ni447x | - |

The filter_string contains the applied filters in order, as well as parameters relevant to the filter. The filters and parameters are separated by ",". 
In the case of multiple applied filters ";" is used to separate the different filters. The read_filter() method reads the filter string into the Trace object.
```python
filter_string = "filter1,parameter1;filter2,parameter2"
example_trace.read_filter(filter_string)
```

The correction factors are stored in the Trace class attributes
* beta_dagger1
* beta_dagger2
* k_dagger1
* k_dagger2

## Example
The correction function performs loading of data from a  .csv or .xlsx file, setting the attributes of a trace object and running the correct() method. The corrected trace object is then returned.

```python
from CTCF import correction
corrected_trace = correction(filename, k1, k2, [filter_string, sheet_name])
```
The correction function assumes the data to be stored in the following order
```
Force | Stdev | Distance | Force_mob | Force_fix | Stdev_mob | Stdev_fix
```

## Contact
Please direct any inquiries about the code to Dieter Kamp (kamp@genzentrum.lmu.de).
