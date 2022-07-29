# Atmospheric Models for Solar Radiation Pressure

## NRLMSISE-00

NRLMSISE-00[^1] is an empirical, global reference atmospheric model of the Earth 
from ground to space. It models the temperatures and densities of the atmosphere's 
components. A primary use of this model is to aid predictions of satellite 
orbital decay due to atmospheric drag. The model also predicts anomalous oxygen. 

[^1]: J. M. Picone, A. E. Hedin, D. P. Drob, NRLMSISE-00 empirical model of the 
atmosphere:Statistical comparisons and scientific issues, JOURNAL OF GEOPHYSICAL 
RESEARCH, VOL. 107, NO. A12, 2002

### Implementation (C++)

The C++ implementation here is designed in an OOP way. It is clear that the 
original FORTRAN code has some kind of "state" anyway, so this is moved to an 
encapsulated design. All of the original switches/options are kept "as-is" 
but basically the only interesting thing for the project is the units 
transformation (aka results in m<sup>3</sup> and kg/m<sup>3</sup> instead of 
cm<sup>3</sup> and gm/cm<sup>3</sup>). To get the output in meters-related units, 
set the corresponding variable `meters_` to `true` in the corresponding 
`dso::nrlmsise00::InParams` instance (call `dso::nrlmsise00::InParams::meters_on()`).

### (Unit) Tests

Using the FORTRAN implementation, i have stacked a series of results based on 
random test cases, including the default test cases used in the original 
FORTRAN distribution. Source code is in `unit_tests/test_msise00.cpp`; the 
program should **NOT** fail.

### Input Data

Input data can be:
  
  * ap - magnetic index(daily), triggered when the input switch sw[8]=-1; the
    function(s) will expect relevant indexes in the `nrlmsise00::ApArray` 
    instance of the `dso::nrlmsise00::InParameters` used to request results. 
    I don't really know wtf this is.

  * F107A - 81 day average of F10.7 flux (centered on day), and

  * F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY

The last two are passed in via the relevant variables in the 
`dso::nrlmsise00::InParameters`. F107 and F107A values used to generate the 
model correspond to the 10.7 cm radio flux at the actual distance of the Earth
from the Sun rather than the radio flux at 1 AU.

Input data can be retrived from the following sources:

  * [NWRA](https://spawx.nwra.com/spawx/env_latest.html)

  * [CelesTrack](https://celestrak.org/SpaceData/) Parsing of CSV files 
   containing flux data from this repository is supported, via the 
   `dso::get_CelesTrack_flux_data` function. First, users need to download 
   the respective CSV file, and the the function can be used to parse it and 
   return the relevant quantities. See also: 
   `dso::utils::celestrack::details::CelestTrackSWFlux` which is a class 
   meant to hold flux data from this resource.

### Usage

For the most part, we will want to use:
  
  * an `dso::nrlmsise00::InParameters` instance, with all switches set to 1. 
    This is by default performed by the constructor of the class.

  * If the densities are requested in m<sup>3</sup> and kg/m<sup>3</sup> 
    (usually the case), call `dso::nrlmsise00::InParameters::meters_on()` on 
    the instance.
  
  * an `dso::nrlmsise00::OutParameters` to hold results.

  * Call the function `dso::Nrlmsise00::gtd7d()` on an already constructed 
    `dso::Nrlmsise00` (defualt constructed) instance.

Note that `UT`, `Local Time`, and `Longitude` are used independently in the
model and are not of equal importance for every situation. For the most 
physically realistic calculation these three variables should be consistent 
(STL=SEC/3600+GLONG/15). The Equation of Time departures from the above 
formulafor apparent local time can be included if available but are of minor 
importance.

### See also

  * https://github.com/graziano-giuliani/Meteostuff/tree/master/NRLMSIS_F90
  
  * https://github.com/alesmorse/nrlmsise-00

  * https://www.brodo.de/space/nrlmsise/
