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
