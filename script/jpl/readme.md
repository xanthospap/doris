## Notes on the NAIF/SCPICE JPL library, see https://naif.jpl.nasa.gov/naif/

  * this library is used (its C version that is) to parse kernel files (aka 
    planetary ephemeris and leap-second files) to compute planet position 
    vectors, namely for Sun and Moon.

  * Download the C version (toolkit) of the library from 
    https://naif.jpl.nasa.gov/naif/toolkit.html

  * Before building, we need to wrap the header files in preprocessor macros, 
    so that we can easily call the lib within C++. Use the c2cpp_header.py
    script file to transform **ALL** header files

  * [OPTIONAL] Change the makefiles to use a C++ compiler (instead of the defined C 
    compiler). To do this, change the script named 'mkprodct.csh' in **ALL** 
    folders within the src/ directory to use the default (or any preferred)
    C++ compiler. Aka, find the line:
    set TKCOMPILER  =  "gcc", and change it to
    set TKCOMPILER  =  "g++"

  * Call the top-level directory 'makeall.csh' to build the library

  * Use the script 'install.sh' to install the library and header files for
    all users, aka system-wide. The header files will be installed in the 
    directory: '/usr/local/include/cppspice'.

Documentation on the (C version) of the library, can be found here:
https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/index.html

### Kernels

Required reads: 
[1](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/Tutorials/pdf/individual_docs/34_naif_server.pdf),
[2](https://naif.jpl.nasa.gov/naif/using_spice_data) and
[3](https://naif.jpl.nasa.gov/naif/data_generic.html)

The library uses the so-called CSPICE ``kernels'' to perform a bunch of operations, 
including time-scale transformations[^1], extract state/position of planets and 
relevant constants. These are:

* (SPK) Ephemerides for planets ([spk data](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/))

* (PCK) Planetary constants kernels ([pck data](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/)), and

* (LSK) Current leapseconds kernel ([lsk data](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/))

The above data repositories contain ``readme'' files describing folder trees and 
data file contents; use them!

[^1]: Within the CSPICE library only; in general, time-scale transfromations and 
datetime manipulation is performed via the [datetime](https://github.com/xanthospap/ggdatetime) 
library.