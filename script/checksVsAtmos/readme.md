# Checking atmospheric models implemenation agains [ATMOS](https://github.com/lcx366/ATMOS)

[ATMOS](https://github.com/lcx366/ATMOS) is a Python library, containing 
implementations for both the NRLMSISE00 and the JB2008 models. It also 
supports downloading and parsing Space Weather data for retrieving the 
relevant quantities needed by the two models.

Usage is shown at the project's homepage. Note that ATMOS uses the 
[Numba](https://numba.readthedocs.io/en/stable/index.html) jit. Normally, 
`pip install numba` will do.
