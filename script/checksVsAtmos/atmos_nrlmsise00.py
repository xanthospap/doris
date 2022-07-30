#! /usr/bin/python

from pyatmos import download_sw_nrlmsise00,read_sw_nrlmsise00
from pyatmos import nrlmsise00

# Download or update the space weather file from www.celestrak.com
swfile = download_sw_nrlmsise00() 

# Read the space weather data
swdata = read_sw_nrlmsise00(swfile) 

t = '2014-07-22 22:18:45' # time(UTC) 
lat,lon,alt = 25,102,600 # latitude, longitude in [degree], and altitude in [km]
nrl00 = nrlmsise00(t,(lat,lon,alt),swdata)
print(nrl00.rho) # [kg/m^3]
print(nrl00.T) # [K]
print(nrl00.nd) # composition in [1/m^3]
