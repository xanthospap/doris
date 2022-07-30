#! /usr/bin/python

from pyatmos import download_sw_jb2008,read_sw_jb2008

# Download or update the space weather file from https://sol.spacenvironment.net
swfile = download_sw_jb2008()
# Read the space weather data
swdata = read_sw_jb2008(swfile)

from pyatmos import jb2008
# Set a specific time and location
t = '2014-07-22 22:18:45' # time(UTC)
lat,lon,alt = 25,102,600 # latitude, longitude in [degree], and altitude in [km]
jb08 = jb2008(t,(lat,lon,alt),swdata)
print(jb08.rho) # [kg/m^3]
print(jb08.T) # [K]
