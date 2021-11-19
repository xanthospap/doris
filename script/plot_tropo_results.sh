#! /bin/bash

# exit when any command fails
set -e

fn=../data/sentinel6a_004_tls_var
sat=Sentinel-6A
rinex="s6arx21004.001"
#fn=../data/cryosat_005_tls_var
#sat=Cryosat
#rinex="cs2rx21005.001"
#fn=../data/saral_001_tls_var
#sat=Saral
#rinex="srlrx21001.001"
  
dt1=$(head -1 $fn | awk '{printf "%.3f", $2}')
dt2=$(tail -1 $fn | awk '{printf "%.3f", $2}')
dt2="59218.4"
#dt2="59215.3"

gnuplot -persist <<-EOFMarker
set term postscript landscape enhanced color
set output "${rinex}-tropo-correction.ps"

set ylabel "Elevation Angle (degrees)"
set xlabel "Modified Julian Day"
set xrange [${dt1}:${dt2}]

set y2tics -2.5, 0.25
set ytics nomirror

set title "Tropospheric Correction for ${rinex}:TLS"
plot "${fn}" u 2:(90-\$9) w lp lw 2 pointtype 5 pointsize 0.7 lt rgb "#000000" title "elevation angle" axis x1y1, \
  "${fn}" u 2:15 w lp lw 2 pointtype 6 pointsize 0.5 lt rgb "#0000FF" title "total tropo delay" axis x1y2
EOFMarker
