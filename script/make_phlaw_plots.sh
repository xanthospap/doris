#! /bin/bash

# exit when any command fails
set -e

../test/testAntennaPcv > phase_law

## Alcatel 2GHz
gnuplot -persist <<-EOFMarker
set term postscript landscape enhanced color
set output "alcatel-phlaw.ps"

set xlabel "Zenith angle (degrees)"
set ylabel "Phase correction for 2GHz in mm"
set xrange [0:90]

set title "Alcatel Ground Antenna/Beacon Phase Law"
plot "phase_law" u 1:2 w lp lw 3 pointtype 7 pointsize 0.9 title "doris\\\_phase\\\_law\\\_antex\\\_alcatel"
EOFMarker

## Starec B/C 2GHz
gnuplot -persist <<-EOFMarker
set term postscript landscape enhanced color
set output "starecbc-phlaw.ps"

set xlabel "Zenith angle (degrees)"
set ylabel "Phase correction for 2GHz in mm"
set xrange [0:90]

set title "Starec B/C Ground Antenna/Beacon Phase Law"
plot "phase_law" u 1:4 w lp lw 3 pointtype 7 pointsize 0.9 title "doris\\\_phase\\\_law\\\_antex\\\_starecBC"
EOFMarker
