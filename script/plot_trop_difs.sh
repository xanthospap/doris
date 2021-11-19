#! /bin/bash

# exit when any command fails
set -e

fn=../data/all_beacons.saral001
sat=Saral

gnuplot -persist <<-EOFMarker
set term postscript landscape enhanced color
set output "${sat}-allbeacons-pressure.ps"
set ylabel "Pressure hPa (RINEX-GPT3)"
set xlabel "Modified Julian Day"
set title "Pressure Values from RINEX vs GPT3 for $sat"
plot "${fn}" u 2:(\$3-\$4) w lp lw 2 pointtype 7 pointsize 0.7
EOFMarker

gnuplot -persist <<-EOFMarker
set term postscript landscape enhanced color
set output "${sat}-allbeacons-ztdh.ps"
set ylabel "Zenith Delay Hydrostatic in meters (RINEX-GPT3)"
set xlabel "Modified Julian Day"
set title "ZTD Hydrostatic from RINEX vs GPT3 for $sat"
plot "${fn}" u 2:(\$5-\$6) w lp lw 2 pointtype 7 pointsize 0.7
EOFMarker

fn="../data/srl_tropo_tlsb"
dt1=$(head -1 $fn | awk '{printf "%.3f", $2}')
dt2=$(tail -1 $fn | awk '{printf "%.3f", $2}')
minp=$(awk 'BEGIN{a=1000}{if ($3<0+a) a=$3; if ($4<0+a) a=$4;} END{printf "%.3f",a-1}' $fn)
echo $minp
maxp=$(awk 'BEGIN{a=1000}{if ($3>0+a) a=$3; if ($4>0+a) a=$4;} END{printf "%.3f",a+1}' $fn)
echo $maxp
#dy1=$(head -1 $dfn | awk '{printf "%.3f", $2}')
#dy2=$(tail -1 $dfn | awk '{printf "%.3f", $2}')

gnuplot -persist <<-EOFMarker
set term postscript landscape enhanced color
set output "${sat}-TLSB-pressure.ps"
set xrange [${dt1}:${dt2}]
set ylabel "Pressure hPa (RINEX-GPT3)"
set xlabel "Modified Julian Day"
set title "Pressure Values from RINEX vs GPT3 for $sat"
plot "${fn}" u 2:(\$3-\$4) w lp lw 2 pointtype 7 pointsize 0.7
EOFMarker

gnuplot -persist <<-EOFMarker
set term postscript landscape enhanced color
set output "${sat}-TLSB-pressure12.ps"
set xrange [${dt1}:${dt2}]
set yrange [${minp}:${maxp}]
set ylabel "Pressure hPa"
set xlabel "Modified Julian Day"
set title "Pressure Values from RINEX vs GPT3 for $sat"
plot "${fn}" u 2:3 w lp lw 2 pointtype 7 pointsize 0.7 title "RINEX value", \
"${fn}" u 2:4 w lp lw 2 pointtype 7 pointsize 0.7 title "GPT3 value"
EOFMarker

gnuplot -persist <<-EOFMarker
set term postscript landscape enhanced color
set output "${sat}-TLSB-ztdh.ps"
set xrange [${dt1}:${dt2}]
set ylabel "Zenith Delay Hydrostatic in meters (RINEX-GPT3)"
set xlabel "Modified Julian Day"
set title "ZTD Hydrostatic from RINEX vs GPT3 for $sat"
plot "${fn}" u 2:(\$5-\$6) w lp lw 2 pointtype 7 pointsize 0.7
EOFMarker
