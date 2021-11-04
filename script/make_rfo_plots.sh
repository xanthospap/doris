#! /bin/bash

# exit when any command fails
set -e

cd ../data
for sat in Cryosat-2 HY-2C Jason-3  Saral  Sentinel-3A  Sentinel-3B  Sentinel-6A ; do
  cd $sat
  echo "Fitting data for satellite $sat"
  files=$(ls *2100[0-9].001)
  ../../test/testRinexRfoFit ${files} | grep -v "Fit" | awk '/^[0-9]/ {print}' > ../rfo.$sat
  ../../test/testRinexRfoFit ${files} | grep "Fit" | awk '{print $2,$3}' > ../fit.$sat
  cd ../
done

for sat in Cryosat-2 HY-2C Jason-3  Saral  Sentinel-3A  Sentinel-3B  Sentinel-6A ; do
  echo "Plotting data for satellite $sat"
  dfn=rfo.${sat}
  ffn=fit.${sat}
  vel=$(tail -1 $dfn | awk '{printf "%.3f", $6}')
  dt1=$(head -1 $dfn | awk '{printf "%.3f", $1}')
  dt2=$(tail -1 $dfn | awk '{printf "%.3f", $1}')
  lx=$(python -c "print('{:.3f}'.format(${dt1} + 0.75*(${dt2} - ${dt1})))")
  dy1=$(head -1 $dfn | awk '{printf "%.3f", $2}')
  dy2=$(tail -1 $dfn | awk '{printf "%.3f", $2}')
  ly=$(python -c "print('{:.3f}'.format(${dy1} + 0.25*(${dy2} - ${dy1})))")
gnuplot -persist <<-EOFMarker
set term postscript landscape enhanced color
set output "${sat}-RinexRfo.ps"

set ylabel "(f-f0)/f0 in 1e-11"
set xlabel "Modified Julian Day"
set xrange [${dt1}:${dt2}]

set title "Relative Frequency Offset from RINEX for $sat"
set label 1 "velocity: ${vel} / day" at ${lx},${ly} center
plot "${dfn}" u 1:2 w lp lw 2 pointtype 7 pointsize 0.7 title "RINEX value", \
  "${ffn}" u 1:2 w l lw 3 title "linear fit"
EOFMarker
done
