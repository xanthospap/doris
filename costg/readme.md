# Results against the [COST-G](https://cost-g.org/) Benchmark Test

For information on the benchmark test, see the paper 
[Benchmark data for verifying background model implementations in orbit and gravity field determination software](https://adgeo.copernicus.org/articles/55/1/2020/)

## 01 earth rotation

### 01earthRotation_interpolatedEOP
```
$>costg/check-01eop.out \
 costG/models/eopc04_14_IAU2000.62-now \
 costG/satellite/01earthRotation_interpolatedEOP.txt \
 | costg/plot_costg_01earthRotation_interpolatedEOP.py
```
![alt text](figures/01earthRotation_interpolatedEOP.png)

### 01earthRotation_rotaryMatrix
```
  $>costg/check-01rot.out \
  costG/models/eopc04_14_IAU2000.62-now \
  costG/satellite/01earthRotation_rotaryMatrix.txt \
  | costg/plot_costg_01earthRotation_rotarymatrix.py
```
![alt text](figures/01earthRotation_rotaryMatrix.png)

## Accelerations

### 02gravityfield_itrf
```
$> costg/check-gravity-field.out \
    costG/models/EIGEN6-C4.gfc \
    costG/satellite/02gravityfield_itrf.txt \
    costG/satellite/00orbit_itrf.txt \
    | costg/plot_costg_02gravityfield_itrf.py
```
![alt text](figures/02gravityfield_itrf.pdf)
Note that $GM$ and $Re$ constants (for Earth) are extracted from the `EIGEN6-C4.gfc` 
data file (header).

### 03directTide[Sun|Moon]_icrf
```
$> costg/check-third-body.out \
    costG/satellite/03directTideMoon_icrf.txt \
    costG/satellite/03directTideSun_icrf.txt \
    costG/satellite/00orbit_icrf.txt \
    data/jpl/de421.bsp data/jpl/naif0012.tls \
    | costg/plot_costg_03directTides.py
```
![alt text](figures/03directTideMoon_icrf.pdf)
![alt text](figures/03directTideSun_icrf.pdf)
Note that $GM$ for Sun and Moon are extracted from the benchmark documentation, i.e. 
`00README_simulation.txt`. Postion are extracted using 
[c-spice](https://naif.jpl.nasa.gov/naif/toolkit.html)

## Other Tests

### ITRF to ICRF transformation
Here, we get the positions (cartesian, ITRF) from the input file `00orbit_itrf.txt`, 
and transform each position vector to ICRF (cartesian). We compare the results 
obtained, with the file `00orbit_icrf.txt`.

```
  $>costg/check-itrf2icrf.out \
  costG/eopc04_14_IAU2000.62-now \
  costG/00orbit_itrf.txt costG/00orbit_icrf.txt \
  | costg/plot_costg_itrf2icrf.py
```
![alt text](figures/00itrf2icrf.png)
