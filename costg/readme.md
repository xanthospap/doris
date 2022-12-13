# Results against the [COST-G](https://cost-g.org/) Benchmark Test

For information on the benchmark test, see the paper 
[Benchmark data for verifying background model implementations in orbit and gravity field determination software](https://adgeo.copernicus.org/articles/55/1/2020/)

## 01 earth rotation

### 01earthRotation_interpolatedEOP
```
$>costg/check-01eop.out \
 costG/eopc04_14_IAU2000.62-now \
 costG/01earthRotation_interpolatedEOP.txt \
 | costg/plot_costg_01earthRotation_interpolatedEOP.py
```
![alt text](figures/01earthRotation_interpolatedEOP.png)

### 01earthRotation_rotaryMatrix
```
  $>costg/check-01rot.out \
  data/EOP_14_C04_IAU2000A_one_file_1962-now.txt \
  costG/01earthRotation_rotaryMatrix.txt \
  | costg/plot_costg_01earthRotation_rotarymatrix.py
```
![alt text](figures/01earthRotation_rotaryMatrix.png)

## Accelerations

### 02gravityfield_itrf
![alt text](figures/02gravityfield_itrf.png)

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
