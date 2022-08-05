# Handling Earth Orientation Parameters

IERS routinely distributes EOP data files [^1] (e.g. in 
[Bulletin B](https://datacenter.iers.org/productMetadata.php?id=207)/C04 files)
containing (among other things) daily final values at 0:00 UT of x, y, UT1-UTC, 
dX, dY, and their uncertainties. The time span is one month with final values, 
one month with preliminary values.

To use these data, we first transform it to a different format, the so called 
`EopFile`s; these are practically the same as the original IERS Bulletin B 
files, but containing only the Sextion 1 data.

To construct such a file for a requested date, use the script 
[fetch_iers_c04.py](https://github.com/xanthospap/doris/blob/main/script/fetch_iers_c04.py) 
which will download the files needed and extract relevant data. The result 
file can then be used to get/handle EOP data.

`EopFile`s can be used to extract EOP data and fill a respective `EopLookUpTable` 
instance. EOP data handling is performed via this instance.

The most useful operation is the interpolation and correction of EOP values for 
a requested date in order e.g. to use them for constructing ICRF top GCRF 
transformation matrix. According to [^2] (see chapter 5.5), eop data to be used 
should be:

> those published by the IERS with additional components to account for the 
> effect of ocean tides and for forced terms with periods less than two days 
> in space
the latter are called *libration effects*.

* ocean tides effects are the diurnal and semi-diurnal variations in pole 
coordinates caused by ocean tides

* libration effects are the variations in pole coordinates corresponding to 
motions with periods less than two days in space that are not part of the 
IAU 2000 nutation model.

Hence, the following should be performed before using EOP:

  1. **Account of ocean tidal and libration effects in pole coordinates xp and yp** 
  (via `iers::iers2010::interp_pole`). This will correct for subdaily variations 
  via interpolation and adding tidal terms and diurnal components 

  2. **Variations in polar motion due to ocean tides** aka diurnal and 
  semi-diurnal variations in pole coordinates caused by ocean tides (via 
  `iers2010::ortho_eop`)

  3. **Variations in polar motion due to libration**; diurnal and semi-diurnal 
  nutations, designated here as *libration*, originate from the direct effect 
  of the external (mainly luni-solar) torque on the non-axisymmetric part of 
  the Earth as expressed by the non-zonal terms of the geopotential. The diurnal 
  components are computed via `iers2010::pmsdnut2` 

Once a `EopLookUpTable` is initialized, these corrections are applied with a 
call to `EopLookUpTable::interpolate`.

[^1]: https://www.iers.org/IERS/EN/Publications/Bulletins/bulletins.html
[^2]: IERS Conventions (2010). Gérard Petit and Brian Luzum (eds.). (IERS Technical Note ; 36) Frankfurt am Main: Verlag des Bundesamts für Kartographie und Geodäsie, 2010. 179 pp., ISBN 3-89888-989-6
