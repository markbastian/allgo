# Cross-check values from Orekit

Values computed with Orekit 12.2 (Apache License 2.0;
https://www.orekit.org), a Java space-dynamics library written separately
from this one, as an independent check on models that are otherwise
tested only against their own source's output. Only these output numbers
are used -- no Orekit code was read or copied. Orekit was run with no data
files: times in TT, and vectors in a stand-in frame, which leaves each
model's raw output untouched.

| File | Orekit | Columns |
|---|---|---|
| `knocke.txt` | `KnockeRediffusedForceModel.computeAlbedo`, `computeEmissivity` | MJD, latitude (deg), albedo, emissivity |
| `sgp4.txt` | `TLEPropagator` over `test/data/sgp4/SGP4-VER.TLE`, each case's start, stop and step | case (index:satnum), minutes, TEME position (km), velocity (km/s) |
| `nrlmsise00.txt` | `NRLMSISE00.getDensity`, a spherical Earth of radius 6378.137 km, the Sun at the longitude UT puts it | MJD, latitude, longitude (deg), altitude (km), F10.7 average, daily F10.7, Ap, Sun longitude, switch 9, density (kg/m^3, GTD7D's, with anomalous oxygen) |
| `harris-priester.txt` | `HarrisPriester.getDensity(sun, position)`, the WGS-84 ellipsoid | position (km), Sun direction, density (kg/m^3) |

Orekit's SGP4 runs in the improved ("i") operation mode. Past a case's
first error code it still returns states, which the tests do not use; and
on three of the verification set's hard cases -- 33333 (e = 0.995), 33335
(near-equatorial geostationary) and 20413's second run, years out -- it
parts from the reference output in `test/data/sgp4/`, which
`allgo.astro.sgp4` matches.
