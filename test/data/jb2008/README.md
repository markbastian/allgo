# JB2008 comparison values

`orekit-12.2-grid.txt` holds densities computed by Orekit 12.2 (Apache
License 2.0; https://www.orekit.org), a Java library with its own JB2008
implementation, called through its raw `JB2008.getDensity(dateMJD, sunRA,
sunDecli, satRA, satLat, satAlt, f10, f10B, s10, s10B, xm10, xm10B, y10,
y10B, dstdtc)`. Only these output numbers are used here -- no Orekit code
was read or copied -- as an independent check of `allgo.astro.jb2008`,
which is implemented from the published papers alone.

Each line is one point, whitespace-separated:

    mjd  dec  lat  lst  alt  f10 f10b  s10 s10b  m10 m10b  y10 y10b  dtc  rho

MJD (UT); the Sun's declination, latitude (deg); local solar time
(hours; the Sun's right ascension was 0 and the satellite's 15 (lst - 12)
degrees); altitude (km); the solar indices; the geomagnetic temperature
rise dTc (K); and the density (kg/m^3).

The file keeps 220, 300 and 400 km: 648 points over three dates in 2006,
low, moderate and high solar activity, dTc 0 and 60 K, four latitudes and
three local times. There the two implementations agree to within 0.2% on
average and 7% at worst. Outside that band they diverge -- most at the
southern winter pole in July, where Jacchia 1970's helium bulge dominates
(up to +99% at 900 km), and by a factor of several at 2000 km -- in ways
the published description does not settle; see `allgo.astro.jb2008/limits`.
