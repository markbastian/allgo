# Independent references for Meeus's tables

The series under `allgo.astro` that come from Meeus's *Astronomical
Algorithms* (converted from a port of the book) are tested against the
book's worked examples. These files check them against sources that owe
nothing to Meeus, each within the accuracy Meeus states for his method.
All are works of the US government, retrieved in September 2026.

| File | Source | Contents |
|---|---|---|
| `horizons-pluto.txt` | JPL Horizons (https://ssd.jpl.nasa.gov/horizons/), DE441 | Pluto (999), heliocentric, ecliptic and equinox of J2000, at 20 dates 1890-2099: JD (TT), x y z (au) |
| `horizons-jupiter-moons.txt` | JPL Horizons | Io, Europa, Ganymede, Callisto (501-504) from Jupiter's center, ICRF: body, JD (TT), x y z (km), 20 dates 1950-2050 |
| `horizons-saturn-moons.txt` | JPL Horizons | Mimas to Iapetus (601-608) from Saturn's center, the same way |
| `usno-seasons.txt` | US Naval Observatory, https://aa.usno.navy.mil/api/seasons | equinoxes and solstices 1972-2040: year month day hh:mm (UT) |
| `usno-moon-phases.txt` | US Naval Observatory, https://aa.usno.navy.mil/api/moon/phases/year | every phase in every fourth year 1972-2040: year month day hh:mm (UT) phase |
| `nasa-solar-eclipses.txt`, `nasa-lunar-eclipses.txt` | Eclipse predictions by Fred Espenak, NASA's GSFC: the *Five Millennium Canon* catalogs for 2001-2100, https://eclipse.gsfc.nasa.gov/ | each eclipse's catalog line: number, date, time of greatest eclipse (TD), delta T, lunation, saros, type, gamma, magnitudes and more |

The Moon's series (chapter 47) is checked separately, against ERFA's
eraMoon98, in `test/allgo/meeus_moon_test.clj`.
