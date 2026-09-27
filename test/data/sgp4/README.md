# SGP4 verification data

Copied unchanged from the code package for Vallado, Crawford, Hujsak and
Kelso, *Revisiting Spacetrack Report #3* (AIAA 2006-6753),
https://celestrak.org/publications/AIAA/2006-6753/ (`AIAA-2006-6753.zip`).
The package carries no license; its FAQ asks that users cite it and link to
that page.

| File | From | Made by | opsmode |
|---|---|---|---|
| `SGP4-VER.TLE` | `sgp4/cpp/testsgp4/` | the paper's verification cases, with each case's start, stop and step in minutes after column 69 | |
| `cpp/*.e` | `sgp4/cpp/testsgp4/TestSGP4/` | Vallado's C++ (`TestSGP4.cpp`), STK ephemeris files, time in seconds | a |
| `tmatverDec2015.out` | `sgp4/mat/` | Vallado's MATLAB (`testmat.m`), December 2015 | a |
| `java_sgp4_ver.out` | `sgp4/java/JAVA_SGP4_v2/` | Shawn Gano's 2009 Java conversion of the C++, shipped in the package | i |

`20413` appears twice in `SGP4-VER.TLE`; the C++ writes its `.e` per
satellite number, so `cpp/20413.e` holds the second run only.
