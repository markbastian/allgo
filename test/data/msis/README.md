# NRLMSISE-00 test output

`nrlmsise00_output.txt` is copied unchanged from the NRLMSISE-00 Fortran
distribution in NASA's CCMC ModelWeb archive,
https://git.smce.nasa.gov/ccmc-share/modelwebarchive (`MSIS/NRLMSIS00/`).
The model is by Picone, Hedin and Drob of the Naval Research Laboratory,
a work of the US government.

The file is the output of the distribution's test driver,
`nrlmsise00_driver.for`, which calls GTD7 for seventeen cases:

- fifteen with the daily Ap, the inputs in the summary tables at the end,
  where each output is given to four significant figures;
- two with the 3-hour ap history (all seven values 100, switch 9 set to
  -1), at 400 km and 100 km and otherwise the first case's inputs.

Each case is also printed, to three figures, in the blocks at the top:
the eight number densities, then anomalous oxygen, the exospheric
temperature and the temperature at altitude, then two lines of internal
values (`DL`) not tested here. Densities are in cm^-3 and g/cm^3.

The Fortran is single precision: at 100 km its anomalous oxygen
underflows to zero where double precision gives some 1e-42 cm^-3.
