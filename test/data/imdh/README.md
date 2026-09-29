# NASA Interplanetary Mission Design Handbook: Earth to Mars, 2026-2045

Transcribed by hand from Burke, Falck and McGuire, *Interplanetary
Mission Design Handbook: Earth-to-Mars Mission Opportunities 2026 to
2045*, NASA/TM-2010-216764 (Glenn Research Center, October 2010), a work
of the US government and so in the public domain:
https://ntrs.nasa.gov/citations/20100037210

- `earth-mars-2026-2045.edn`: tables 2-11, the ballistic energy minima of
  each opportunity -- least C3 and least Mars arrival excess speed, for
  type I and type II transfers -- with the launch asymptote's right
  ascension and declination.
- `earth-mars-2005.edn`: appendix A's tables 12 and 13, the 2005 minima
  as the handbook's own program (MIDAS) and its reference 6 (JPL) give
  them.

The rows are as printed. The handbook computes each point on a two-day
grid of dates, at 0h, and prints two decimals of C3; its asymptote
angles are referred to the mean equator and equinox of the launch date.
`mission_handbook_test` checks every row against
`allgo.astro.interplanetary/transfer`, and in doing so found these slips
in the tables, which it corrects and says why:

- 2035, 2037 and 2039 (tables 6-8): every date is early, by 63, 81 and 75
  days. Shifted, all twelve rows agree in all four quantities, and the
  shifted minima are where the handbook's own contour plots put them.
- 2028 (table 3), the type I least-C3 row: its declination, 1.581, is the
  row below's.
- 2041 (table 9), the type I least-arrival-speed row: a type II transfer.
- 2043 (table 10), the type I least-arrival-speed row: its right
  ascension and declination are the row above's.
