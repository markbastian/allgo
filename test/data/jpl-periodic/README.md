# JPL Three-Body Periodic Orbits: Earth-Moon sample

A sample of the periodic orbits JPL's Solar System Dynamics group
catalogs in the circular restricted three-body problem, served by its
Three-Body Periodic Orbits API
(https://ssd-api.jpl.nasa.gov/periodic_orbits.api; documentation at
https://ssd-api.jpl.nasa.gov/doc/periodic_orbits.html), retrieved
September 2026:

- `family=halo&libr=1&branch=N` and `family=halo&libr=2&branch=N`: four
  L1 and three L2 northern halo orbits, at heights z of about 0.02, 0.05,
  0.1 and 0.15, those starting on the family's side of the point;
- `family=lyapunov&libr=1` and `family=lyapunov&libr=2`: three planar
  Lyapunov orbits each, small to large.

`earth-moon.edn` keeps each orbit's start in the x-z plane (x, z, vy; y,
vx and vz are zero to the catalog's precision), Jacobi constant, period
and stability index as served, with the system's mass ratio, collinear
points and units. `halo_test` checks `allgo.astro.halo` against them.
