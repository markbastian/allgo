# TODOs

Open work only. The reasoning and measurements behind finished work are in
[NOTES.md](NOTES.md).

## Algorithms to add

- **Thermal erosion.** Talus slippage, to pair with the droplet model in
  `allgo.procedural.erosion` and put scree under the cliffs it leaves bare.
- **Tectonic plates.** Voronoi plates driving mountain ranges and
  coastlines, on `allgo.geometry.dual-mesh` or the planet sphere.
- **States, religions, population and history** (the rest of Azgaar). This
  is simulation on top of a finished map rather than geometry, so it is
  arguably a separate project.

## Physics performance

On a settled 16x16 wall (JVM) a step is 14.7 / 19.6 / 32.9 ms for
sequential impulse / TGS / XPBD, and the solve is the smaller part. In order
of expected return:

1. **XPBD velocity pass.** `solve-xpbd-velocities!` is about 12 ms of
   XPBD's 33 ms step and still allocates a vec3 per contact per use. Move it
   onto primitive doubles, as the position pass was.
2. **Narrow phase on primitive doubles.** The separating axis test in
   `allgo.physics.contact` allocates a vector per dot and cross product.
3. **Persistent manifolds.** Keep contacts across frames: re-run the narrow
   phase only for pairs that moved, and warm start by identity instead of a
   map lookup.
4. **SIMD and graph coloring** (Box2D v3). Last: it speeds up only the solve.

## Stacking

- **Tall columns still topple at the default settings.** Real blocks
  would stand at forty courses; TGS holds about ten. Forty stands with
  stiff contacts, a 1 mm slop and body damping at sixteen substeps --
  four times the cost of a step -- or by starting asleep. Load-scaled
  stiffness and shock propagation were tried and did not pay (NOTES,
  "Tried and not kept"). What is left:
  - **An implicit contact solver.** Contacts as a convex problem solved
    by Newton's method, as MuJoCo and Drake's SAP do, take stiff
    contacts at large steps without the quarter-step-rate cap that
    limits this one. The principled fix, and a solver of its own.
  - XPBD: re-measure penetration against the geometry each substep
    instead of `depth0` plus anchor drift.
  - Settle faster, so sleeping can start before bricks are shed.
- **Revisit manifold reduction together with per-pair stiffness.** Reducing
  to four points is correct and cheap, but it costs stack depth at the
  stiffness that keeps a stack still.

## Mixed scenes (`allgo.physics.world`)

- Only sequential impulse runs here. Missing: substepping, XPBD, sleeping,
  and warm starting across steps. `allgo.physics.solver` has all four.

## Articulated bodies

- A single model cannot close a loop (it is a tree). Pinning two models
  together through `allgo.physics.world` is the workaround the motorcycle
  rider uses.

## Solar system demo

- **The Earth and Moon's shadows.** A lunar eclipse tints the Moon from
  the true geometry, but nothing shows a solar eclipse's shadow on the
  Earth, and Saturn's rings neither cast a shadow on the planet nor fall
  into its shadow. The traced shadows that Jupiter's moons use would do
  for the rings; the Earth and Moon would need drawing to one scale near
  an eclipse.
- **The sky view** shows the planets as points; their phases and disks
  (Venus's crescent, Saturn's rings) could be drawn when zoomed in, as
  the Moon is.
- **More events.** The Events menu lacks Venus's elongations and
  stations, Mercury's and Mars's other stations, and the conjunctions of
  Mars, Jupiter, Uranus and Neptune: the rows of Meeus's table 36.B that
  could not be checked against a second source (NOTES).

## Astronomical and Astrodynamic Calculations

- **Vallado, *Fundamentals of Astrodynamics and Applications*.** Worked
  through in chunks, each built on what `allgo.astro` already has from
  Montenbruck & Gill and Meeus, implemented from the book's equations and
  algorithms. The tests check each result independently: by flying it
  with the two-body propagator, by the identities it must satisfy, and
  against ERFA and published standards where those cover it. SGP4 is
  transcribed from the code that accompanies *Revisiting Spacetrack Report
  #3* and checked against that package's own verification output.
  - Add the book's inline worked examples as unit tests, with their
    printed values, alongside the independent checks.
  1. ~~Two-body core (ch. 2)~~ -- `allgo.astro.universal`.
  2. ~~Maneuvers (ch. 6)~~ -- `allgo.astro.maneuvers`.
  3. ~~Initial orbit determination (ch. 7)~~ -- `allgo.astro.iod`, with
     three independent Lambert solvers: universal variables, Lagrange's
     equation and Battin's.
  4. ~~Coordinates and frames (ch. 3-4)~~ -- `allgo.astro.reduction`,
     `allgo.astro.states` and `allgo.astro.cio`.
  5. ~~SGP4/SDP4 and TLEs~~ -- `allgo.astro.sgp4`.
  6. ~~Celestial (ch. 5)~~ -- `allgo.astro.visibility` for sight, shadow,
     eclipses and naked-eye visibility, and twilight in `allgo.astro.rise`;
     Sun and Moon positions and rise and set were already there.
  7. Perturbations (ch. 8-9): done -- `allgo.astro.perturbations` (Gauss's
     variational equations, orbit-averaged rates, J2's secular rates,
     King-Hele's drag decay, radiation pressure's secular rates) and
     `allgo.astro.gravity` (Pines and the spherical partials) and J2's
     short-period terms as Spacetrack Report No. 3 gives them, and Lear's
     and Gottlieb's normalized gravity algorithms from NASA/TP-2016-218604
     -- and Brouwer's short-period terms in full, with the J2 e terms the
     report drops, his second-order secular rates and his long-period
     terms, in `allgo.astro.brouwer`.
  8. ~~Mission geometry (ch. 11)~~ -- `allgo.astro.mission`, with the
     sun-synchronous inclination in `allgo.astro.perturbations`.
  9. ~~Covariance transformations~~ -- `allgo.astro.covariance`: the
     classical and equinoctial elements' Jacobians analytic, the flight
     elements' numerical.
  10. ~~Interplanetary (ch. 12)~~ -- `allgo.astro.interplanetary`: spheres
      of influence, departure, capture and flyby hyperbolas, Hohmann and
      Lambert transfers between planets.
  11. ~~Encke's method (ch. 8)~~ -- `allgo.astro.encke`.
  12. ~~Satellite pass prediction (ch. 11)~~ -- `allgo.astro.passes`, with
      naked-eye visibility from `allgo.astro.visibility`.
  13. ~~The circular restricted three-body problem (ch. 2)~~ --
      `allgo.astro.cr3bp`.
  14. ~~IAU 2006/2000A CIO-based reduction~~ -- `allgo.astro.cio`.
  15. ~~Move what is plain mathematics out of `allgo.astro`~~ --
      `allgo.numerics` differentiation, roots, quadrature, special and
      linear's small solves; `allgo.math` fmod and frac; `allgo.geometry`
      rotation, sphere and disk, and `vec3/angle`.
  16. ~~Gauss-Jackson integration (ch. 8)~~ -- `allgo.numerics.gauss-jackson`.
  17. ~~Element sets by differential correction~~ -- `allgo.astro.sgp4-fit`.
  18. ~~Atmosphere models (ch. 8)~~ -- NRLMSISE-00 (`allgo.astro.msis`),
      Jacchia 1971 and 1970 (`allgo.astro.jacchia`) and Jacchia-Roberts
      (`allgo.astro.jacchia-roberts`), with the 1976 standard atmosphere
      and exponential model (`allgo.astro.us76`) and Harris-Priester
      (`allgo.astro.drag`).
  19. ~~The rest of what the book names, from primary sources~~ --
      polynomial roots (`allgo.numerics.polynomial`), Gooding's angles-only
      IOD (`iod/gooding`), the B-plane (`allgo.astro.bplane`), fixed-delta-v
      maneuvers, the sequential batch and extended and unscented Kalman
      filters (`allgo.astro.od`), Kaula's functions (`allgo.astro.kaula`),
      Earth radiation pressure (`allgo.astro.earth-radiation`), Walker
      constellations and coverage (`allgo.astro.constellation`), and a
      DSST-style semi-analytic theory (`allgo.astro.semianalytic`).
      JB2008 (`allgo.astro.jb2008`) from its published papers alone, its
      distributors' code unused; within 0.5% of Orekit's on average from
      220 to 400 km, diverging above 500 km at the winter pole in ways
      the papers do not settle (`jb2008/limits`).

## Time

- All calculations from Calendrical Calculations by Edward M. Reingold and Nachum Dershowitz
  - Book examples are in Lisp, which should port well to Clojure.
