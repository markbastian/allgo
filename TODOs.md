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
  3. Initial orbit determination (ch. 7): done except Battin's Lambert
     solver (a second, independent one to check the universal solver
     against).
  4. Coordinates and frames (ch. 3-4): done -- `allgo.astro.reduction` and
     `allgo.astro.states` -- except the IAU 2006/2000 CIO-based reduction.
  5. ~~SGP4/SDP4 and TLEs~~ -- `allgo.astro.sgp4`.
  6. ~~Celestial (ch. 5)~~ -- `allgo.astro.visibility` for sight, shadow,
     eclipses and naked-eye visibility, and twilight in `allgo.astro.rise`;
     Sun and Moon positions and rise and set were already there.
  7. Perturbations (ch. 8-9): J2 secular and periodic rates, gravity
     field algorithms (Pines, Gottlieb, Lear), drag and SRP analytic.
  8. Mission geometry (ch. 11): repeat ground tracks, sun-synchronous and
     frozen orbits, field of view, range and azimuth between sites.
  9. Covariance transformations between Cartesian, classical, equinoctial,
     flight and RSW/NTW.

## Time

- All calculations from Calendrical Calculations by Edward M. Reingold and Nachum Dershowitz
  - Book examples are in Lisp, which should port well to Clojure.
