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
