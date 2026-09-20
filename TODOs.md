# Done

Removed from the list below, with where each one landed:

- **Perlin noise and fBm** — `allgo.procedural.noise` (value, gradient, 4D
  gradient, sparse convolution, domain warps) and `allgo.procedural.fractal`
  (fBm, turbulence, and the three multifractals).
- **Midpoint displacement** — `allgo.procedural.terrain`, as diamond-square.
- **Voronoi diagrams** — `allgo.geometry.voronoi3d` for the 3D cells, and
  `allgo.procedural.cellular` for Worley's F1 / F2 / F2-F1 bases and the
  distance metrics. The primitives for tectonic plates and cellular biomes
  are there; nothing yet *uses* them for that.
- **Stable fluid modelling on a sphere** — `allgo.physics.sphere-fluid`,
  with `allgo.numerics.fft` and `allgo.numerics.tridiagonal` under it and
  a demo in `allgo.demo.sphere-fluid`. Solved in vorticity and
  streamfunction rather than velocity and pressure, which is what makes
  incompressibility structural instead of enforced, and what makes the
  polar singularity a property of the coordinates rather than of the
  solver. The Poisson solve is direct: FFT along longitude, one
  tridiagonal system per zonal wavenumber.
- **Hydraulic erosion** — `allgo.procedural.erosion`, as the droplet
  model rather than the pipe model: a droplet is a particle that walks
  downhill carrying sediment, and the drainage network is accumulated out
  of many independent paths rather than solved for. Runs on anything with
  a flat `:heights` array and a `:dim`, which is what
  `allgo.procedural.terrain` already hands back. Also returns `:flux`,
  the drainage accumulation, which is what to threshold to place rivers.

  Two things worth knowing before using it. `:min-slope` is scale
  dependent and the usual published default of 0.01 is far too high for
  this repository's grids — it sits above 94% of the slopes in a
  diamond-square map and turns erosion into a blur; the default here is
  0.0005. And erosion is tapered to zero at the map border, which is not
  cosmetic: without it the cells beside the edge are mined by every
  droplet that leaves and refilled by none, and the map grows a trench
  that deepens without bound.

  *Thermal erosion is not implemented* — Musgrave, Kolb and Mace pair
  the two, and talus slippage is what would put screes under the cliffs
  this leaves bare.
- **Simplex noise** — `allgo.procedural.noise/simplex-basis` and
  `simplex-basis-4d`, sharing the permutation table, hashing and gradient
  sets with the Perlin bases beside them. Four corners a sample instead of
  eight, five instead of sixteen in 4D.

  The corner counts oversell the speed: measured, it is about a sixth
  faster in 3D and a fifth in 4D, not half and not three times. Finding
  which simplex you are in is work gradient noise does not do, and in 4D
  both bases are held back by the same thing — `hash4` and `grad4` take
  five arguments, and Clojure's primitive interfaces stop at four, so
  every corner boxes. *Measured, that explanation is wrong*: a
  four-argument hasher with the dot product inlined buys 12%, not a
  factor — the JIT is already removing most of those boxes. Not worth
  the arithmetic it would inline into the two 4D bases, and the var's
  docstring records the numbers so nobody tries it twice.

  The real reason to reach for it is isotropy. Perlin noise fades along
  each axis and carries a faint squareness; simplex sums radial bumps and
  has no axis to align to. Selectable in the planet demo, where it reads
  as coastlines that stop preferring the compass points.

  Two things worth remembering. It is louder — standard deviation near
  0.39 against gradient noise's 0.27 — so a multifractal tuned against
  Perlin is not tuned against this. And there is no `:period`: tiling
  comes from folding the integer lattice, and the simplex tiling is that
  lattice sheared, so an axis-aligned repeat is not a symmetry of it.

  On the patent that prompted the note here: it covered simplex in three
  dimensions and up, and it has expired, so this is the classic
  construction rather than OpenSimplex — which existed to route around
  it. Worth confirming independently if it ever matters commercially.
- **Island generation** — `allgo.geometry.dual-mesh` and
  `allgo.procedural.island`, following redblobgames' mapgen2, with a demo
  in `allgo.demo.island`.

  Split in two on purpose. The dual mesh is the general structure — the
  Delaunay triangulation and its Voronoi dual as one graph, with centers,
  corners and edges that reference each other by id — and it is what the
  Voronoi entry above meant by primitives nothing was using yet. The
  island namespace is the passes on top, and every one of them is a graph
  traversal rather than a formula: ocean by flood fill inwards from the
  border, elevation breadth-first from the coast, rivers by following a
  downhill pointer, moisture breadth-first from fresh water. Only the
  biome step is a lookup, on Whittaker's diagram.

  Elevation and water live on the corners, biomes on the cells. That is
  the part worth remembering: water runs along the boundaries between
  regions rather than through the middles, so a river is a border.

  Watersheds and noisy edges are in too. Watersheds follow the downhill
  pointer to the sea and name the basin by where it arrives, so the
  ridges between basins fall out of the disagreements. Noisy edges draw
  every boundary as a wandering path confined to the quad of the edge's
  two corners and the two cell centres either side — stored on the edge,
  so the cells sharing it still tile exactly.

  **Cost**: `delaunay/triangulate` is Bowyer-Watson against a linear scan,
  so the whole thing is quadratic. Rewriting its inner step as one pass
  over transients roughly halved it, and it is still the ceiling — five
  hundred points and one relaxation round is about a second and a half in
  a browser. *Giving the triangulation a spatial index is what would lift
  it*, and would speed up the Delaunay, TIN and boids-voronoi demos with
  it.
- **Settlements, roads and territories** — `allgo.procedural.settlement`,
  the human layer on top of a finished island, and the part of Azgaar's
  generator that is still geography.

  Three passes, each the same kind of traversal the physical map was
  built from. Towns are scored for fresh water, low flat ground, a coast
  and what the land grows, then taken best-first with a refusal rule —
  without which every town lands in the same river valley, because the
  second best site in it beats anywhere else. Roads are a minimum
  spanning tree over the towns with each link routed by A*, so a road
  bends round a mountain instead of going over it. Territories are one
  multi-source Dijkstra outward from every town at once.

  The territory pass is the reason to do this on a graph. Assigned by
  *distance* it would be a Voronoi diagram of the towns and would know
  nothing about the land; assigned by *travel cost*, a ridge pushes the
  border away from itself because crossing it is dear from both sides.
  Frontiers land on watersheds and mountain chains without anything ever
  looking for one — measured, border cells average an elevation of 0.36
  against 0.21 for interiors, and there is a test pinning it.

  Names and cultures are in too — `allgo.procedural.naming` invents a
  language and then draws names out of it, rather than drawing syllables
  out of one shared bag. A language is a small phonology picked from a
  much larger pool, so the family resemblance within one is a consequence
  rather than an effect: if a language never drew `k`, none of its towns
  have a `k` in them. Cultures are coarser than realms on purpose — a
  realm holds one town, so a language per realm would mean no two places
  on the map ever share one.

  **Not implemented**: states as distinct from realms, religions,
  population, history. See below.
- **A spatial index for the triangulation** —
  `allgo.geometry.delaunay` indexes triangles by the bounding box of
  their circumcircle and looks only in the cell the new point falls in.
  Exact rather than approximate: a circumcircle containing the point has
  a box containing it, so every triangle the old linear scan would have
  found is in that cell. No adjacency, no point location and no argument
  about the cavity being connected — which is what the textbook fix needs
  and where its bugs live.

  Two things had to be measured rather than reasoned about. The threshold
  for holding an oversized triangle aside is bad at both ends: too low
  and ordinary triangles land in a list scanned on every insertion, too
  high and the huge circumcircles near the super-triangle get filed into
  thousands of cells each. And the first version used nested persistent
  maps, which was a large win on the JVM and made a 500-point
  triangulation *slower* in a browser than the scan it replaced. A flat
  array with lazy deletion fixed that.

  Browser, `triangulate` alone: 500 points 1193ms → 26ms, 4000 points
  13.4s → 244ms. There is a test comparing the indexed result against the
  old linear scan, because an index bug would still produce a plausible
  triangulation — just not the Delaunay one.
- **Diffusion-limited aggregation** — `allgo.procedural.dla`, with a
  demo in `allgo.demo.dla`. Release a particle far from a seed, let it
  wander, stick it where it first touches the cluster. The shape that
  grows is dendritic, and nothing in the rule mentions branching: the
  cluster shadows itself, so a walker is far likelier to meet a tip
  sticking out than to find its way into a gap between two, and tips
  outrun gaps.

  The cluster is a *tree*, not a set of cells — every particle remembers
  what it stuck to — and that is what makes it terrain rather than a
  texture. A node carrying half the cluster behind it is a trunk and
  should be high; a node with nothing behind it is a twig. Height is read
  off the subtree size, so the ridge line follows the branching structure
  rather than the geometry, then blurred and refined to give the skeleton
  flanks.

  This is the one generator here whose branching has a *direction*, out
  from the trunk to the tips, which is what `allgo.procedural.fractal`
  cannot produce at any setting.

  Two bounds make it finish rather than wander forever: walkers are born
  on a circle just outside the cluster's reach, and abandoned if they
  stray well past it — by symmetry a walker that far out is as likely to
  return from anywhere else, so restarting loses nothing. And, as with
  the triangulation, flat arrays rather than a set of `[x y]` vectors:
  hashing those to answer "is this cell taken?" a few million times cost
  six times the rest of the algorithm.

# Map Generation

Most of the geography layer is done — see **Settlements, roads and
territories** above. What is left in
[Azgaar](https://azgaar.github.io/Fantasy-Map-Generator/) is the part
that makes it a *fantasy* map rather than a map:

## States, religions and history

Cultures and names are done. What is left in Azgaar is states as a layer
distinct from the realms here (with capitals, diplomacy and borders that
move), religions, population and a potted history run over all of it.

None of that is geometry. It is simulation on top of a map that already
exists, which is a different project from the rest of this repository —
worth saying out loud rather than leaving as an open item that looks
like the others.

Rivers have bodies now. `allgo.procedural.island` traces each one from
its mouth up the fuller branch at every fork, which is the rule that
decides which stream is the same river as the one below, and gives it a
length to be labelled along — the demo does. Catchment is a separate
question and already answered separately, by the watershed pass: the
stem is the river, the watershed is everything draining into it.

## Physics Engines

Sequential impulse, TGS and XPBD are in `allgo.physics.solver`, over
contact manifolds from `allgo.physics.contact`, with the brick demo in
`allgo.demo.bricks`. All three run the same scene so they can be
compared; the differences are in the docstrings and are real — TGS holds
a stack at three iterations that sequential impulse lets sag, because it
moves the bodies between iterations and re-measures.

There are four scenes now, because one was answering one question. A
running-bond **wall** asks how much weight a solver holds up; a
**column** of single bricks asks how straight it can stand something
that has no business standing at all; a **Jenga** tower asks how quiet
it is once nothing is happening; a **keep** inside a semicircular
rampart asks all of that again about rings of boxes, where no two faces
in contact are parallel.

They do not agree, which is the useful part. TGS settles the Jenga
tower in two seconds and XPBD never does; on the keep it is the other
way round. Sequential impulse is last in every one of them. Anything
claiming a solver is simply better than another should have to explain
both of those.

**It is not yet fast, and the reason is not the solvers.** Measured in a
browser on a settled 9x8 wall (72 bricks, 581 contacts), a step is about
142ms: 78ms of collision detection and 64ms of everything else. Turning
the iteration count *down* makes it slower, because a wall that is not
held up spreads out and touches more — which is the clearest evidence
that the solve is not the bottleneck. About twenty bricks runs at fifty
frames a second today.

What would fix it, in order of expected return:

1. **`allgo.physics.contact` on primitive doubles.** The separating axis
   test asks fifteen questions of every touching pair, each a handful of
   dot and cross products, and in JavaScript every one allocates a
   three-element vector. The solve was rewritten this way and went from
   139 microseconds a contact to a fraction of it; the collision
   detection has not been.
2. ~~**A broad phase that is not every pair against every other.**~~
   Done: `allgo.physics.contact/all` sweeps and prunes through
   `allgo.spatial.sweep`, keeping the sorted order in the world between
   steps, which is what makes it cost about `n` rather than `n log n` on
   a scene that barely moves. Detection alone, a settled wall:

       bodies               257    577    901
       all pairs (ms)      8.04  24.98  48.97
       sweep and prune     6.34  15.16  24.61

   The gap widens with the body count, which is the whole argument. It is
   not larger at small sizes because there was already an axis-aligned box
   rejection in front of the exact test -- what sweep removes is the `n^2`
   box comparisons, not the exact tests, and those were already few. Which
   means the narrow phase is now nearly all of it, and item 1 is the
   remaining lever.
3. **Leaner contact preparation.** `prepare` builds a dozen typed arrays
   per step out of lazy sequences, and computes the tangent basis and the
   body-local anchors through persistent vectors.

### Standing still

A wall or a column left alone used to come down by itself, and the three
things that made it are dealt with. What is left is written after them.

**Contacts are soft now, and that is what stopped the energy.** A rigid
contact solved with a Baumgarte bias separates overlapping bodies at
`bias/dt` times the overlap and leaves that velocity in them; gravity
puts the overlap back, and the pair does it again. A stack too deep to
converge in the iterations it has gets pushed hard enough to leave the
ground: a twenty brick column threw its own bricks into the air, lost
the contacts and the impulses remembered against them, fell back deeper,
and pushed harder still — apart in three seconds. `allgo.physics.solver`
now solves the normal as a spring of a stated frequency and damping
(Catto's soft constraints; Box2D v3 and PhysX TGS Soft are both this),
which bounds the pushout by construction and gives a little of the
accumulated impulse back each iteration, so the constraint cannot store
energy. A six course column is now still to every digit — `KE` exactly
zero, positions unchanged to four decimals, for sixty seconds — where it
used to wander a centimetre and jitter at 2 cm/s forever.
`stays-put-test` is that, and it fails without this.

**Contacts have identity now.** Warm starting matched last step's
impulse by rounding the contact point to two centimetres, which is a
good key for a stack already still and a bad one for a stack that is
moving: the points slide across the rounding and the impulse history is
lost exactly when it is most needed. `allgo.physics.contact` stamps each
point with the feature that made it — which corner of which face is
pressed into which face — and the solver matches on that. Two
consequences: the match is exact, and the tangential impulses can come
back too, which they never did. Friction used to start every step from
nothing and be rediscovered in the iterations it had left, which is what
let a stack creep sideways while it stood.

**And the manifold is stable.** Two identical bricks stacked square have
all four corners of the incident face sitting exactly on the reference
face's clip planes. Clipped strictly each of those is a cut, not a
corner, the answer comes back with six points instead of four, and which
six depends on the last bit of the arithmetic — so it was a different
manifold every step and nothing could be matched to anything. Corners
within `clip-eps` of a plane count as inside it now.

Measured on the running-bond wall, no projectile, forty seconds:
sequential impulse settles to `mean |v|` of exactly zero and does not
move again; it used to sit at 0.011 and topple at thirty-four seconds.
TGS stops toppling too, though it keeps its own floor, below.

### What is still wrong

- **A tall column topples, and by now that is mostly honest.** The
  telescoping is gone. A stack used to sink into itself until the boxes
  were more than half overlapped, at which point the separating axis test
  picked a different axis, the normal flipped, and they passed through
  each other. Speculative contacts stopped that — the solver is told
  about a touch before it is a penetration — and small steps hold the
  rest. Courses still standing after twenty seconds, TGS:

      column of          6    10    14    20
      before             6     3     2     2
      after              6    10     3     8

  What remains is a real toppling mode rather than a numerical one. A
  nine metre column half a metre thick is an unstable equilibrium, and
  the lean grows as `e^(0.4 t)` — *slower* than the `e^(1.05 t)` an
  inverted pendulum that tall would manage, so friction is doing its work
  and the solver is not adding to it. What nothing here does is damp the
  perturbation to nothing, and nothing will: the seed is the solve's own
  asymmetry, a millimetre a second or so, and an unstable equilibrium
  amplifies whatever it is given.

  For XPBD specifically, what would raise the ceiling is a penetration
  re-measured against the geometry each substep rather than `depth0` plus
  the anchor drift, which is a linearisation that goes stale as the
  bricks turn.

  Sleeping is in now, and it is the answer for a column that *reaches*
  rest: a six course one is asleep two seconds in and then does not move
  at all, ever, because it is not being integrated. It is not the answer
  for one that never gets there. Twenty courses are still shedding bricks
  at half a second, well inside the half second of stillness sleeping
  asks for, so the clock never starts on the part that matters. What
  would help is the settling itself being quicker — the compression wave
  from twenty bricks landing on each other takes about two seconds to die
  out, and every course above the one that is still moving has to wait
  for it.

- **Sequential impulse does not substep, so it got none of this.** The
  small-steps rework is `step-tgs` alone. Sequential impulse keeps its
  one linearisation per step deliberately — that is what it *is*, and
  comparing the three is the point of having three — so its column limit
  is still about six courses. If the demo wants a default that stacks,
  the default should be TGS, which it already is.

- **TGS's velocity floor is gone, and the relax pass is why.** It used to
  leave about one substep of gravity in every body — 0.05 m/s, downward,
  more of it the higher up the wall, which is the shape of a Gauss-Seidel
  correction that has not reached the top. What fixed it was not more
  sweeps but the right ones: warm starting every substep instead of only
  the first, and a relax pass after the positions move that solves for
  nothing but *stop closing*, taking back the velocity the push put in.
  That pass is the cheap standard alternative to a pseudo velocity, and
  unlike a pseudo velocity it goes through the friction solve. The cost
  is real — a settled 5x4 wall went from 1.5ms a step to 2.3ms.

  For the record, since establishing it cost a day: the floor was never
  Baumgarte. It was flat against the old `bias` — 0.2, 0.1 and 0.05 all
  gave 0.049 — and it tracked the sweep budget exactly:

      sweeps per substep   8     4     2     1
      floor (m/s)          0.011 0.035 0.049 0.058

- **Split impulse was tried and does not pay.** Written and measured:
  the bias solved into a pseudo velocity of its own, cleared every
  substep, added to the real velocity to move the bodies and then
  dropped. It does what it says and makes things *worse*, because the
  separating motion no longer passes through the friction solve, so the
  wall slides apart instead — sequential impulse went from toppling at
  thirty-five seconds to nineteen. Solving friction against the pseudo
  velocity too recovers most of that and still loses to leaving it alone.
  Soft contacts were the right answer to the same problem.

- **The three solvers part company at about eight courses, and that is
  documented rather than outstanding.** All of them hold a four-wide
  running bond up to eight; at twelve only TGS does, and at sixteen it
  holds 63 of 64. The table and the reasoning are in
  `allgo.physics.solver`'s docstring and in the demo's.

  Sequential impulse is not fixable here and should not be fixed: it
  linearises once per step, which is what it *is*, and iterating harder
  does nothing — sixteen iterations and sixty-four leave the same
  nineteen bricks standing. TGS exists because of this.

  XPBD used to *explode* past its limit rather than fall over, reaching
  sixty metres a second on a sixteen course wall from a standing start.
  That was a missing bound, not a depth limit: it reads velocity back off
  the position correction over the substep, so an unbounded correction is
  an unbounded velocity. It is now bounded by `max-push-speed`, the same
  number the impulse solvers have always used, and `pushout-bound-test`
  holds it there. It still cannot stack past eight courses, and the
  reason it cannot is the one below.

### Tried and not kept

- **Manifold reduction to four points.** Standard practice — Box2D and
  Bullet both cap — and it works as advertised: a settled 9x8 wall went
  from 669 contact points to 556, with every pair at four, and the four
  chosen are the right four (the corners of the overlap region, picked by
  spread from a stable anchor; choosing by *depth* instead is the obvious
  thing and is wrong, because across a resting face the depths differ by
  microns while the points differ by centimetres, so noise picks a
  different four every step and warm starting loses its history).

  It still costs stack depth, and the reason is worth knowing: with soft
  contacts the point count is doing double duty. Six points on a pair is
  six constraint rows, which is half again as much solving for that pair
  per sweep as four. Removing the two surplus rows at the same iteration
  count is a real loss of convergence, and a fourteen course column went
  from three bricks standing to two.

  It can be paid for with stiffness — `contact-hertz` 45 instead of 30
  makes the fourteen course column stand *fully*, with or without the
  reduction. But that stiffness is not free either: it breaks
  `stays-put-test` and `restitution-test`, which is to say it buys stack
  depth by spending the stillness. Sitting on the clamp (60, a quarter of
  the substep rate) is worse still — marginal by construction, and it
  showed, with a settled column that would not stop and a suite three
  times slower for the contacts a shivering wall makes.

  So: reduction is correct, cheap, and a wash at the stiffness that keeps
  a stack still. Revisit it together with per-pair stiffness rather than
  on its own.

- **Ordering contacts for Gauss-Seidel.** The claim is that solving a
  stack bottom-up lets one sweep carry the floor's support to the top.
  Measured: forcing bottom-up is *identical* to what we already do, and
  top-down is slightly worse. Sweep and prune already sorts along the
  axis the scene is most spread along, so a column is visited bottom-up
  for free. Nothing there to win.

### Other techniques still missing

- **Featherstone / articulated bodies** — started.
  `allgo.physics.articulated` has the spatial algebra, the recursive
  Newton-Euler algorithm for inverse dynamics and the articulated body
  algorithm for forward dynamics, over a tree of one degree of freedom
  joints on a fixed base.

  Why it exists, measured rather than asserted: two metre-long links off
  a fixed base, the lower heavier than the upper, worst separation
  between a joint's two anchors — anchors that are meant to be the same
  point — over two seconds of `allgo.physics.joint`:

      mass ratio   5 substeps   10      20
               1      0.00475   0.00126   0.00029
              10      0.01633   0.00393   0.00098
             100      0.11241   0.04149   0.01001
            1000      2.75497   0.29095   0.09466

  Substeps help and do not cure; every column climbs with the ratio.
  Reduced coordinates have no such number, because a joint angle has
  nowhere to put a violation.

  It agrees with the closed form for a hinged rod to twelve digits, and
  the fast path agrees with `M^-1 (tau - bias)` built out of the slow one
  to about fourteen. What is left:

  - **Contacts.** The reason the constraint formulation is still the
    right shape for a pile of bricks is that a contact and a joint go
    through the same solver. A reduced-coordinate ragdoll that cannot be
    hit is not much of a ragdoll. This is the big one.
  - **A floating base**, which is a six degree of freedom joint to the
    world. Without it the chain is bolted to the origin.
  - **Multi-degree-of-freedom joints.** A shoulder is a ball joint;
    today it has to be spelled as three hinges with massless links
    between them.
  - **A ragdoll demo**, which is what all of the above is for.
- **Continuous collision detection.** Half done. Speculative contacts
  — the cheap version — are in: the margin out to which a gap still
  counts as a contact is now the fixed 2cm plus how far the pair can
  travel in the step, and the broad phase sweeps its boxes over the
  step so the pair survives to be asked about at all. A ball fired at a
  static slab, highest speed at which it is still stopped:

      sequential impulse    26 m/s → over 2000
      tgs                   26 m/s → over 2000
      xpbd                  26 m/s → 120

  The demo's own speed slider goes to 60, so this was reachable by
  hand rather than theoretical. It costs nothing measurable: a settled
  wall is untouched because a sleeping body contributes no travel, and
  a wall being hit costs the same at 10 m/s as at 120. Firing a ball
  past a brick rather than at it produced no deflection at any speed
  or offset tried, because a pair whose swept boxes do not overlap is
  never handed to the narrow phase however wide the margin is.

  **Conservative advancement is what is left, and XPBD is why.** It
  integrates the substep and only then pushes overlaps apart, so a gap
  it has not yet reached buys it nothing and it steps over a thin slab
  above about 120 m/s. The impulse solvers do not need it. Finding the
  time of impact and advancing to it would fix XPBD and would also
  replace `stopped at the surface some time within the step` with
  `stopped at the instant of arrival`, which is what a bullet wants.
- ~~**Sleeping.**~~ Done. An island whose bodies have all been under the
  speed thresholds for half a second is frozen: velocities zeroed, not
  integrated, not solved, and no contact generated between two bodies
  that are both inert. It wakes when something awake is found touching
  it, in the same step rather than the next, because a step late is a
  ball halfway through a wall.

      settled wall      5x4     9x8
      before (ms)      2.198   8.968
      asleep (ms)      0.198   0.599

  The scene that never settles pays about five percent for asking --
  a 16x16 wall, which topples rather than resting, went from 32.2ms to
  34.0ms. Most of that was bought back with an early-out: an island
  sleeps on the clock of whichever member has been still least long, so
  if no body anywhere has reached the threshold then no island can, and
  the union-find is work with a known answer.
