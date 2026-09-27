# Notes

How the finished parts of this repository work, what was measured along
the way, and what was tried and not kept. Open work is in
[TODOs.md](TODOs.md); how to build and run is in the [README](README.md).

The namespace docstrings are the reference. These notes are the
reasoning that does not fit in one docstring: decisions that cut across
namespaces, and the numbers behind them.

# Procedural generation

## Texturing & Modeling

Seven namespaces under `allgo.procedural` follow Ebert, Musgrave, Peachey,
Perlin and Worley, *Texturing & Modeling: A Procedural Approach*, and are
arranged the way chapter 20 argues they should be -- a small set of parts
that compose freely, rather than a set of finished effects:

| Namespace                     | Chapter | What it is                                                       |
|-------------------------------|---------|------------------------------------------------------------------|
| `allgo.procedural.shaping`    | 2       | `step`, `smoothstep`, `bias`, `gain`, `spline` -- the vocabulary |
| `allgo.procedural.noise`      | 2, 12   | Value, gradient and 4D gradient noise, sparse convolution, warps |
| `allgo.procedural.cellular`   | 4       | Worley's F1, F2, F2-F1 and the metrics that reshape the cells    |
| `allgo.procedural.fractal`    | 14, 16  | fBm, turbulence, and the three multifractals                     |
| `allgo.procedural.qaeb`       | 17      | Error-bounded ray marching of a surface that is only a function  |
| `allgo.procedural.atmosphere` | 18      | Exponential density, Rayleigh and Mie scattering, aerial haze    |
| `allgo.procedural.planet`     | 20      | The lot, assembled into a world                                  |

The contract between them is one line: a *basis* is a function of `[x y z]`
to a number near [-1, 1], a *fractal* takes a basis and returns a basis,
and so a fractal can be the basis of another fractal. Everything else
follows from that.

    (planet/planet
     {:sea-level 0.0
      :terrain (-> (noise/gradient-basis {:seed 7})
                   (noise/distorted (noise/vector-basis {:seed 8}) 0.4)
                   (fractal/ridged-multifractal {:octaves 9 :H 0.9}))})

Nothing is stored: the coastline, the color of the ground and the cover of
the clouds are all evaluated from the 3D direction of the point being
asked about, which is what keeps a planet free of the seam and the polar
smear that a latitude-longitude map cannot avoid.

The other half of making that look right is knowing when to stop.
`fractal/octaves-for` says how many octaves a given sampling rate can
actually carry; past that, octaves do not arrive as detail but as speckle
that moves when the camera does, and they are paid for at full price.
Passing `:sample-spacing` to a planet cuts every layer to what its mesh
can hold -- which is both cleaner and, at a typical resolution, about
three times faster. The detail that no longer fits in the geometry goes
where it can be seen instead: `noise/gradient-basis` takes a `:period`,
so a tile of fine relief can be baked and read per pixel without a seam.

## Perlin noise and fBm

`allgo.procedural.noise` (value, gradient, 4D
gradient, sparse convolution, domain warps) and `allgo.procedural.fractal`
(fBm, turbulence, and the three multifractals).

## Midpoint displacement

`allgo.procedural.terrain`, as diamond-square.

## Voronoi diagrams

`allgo.geometry.voronoi3d` for the 3D cells, and
`allgo.procedural.cellular` for Worley's F1 / F2 / F2-F1 bases and the
distance metrics. `allgo.geometry.dual-mesh`, under island generation
below, is what puts the 2D diagram to work.

## Simplex noise

`allgo.procedural.noise/simplex-basis` and
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

## Hydraulic erosion

`allgo.procedural.erosion`, as the droplet
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

## Diffusion-limited aggregation

`allgo.procedural.dla`, with a
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

## Island generation

`allgo.geometry.dual-mesh` and
`allgo.procedural.island`, following redblobgames' mapgen2, with a demo
in `allgo.demo.island`.

Split in two on purpose. The dual mesh is the general structure — the
Delaunay triangulation and its Voronoi dual as one graph, with centers,
corners and edges that reference each other by id — and it is the
2D Voronoi structure the rest of the map is built on. The
island namespace is the passes on top, and every one of them is a graph
traversal rather than a formula: ocean by flood fill inward from the
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
two corners and the two cell centers either side — stored on the edge,
so the cells sharing it still tile exactly.

**Cost**: `delaunay/triangulate` is Bowyer-Watson against a linear scan,
so the whole thing is quadratic. Rewriting its inner step as one pass
over transients roughly halved it, and it is still the ceiling — five
hundred points and one relaxation round is about a second and a half in
a browser. That ceiling is gone now: see the spatial index for the
triangulation, below.

Rivers have bodies now. `allgo.procedural.island` traces each one from
its mouth up the fuller branch at every fork, which is the rule that
decides which stream is the same river as the one below, and gives it a
length to be labeled along — the demo does. Catchment is a separate
question and already answered separately, by the watershed pass: the
stem is the river, the watershed is everything draining into it.

## Settlements, roads and territories

`allgo.procedural.settlement`,
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
population, history. See [TODOs.md](TODOs.md).

## A spatial index for the triangulation

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

# Flow on a sphere

`allgo.physics.sphere-fluid` is the same physics as the flat
`allgo.physics.fluid` and a deliberately different answer to it. The flat
one tracks velocity and cancels the divergence with a pressure solve; this
one never represents divergence at all.

A sphere is closed and has genus zero, so the Hodge decomposition leaves
no harmonic part and *every* divergence-free field on it is the rotated
gradient of a streamfunction, with nothing left over. The whole state is
then one scalar -- the vorticity -- and incompressibility stops being a
constraint to enforce and becomes a property of the representation. There
is no pressure, no projection and no residual; the discrete velocity is
divergence-free to the last bit, poles included, and there is a test that
says so.

The classical pole problem largely evaporates with it. The singularity in
latitude and longitude is a singularity of the *frame*, not of the
physics, and a solver that advects only scalars -- by rotating points
along great circles in 3D, where no frame appears -- never meets it. What
the axis buys in exchange is a fast direct Poisson solve: the operator is
separable and longitude is periodic, so an FFT along it leaves one
tridiagonal system per zonal wavenumber. Transform, solve, transform
back, with no iteration and no tolerance.

The same machinery takes three operators, and the third is what lets the
planet spin. `laplacian` is Poisson, for the streamfunction;
`laplacian - alpha` is Helmholtz, for implicit diffusion; and
`laplacian + c d/dlam` is the Coriolis term solved implicitly rather than
stepped. That last one matters more than it looks. Rotation makes the
vorticity-streamfunction loop oscillate -- those are Rossby waves, and the
fastest are the largest-scale ones, with frequency near the rotation rate
itself. Stepping them from the start of the interval is forward Euler on
an oscillation, which grows at every step size there is: energy came out
ten thousand times too large at a rotation of 80, however small the step.
Solving the term implicitly removes the restriction and costs nothing,
because `d/dlam` is diagonal in the same Fourier basis -- the tridiagonal
system merely becomes complex.

The trapezoidal rule, specifically, and not the obvious implicit one:
backward Euler is stable *because* it damps, and what it would damp here
are precisely the waves that organize the flow into bands.

Bands do not emerge from scattered vortices in any reasonable time -- on a
real planet they are held up by forcing from below over geological time --
so `zonal-jets!` puts them in, and the simulation shows what the rotation
does to them. They shred without it and hold with it. The vortex chains
along each band appear only while the jets are still barotropically
unstable, which is `rotation < amplitude * jets / 2`; past that the bands
go glassy and the whorls disappear.

# Physics

## Three contact solvers

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

## Where a step goes

**It is not yet fast, and the reason is not the solvers.** Turning the
iteration count *down* makes it slower, because a wall that is not held
up spreads out and touches more — which is the clearest evidence that
the solve is not the bottleneck.

On a settled 16x16 wall (257 bodies, ~2,400 contacts) on the JVM, a step
and where it goes:

                          sequential impulse    tgs    xpbd
    total (ms)                         14.7   19.6    32.9
    collision detection                 6.8    6.8     7.5
    the constraint solve                2.2    3.0    16.4
    everything else                     5.7    9.8     9.0

The solve is a sixth of a sequential impulse step and a seventh of a TGS
one. What the other five sixths were, until recently, was the cost of
carrying bodies and contacts through persistent maps and vectors: the
XPBD position pass rebuilt a 257-element vector of body maps twice per
contact, `refresh-anchors!` read four keyword lookups and allocated
three vectors per contact twice a substep, and `prepare` built a dozen
typed arrays out of lazy sequences of boxed doubles. Moving those onto
`body-arrays` took the three solvers from 19.9 / 32.9 / 91.0 ms.

The ranking of what is left is in [TODOs.md](TODOs.md). Collision
detection is at the top of it for the first time honestly rather than by
default.


### Sweep and prune

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
 means the narrow phase is now nearly all of it, and item 2 is the
 remaining lever.

### Leaner contact preparation

`prepare` is one pass
 writing the arrays directly, and takes the poses off `body-arrays`
 rather than the body maps; the world lever arm and the body-local
 anchor now come out of one subtraction instead of four vector
 operations and two `nth`s. 5.4ms to 2.3ms.

## Standing still

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
used to wander a centimeter and jitter at 2 cm/s forever.
`stays-put-test` is that, and it fails without this.

**Contacts have identity now.** Warm starting matched last step's
impulse by rounding the contact point to two centimeters, which is a
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
TGS stops toppling too; its velocity floor is below.

### Tall columns

The telescoping is gone. A stack used to sink into itself until the boxes
were more than half overlapped, at which point the separating axis test
picked a different axis, the normal flipped, and they passed through
each other. Speculative contacts stopped that — the solver is told
about a touch before it is a penetration — and small steps hold the
rest. Courses still standing after twenty seconds, TGS:

    column of          6    10    14    20
    before             6     3     2     2
    after              6    10     3     8

What remains is a real toppling mode rather than a numerical one. A
nine meter column half a meter thick is an unstable equilibrium, and
the lean grows as `e^(0.4 t)` — *slower* than the `e^(1.05 t)` an
inverted pendulum that tall would manage, so friction is doing its work
and the solver is not adding to it. What nothing here does is damp the
perturbation to nothing, and nothing will: the seed is the solve's own
asymmetry, a millimeter a second or so, and an unstable equilibrium
amplifies whatever it is given.

For XPBD specifically, what would raise the ceiling is a penetration
re-measured against the geometry each substep rather than `depth0` plus
the anchor drift, which is a linearization that goes stale as the
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

### Sequential impulse does not substep

The small-steps rework is `step-tgs` alone. Sequential impulse keeps its
one linearization per step deliberately — that is what it *is*, and
comparing the three is the point of having three — so its column limit
is still about six courses. If the demo wants a default that stacks,
the default should be TGS, which it already is.

### TGS's velocity floor

The relax pass removed it. It used to leave about one substep of gravity
in every body — 0.05 m/s, downward, more of it the higher up the wall,
which is the shape of a Gauss-Seidel correction that has not reached the
top. What fixed it was not more sweeps but the right ones: warm starting
every substep instead of only the first, and a relax pass after the
positions move that solves for nothing but *stop closing*, taking back
the velocity the push put in. That pass is the cheap standard
alternative to a pseudo velocity, and unlike a pseudo velocity it goes
through the friction solve. The cost is real — a settled 5x4 wall went
from 1.5ms a step to 2.3ms.

For the record, since establishing it cost a day: the floor was never
Baumgarte. It was flat against the old `bias` — 0.2, 0.1 and 0.05 all
gave 0.049 — and it tracked the sweep budget exactly:

    sweeps per substep   8     4     2     1
    floor (m/s)          0.011 0.035 0.049 0.058

### The three solvers part company at about eight courses

All of them hold a four-wide running bond up to eight; at twelve only
TGS does, and at sixteen it holds 63 of 64. The table and the reasoning
are in `allgo.physics.solver`'s docstring and in the demo's.

Sequential impulse is not fixable here and should not be fixed: it
linearizes once per step, which is what it *is*, and iterating harder
does nothing — sixteen iterations and sixty-four leave the same
nineteen bricks standing. TGS exists because of this.

XPBD used to *explode* past its limit rather than fall over, reaching
sixty meters a second on a sixteen course wall from a standing start.
That was a missing bound, not a depth limit: it reads velocity back off
the position correction over the substep, so an unbounded correction is
an unbounded velocity. It is now bounded by `max-push-speed`, the same
number the impulse solvers have always used, and `pushout-bound-test`
holds it there. It still cannot stack past eight courses; why is the
tall-column entry in [TODOs.md](TODOs.md).

## Tried and not kept

- **Manifold reduction to four points.** Standard practice — Box2D and
  Bullet both cap — and it works as advertised: a settled 9x8 wall went
  from 669 contact points to 556, with every pair at four, and the four
  chosen are the right four (the corners of the overlap region, picked by
  spread from a stable anchor; choosing by *depth* instead is the obvious
  thing and is wrong, because across a resting face the depths differ by
  microns while the points differ by centimeters, so noise picks a
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

- **Split impulse was tried and does not pay.** Written and measured:
  the bias solved into a pseudo velocity of its own, cleared every
  substep, added to the real velocity to move the bodies and then
  dropped. It does what it says and makes things *worse*, because the
  separating motion no longer passes through the friction solve, so the
  wall slides apart instead — sequential impulse went from toppling at
  thirty-five seconds to nineteen. Solving friction against the pseudo
  velocity too recovers most of that and still loses to leaving it alone.
  Soft contacts were the right answer to the same problem.

- **Stiffening a contact by the load it carries.** The idea: a real
  column of perfect blocks stands at forty courses, and this one buckles
  at about eleven because its soft contacts are springy joints. Treat
  each interface as a torsion spring and the column as one buckling
  under its own weight, and the critical height goes as `hertz^(2/3)`.
  Measured, that holds (demo column, TGS, 20 s):

      stiffness      30 Hz   60 Hz   90 Hz   240 Hz
      stands up to      10      16   18-20       20
      predicted         10      16      21       40

  So the plan was to scale each contact's frequency by the square root
  of the weights it holds up -- last step's normal impulse over
  `m_eff g` -- so every interface sags the same and the top stays soft
  enough to be still. It is worse than uniform stiffness, at every
  setting tried. The frequency cap of a quarter of the substep rate is
  what binds: past it a contact is a rigid one with a full position
  correction and pumps energy into the stack, which is the pogo soft
  contacts were brought in to stop. Uncapped, columns fell *sooner* than
  at the default; capped, it is no better than raising the stiffness
  everywhere to the cap.

  The top row stops following the prediction at 240 Hz, and a trace says
  why: stiff enough, a column no longer buckles but sways, the top
  swinging ±15 cm with a period of about two and a half seconds and
  almost no damping, until a swing carries it past the footprint. The
  5 mm slop is part of that -- an interface can rock inside it for
  nothing -- and cutting it to 1 mm at 240 Hz with sixteen substeps
  stands thirty courses. Forty still falls, at four times the cost of a
  default step.

  Damping the sway is what finishes it. Angular damping of 0.5 on each
  brick (linear 0.25), with the same 240 Hz, 1 mm slop and sixteen
  substeps, leaves all forty courses in place after twenty seconds. It
  does nothing for the default column, which is the diagnosis
  confirming itself: buckling is a static instability, and damping only
  slows it. So forty courses is reachable honestly, at four times the
  price of a step, and not at the default settings.

- **Shock propagation.** Guendelman, Bridson and Fedkiw (2003): order
  the contact graph from the ground up and, in a final pass, treat each
  lower body as immovable, so support reaches the top in one sweep.
  Tried as a velocity pass inside TGS -- the last solve iteration, the
  relax pass, or both. One sweep made things worse, and fast: a pair's
  four corners visited in the same order, with only the upper body free
  to answer, leave it turned a little the same way every time, and a
  twenty course column walked twenty-five centimeters in a second.
  Sweeping each level eight times before the next fixes that, and in the
  relax pass raises the default column from ten courses to sixteen.
  That is the whole gain. Combined with stiffer contacts it is worse
  than the stiff contacts alone. The paper applies it to rigid contacts
  at the position level; a velocity pass over soft ones is a weaker
  version of it, and this is a result about that version.

## Articulated bodies

`allgo.physics.articulated` has the spatial algebra, the recursive
Newton-Euler algorithm for inverse dynamics and the articulated body
algorithm for forward dynamics, over a tree of hinges, sliders and
ball joints on a base that is bolted or free, with joint limits,
contacts against static geometry and against the model's own limbs.
`allgo.demo.ragdoll` drops a figure on the floor with it.

What it cannot do is close a loop, which is what a tree means. Two
models pinned together through `allgo.physics.world` is the way round
that: the motorcycle rider is a second model held to the bike at the
seat, the grips and the pegs.

Why it exists, measured rather than asserted: two meter-long links off
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

The root can be free as well as bolted: `{:base {...} :links [...]}`
gives it a body whose pose and spatial velocity live in the state
rather than in `q`, so nothing has to parameterize a rotation with
three numbers. It costs one 6x6 solve a step and changes nothing about
the joints.

It agrees with the closed form for a hinged rod to twelve digits, the
fast path agrees with `M^-1 (tau - bias)` built out of the slow one to
about fourteen, and a free chain under no gravity conserves both its
momenta. That last is the one that earns its keep: it caught gravity
being applied in the world's frame rather than the root's, which is
invisible until the root both falls and spins. Added since:

- **Contacts,** against static geometry. A part with a
  `:shape` is handed to `allgo.physics.contact` where it is, and
  `solve-contacts` runs the same sequential impulse loop the rigid
  solvers do -- it can, because a contact solver only ever asks how
  much velocity a push buys and then pushes, and both answers come out
  of `H^-1`. A box comes to rest at half its own height and stays
  there, a slope is slid down when slippery and not when gripping, and
  a bounce returns to `0.2 + 0.8 e^2` to within two centimeters.

  Contact against things that can move back is `allgo.physics.world`,
  below.
- **Multi-degree-of-freedom joints.** `:spherical` is three
  numbers of velocity and a quaternion of configuration, for the same
  reason the base is. Three stacked hinges would be the alternative
  and they gimbal lock, which a tumbling ragdoll finds. Checked
  against a hinge: the same rod swinging about one axis gives the same
  angle to 1e-6 after two seconds, and the same effective mass to
  1e-9.
- **Speed.** Thirty times, and no longer the thing in the way. A
  chain lying on the floor, JVM, per step:

      links         0      2      4      8     16
      contacts      4      6      8     16     52
      at first   0.71   4.92  12.19  29.85  154.1
      structure  0.67   2.03   4.24   6.87   31.35
      flat       0.45   0.86   1.36   1.78    5.21

  The structural half: a contact used to cost a whole inverse inertia
  matrix and is now three O(n) impulse responses, with the
  per-configuration data built once a step. The other half: the inner
  loops run on flat `double` arrays rather than vectors of vectors,
  and the contact sweep keeps one mutable generalized velocity.

  A ragdoll is around eleven links, so about 1.3ms a step on the JVM
  settled. Worth knowing: the first flat-array attempt was *slower*,
  from three missing type hints -- reflective `aget` costs two orders
  of magnitude, and `make reflect` is what catches it.

  The recursion over that arithmetic went the same way afterward:
  `deftype`s instead of maps, so a link's transform is a field read
  rather than a hash lookup, and the total degrees of freedom read off
  the last link rather than walked for on every one of the seventy
  calls a step makes. Worth 35% in a browser and lost in the noise on
  a JVM, which is the expected shape.

  **How to measure this from a browser, because it was nearly got
  wrong twice.** An automated tab is backgrounded and Chrome gives it
  less of a processor: ten million square roots, no Clojure in them,
  take 6.3ms on this JVM and 110ms in that tab. Absolute figures from
  there mean nothing, and one was published here before anybody
  checked. Ratios do -- run a calibration loop beside the thing being
  measured, in the same tab, and divide. That is how the 35% was
  measured, and it is also what showed the ragdoll to be about twice
  the JVM's cost rather than thirty times: the ordinary price of
  Clojure in JavaScript, and not a bug. A day was spent hunting it.

- **Joint limits,** needed sooner than expected: without
  them a body settles with its head folded back on itself and its
  knees bent the wrong way. `:limit [lo hi]` for a hinge, `:cone` and
  `:twist` for a ball joint, solved as one-sided constraints in the
  same sweep as the contacts -- a limit is a push between a link and
  its own parent where a contact is a push between a link and the
  floor, and the impulse response does not need to know which. Held
  to within half a degree of the limit in every direction tested.

- **Self-collision between links.** `:self-collide?` adds
  the model's own parts against each other; a link and its parent are
  never tested, since they overlap at the joint by construction, and
  `:no-collide` names any other pair to leave alone. The settled
  ragdoll goes from thirty-two contact points' worth of limb inside
  other limb to two. An internal impulse leaves both momenta exactly
  where they were, to 1e-12, and the mass it meets is a reduced mass
  -- `1/m_pair = 1/m_i + 1/m_j` to four digits when the two limbs are
  made independent.

- **A ragdoll demo.** `allgo.demo.ragdoll`: eleven boxes, ten
  joints, twenty-eight degrees of freedom, a floor, a tilted slab to
  fall off and a wall of loose bricks to land in -- dropped on four
  courses it knocks all sixteen askew, which is `allgo.physics.world`
  doing the thing neither half of this library could do alone. On the
  JVM a step is 0.17ms in flight and 1.3ms once
  it has landed on eighteen contacts, so one substep a frame, which
  is stable. The browser figure that used to be quoted here has been
  withdrawn rather than corrected -- see the note on measuring from a
  backgrounded tab, above. The `joint limits` toggle is the one worth
  trying: off, the figure settles into a pile of sticks.

## Mixed scenes

`allgo.physics.world`: a scene of loose rigid bodies and jointed
models that collide with each other. A ragdoll can knock a brick off
a wall now, which is the first thing either half of this library has
done that the other could not.

The shape it took is the part worth keeping. It is not a connector
between the two solvers, because `allgo.physics.solver`'s sweep is
written around flat per-body arrays of inverse mass and inverse
inertia and an articulated model has no `1/m` to put in one. The
split that works is per *contact*: a contact knows what is on each
side of it, and every kind of body can answer the same three
questions -- how fast is the surface here, what does a unit impulse
buy, take this impulse. A rigid body answers from its inverse
inertia, an articulated link from a walk of its tree, a static body
answers nothing. One Gauss-Seidel sweep covers all of them together,
which is what makes a brick pushed by both a hand and the brick below
it see both pushes in the same iteration.

Measured: a limb swung at 4 rad/s into a loose brick is slowed to
0.67 and sends the brick off at 1.3 m/s; with nothing there it keeps
4.0 exactly. A model alone in a world moves identically to the same
model under `allgo.physics.articulated/step`, to 1e-12 -- the world
adds nothing when there is nothing to add. And the exchange loses
momentum only as fast as the integrator does, which is first order:
halve the step and the loss halves.

What it is not: sequential impulse only, with no substepping, no
XPBD, no sleeping, and warm starting that lasts a step rather than
crossing between them. `allgo.physics.solver` keeps all of those and
is still the thing for a scene that is only bricks. Two articulated
models touch each other through the same sweep and can be pinned
together, which is how the motorcycle carries its rider.

## A rider who stays on, and a crash that stays in one piece

The motorcycle's rider is a second articulated model held to the bike by
five pins -- seat, grips, pegs -- and two things went wrong with that,
both found by riding the three courses flat out, at three-quarter
throttle, leaned hard over and without levelling in the air.

**Riders came off upright bikes.** Five rides in twelve threw the rider
while the bike stayed up. None of them was the rider being overloaded: a
pin's force sat well under its limit and then passed it in one 1/240 s
step -- the seat went from a steady hundred newtons to over twenty
thousand. A pin is rigid, so when a wheel strikes something and the
bike's velocity changes in a step, the pin changes the rider's in the
same step, and eighty kilograms jolted by a meter a second in 1/240 s is
nineteen kilonewtons. A person's body takes that over tens of
milliseconds. So a pin can now carry a `:break-time`, and breaks on its
force averaged over that long (`allgo.physics.world`); the rider's are a
twentieth of a second. And a foot knocked off its peg now frees that
foot rather than throwing the rider -- people ride on after clipping a
rock. Only the seat or a grip letting go, or the bike going over, throws
them. After both, none of the twelve rides throws a rider off an upright
bike.

**Crashes took the solver apart.** Two separate mechanisms, with the
same look -- limbs spinning at thousands of radians a second, then NaN.

The first was the rider's muscles, while still pinned on. They are
springs and dampers applied as torques, which is explicit, and an
explicit damper overshoots once `kd dt` passes about twice the inertia
it acts on. A forearm turning about its own length has almost none; in a
hard enough lean the forearm's damper was around four times past that,
and the rider's joint rates went 13, 30, 76, 212, 560, 1457 rad/s on
successive steps with nothing touching them. The fix is armature
(`:armature` in `allgo.physics.articulated`): inertia a joint has of its
own, along its own axes, which is what makes a damper of `kd dt` of it
implicit near enough. Each rider joint gets `kd dt + kp dt^2`. It is the
joint's inertia to forward dynamics, inverse dynamics and the contacts
alike, and there is a test that the fast path and the mass matrix still
agree to 1e-9 with it on.

The second was the thrown ragdoll, which has no muscles at all. Stepped
alone -- no gravity, no contacts, no limits -- a rider thrown from a
crash gains energy it cannot have, and the gain is the integrator's:

    step          energy gained in 0.25 s
    1/240 s       +15.5%
    1/480 s       +4.8%
    1/960 s       +2.1%
    1/1920 s      +1.0%

First order in the step, which is semi-implicit Euler on light links
spinning at seventy radians a second. At 1/240 s it compounds. Halving
the step, or giving the limbs some damping, each fell short alone --
four times the steps with no damping still blew up, and damping up to
half the muscles' at the normal step did too -- but together they hold:
once thrown the rider is stepped in two halves and keeps a tenth of its
muscles' damping (`limp` in `allgo.simulation.motorcycle`). Riding pays
for neither. Every one of the twelve crashes now falls, lands and comes
to rest, including a trials run at a hundred kilometers an hour.

Tried on the way and not enough on its own: a cap on the thrown rider's
joint speed. At 40 rad/s it still blew up two seconds in, through the
bodies' own motion rather than the joints'.

## Continuous collision detection

**Speculative contacts**, the cheap half: the margin out to which a
gap still counts as a contact is the fixed 2cm plus how far the pair
can travel in the step, and the broad phase sweeps its boxes over the
step so the pair survives to be asked about at all. That fixed the
impulse solvers, which are told about a gap and refuse to close it
faster than it is wide.

**Conservative advancement**, the thorough half, in
`allgo.physics.toi`. Take the distance between two shapes and the
fastest either could be closing it; their quotient is a length of
time in which nothing can happen, so advance by it and ask again. The
sequence walks down onto the moment of first contact without ever
stepping past it. Distance comes from `allgo.geometry.gjk`, so a box
and a ball are the same problem, and rotation is bounded by
`|omega| r` rather than simulated -- which over-states how fast a
spinning body closes a gap, and is the right way to be wrong.

XPBD is why it was needed. It integrates a substep and only then
pushes overlaps apart, so a gap it has not reached buys it nothing.
Its substeps are now variable: as long as the nominal one ordinarily,
and no longer than the time of impact when something is about to be
stepped over.

    a ball fired at a static slab, highest speed still stopped

    sequential impulse    26 m/s  ->  over 20000
    tgs                   26 m/s  ->  over 20000
    xpbd                  26 m/s  ->  120  ->  10000

The catch worth writing down: advancing *onto* the moment of impact
is not enough. A position solver corrects `overlap - slop` and has
nothing to do at an overlap of zero, so the ball arrived at the
surface still traveling and the next slice carried it through. The
advance goes to the impact plus the time to reach twice the slop,
which is the shallowest overlap certainly worth solving.

It costs nothing when nothing is moving fast, which is the point of
asking first: a settled 9x8 wall under XPBD is 0.544ms a step against
TGS's 0.611, unchanged by any of this, because no body crosses half
its own thinnest dimension in a substep and the time of impact is
never asked for. XPBD's ceiling is now the cap on how finely a step
may be cut -- 64 slices of a sixtieth of a second is a quarter of a
millisecond, and 20 km/s crosses the slab in less than that.

## Sleeping

An island whose bodies have all been under the
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

A body can also *start* asleep: `rigid/asleep`. That is what games do
with towers, and in the bricks demo it is the `start asleep` toggle. A
forty course column placed that way stands under every solver, where
placed awake no solver here holds more than about eleven -- because it
is not being simulated, which is the honest description of it.

Waking has to take the whole pile. Nothing generates contacts between
two sleeping bodies, so a sleeping stack woken only where an awake body
touches it wakes a body a step from the bottom up -- and a shot that
clears the bottom of a column faster than that left the upper courses
asleep in mid air. When anything wakes now, so does every sleeping
body touching it, directly or through a chain, found by bounding box at
that moment: the island an engine would wake, without keeping the
islands of sleeping bodies between steps.

The opposite problem is a stack that starts awake and exactly touching.
Its contacts are springs that push only once compressed, with no
remembered force, so a reset opens with the whole stack squashing into
place: an eight course column sinks 11 cm under sequential impulse, 6 cm
under TGS and 4 cm under XPBD, and the first and last of those bounce.
Most of it is the 5 mm slop, once per interface; the rest is the
springs taking up the weight. `solver/settled` runs the solver out of
sight for half a second with every velocity cut to seven tenths after
each step -- dynamic relaxation -- and the scene is shown with the
squash done: a wall or column then sinks a few millimeters more under
TGS and none under XPBD, against three to fifteen centimeters placed.

The damping is for scenes that start asleep, which are settled and then
put to sleep. Settled undamped, anything about to tip had already
started -- the half-overhanging ends of a running bond, the top of a
wall too tall for the solver -- and sleep froze it there, up to twelve
degrees over. Damped, nothing gathers momentum and the worst tilt is
under a degree. Zeroing velocities outright creeps, and takes twice as
long to get as far.

A fixed half second rather than "until at rest", which was tried and
cost seconds a reset: TGS creeps at a centimeter a second for longer
than anyone would wait, and by its own sleep speeds a five by four wall
never settled at all. It is also the wrong question for a scene that is
going to fail, which would be shown already fallen.

Running the solver is the point. The first version placed the bodies
quasi-statically -- stepping with every velocity zeroed after each step
-- and found a pose that was at rest but not the rest the solver keeps:
an eight course TGS column sat still for half a second and then rattled.
The solver's own rest is the only one it will hold.

# Astronomy

## Meeus, *Astronomical Algorithms*

Every chapter of the second edition that is an algorithm is implemented,
spread over the package by topic rather than gathered in one namespace,
because half of them already had a home. Where Montenbruck & Gill had
already supplied something -- Julian dates, sidereal time, nutation,
precession to J2000, Kepler's equation, vis-viva -- the Meeus chapter
extends that code instead of repeating it.

| Chapter | Topic                                          | Where                                                     |
|---------|------------------------------------------------|-----------------------------------------------------------|
| 3       | Interpolation                                  | `allgo.numerics.interpolation`                            |
| 4, 5    | Curve fitting, iteration                       | `allgo.numerics.fit`                                      |
| 7       | Julian Day, day of week and year               | `allgo.astro.time`, `allgo.astro.calendar`                |
| 8, 9    | Easter; Jewish and Moslem calendars            | `allgo.astro.calendar`                                    |
| 10      | ΔT                                             | `allgo.astro.time/delta-t`                                |
| 11      | The Earth's globe                              | `allgo.astro.geodesy`                                     |
| 12      | Sidereal time                                  | `allgo.astro.time/gmst`, `allgo.astro.frames/gast` (existing) |
| 13, 14  | Coordinate transformations, parallactic angle  | `allgo.astro.coordinates`                                 |
| 15      | Rising, transit and setting                    | `allgo.astro.rise`                                        |
| 16-20   | Refraction, separation, conjunctions, alignments, smallest circle | `allgo.astro.coordinates`              |
| 21, 24  | Precession, reduction of elements              | `allgo.astro.precession`                                  |
| 22      | Nutation and obliquity                         | `allgo.astro.frames`                                      |
| 23      | Apparent place of a star                       | `allgo.astro.apparent`                                    |
| 25-29   | The Sun, seasons, equation of time, solar disk | `allgo.astro.solar`                                       |
| 30, 34, 35 | Kepler's equation, parabolic and near-parabolic orbits, ellipse length | `allgo.astro.kepler`           |
| 31, 37-39 | Mean elements, Pluto, apsides, nodes         | `allgo.astro.planet-orbits`                               |
| 32      | VSOP87                                         | `allgo.astro.vsop87`, `allgo.astro.vsop87-data`           |
| 33      | Elliptic motion                                | `allgo.astro.elliptic`                                    |
| 36      | Planetary phenomena                            | `allgo.astro.phenomena`                                   |
| 40      | Parallax                                       | `allgo.astro.coordinates`                                 |
| 41, 55  | Illumination, magnitudes, semidiameters        | `allgo.astro.illumination`                                |
| 42, 43, 45 | Mars, Jupiter, Saturn's ring                | `allgo.astro.physical`                                    |
| 44, 46  | Satellites of Jupiter and Saturn               | `allgo.astro.jupiter-moons`, `allgo.astro.saturn-moons`   |
| 47, 48, 53 | The Moon's position, phase and libration    | `allgo.astro.moon`                                        |
| 49-52   | Lunar phases, apsides, nodes, declinations     | `allgo.astro.lunar-events`                                |
| 54      | Eclipses                                       | `allgo.astro.eclipse`                                     |
| 56, 57  | Stellar magnitudes, binary stars               | `allgo.astro.stars`                                       |
| 58      | Sundials                                       | `allgo.astro.sundial`                                     |

Chapters 1, 2 and 6 (hints, accuracy, sorting) are advice rather than
algorithms, and have no code.

**Conventions.** The package's, not the book's: radians, MJD (TT unless
a name says otherwise), geographic longitude positive east, azimuth
from north. Meeus counts longitude west and azimuth from the south, so
his hour angle θ0 − L − α is θ0 + L − α here and his azimuths are 180°
off. VSOP87 and the planets' orbits are in AU; the Moon's distance is in
kilometers.

**Where the data comes from.** Nothing was typed in from memory. The
coefficient tables -- nutation, the Moon, the lunar events, Pluto, E5 for
Jupiter's moons, Saturn's satellites, the phenomena -- were converted
by script from Sonia Keys's MIT-licensed Go port of the book
([soniakeys/meeus](https://github.com/soniakeys/meeus)), whose tests
also supplied most of the worked examples the `meeus-*-test` namespaces
check against. VSOP87 came from the original files instead (CDS
catalog VI/81), since Meeus's appendix is a truncation that no open
source reproduces term for term.

**The VSOP87 truncation.** Every term of at least 2e-7 rad in longitude
and latitude, and 2e-7 times the mean distance in radius: 2802 terms,
about the size of Meeus's appendix. Against the full series at forty
random epochs from 2000 BC to AD 3000 the worst errors were 1.7″ in
longitude (Jupiter), 1.8″ in latitude (Saturn) and 1.5″ in the radius
as seen from the Sun. 1e-7 would halve those at 3755 terms; 5e-7 doubles
them at 1895. The one place it shows in the book's examples is an apsis:
the distance is flat there, and Saturn's 1944 perihelion comes out 1944
September 7 at 22h rather than on the 8th.

**What the book had and this does not.** Table 36.B's rows could only be
checked for the phenomena Keys transcribed, so Venus's superior
conjunction, elongations and stations, the other stations of Mercury and
Mars, and the conjunctions of Mars, Jupiter, Uranus and Neptune are
missing rather than entered unverified. Table 31.B, the mean elements
referred to J2000, is derived from 31.A by precessing the orbit
(agreeing with chapter 24's formulae to 1e-12), not copied.

**Changes to what was there.**

- `time/calendar->mjd` put every Julian-calendar date -- anything before
  1582 October 15 -- two days early: it used B = -2 where the formula
  has B = 0. The existing tests only looked at dates after the reform.
- `frames/nutation-angles` now sums the 63 terms of Meeus's table 22.A
  rather than 20, which takes the truncation from about 0.01″ to
  0.0003″. The frame cache makes the extra terms free for the orbit
  propagator.
- One transcription error in the source was corrected: Venus's 1984
  Almanac magnitude has +0.000239 i² − 0.00000065 i³, not the reverse,
  which would have Venus brightening as it thinned.

## The solar system demo on Meeus

The demo moved onto the Meeus work in five steps; what was decided and
measured along the way:

- **Sources.** Planets from VSOP87 and the Moon from ELP, carried into
  EME2000 (`vsop87/equatorial-J2000`, `moon/geocentric-J2000`). The panel
  switches back to Standish's elements and Montenbruck & Gill's lunar
  series, and the readout gives the angle between the two for the body in
  focus -- 3 to 6 arcminutes for the giants, about 1 for the Moon.
  Timed in an unoptimized ClojureScript build under Node, one frame's
  worth costs: all eight VSOP87 planets 0.58 ms, the Moon 0.09 ms,
  Jupiter's moons 0.24 ms, Saturn's 0.19 ms, Pluto 0.02 ms. The Moon's
  month-long orbit line is 130 ELP evaluations, so it is resampled only
  when the Moon has moved a fifth of a day.
- **Pluto** is drawn only over 1885-2099, the span of Meeus's fit, and
  its orbit only along that arc. Its rotation model came from
  `pck00011.tpc`, like the others.
- **Moons in three dimensions.** Meeus's satellite theories end in a
  projection onto the sky as seen from the Earth; `positions-3d` stops
  before it. Tested by the moons lying in their planets' equators as the
  IAU poles have them (Galileans within 0.6 degrees, Saturn's inner seven
  within 2) -- two independent sources agreeing. The moons are drawn at
  true scale against their planet, which is itself drawn hugely
  oversized, so they appear only for the planet in focus.
- **Shadows** are traced, not mapped: a shadow map covering the solar
  system would put Io's shadow across one texel. Each fragment of a
  giant or its moons tests its line to the Sun against the others'
  spheres, with a penumbra as wide as the Sun's disk seen from there. Io's
  shadow on 2024 December 2 at 02:40 TT falls beside the Red Spot. The
  Earth and Moon are not drawn to scale with each other, so a lunar
  eclipse is shown by tinting the Moon from the true geometry (Danjon's
  shadow radii) rather than by tracing.
- **Magnitudes.** The readout uses the 1984 Astronomical Almanac formulas.
  Mueller's, which Meeus prints first, gave Venus -3.7 on 2024 December
  2 against the -4.2 observed; the Almanac's give -4.2.
- **The sky from the Earth** turns the star cloud by one matrix per frame
  -- J2000 to true of date to the horizon -- and refracts in the vertex
  shader by the formula `coordinates/refraction` uses for the bodies.
  Checked against London on 2024 December 2: sunrise 07:45 and sunset
  15:54 UT, Jupiter rising 16:08 five days before opposition; Orion
  upright on the meridian at sidereal time 5h18m; and the first quarter
  of December 8 lit toward the setting Sun with Mare Crisium on the
  eastern limb.

## SGP4, from *Revisiting Spacetrack Report #3*

`allgo.astro.sgp4` is transcribed from the C++ in the code package for
AIAA 2006-6753 (SGP4.cpp, 12 March 2020), not from the paper's equations
or any other port: where anything disagrees, the reports and their code
are authoritative. Each expression keeps the C++'s names and order of
operations, because the aim is the C++'s numbers, not just its theory.
Against the package's published output -- the C++ `.e` files and the
MATLAB file in opsmode a, the shipped Java conversion in opsmode i -- every
printed digit is reproduced, on the JVM and in ClojureScript alike.
Against a build of the package's C++ at full precision, 758 of 1338 states
agree to the last bit and the rest to 3e-8 km, the platform's sin, cos
and pow accounting for what is left.

What it took to get there, none of it visible in the theory:

- **`rem` is not `fmod`.** Clojure's `rem` on doubles computes
  x - trunc(x/y) y, which rounds; ClojureScript's does the same. C's
  `fmod` is exact. It matters most for the sidereal time, a large multiple
  of 2 pi reduced, which came out 1e-12 off. The namespace has its own
  exact `fmod`: `IEEE-remainder` shifted by y where its sign disagrees
  with x on the JVM, `js-mod` in JavaScript. C's `fmod` also keeps the
  sign of x, which Clojure's `mod` does not; several of SGP4's angles go
  negative.
- **Kepler's equation's sine and cosine lag a step.** The C++ loop takes
  sin and cos of `eo1` at the top of each pass, so after the loop they
  belong to `eo1` before its last correction. A tidier loop that returns
  sin and cos of the final `eo1` does not match.
- **`gsto` is `gstime` in both opsmodes.** The C++ still computes AFSPC's
  1970-based sidereal time in `initl`, then discards it. Ports of older
  versions of the code use it in opsmode a.
- **The epoch is split, and years divisible by four are leap years.**
  `days2mdhms` and `jday` give a Julian day and fraction, and `sgp4init`
  gets their sum less 2433281.5.
- **Grouping matters to the last bit.** `0.375*j2*tsi/psisq*con41` is not
  `0.375*j2*(tsi/psisq)*con41`, and `pow(x, 4.0)` is not `x*x*x*x`;
  dspace's `nm = no + (nm - no)` round trip is kept too.
- **The integrator's memory is a cache.** The C++ keeps the resonance
  integrator's last step in the element record and resumes from it. The
  steps fall on the same 720-minute grid from the epoch either way, so
  `sgp4` here integrates from the epoch each call and stays pure.

The published files need two allowances. For a case that fails at the
epoch (33334), the drivers still print a row of whatever vectors the
previous case left. The MATLAB file ends with a stray copy of 08195's
120-minute row. The tests drop both and say so. The package's FAQ
recommends opsmode a to match the US Air Force, and it is the default.
