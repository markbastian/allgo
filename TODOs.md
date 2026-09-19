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

# Terrain Generation

## Diffusion-Limited Aggregation (DLA)

Particles move randomly until they stick to existing structures, naturally
growing branching patterns that resemble organic mountain ridges.
[1](https://www.youtube.com/watch?v=gsJHzBTPG0Y&t=568)

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

- We should implement Sequential Impuse (SI) and Temporal Gauss-Seidel (TGS)
- What other physical simulation techniques are we missing?
- I'd like a demo where I get to shoot a 3D spherical projectile at a wall of bricks and have the physics be fast
  - Would it make sense to use all three techniques (the above plus XPBD) to show differences?
