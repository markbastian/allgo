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
  every corner boxes. Giving those an unboxed path would widen the gap
  and speed up `gradient-basis-4d` too. *That is the obvious next
  optimisation here.*

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

# Terrain Generation

## Diffusion-Limited Aggregation (DLA)

Particles move randomly until they stick to existing structures, naturally
growing branching patterns that resemble organic mountain ridges.
[1](https://www.youtube.com/watch?v=gsJHzBTPG0Y&t=568)

# Island Generation

- https://www.redblobgames.com/maps/mapgen2/
- https://www.redblobgames.com/maps/mapgen4/
- https://www.redblobgames.com/maps/terrain-from-noise/islands.html

Distinct from what `allgo.procedural.planet` does: these build a *polygon*
map -- Voronoi cells with elevation, moisture and biome assigned by graph
traversal from the coast -- rather than sampling a continuous function.
The polygon map is what makes rivers, watersheds and named regions
tractable, which a pure noise field does not give you.

# Map Generation

- https://azgaar.github.io/Fantasy-Map-Generator/
