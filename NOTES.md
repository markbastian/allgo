# Notes

How the finished parts of this repository work, what was measured along
the way, and what was tried and not kept. How to build and run is in the
[README](README.md).

The namespace docstrings are the reference. These notes are the
reasoning that does not fit in one docstring: decisions that cut across
namespaces, and the numbers behind them.

# Procedural generation

## Texturing & Modeling

Fifteen namespaces under `allgo.procedural` follow Ebert, Hart,
Musgrave, Peachey, Perlin and Worley, *Texturing & Modeling: A
Procedural Approach*, and are arranged the way chapter 20 argues they
should be -- a small set of parts that compose freely, rather than a
set of finished effects:

| Namespace                     | Chapter | What it is                                                       |
|-------------------------------|---------|------------------------------------------------------------------|
| `allgo.procedural.shaping`    | 2       | `step`, `smoothstep`, `bias`, `gain`, `spline` -- the vocabulary |
| `allgo.procedural.noise`      | 2, 12   | Value, gradient and 4D gradient noise, sparse convolution, warps |
| `allgo.procedural.surflet`    | 12      | Noise as a sum of surflets: separable, radial, editable by hand  |
| `allgo.procedural.cellular`   | 4       | Worley's F1, F2, F2-F1 and the metrics that reshape the cells    |
| `allgo.procedural.texture`    | 2-4, 15 | Brick, stone, metal; clouds, water, earth, gas giants; bump maps |
| `allgo.procedural.antialias`  | Worley  | Box-filtered patterns, octave clamping, and three supersamplers  |
| `allgo.procedural.volume`     | 6-8, 12 | Hypertexture, fire, steam, cumulus, cirrus, stratus, ray marched |
| `allgo.procedural.gesture`    | 12      | Procedural gesture: sines and noise, layered and blended         |
| `allgo.procedural.instancing` | Hart    | Recursive, parameterized instancing: trees, forests, fractals    |
| `allgo.procedural.lsystem`    | Hart    | L-systems: rewriting, and a turtle with tropism                  |
| `allgo.procedural.fractal`    | 14, 16  | fBm, turbulence, and the three multifractals                     |
| `allgo.procedural.qaeb`       | 17      | Error-bounded ray marching of a surface that is only a function  |
| `allgo.procedural.atmosphere` | 18      | Exponential density, Rayleigh and Mie scattering, aerial haze    |
| `allgo.procedural.genetic`    | 19      | Textures as expression trees, bred by mutation and crossover     |
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

Geometry is decided the same way. Hart's "Procedural Synthesis of
Geometry" (its number in the third edition is not checked here, hence
the author in the table, as for Worley's "Advanced Antialiasing") describes a model as rules that place copies of
each other, and `instancing/expand` measures each copy before expanding
it, keeping a proxy -- a crown for a bough, a cube for a sponge -- for
anything smaller than the cutoff. It expands the largest copies first,
so when its budget of primitives runs out what is left undone is the
finest level, everywhere; depth first spent the whole budget of a
Menger sponge on one corner and drew the rest as eight big cubes. And
because every random choice is drawn from the instance's own seed, a
branch is the same branch at any detail, which the tests check.

The antialiasing demo is the argument for integrating where possible.
On the checkerboard at its defaults, the analytic filter costs one
sample a pixel and leaves the horizon a smooth gray; jittered 4x4
supersampling costs sixteen and leaves it noisy; adaptive sampling cost
28 a pixel there and still drew moiré at the horizon, because when all
five of its first samples land on squares of one color it believes the
pixel is flat. Its corners are not shared between pixels, which a real
renderer would do, and which would bring that count down -- but not fix
the moiré.

Perlin's surflets are the same noise read another way, and
`surflet/basis` proves it: summing one hump-weighted ramp per lattice
point reproduces `noise/gradient-basis` to 1e-15, because a corner's
blend weight along an axis is one minus the fade, which is the hump.
What the reading buys is locality -- an edited gradient changes the
noise within a cell of its point and nowhere else, which the tests
check point by point -- and the freedom to make the hump round, at the
price of reaching 64 lattice points in 3D instead of 8.

`gesture` is the other half of the same chapter: noise used for time
rather than space. Its actions are ranges of joint angles driven by
sines and noise, in three layers that blend; a change in the middle
of a change freezes the half-finished blend as its starting point, so
the angles never jump, and the legs are solved by two-bone IK rather
than animated, so the feet stay planted whatever the hips do.

The volumes are animated Ebert's way, by moving the noise through a
gas that stays put, and every one but the two solid hypertextures takes
`:time`; the catalog marks which (`:moves?`), and the tests hold it to
that. Playing them is the problem: a frame is hundreds of thousands of
ray marches on the CPU, so the demo renders a loop once, small -- a
frame of 240 pixels takes about two seconds -- at eight frames a second,
and plays it back fading each frame into the next, which gas forgives.
Noise in time never repeats, so a quarter as many frames again are
rendered past the end and faded into the start; measured frame to frame
across the seam, the picture changes no more there than anywhere else.

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

Thermal erosion is `erosion/slump!`, the talus slippage Musgrave, Kolb
and Mace pair with the water: material slides wherever neighboring cells
are steeper than `:talus`, and run after the rain it puts screes under
the cliffs the droplets cut. Two things changed from the paper, and both
were measured rather than assumed. Every move in a pass is computed from
the same heights before any is applied, since sweeping in place drags
material along the sweep. And material moves by pairs of cells, at most
an eighth of each excess a pass, instead of the paper's share of a
cell's largest excess split among its lower neighbors: applied all at
once that rule took a one-cell pit to nearly four times its rim in one
pass and left it a mound, and by pairs no pass can make a new peak or
pit. Mass is conserved exactly -- nothing slides off the map -- and the
default talus of 0.01 touches only the steepest few percent of an
eroded diamond-square map.

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
population, history.

## Cities

Six namespaces under `allgo.procedural.city`, from three papers:
Parish and Müller, "Procedural Modeling of Cities" (SIGGRAPH 2001),
for the pipeline; Lechner, Watson, Wilensky and Felsen, "Procedural
City Modeling", for the agents; and Kelly and McCabe's survey of the
field, for Greuter's buildings. Plane geometry is in
`allgo.geometry.polygon`.

| Namespace   | From                          | What it is                                                                        |
|-------------|-------------------------------|-----------------------------------------------------------------------------------|
| `land`      | Parish & Müller §2            | The input maps, as functions: elevation, water, population, street pattern, zones |
| `roads`     | Parish & Müller §3            | Highways and streets by the extended, self-sensitive L-system                     |
| `lots`      | Parish & Müller §4.1          | Blocks as faces of the road graph, cut into lots                                  |
| `buildings` | Parish & Müller §4.2, Greuter | A building per lot by a stochastic L-system from its bounding box                 |
| `facade`    | Parish & Müller §4.3          | Layered grids: interval groups, layers, stacks                                    |
| `agents`    | Lechner et al.                | Extenders, connectors, developers and speculators on a grid of patches            |

**The L-system is a queue.** The paper's road productions only decide
*when* each proposed segment is considered: a branch waits out its delay,
a road module proposes its successor, and everything else -- where a road
goes and whether it may -- is in the global goals and local constraints
they call. So the productions are run as what they come to, a queue of
proposals ordered by when each is due, and popping one and checking it is
the insertion query `?I`. The paper's own rules are all there: highways
cast a fan of rays and turn toward the population they sum; streets keep
to the pattern of their district and stop where nobody lives; patterns
blend by their weights; a segment into water is pruned, turned or (a
highway) bridged; one that crosses a road stops there, one ending near a
junction is pulled to it, one ending short of a road is stretched to it.
One rule is added: a highway with nobody ahead of it goes straight,
because the ray casting otherwise turns it back toward the people behind
it, and the first highways all closed into rings.

**Blocks fall out of the junction rule.** Every crossing is made by
cutting both roads at one point, so the network is a plane graph with
exactly shared junctions, and its faces are found by walking it --
sharpest left at every junction -- after pruning dead ends. A block is
its face shrunk by half of each surrounding road's width. Lots are the
paper's recursive rule, a cut across the longest side until under the
zone's area, and as there, lots that are concave, too small or without a
street side are thrown away: for seed 2, 3714 lots on 68% of the block
area, the rest being the middles of deep blocks.

**A building's rounds are its levels of detail.** The axiom is the lot's
bounding box at its height, every production refines a box into what is
inside it, and the tests hold every round inside the axiom. So round 0 is
the city as boxes and the last is the city, and the demo's `detail`
slider walks between them -- the paper notes this and does not show it.
Downtown towers step back as they rise, offices are Greuter's plans (each
section down adds a primitive), low commercial blocks carry plant, and
houses get gable or hipped roofs. For seed 2 nothing needs more than 7
rounds; the tallest building is 226 m.

**Facades.** A style is a stack of layers, each two interval groups and a
function of their pulses, with the wall at the bottom. The house style has
the paper's own example of a function between layers -- lintels of darker
brick that exist only over an opening -- and the ground floor is always a
layer of its own, which is why the 3D demo bakes it as a separate texture
that repeats only along the street. Window lights are hashed through
Perlin's permutation table: a seeded stream folded over neighboring cells
comes out a constant step apart, and lit whole runs of windows at once.

**Agents.** Lechner et al. ran in NetLogo; this is arrays mutated in
place, so the browser steps it every frame. Their grid constraint is
taken as "every patch of a road within the deviation of one grid line":
at 0 the town is a strict grid, at 1 a loose one, from half the spacing
up organic -- their Figure 5's three towns. A test checks that commerce
sits on thicker road than housing, which is the pattern their Chicago
comparison is about. Two choices are mine where the paper gives none: the
weights in the land value, and that a road laid over a building clears
it.

**Speculators.** The paper's future work gives them one paragraph:
"pseudo-developers" reserving large tracts to sell to developers later,
for residential districts, universities and parks. Everything past that
is this code's. A speculator roams the edge of town and, on land the roads
do not serve yet, reserves the free patches joined to it within a disc of
100 patches (it keeps the tract if it gets 60 of them), holding one tract
at a time. Rough or waterside land it holds for a park half the time;
flat land for a campus, one per town; the rest for housing. It sells
once the roads serve a quarter of the tract and it has been held 120
ticks, and gives up a tract held five times that. Sold, a residential
district takes houses only, a park takes nothing -- not even a road --
and counts for housing as water does, and a campus is built by an
institution one hall at a time, never beside another. Commerce is what
they cost: on seed 4 it falls from 829 patches to 412 in 1500 ticks,
shut out of the districts and parks, while the speculators add nothing
measurable to the time a tick takes.

**Same seed, same city -- on one platform.** The land and the RNG are
exact on both, but the road rules take `sin`, `cos` and `exp`, which the
JVM and V8 may round differently in the last place, and over a thousand
decisions one near-tie flips: seed 2 grows 1211 segments on the JVM and
1187 in the browser. Each is repeatable on its own platform. A related
trap found here: ClojureScript's `(int-array n)` is `n` nils, not zeros,
so the agents' arrays are all given their zero explicitly.

**Not implemented**: the paper's smoothing of sharp bends, its streets
lowering the population around them, the zone and height *maps* as input
images (they are functions of the land here), Sun et al.'s template road
networks and Wonka's split grammars (the survey describes both too
briefly to build from), and the agents' multi-patch plots, which the
paper lists as future work beside the speculators.

## Search, after Russell & Norvig

Chapters 3 and 4 of *Artificial Intelligence: A Modern Approach* (4th
ed.), from the book's pseudocode (aimacode/aima-pseudocode, CC BY 4.0).
There is no copy of the book in `nocommit/`; the figures it gives for its
worked problems are in the tests, with the problems themselves.

| Namespace                        | Chapter | What it is                                                                 |
|----------------------------------|---------|----------------------------------------------------------------------------|
| `allgo.search`                   | 3       | Best-first (BFS, DFS, uniform cost, greedy, A*, weighted A*, beam), deepening (depth-limited, iterative, IDA*), RBFS, bidirectional |
| `allgo.search.local`             | 4       | Hill climbing and its variants, random restarts, simulated annealing, local beam, the genetic algorithm |
| `allgo.search.nondeterministic`  | 4       | AND-OR search for conditional plans; sensorless problems as belief-state search |
| `allgo.search.online`            | 4       | Online DFS and LRTA* agents                                                 |
| `allgo.search.problems`          |         | n-queens, a traveling salesperson, a maze, sudoku and map coloring, for the demos |

**Every search is a step.** `allgo.search` began as breadth-first,
depth-first, Dijkstra, greedy and A* written as one step function over
one search state, run by `iterate`, the algorithms differing only in
their frontier and its priority. Everything added keeps that shape:

- Weighted A* is A*'s step with another priority; beam search is A*'s
  step on a priority frontier with a size limit, which drops its worst.
- Depth-limited search, iterative deepening and IDA* are one step over an
  explicit stack of nodes, differing only in what they bound -- depth or
  g + h -- and the deepening ones restart from the root, inside the step,
  with the bound raised to the smallest that cut something off.
- RBFS and AND-OR search are recursive, and run by `call-seq`: recursion
  reified as a stack of frames, where each step is a call (push a frame)
  or a return (pop one and hand its value to the frame below). RBFS's
  backed-up f-values and AND-OR's mutually recursive OR and AND searches
  are each one frame function.
- Local search steps a move, or a generation; random restarts are a
  composition of any of them, concatenated seqs from fresh states.
- An online agent is `(fn [memory percept]) -> [action memory']`, the
  book's persistent tables passed back in, and `run-seq` iterates it in
  its environment.
- Sensorless search is no new algorithm: a function from a problem to
  a belief-state problem that breadth-first search then solves.

So every algorithm can be watched, stopped early, or counted, the same
way; the search demo is one loop over whichever seq it is given.

**A bug found on the way.** Dijkstra and A* stopped as soon as the goal
was *reached*, so on a graph where a costly road reaches the goal before
a cheap one does, they returned the costly one. The book's uniform-cost
search keeps going while the frontier holds anything cheaper than the
best solution; here they stop when the goal comes off the frontier,
which is the same. Settlements and the hex demo use A*; the hex tests
still pass, their moves all costing the same.

**The book's examples, as tests.** On Romania (roads as aima-python and
aima-java, both MIT, have them; aima-java alone has the book's
straight-line distances) uniform cost, A*, IDA*, RBFS, bidirectional and
beam search all find Arad-Sibiu-Rimnicu Vilcea-Pitesti-Bucharest, 418 km;
greedy goes by Fagaras, 450. AND-OR search on the erratic vacuum world
returns the book's plan -- Suck, then if the right square is still dirty
Right, Suck -- and the sensorless vacuum world is cleaned blind in four
actions (Left, Suck, Right, Suck; the book's is the mirror image).

**Measured.** 8-queens by steepest ascent, 2000 random boards: 13.6%
solved, in 4.0 moves when solved and 3.1 when stuck; with up to 100
sideways moves, 94.6%, in 19.6 and 64.6. These are the figures I recall
the book giving (14%, 4 and 3; 94%, 21 and 64) -- worth checking against
it. Random restarts average 6.3 restarts, as a 14% success rate implies.
On 20 8-puzzles 16 to 26 moves deep (mean 22), the steps to solve:

| Search                         | Steps (mean) |
|--------------------------------|--------------|
| Breadth-first                  | 68394        |
| A*, h1 (misplaced tiles)       | 11512        |
| A*, h2 (Manhattan distance)    | 1142         |
| IDA*, h2                       | 3253         |
| RBFS, h2                       | 3347         |
| Bidirectional uniform cost     | 1630         |

IDA* and RBFS do about three times A*'s work for memory linear in the
depth. Iterative deepening, a tree search, takes 1.5 s at 20 moves and
about three times longer a move deeper, so the tests run it only on
shallow boards; on a grid full of loops it, IDA* and RBFS revisit the same
cells so often that the demo gives up at a million steps -- which is what
it is there to show. A measuring loop that kept every state of a 200,000
step run exhausted a 17 GB heap; counting with `take-while` without
holding the seq does not.

**Choices of my own.** The tree searches -- deepening, RBFS -- skip a
child whose state is already on its own path, as the book's 4th-edition
IS-CYCLE does. The deepening step checks the bound before the goal, so
that depth-limited search to depth 2 does not report a three-move
solution. Online DFS does not record a move into a wall as a way back,
which the figure would, and then walk into the wall to backtrack.
Bidirectional search stops when the best meeting is no worse than the
lower bound -- the two frontiers' g's summed without a heuristic, the
larger f with one -- and when one frontier is empty, on the other's
alone; without that it stopped early on a directed graph whose backward
search had run out.

**Cyclic plans and percepts, added later.** AND-OR search with
`:cyclic?` returns `[:loop state]` where it meets a state already on its
path, rather than failing; an OR node takes a plan only if a goal is
reachable from every point of it, so a plan that could only loop is
refused. On the slippery vacuum world it finds the book's plan -- Suck,
then try Right until it works, then Suck. `partially-observable` makes a
problem with percepts an AND-OR problem over belief states, an action
predicting the states it may lead to and the percepts splitting them;
on the local-sensing vacuum world it plans the book's Suck, Right, then
Suck only if the square is dirty. `execute` follows loops, with a cap.

**Not implemented here:** continuous local search is
`allgo.numerics.optimize`, minimizing where the book maximizes; and the
genetic algorithm on tours, for which one-point crossover makes no tour.

## Games, after Russell & Norvig

Chapter 5, `allgo.search.adversarial`, with tic-tac-toe and Connect Four
in `allgo.search.games` for the demo; the book's two-ply game and a
chance game are in the tests.

**The same shape again.** Minimax, alpha-beta and expectiminimax are one
frame function for `call-seq` -- a decision node tries its moves as
child frames, a chance node sums its outcomes -- with pruning, a depth
cutoff, move ordering and a transposition table as options to it. To
count nodes and share a table across frames, `call-seq` gained an
optional context: state beyond the frames, which a frame may replace as
it calls or returns. The root's best move and each move's score live
there too, so a search part way through can be shown. Iterative
deepening is a composition, alpha-beta passes concatenated, each trying
the last pass's choice first and inheriting its table; Monte Carlo tree
search steps one playout at a time. The demo spreads every search across
animation frames -- a few thousand steps each -- so a full minimax of
tic-tac-toe runs in the browser without the page stopping.

**Checked against the book.** The chapter's two-ply game: minimax values MIN's
three nodes at 3, 2 and 2 and plays a1 for 3; alpha-beta finds the same
and never looks at the two leaves the book shows it pruning. Tic-tac-toe is a
draw by every search; full minimax visits 549,946 nodes (the commonly
quoted size of its tree), alpha-beta 20,866, with the center-corners-
edges ordering 7,865 and with a transposition table as well 2,178. On
Connect Four cut off at 4 plies alpha-beta visits 858 nodes to
minimax's 2,801 for the same value; at 8, ordered and with a table,
13,734 in 0.2 s on the JVM. A test checks, on random positions, that
alpha-beta's value is minimax's at the same depth with and without
ordering and the table, and that every search takes a win in one and
blocks a loss in one. Expectiminimax is tested on a chance game where an
order-preserving change to the utilities changes the move, which it
cannot do to minimax -- the book's point about games of chance.

**Choices of my own.** A win is worth 1 plus a hundredth for each empty
cell left, so a search that sees two wins takes the nearer; Connect
Four's evaluation (open lines of four weighted 1, 4, 16 by the cells
held, through tanh) stays inside 0.9, so a real win outranks any guess.
Under alpha-beta a move cut short at the root has only an upper bound
for a score, and the context marks it (the demo writes it with a ≤).
Monte Carlo picks the move with the most playouts, where figure 5.11's
pseudocode says the highest average, which a lucky move with few
playouts can top; UCB1's constant is sqrt 2. Its tree is keyed by state,
so transpositions share statistics.

**Not implemented here:** multiplayer games' vectors of utilities, and
forward pruning (beam search on moves, ProbCut), which the chapter
describes without pseudocode.

## Constraint satisfaction, after Russell & Norvig

Chapter 6 of *Artificial Intelligence: A Modern Approach* (4th ed.):
`allgo.search.csp`, with sudoku, n-queens and a random map to color in
`allgo.search.problems`; Australia and TWO + TWO = FOUR are in the tests.

| Algorithm            | Shape                                                                    |
|----------------------|--------------------------------------------------------------------------|
| AC-3                 | One arc revised a step, from a queue                                     |
| Backtracking         | `call-seq` frames, each holding its assignment and pruned domains; MRV with the degree heuristic, LCV, forward checking, MAC as options |
| Min-conflicts        | One conflicted variable moved a step                                     |
| Tree solver          | Arcs made consistent leaves up, then values assigned root down, a step each |
| Cutset conditioning  | Tree-solver seqs concatenated, one per consistent assignment of the cutset |

**Checked.** Australia has 18 three-colorings -- six for the mainland, three
for Tasmania -- and the cycle cutset found is SA alone, the book's own
example; with two colors backtracking finds none but AC-3 finds nothing
wrong, the case the book makes for path consistency. Forward checking
after WA = red and Q = green leaves NT and SA only blue, and V = blue
empties SA, as the book walks through. Every solution of TWO + TWO =
FOUR backtracking finds is checked against brute force, and it finds them
all. The tree solver agrees with backtracking on random tree CSPs, solvable
or not.

**Measured.** On 12-queens, plain backtracking makes 3,066 assignments;
MRV alone, the same, since without inference no domain shrinks; with
forward checking, 153; with MRV, LCV and MAC, 17. Min-conflicts solves
100 queens in 47 steps. Plain backtracking ran past ten minutes on a
generated sudoku with 55 blanks, which with MRV and MAC takes a second.
Cutset conditioning is exponential in the cutset, as the book says, and
the greedy cutset of a random map grows fast:

| Regions | Cutset | Steps        |
|---------|--------|--------------|
| 20      | 10     | 22           |
| 30      | 11     | 40           |
| 45      | 22     | 86,680       |
| 60      | 26     | over 300,000 |

So the demo offers it only on the map, at 30 regions by default.

**Choices of my own.** The cycle cutset is greedy -- most constrained
first -- not minimal. Higher-order constraints are checked once their
scope is assigned, so AC-3 and the tree solver see only the binary ones.
The sudoku generator fills a grid by backtracking with shuffled digits and
blanks cells while the puzzle keeps exactly one solution; the demo's four
puzzles are its output for seeds 11 to 14.

**Not implemented here:** path consistency and k-consistency, bounds
propagation, conflict-directed backjumping and constraint learning, and
tree decomposition, which the chapter describes without pseudocode.

## Logical agents, after Russell & Norvig

Chapter 7 of *Artificial Intelligence: A Modern Approach* (4th ed.), in a
package of its own:

| Namespace                     | What it is                                                                  |
|-------------------------------|-----------------------------------------------------------------------------|
| `allgo.logic.propositional`   | Sentences and a parser, CNF, TT-ENTAILS?, PL-RESOLUTION, forward and backward chaining, DPLL, WalkSAT |
| `allgo.logic.wumpus`          | The wumpus world and the hybrid agent                                       |
| `allgo.logic.satplan`         | SATPlan: a STRIPS-like problem to CNF for each plan length, solved by DPLL  |

Its tests are under `test/allgo/logic/`.

**Steps again.** Model checking steps a model, resolution a pair of
clauses, forward chaining a symbol off the agenda, WalkSAT a flip;
backward chaining and DPLL are recursive and run on `call-seq`; SATPlan
concatenates DPLL's seqs, one per plan length; the agent is a function
of its memory and a percept, iterated with the world.

**Checked against the book.** On the chapter's knowledge base after a
breeze in [2,1] -- R1 to R5 -- model checking (all 128 models of its seven
symbols), resolution and DPLL all entail no pit in [1,2], and none of them
decides [2,2]. Forward chaining on the chapter's definite clauses infers
A, B, L, M, P, Q in that order; backward chaining agrees. SATPlan solves
the 'have cake and eat cake too' problem with Eat, then Bake. DPLL is
checked against trying every model on 200 random 3-SAT problems.

**DPLL keeps its clauses simplified.** The figure tests every clause
against the partial model on every call; here each call carries the
clauses left under its model -- satisfied ones dropped, false literals
removed -- so a unit clause is one of length one and a conflict an empty
one. Same algorithm, same answers (checked against trying every model on
600 random problems); a wumpus safety query went from 12.6 ms to 1.5, and
200 games from 172 s to 22. It also splits on a symbol of the first
clause left rather than the next in a fixed order.

**A bug the agent found.** At first the world sensed a stench in the
wumpus's own square, as well as next to it. The physics says next to it
only, so once the agent stood on the dead wumpus its knowledge base
became contradictory -- and a contradiction entails everything, so every
square was proved safe and the agent walked into pits it had 'proved'
empty. A test now checks the physics against every cave's percepts, and
that the agent never dies in a square it proved safe.

**How it plays.** Over 200 random caves the agent brings out the gold in
101, climbs out empty-handed in 7, and dies in 92, for a mean score of
+30; it shoots in 148 games and kills the wumpus in 46. The deaths are
the book's policy: with no square left proved safe, it takes a risk on
one not proved deadly rather than give up.

**Choices of my own.** The agent knows where it is and which way it
faces by counting, so its knowledge base holds only the atemporal physics
and the percepts -- the book keeps location in the logic too, with
successor-state axioms. A shot without a scream is told directly as 'no
wumpus on the arrow's line', which the temporal axioms would infer. Only
squares next to visited ones are asked about each turn, and a proof once
made is kept, since the knowledge base only grows. Scores are
aima-python's (-1 an action, +1000 for the gold out, -1000 for dying);
the book's extra charge for the arrow is `:arrow-cost`, default -1 like
any action, since neither aima-python nor aima-java confirms it.

**Not implemented here:** the generic KB-AGENT (figure 7.1) as anything
but the hybrid agent's shape, and successor-state axioms for the agent's
own location.

## First-order logic, after Russell & Norvig

Chapters 8 and 9 of *Artificial Intelligence: A Modern Approach* (4th
ed.): chapter 8 is the language, without algorithms, and chapter 9 the
inference, in `allgo.logic.first-order`. Variables are lowercase symbols
and constants, predicates and functions capitalized keywords, so the
book's notation reads straight in: `parse` takes `Owns(Nono, M1)`, `=>`,
`forall x, y: ...`.

**Steps again.** Forward chaining steps a rule at a time through rounds;
backward chaining -- the book's generators FOL-BC-OR and FOL-BC-AND --
runs as one explicit stack of goal lists, a step resolving the top list's
first goal against every rule, which is the same depth-first order the
generators take; resolution steps a given clause at a time, the negated
query as its set of support, each clause keeping its parents so that a
refutation reads out as a proof.

**Checked.** The chapter's crime -- is Colonel West a criminal? -- as
aima-python (MIT) states it: forward chaining, backward chaining and
resolution all prove it, forward chaining in two rounds; resolution in
68 steps from the quantified sentences, an existential Skolemized. Every
answer aima-python's tests expect of its crime and 'Peter Rabbit'
knowledge bases comes out the same by both chainings, and its unification
tests pass, compared by what each substitution does rather than its
exact form. Append(x, y, [1, 2]) gives its three answers in depth-first
order.

**A bug found on the way.** Standardizing a rule apart renamed its
variables with SUBST, which passed the rule -- a map -- through
untouched, so a rule's `x` was the query's `x`: backward chaining could
not prove Criminal(x) though it proved every subgoal, and Append found
one answer of three. Forward chaining hid it, its facts having no
variables.

**Choices of my own.** Forward chaining keeps going past the first answer
and gathers them all, as aima-python's does; the book's stops at the
first. Backward chaining cuts proofs off at `:max-depth` (default 50),
since the wrong order of conjuncts sends depth-first search down for
ever -- aima-python's knowledge base carries a comment about exactly
that rule. Resolution orders its queue shortest clause first, and is
bounded by the caller: first-order entailment is semidecidable.

**Not implemented here:** equality and paramodulation, indexing and
retrieval beyond a linear scan, the Rete algorithm, and constraint logic
programming.

## The rest of the book, after Russell & Norvig

Chapters 2 and 11 to 26 of *Artificial Intelligence: A Modern Approach*
(4th ed.): every algorithm the book gives as pseudocode, and the main
ones it teaches in its text. Problems and expected answers come from
aima-python (MIT) and aima-pseudocode (CC BY 4.0); the restaurant
examples of chapter 19 are inlined in their test.

| Namespace                          | Chapters | What it is                                                                 |
|------------------------------------|----------|----------------------------------------------------------------------------|
| `allgo.agents`                     | 2        | The vacuum world; table-driven, reflex, simple reflex and model-based agents |
| `allgo.planning`                   | 11       | Classical planning by A* over ground actions; air cargo, spare tire, blocks world; hierarchical and angelic search; job-shop scheduling |
| `allgo.decision.utility`           | 12, 16   | Expected utility, the decision-theoretic agent, the value of information, the information-gathering agent |
| `allgo.probability.bayes`          | 13       | Bayesian networks: enumeration, variable elimination, prior, rejection and likelihood-weighted sampling, Gibbs |
| `allgo.probability.temporal`       | 14       | HMMs: filtering, forward-backward, fixed-lag smoothing, Viterbi, particle filtering; the Kalman filter |
| `allgo.probability.programs`       | 15       | Probabilistic programs: `sample!`, `observe!`, likelihood weighting       |
| `allgo.decision.mdp`               | 17       | Value and policy iteration on grid MDPs; POMDP value iteration            |
| `allgo.learning.trees`             | 19       | Decision trees, cross-validation and model selection, decision lists, AdaBoost, random forests |
| `allgo.learning.models`            | 19, 20   | Linear and logistic regression, the perceptron, k nearest neighbors, SVMs, naive Bayes |
| `allgo.learning.unsupervised`      | 20       | k-means, EM for Gaussian mixtures, Baum-Welch                             |
| `allgo.learning.probabilistic`     | 20       | Bayesian, MAP and ML learning; parameters by counting and by beta and Dirichlet priors; Gaussian and linear Gaussian ML; structure by hill climbing with restarts; k-NN and kernel densities; EM for a naive Bayes mixture |
| `allgo.learning.neural`            | 21       | Dense, convolutional, pooling, residual and batch-norm layers; back-propagation, minibatch SGD, weight decay, dropout |
| `allgo.learning.generative`        | 21       | Probabilistic PCA by EM; linear and deep autoencoders; a variational autoencoder; Yule-Walker and deep autoregressive models; a GAN |
| `allgo.decision.reinforcement`     | 22       | Direct utility estimation, passive ADP and TD, Q-learning, SARSA, REINFORCE; safe exploration's Bayesian and robust policies; reward shaping; linear and deep function approximation, experience replay |
| `allgo.decision.hierarchical`      | 22       | Hierarchical RL: partial programs with choose, learned by Q-learning over the joint state space |
| `allgo.decision.apprenticeship`    | 22       | Imitation learning; inverse RL by feature matching; Boltzmann-rational likelihoods |
| `allgo.language.parsing`           | 23       | CYK, A* and shift-reduce parsing with a probabilistic grammar; dependency parses from heads and back; treebank trees; a grammar learned from a treebank, and by inside-outside |
| `allgo.language.augmented`         | 23       | Augmented grammars by unification: agreement, lexicalized probabilities, compositional semantics with lambda-expressions |
| `allgo.language.ngrams`            | 23       | N-gram models: Laplace and interpolated smoothing, EM for the weights, perplexity, <UNK>, sampling, language identification |
| `allgo.language.tagging`           | 23       | Part-of-speech tagging with an HMM decoded by Viterbi, and by logistic regression decoded greedily or by beam search |
| `allgo.robotics.localization`      | 26       | Monte Carlo localization with a range scanner                             |

The tests are under `test/allgo/decision/`, `probability/`,
`learning/` and `robotics/`.

**Steps again.** Every iterative algorithm is a seq of its states: value
iteration a sweep, policy iteration an evaluate-and-improve, POMDP value
iteration a depth of plans; each sampler an estimate after each sample;
filtering a message per observation, Kalman a step; AdaBoost a stump,
k-means and EM an iteration, Baum-Welch a re-estimation, a network an
epoch, an SVM a pass of SMO; CYK a span length; Monte Carlo localization
a move, a scan and a resampling. Hierarchical search runs on the same
frontier as chapter 3's. Agents -- the vacuum world's and the
reinforcement learners -- are a function of memory and percept, run by
`run-seq` or `trial-seq` against the world.

**Checked against the book's answers.**

- Value iteration on the 4 x 3 world gives aima-python's utilities exactly,
  and policy iteration its policy.
- P(Burglary | JohnCalls, MaryCalls) = 0.284 by enumeration and by
  elimination, which agree on every query of the sprinkler network; every
  sampler comes near the exact answer, and Gibbs gives P(Cloudy | Rain) = 0.8.
- The umbrella world filters to 0.818 and 0.883. Fixed-lag smoothing
  gives aima-python's [0.2259 0.7741], and Viterbi its paths and messages.
  The Kalman filter's closed form gives mean 1.875 and variance 0.75.
- On the restaurant data, Patrons has the greatest information gain,
  0.541 bits, and the tree tests Patrons, Hungry, Type and Fri/Sat.
- The survey in the oil-block example is worth C/n.
- The blocks world plan is the figure's three moves.
- The car assembly's critical path is 85.
- "The robot is good" parses with probability 0.015.
- Back-propagated gradients match numerical ones to 1e-7, through a
  convolution and two dense layers.

**Dropout, section 21.5.4, checked against the 4th edition's text.**
Each training step goes through a thinned network: every hidden unit
(`:keep p`) and input unit (`:keep-input p`) is kept with probability p
and its output multiplied by 1/p, or fixed at zero. Back-propagation
runs through the same masks, and at test time the network runs with no
dropout. The book's p is the probability of *keeping* a unit, 0.5 for
hidden units and 0.8 for inputs; `with-dropout` takes them in that
order. The tests check:
- the masks' values and kept fractions, and that the output layer is
  never masked;
- the thinned network's gradients against numerical ones, and that a
  dropped unit learns nothing on that step;
- the book's remark that the ensemble claim holds analytically for
  linear models: over 20,000 masks, a linear network's thinned outputs
  average to the whole network's;
- that p = 1 trains exactly as without dropout;
- the book's "harder to fit the training set": on 50 noisy points the
  training loss is higher with dropout.

The other half of the claim, a lower test error, did not hold reliably
at this library's sizes. On noisy circles with 20 hidden units (five
seeds), dropout lowered the test loss on four and raised it on one,
where keeping half of twenty units left too little to learn the
boundary. On 20-dimensional data with 60 examples it lowered the test
loss on four of six seeds, and the accuracy barely moved. On noisy
curve fitting with an identity output it made things worse. The book
adds that dropout usually needs a larger model trained longer, and so
the tests do not assert it. The classifiers demo shows it instead: a
"wide network" (two hidden layers of 32) beside the same network with
dropout, with held-out accuracy now in the caption.

**The rest of chapter 21, from the 4th edition's text.** Each claim the
chapter makes is a test where it could be checked.

- *Pooling* (21.3.1): max-pooling the chapter's [5, 9, 4] gives 9.
  Average-pooling equals a convolution with the uniform kernel. A
  stride of l downsamples by l.
- *Residual layers* (21.3.3, equation 21.10): with V = 0, relus pass
  their input through unchanged, the chapter's derivation. Thirty
  layers deep, the first layer's gradient is about 1e-27 through plain
  sigmoid layers. Through plain relu layers it is 5e-13, or exactly
  zero for another seed. Through residual layers it is about 0.8.
- *Batch normalization* (21.4.2):
  - In training, each unit has mean beta and standard deviation gamma
    over the minibatch.
  - Its gradient runs through every example's dependence on the
    others, and is checked numerically.
  - With weights too small, the spread of a deep layer's values decays
    to near zero without it (4e-7 by the eighth layer here), but not
    with it.
  - Six small-weighted tanh layers, trained 40 epochs, reach a loss of
    0.13 to 0.33 and 92% to 98% accuracy with it, over the two seeds
    tried. Without it they are stuck at 0.69, which is chance.
  - The chapter does not say how to normalize at test time, when there
    is no minibatch. We use running averages of the minibatches'
    statistics, the usual way.
- Training now takes minibatches (`:batch`). A minibatch of one gives
  exactly the old per-example gradients, so every earlier test and
  demo is unchanged.
- *Weight decay* (21.5.3): on one linear unit, gradient descent with
  the penalty lambda w^2 finds ridge regression's cov(x, y) / (var x +
  lambda) to 1e-9. The penalty applies to the weights and kernels, not
  the biases.
- *Probabilistic PCA* (21.7.1, equation 21.16):
  - The chapter names EM but gives no update, so the step is Tipping
    and Bishop's.
  - The check does not rest on that paper: EM never lowers the
    likelihood, and it converges to the closed-form maximum. That is,
    sigma^2 is the mean of the discarded eigenvalues, and W W^T +
    sigma^2 I has the top eigenvalues and eigenvectors.
  - The chapter's model has mean zero; ours is fitted about the data's
    mean.
- *Linear autoencoder*: trained by Adam, W^T W becomes the projection
  onto the top two principal components (to 1e-9), as the chapter says.
- *Deep autoencoder*: on an arc, a one-number deep autoencoder
  reconstructs 70 times better than a linear one.
- *Variational autoencoder*:
  - The bound is computed in the chapter's own decomposition, H(Q) +
    E log P(z, x). The sample's gradient reaches the encoder by Kingma
    and Welling's reparameterization, which the chapter does not
    describe.
  - With a linear decoder the exact log P(x) is known, and the bound
    stays below it while rising toward it.
- *AR models*: the Yule-Walker fit recovers an AR(2) model. On the
  logistic map, the linear model's rms error is 0.24 and the deep
  one's 0.01.
- *GAN*:
  - On N(4, 0.5^2), the generator's mean comes out right, and the
    discriminator ends near 0.5 on both real and generated data.
  - The spread falls short, though: about 0.3 for the data's 0.5, a
    known tendency of GANs. A tanh discriminator did worse, about 0.15,
    so the discriminator uses relu.
  - The test asserts only the looser range.
- In 2-D, the GAN fits a ring well (radius 1.96 for 2), as the demo
  shows.

The new demo, *Deep Generative Models*, shows all of these on 2-D data
and on the logistic map. The classifiers demo gains a wide network with
weight decay, a deep network with batch normalization, and a residual
network. Left out: momentum and learning-rate schedules, which 21.4
mentions only as considerations; 21.7.2's transfer learning; and
unsupervised translation, which points ahead to chapter 25.

**The rest of chapter 22, from the 4th edition's text.** As with chapter
21, each claim that could be checked is a test.

- *Safe exploration* (22.3.2): the Bayesian and robust policies, each
  chosen among given candidates. The test has two hypotheses about the
  4 x 3 world, the pit costing 1 or 20. The Bayesian policy takes the
  risk while the costly pit has probability 0.05, and not at 0.3. The
  robust policy never takes it.
- *Widrow-Hoff* (equation 22.10): the chapter's worked example passes.
  With theta = (0.5, 0.2, 0.1), U(1,1) = 0.8, and an observed 0.4 lowers
  each parameter by 0.4 alpha.
- *The 10 x 10 claims* (22.4.1), checked by least squares on the trials'
  examples, which the chapter also names:
  - with +1 at (10,10), (1, x, y) fits the true utility to 0.015 rms;
  - with +1 at (5,5), the pyramid, it cannot do better than 0.11;
  - adding the distance to the goal brings that to 0.02.
- *Online fitting* is noisier. The Widrow-Hoff rule and TD (equation
  22.11) come within 0.06 and 0.08 rms over 2000 trials.
  - The step must be small: the coordinates run to 10, and with alpha =
    0.01 TD diverged, which is a step-size instability, not the
    divergence the chapter warns of.
  - Five trials from (1,1) leave most of the table empty, but the linear
    estimate is good there too: the generalization the chapter
    describes.
- *Linear Q-learning* (equation 22.12) with Q linear in x and y for each
  action comes within 0.02 of optimal on the (10,10) world.
- *Deep Q-learning* (22.4.3): on a 5 x 5 world with +1 in the middle,
  where no linear Q can do it, a network on (x, y) learns a nearly
  optimal policy.
  - It needed experience replay, centered inputs and 32 hidden units.
    Earlier tries -- no replay, 16 units or inputs from 0 to 1 -- failed
    in up to 600 trials: the chapter's
    "delicate art".
- *Reward shaping* (22.4.4): exact Q iteration shows that shaping by a
  potential changes Q by exactly -Phi(s), and the optimal policy not at
  all.
  - The potential matters for learning speed. Shaping a sparse 10 x 10
    world with Phi = -0.1 x distance made Q-learning far slower.
  - Shaping by a potential works like starting Q at that potential, and
    a negative one makes bumping into a wall, whose value settles toward
    0, look better than moving on.
  - With Phi = gamma^distance, an estimate of the utility, it cut
    epsilon-greedy Q-learning's moves by about two thirds.
  - A bonus that is not a potential difference changes what is optimal:
    beside the goal, it makes the agent stay rather than finish.
- *Hierarchical RL* (22.4.5): partial programs that choose, act or
  compute, with Q-learning over the choice states of the joint state
  space.
  - The rewards between two choices are discounted by gamma^k, the
    bookkeeping the chapter's footnote mentions.
  - The trivial program, while true do choose(A(s)), is ordinary
    Q-learning.
  - In two rooms, a program choosing door or goal, with shaping by the
    distance to the routine's target, learns in under half the moves
    flat Q-learning takes. That target is no part of the physical
    state.
  - Flat Q-learning shaped by the distance to the goal was faster
    still, since that potential is optimistic and drives exploration.
    The chapter claims only that such advice cannot be given to a flat
    learner, not that it beats everything.
- *Inverse RL by feature matching* (22.6):
  - Abbeel and Ng's projection method, in place of the chapter's
    max-margin step, which needs a quadratic program.
  - U(s) = theta . mu(pi) holds exactly.
  - From 20 of the expert's trajectories, and never its reward, it finds
    a policy as good as the expert's from the start and from every
    state.
  - Imitation learning from the same trajectories does worse, and worst
    off the expert's path.
  - Under Boltzmann rationality, R = 0 gives each choice probability 1/4
    and the true reward explains the trajectories far better: the
    chapter's argument against R = 0.

The grid-world demo gains three learners: TD with a utility linear in x
and y, Q-learning with Q linear in x and y for each action, and
Q-learning shaped toward the +1 exit.

**The rest of chapter 23, from the 4th edition's text and figures.**

- *Logistic-regression tagging* (23.1.6): one logistic model per tag,
  over features from the chapter's table. Those are the word, its
  neighbors, the tag before, the word two back with the tag before, and
  the spelling (suffixes, prefix, hyphen, digit, capitals).
  - It is trained by SGD with the true previous tags, and decoded by a
    beam search, greedy being b = 1.
  - It tags the ambiguous grammar's held-out sentences perfectly.
  - Restricted to the word and the tag before, greedy search gets the
    look-ahead sentences ("duck swim today" / "duck stones today") 84%
    right and a beam of 2 gets 93%. With the next word as a feature, as
    the chapter allows, greedy gets them all.
  - Words it never saw that end in "-ous" are tagged adjectives: the
    chapter's example of what an HMM cannot be told.
- *Dependency parses* (23.3.1, figure 23.7):
  - Head rules read off the figure turn its phrase structure tree into
    exactly its dependency tree: I and wumpus depend on detect; the,
    smelly and me on wumpus; near on me.
  - The figure's dependency tree makes the noun phrase, not the
    preposition, the head of "near me", and the head rules follow it.
  - Going back, the made-up categories the chapter mentions give a tree
    whose dependencies are the same again.
- *A\* parsing* (23.3): states are lists of items, the cost is -log of
  each rule's probability, and the heuristic is each unparsed word's
  cheapest lexical rule, which is admissible. Its first parse is CYK's,
  tree and probability, on 30 sampled sentences.
- *Shift-reduce parsing* as a beam search (23.3):
  - My first version let shift push a bare word at no cost. The beam
    then filled with configurations that shifted everything and could
    never reduce. So shift now takes the word as one of its parts of
    speech, and every parse is 2n - 1 actions.
  - With b = 1, the deterministic parser, only 4 of 30 sentences parse.
    With b = 4, 23 do, and with b = 16 every one gets the most probable
    parse.
- *Augmented grammars* (23.4): rules over terms with variables, a chart
  parser that unifies (`allgo.logic.first-order`'s unifier), and rule
  probabilities that may be functions of the bindings.
  - Figure 23.10's arithmetic grammar parses 3 + (4 ÷ 2) as Exp(5), with
    figure 23.11's value at every node.
  - Figure 23.12's grammar parses "Ali loves Bo" as S(Loves(Ali, Bo)),
    the VP being lambda x Loves(x, Bo) after beta reduction, and parses
    the figure's three other sentences.
  - Agreement follows figure 23.9, except that the figure's VP rule
    gives the object the verb's person and number. Here the object has
    its own. "I saw she", "I sees her" and "me see them" are rejected.
  - For the lexicalized PCFG the chapter gives P1 and P4 no values, so
    the test's are ours. With them, "ate a banana" beats "ate a
    bandanna", and head words choose [[spaghetti and meatballs] or
    lasagna] but [spaghetti and [pie or cake]].

The language demo gains parsing (a sentence's tree, its dependency arcs,
A*'s expansions, the shift-reduce beam needed) and semantics (figure
23.11's and 23.12's trees with each node's meaning). Left out:
learning semantic grammars (23.4.2), which the chapter only surveys.

**Chapter 20, checked against the book's numbers.** Written from the
3rd edition's chapter 20, then checked against the 4th's: the same
worked examples and numbers, renumbered in places (the kernels are
Figure 20.10, the nearest neighbors 20.9, EM's first iteration equation
20.11). They are the tests:
- *Bayesian learning:* the candy bags of Figure 20.1. h3 is the most
  probable after none and one lime, h4 after two, h5 after three or
  more, and P(h5) is 0.8956 after ten. After three limes the MAP learner
  says lime for certain where the Bayesian prediction is 0.796, "about
  0.8".
- *Beta and Dirichlet:* conjugate updates as virtual counts, and
  beta[3,1], [6,2], [30,10] narrowing round 0.75 as Figure 20.5 draws
  them.
- *Bayesian linear regression* (section 20.2.6, new in the 4th
  edition): a slope through the origin with a Gaussian prior and known
  noise. The closed form matches the posterior computed on a grid, a
  point at a time ends where all at once does, and the predictive
  variance sigma^2 + sigma_N^2 x^2 widens away from the origin.
  Figure 20.7 gives only sigma^2 = 0.2 and the slope variances 0.3861
  and 0.0286, not its data; taking the first as the prior and adding
  the two points the figure draws at about x = +-1.8 gives 0.0286.
- *EM on two bags mixed* (section 20.3.2): from the book's table of a
  thousand candies and its starting values, one iteration gives
  equation 20.11's parameters to four decimals -- 0.6124, 0.6684, 0.6483,
  0.6558, 0.3887, 0.3817, 0.3827 -- and the log likelihood goes from
  about -2044 to -2021. By the tenth iteration it passes the generating
  model's -1982.214.
- *Structure learning:* hill climbing over links added, deleted or
  reversed, scored by log likelihood less BIC's penalty. A single climb
  from no links found the sprinkler network from 3000 samples on only
  five of ten samples: the others stopped at local maxima that the
  true network outscores. Ten random restarts (chapter 4's remedy; the
  section itself only says to start from no links or from a guess)
  found it -- its links, and the two into WetGrass directed -- on all
  ten, at about 1.5 s each. With no penalty the search adds links the
  network does not have, as the chapter says it must.
- *Density estimation* on three Gaussian blobs: Gaussian kernels of
  width 0.07 come nearer the true density than 0.02 or 0.2, and k = 10
  neighbors nearer than 3 or 40, as Figures 20.9 and 20.10 show. Held-out
  likelihood is the wrong judge for nearest neighbors -- their estimate
  is not a density, its integral diverging -- so the tests measure the
  distance from the truth instead.

The demo, Learning Probabilistic Models, has five modes:
- *Candy:* Figure 20.1 drawn as candies are unwrapped.
- *Beta:* the posterior narrowing.
- *EM:* the mixed bags, iteration by iteration, against the true model.
- *Density:* the estimate beside the truth, for each width or k.
- *Structure:* the search, move by move and restart by restart, walked
  lazily so the page never stalls.

**Value iteration returns the book's U.** The book's loop updates U' from
U and returns U when δ is small enough, so the answer is the utilities
before the last sweep. Returning U' instead was 6e-4 off aima-python's
values; returning U, as the book does, matches them exactly.

**POMDP plans pruned on a lattice of beliefs.** The book prunes dominated
conditional plans; doing that exactly takes a linear program per plan,
and ours was slow and failed on degenerate ones. Instead a plan is kept
only if it is best at some belief of a lattice (200 points a side for two
states), which keeps every plan that matters except on slivers between
lattice points (a coarser lattice past three states). Each depth's plans are worth what expectimax lookahead
over beliefs says, to 1e-9 at every belief tried.

**Reinforcement learners that finish.** The learners follow the book --
α = 60/(59 + n), exploration f(u, n) with Ne = 5 and R+ = 2 -- with three
additions. Ties between actions are broken at random: with a fixed order
SARSA kept bumping into a wall at a Q of -0.4 for ever. A trial stops
after `max-moves` (default 10000). And SARSA takes an `:epsilon`, the
chance of a random action, since on-policy learning with only the
book's optimism stops exploring once the counts reach Ne; the demo uses
0.1. Q-learning's test uses Ne = 30, as with 5 some seeds settle too
soon on a worse route. REINFORCE learns a softmax policy over a preference for each state
and action, from the returns of whole trials.

**One linear regression.** The chapter teaches gradient descent on the
squared loss as well as the closed form. Linear regression here is the
closed form only (`allgo.numerics.regression`), being exact and faster;
`gradient-descent-seq` is the one optimizer, used by the models that have
no closed form, and its test checks it reaches the closed form's
weights.

**Baum-Welch kept monotone.** The first version let the likelihood fall:
it re-estimated the prior from the first observation's smoothed state
rather than the state before it, and left the transition from that state
to the first out of the counts. With both fixed the likelihood never
falls, as EM promises.

**Scheduling: minimum slack is not optimal.** The book's minimum-slack
heuristic, with one engine hoist, schedules engine 2 first and takes 130
minutes; `optimal-schedule`, a branch and bound over the order, finds
115 by doing engine 1 first. The test asserts both.

**Angelic search needs its progress check.** Without the figure's
MAKING-PROGRESS test, a high-level action that refines to itself
recursed until the stack overflowed.

**Localization that recovers.** Monte Carlo localization in a 10 x 6
room with an L-shaped wall and a pillar, eight range rays a scan. The
plain algorithm collapsed onto a wrong place now and then, and stayed
there. What fixed it: a few percent of particles thrown anywhere each
step (the book's remedy for the kidnapped robot), motion and sensor noise
broader than the truth's, a route that runs past what tells look-alike places apart, and an estimate taken from the densest cluster rather than the
mean of every particle. From not knowing where it starts, the estimate
is within a meter at the end of the route on every seed tested, and
within half a meter on most.

**Chapter 15 from my own examples.** The book's programs for the chapter
-- citations, NET-VISA, radar tracking, image generation, Markov letters --
are listed in aima-pseudocode's index, but its files for them are missing
(404). So `allgo.probability.programs` gives the idea -- a model as a
program that samples and observes, inferred by likelihood weighting --
on two examples of my own, each checked against an exact answer: a
Markov chain of letters read through noise, against HMM smoothing of
the same model, and an open-universe count of aircraft on a radar that
misses some and invents others, against summing over the detections.

**Demos.** Grid World runs every chapter 17 and 22 method on the 4 x 3
world; Classifiers runs ten learners on five two-dimensional datasets,
drawing their regions; Localization runs MCL with the scan and the
particles.

**Not implemented here:** the book's chapter 15 programs themselves;
exact POMDP pruning by linear programming; GraphPlan, which the 4th
edition dropped; partial-order planning; and the deep learning chapter's
recurrent networks, attention and anything that needs a GPU.

## Robotics and multiagent decisions, after Russell & Norvig

The rest of chapter 26 and all of chapter 18 of *Artificial Intelligence:
A Modern Approach* (4th ed.). Neither gives figures of pseudocode; the
algorithms are the ones their text teaches.

| Namespace                       | What it is                                                                 |
|---------------------------------|----------------------------------------------------------------------------|
| `allgo.robotics.motion`         | Configuration space for a disk robot among polygons, and a translating polygon's by Minkowski sum; the visibility graph, the Voronoi graph, cell decomposition and its refinement by subdivision, PRM, RRT, bidirectional RRT and RRT*, hybrid A*, trajectory optimization, shortcutting, potential fields; a two-link arm's configuration space |
| `allgo.robotics.localization`   | Monte Carlo localization, its low-variance resampling, and now EKF localization among landmarks |
| `allgo.robotics.slam`           | EKF SLAM: the pose and every landmark in one Gaussian                       |
| `allgo.robotics.fastslam`       | FastSLAM: a particle per path, a little Kalman filter per landmark in each  |
| `allgo.robotics.graph-slam`     | Graph SLAM: the whole run's poses and landmarks by Gauss-Newton least squares |
| `allgo.robotics.control`        | PID as a step; LQR by the Riccati recursion; iterative LQR; model predictive control |
| `allgo.robotics.tracking`       | Retiming by a trapezoidal profile; a two-link arm's dynamics; feedforward, PD and computed torque control |
| `allgo.robotics.guarded`        | Guarded moves under motion uncertainty: a peg into a hole                  |
| `allgo.robotics.humans`         | Predicting a person's goal from noisily rational steps; behavioral cloning and DAgger |
| `allgo.robotics.terrain`        | Self-supervised learning of drivable ground, a laser teaching a camera     |
| `allgo.robotics.reactive`       | A hexapod's leg machine and tripod gait; subsumption                       |
| `allgo.robotics.kinematics`     | Craig's chains and the IK solvers, moved here from `allgo.kinematics`      |
| `allgo.decision.contract-net`   | The contract net: announcement, marginal-cost bids, award, subcontracting; recontracting by Sandholm's O-, C-, S-, M- and OCSM-contracts; routes by Held and Karp |
| `allgo.decision.games`          | Normal-form games, dominance and its iterated elimination, Nash equilibria pure and mixed, zero-sum games by linear programming, repeated games, bargaining, the Shapley value, voting, auctions and VCG |

**Reused, not rewritten.** The arm is an `allgo.robotics.kinematics`
chain -- Craig's parameters, two revolute joints -- rather than its own
trigonometry. The 2D segment-to-segment distance went into
`allgo.geometry.polygon` beside the segment helpers already there; the
3D ones in `physics.contact` return points, and GJK works on general
convex shapes. The filters use the scalar Joseph-form update that orbit
determination had, now moved to `allgo.probability.temporal` with the
chapter 14 Kalman filter, so that robotics does not depend on astro.
`allgo.kinematics` moved under `allgo.robotics` for the same reason of
grouping: the arms, the planners and the filters are one package.

**Steps again.** PRM grows a sample a step; RRT and RRT* a node; the
cell planner is `allgo.search`'s A*; the potential field a move down the
gradient; shortcutting a try; EKF localization and SLAM a control and
its observations; PID a control; the Riccati recursion a step further
from the horizon; iterated elimination a strategy removed; instant
runoff a round; the English auction a raise; bargaining a round of
backward induction.

**Checked.**

- The visibility graph's path round a square, for radii 0, 0.2 and 0.5,
  is the exact shortest, 2 sqrt((1 - r)^2 + (1 + r)^2) + 2 + 2r, to 1e-9.
  Every other planner's path is collision free and no shorter; RRT*'s
  shortens as it grows (5.69, 5.44, 5.34, 5.30 m against 5.284 at 200 to
  2000 nodes).
- A passage narrower than a cell is lost to cell decomposition and found
  by the visibility graph -- the incompleteness the book warns of. The
  potential field reaches the goal past a post and stalls in a cup.
- The arm's configuration space agrees with a direct collision test at
  every cell, and its plans never collide.
- EKF localization stays within 15 cm over a lap, against dead
  reckoning's metre or more, with a normalized position error averaging
  near 2 -- the filter's covariance is honest. EKF SLAM maps all twelve
  landmarks to within 0.3 m and ends within 10 cm.
- PID on a mass: P alone oscillates for ever between 0 and twice the
  target; PD settles; a constant push of d leaves PD d / kp off; the
  integral removes it. LQR for x' = x + u with unit costs gives P = the
  golden ratio and K its reciprocal.
- Two-finger Morra's value is -1/12, the even player showing one finger
  7/12 of the time, from the linear program and from support
  enumeration; battle of the sexes has its two pure equilibria and the
  mixed (2/3, 1/3); the glove game's Shapley values are 2/3, 1/6, 1/6;
  bargaining converges to Rubinstein's split; truthful bidding is never
  beaten in a second-price auction.

**The Joseph update made quadratic.** The scalar update formed I - K H
as a matrix and multiplied it out, cubic in the state's size: fine for an
orbit's six, 24 s for SLAM's 27 over a lap. For one row H the Joseph form
is outer products -- A = P - K (P H)^T, then A - (A H) K^T + R K K^T --
quadratic, the same to 4e-15, and symmetrized; SLAM's lap now takes 5 s.

**The demos.** Motion Planning runs each planner in one room, the arm
beside its configuration space; click to move the goal. SLAM drives laps
of a field of landmarks, steering by its own estimate -- so under dead
reckoning the estimate traces a perfect rectangle while the robot
wanders off it, and under the filters it keeps to the route. Both demos
and Robot Localization are in the gallery's Robotics group, with the arms.

**The rest of chapter 26, later.** Three more of the chapter's methods,
each reusing what was here:

- *The Voronoi graph* -- the points equidistant from two or more
  obstacles -- is made from `allgo.geometry.delaunay`: points every `h`
  along each obstacle's edges and each wall, triangulated, and of the
  dual Voronoi diagram's edges those that part points of different
  obstacles, where the robot fits. Points of one obstacle part only from
  each other, so dropping their edges leaves the graph, to within about
  `h`. The plan is the book's: a straight move onto it, along it by
  A*, and off. Round the test square in its 6 m room every node is
  equidistant from its two nearest obstacles to within `h`, and the path
  keeps a metre from the square where the shortest grazes it at 0.2 m --
  and is 7.6 m to its 5.3.
- *FastSLAM*, the Rao-Blackwellized particle filter the chapter's notes
  credit: given the path, the landmarks are independent, so each
  particle is a pose and a 2 by 2 Kalman filter per landmark, weighted by
  the likelihood of its sightings and resampled -- by MCL's low-variance
  resampling, pulled out of `mcl-step` into `localization/resample` for
  both. On EKF SLAM's test lap, with 50 particles, all twelve landmarks
  within 0.17 m and the pose within 0.18 m all the way, in a quarter of
  a second.
- *Graph SLAM*: poses kept every ten steps and the landmarks as one
  least-squares problem. Odometry constraints are the controls between
  two kept poses composed by `localization/ekf-predict` from a pose known
  exactly, which gives the relative move and its covariance at once;
  sightings are ranges and bearings with `range-bearing-jacobian`'s
  partials. Gauss-Newton, an iteration a step, solved by
  `numerics.linear`'s Cholesky with the first pose held. On the same lap
  it converges in seven iterations, the weighted squared error falling
  from 13718 to 288 -- about the measurements' degrees of freedom, as
  it should for a consistent model -- and the path error from dead
  reckoning's 0.81 m to 0.18, the start of the run included.

The demos have them: Motion Planning's Voronoi graph, and SLAM's
FastSLAM, its particles and the best one's map, and graph SLAM, which
records a lap and then solves it an iteration at a time, dead
reckoning's drifted path snapping onto the true one. Graph SLAM's lap
is steered by the truth, since it estimates nothing until the lap is
over.

**The rest of chapter 26, read against the 4th edition.**

- *A polygon's configuration space* (`motion/cspace-obstacle`): for a
  robot that translates without turning, the obstacle's Minkowski sum
  with the robot turned half round, by `geometry.polygon`'s new
  `minkowski-sum` (edges merged in order of direction). Section
  26.5.1's right triangle against a rectangle gives a five-sided
  C-space obstacle, as the section says: the triangle's legs run
  parallel to two of the rectangle's sides, and parallel edges join.
- *Subdivision* (`motion/subdivision-seq`), the section's way to make
  cell decomposition complete: cells free, blocked or mixed by the
  signed clearance at their centers against half their diagonal, and
  the mixed ones split in four a step. The passage that metre-free cells
  miss (the existing test's 0.5 m gap for a 0.1 m robot) is found after
  four splits.
- *The bidirectional RRT* of the section's figure: two trees, from each
  end, joined when a milestone is reachable from both.
- *Hybrid A** (`motion/hybrid-a-star-seq`): A* over cells and heading
  bins, but each cell keeps the continuous pose that reached it, and a
  move drives a car's exact arc -- hard left, straight, hard right --
  from that pose. The heading changes by exactly step / radius between
  poses, never a corner. Forward only, as the book describes it: a car
  started a metre from the square, facing it, with a one-metre turning
  radius, has no path, and none is found. Hybrid A* as usually built
  also reverses.
- *Trajectory optimization* (`motion/trajectory-optimization-seq`),
  after CHOMP: J = J_obs + lambda J_eff, the cost field from the signed
  distance integrated by arc length, and gradient descent from the
  straight line. Through the two posts of the figure, the line comes
  clear by 0.4-0.7 m; with nothing in the way, a wiggly path straightens
  onto its chord, the Euler-Lagrange optimum. The gradient was at first
  central differences of the whole J, 33 ms a step; a waypoint is in
  only its two legs' terms, so differencing those is the same gradient
  at a 25th of the cost. Straight through the middle of a large
  obstacle it stays stuck -- halves of the path pushed out opposite
  sides, a leg left across it -- the local minimum the section warns of.
- *Tracking control* (`allgo.robotics.tracking`): retiming by the
  trapezoidal profile (D / v + v / a, or 2 sqrt(D / a) when too short),
  a planar two-link arm's dynamics with point masses (as in Craig's
  *Introduction to Robotics*; energy conserved to 1e-6 under no torque,
  inverse undoing forward to 1e-9), and the three control laws. Moving
  between two poses in 2.3 s, with a model a fifth too light: open
  loop is 2.5 rad off, PD alone 0.38, computed torque 0.04. Even with
  the exact model, open loop drifts 0.16 rad, its errors unchecked;
  computed torque with it, 0.0002.
- *iLQR and MPC* (`allgo.robotics.control`): iterative LQR with the
  dynamics linearized by differences, the backward Riccati pass for a
  feedforward and a feedback gain at each step, and a halving line
  search. On linear dynamics a second iteration changes the cost by
  under 1e-9 of it. A pendulum swung up from hanging, 60 steps, ends
  upright to 3e-5 rad. Run open loop on one 30% heavier, the plan ends
  at 0.8 rad; MPC with the same wrong model, replanning 40 steps ahead
  every step, holds it upright to 5e-4.
- *Guarded moves* (`allgo.robotics.guarded`), Figures 26.23-26.25: a
  peg's motion errs up to 20 degrees. Straight down, stopping at
  contact, it entered the hole in 57 of 200 trials; down and left to the
  floor, then down and right sliding along it, it reached the bottom in
  all 200.
- *Predicting people* (`allgo.robotics.humans`), equations 26.8-26.9:
  a person's steps scored by a softmax of -beta Q, Q the step plus the
  straight-line distance left, and the belief over goals updated by
  each. Walking right from the middle of a room, the belief in the
  window rises 0.33, 0.67, 0.89, 0.97 ... to 1.0, the kitchen behind
  falling to nothing.
- *DAgger*: lane keeping, a nearest-neighbor learner, two
  demonstrations. Without noise the cloned policy retraces the
  demonstrations exactly and DAgger has nothing to add, so the car is
  buffeted: behavioral cloning's mean offset is then 1.1-3.4 times the
  expert's on four seeds, and four rounds of DAgger bring it within 1.3.
- *Self-supervised terrain* (`allgo.robotics.terrain`), section
  26.4.3's adaptive vision: the laser labels a patch just ahead,
  `learning.unsupervised`'s EM fits a mixture over those colors, and the
  camera classifies everything by it. Taught on a road, 99% of road and
  no grass is drivable; after the robot turns onto grass, 98% of grass
  too; rocks never. The laser first swept whole rows, roadsides
  included, and learned grass at once; it sees a patch, not a row.
- *Reactive control* (`allgo.robotics.reactive`), section 26.9: the
  leg's machine lifts, swings, sets down and stands, and a blocked
  swing retracts and lifts 0.05 higher. A 0.27 rock costs four retries.
  The tripod gait waits for all three legs of a group before the other
  lifts, so three always stand. The first version counted a group done
  only when its legs landed on the same tick, and a leg finished early
  stepped again while a stuck one retried; each group now remembers
  which legs have stepped. Subsumption is the highest layer with
  something to say.

**The demos.** Motion Planning gains subdivision, the bidirectional
RRT, hybrid A* and trajectory optimization. A new Robot Control demo
has the arm's three controllers side by side, the pendulums open loop
and under MPC, the peg's trials, the hexapod, a person's goal being
guessed, and the laser teaching the camera.

**Left out:** coastal navigation and information gain in the cost;
domain randomization; the game-theoretic coordination of Figure 26.28;
learning cost functions from demonstrations, which
`decision.apprenticeship`'s inverse RL already does for MDPs.

**The contract net, later.** `allgo.decision.contract-net` allocates
tasks among cooperating agents the way chapter 18 describes. It follows
Smith's protocol, with Sandholm's marginal-cost bids. A manager
announces each task; every agent replies with a bid, the length the
task would add to its route at the cheapest place in it, or declines
because it lacks the skill or is full; once every reply is in, the
cheapest bid wins and the other bidders are told they lost. A task with
parts makes its contractor the manager of those, so the net is
recursive. A message is a step. Then recontracting, after the contract
types of Sandholm's "Contract Types for Satisficing Task Allocation: I
Theoretical Results" (AAAI Spring Symposium, 1998): each agent first
replans its route, the shortest through its tasks (exact, by Held and
Karp's dynamic program, up to a dozen stops), so that what its work
costs it depends on which tasks it holds and not on the order they came
in. Then agents make offers: one task (an O-contract, the contract
net's own), a task in exchange for one of the taker's (an S-contract, a
swap), two or more tasks together (a C-contract, a cluster), a task
handed on as part of a contract among three or more agents, each
handing on at most one (an M-contract), or a set of tasks, any of
anyone's, put up to be reallocated among all (an OCSM-contract, which
can be any of the others and any mix of them). Of the contracts an
offer admits, the one that lowers the summed cost of those party to it
most is made, if any does -- with the gain shared out by payments, it
is individually rational for each -- and they replan. Passes of every
offer repeat until one moves nothing. Total cost only falls. A step
weighs at most a few hundred contracts, so an M- or OCSM-offer, which
admits thousands, takes several steps; one-task, swap and cluster
contracts are the default, M and OCSM asked for.

The 4th edition's text of the chapter was not at hand -- the 3rd does
not mention the contract net -- so the protocol is Smith's, as he
gives it.

Checked: every task placed with an agent able to do it, and a
subtask announced by its parent's contractor; per task, an
announcement and a reply from every agent, one award and a rejection
for every other bidder; the award to the cheapest bid, and the winning
bids summing to the total route length exactly. An agent at capacity
declines, and a task no one can take is left unplaced. Against the best
allocation found by trying every one on three agents and six tasks,
neither the contract net nor recontracting is ever better, and
recontracting is never worse, with any of the types. With one-task,
swap and cluster contracts it ends where no task, cluster or swap lowers
the cost, checked by trying them all. On twelve instances:

| Contracts                     | Optimal on |
|-------------------------------|-----------:|
| one-task                      |    6 of 12 |
| M alone                       |    5 of 12 |
| one-task, swap, cluster       |   10 of 12 |
| ... and M                     |   10 of 12 |
| OCSM, up to two tasks         |    8 of 12 |
| OCSM, up to three tasks       |   11 of 12 |
| OCSM, any number              |   12 of 12 |

On four, one-task contracts stop at least half a meter short and with
swaps and clusters reach the optimum exactly. Small cases need each
type: a swap between two full agents, neither able to take the other's
task; a cluster of two tasks neither of which is worth moving alone;
an M-contract among three full agents with a ring of skills -- each
task doable by two of them -- where no one can take a task or trade one
but the three can pass them round. One instance is a local optimum of
one-task, swap, cluster and M-contracts together: the best allocation
needs b to hand three tasks to c as c hands one to a, which is none of
them, and an OCSM-contract makes exactly that. That is Sandholm's
result: each of O, C, S and M alone can stop short of the best
allocation, and OCSM-contracts reach it from any start. They do
because one may be the whole problem -- with every task put up, its
contracts are every allocation, n^k of them for k tasks among n agents
-- so they must be bounded past a handful, and bounded they can stop
short too. Held and Karp's routes are checked against every order of
up to seven stops.

The demo, Contract Net, has a dispatcher and four vehicles, some with
cranes, placing a day's jobs a message at a time; the message flies
between sender and receiver, labeled with its type and bid. Then
recontracting: the vehicles replan, and each contract made is shown
with what moved and what it saved, then a new day. A Contracts control
chooses one-task contracts; one-task, swap and cluster; those and M;
or OCSM up to three jobs. On the first day the first takes total route
length from 27.24 m to 23.90 m and the others to 22.99 m. The offers
are weighed a few milliseconds a frame: with twenty jobs an
M-contract's offer admits some fourteen thousand contracts.

**Not implemented here:** kinodynamic planning and trajectory
optimization; data association for landmarks that cannot be told apart;
Lemke-Howson, for games larger than support enumeration can handle.

## Computer vision, after Russell & Norvig

Chapter 25 of *Artificial Intelligence: A Modern Approach* (4th ed.), in
a package of its own. The chapter gives no figures of pseudocode; these
are the algorithms its text teaches.

| Namespace                     | What it is                                                                 |
|-------------------------------|----------------------------------------------------------------------------|
| `allgo.vision.image`          | Images as values, separable convolution, Gaussian smoothing, gradients, bilinear sampling; scenes drawn from shapes |
| `allgo.vision.camera`         | Perspective and scaled orthographic projection, Lambertian shading, triangulation from two views, stereo depth |
| `allgo.vision.edges`          | Canny's detector a stage at a time; orientation histograms for texture     |
| `allgo.vision.flow`           | Optical flow by block matching and by Lucas-Kanade                         |
| `allgo.vision.segmentation`   | Felzenszwalb-Huttenlocher graph-based segmentation, an edge a step; k-means on color |
| `allgo.vision.stereo`         | Disparity by matching along rows, a row a step; depth from disparity       |
| `allgo.vision.detection`      | A convolutional classifier with augmented data; sliding windows, IoU, non-maximum suppression; hard-negative mining |
| `allgo.vision.shape`          | Shape from shading (Ikeuchi and Horn's relaxation, the occluding contour fixing it), texture (each texel's slant and tilt; a plane from their sizes) and contour (the ground plane; T-junctions and the depth order they give) |
| `allgo.vision.structure`      | Structure from motion: Tomasi and Kanade's factorization, a stage a step |
| `allgo.vision.proposals`      | Selective search, a merge a step; anchors scored for objectness; R-CNN's box regression, an example a step |
| `allgo.vision.boundaries`     | Pb: half-disk histograms of brightness and texture compared, a boundary's probability learned |

**Reused, not rewritten.** Segmentation keeps its regions in
`allgo.graph`'s disjoint sets, joining smaller into larger so the trees
stay shallow; color clustering is `allgo.learning.unsupervised`'s
k-means, with restarts; the classifier is `allgo.learning.neural`'s conv
and dense layers; subpixel sampling is `allgo.numerics.splines`'
bilinear interpolation; cameras turn by `allgo.geometry.quaternion` and
the triangulation's two unknowns are a small `allgo.numerics.linear`
solve. The image code is new: nothing here had images, and `dla`'s blur
is private to its grid.

**Images are values.** `{:w :h :px}`, the pixels a vector of numbers or
of `[r g b]`; every operation returns a new image, and a color image
filters a channel at a time. At the sizes a demo or a test uses -- a
hundred pixels across -- that costs nothing worth trading clarity for.

**Steps again.** Canny is a stage a step and its hysteresis a ring of
growth; graph-based segmentation an edge; k-means an iteration;
Lucas-Kanade a warp and correction; disparity a row; non-maximum
suppression a window; the classifier an epoch, or in the demo a few
examples a frame (`neural/sgd-step`).

**Checked against drawn scenes, whose answers are known.**

- Canny on a noisy square: 84 edge pixels against a perimeter of 88, every
  one within half a pixel of the true boundary, one pixel thick; noise
  alone gives none.
- Projection: a unit up at depth 5 lands f / 5 above the center, twice as
  far half as high; parallel lines converge on the vanishing point;
  triangulation recovers a point from its two images to 1e-9.
- Flow: block matching finds a (2, 1) move exactly at every grid point;
  Lucas-Kanade finds it to 1e-6, and a subpixel (0.6, -0.3) to a few
  hundredths; a featureless window gives nothing -- the aperture problem.
- Segmentation of three noisy colored regions: graph-based 97.5% pure in
  at most six regions, k-means with restarts exactly the three.
- Stereo: the exact disparity at over nine pixels in ten, the misses at
  the occluding box's borders.
- Detection: every disk found among squares in new scenes, and the false
  alarms down to none after mining.

**Two remedies the chapter's methods need.** k-means from random centers
can put two in the background and merge the shapes: 85% pure on the test
scene, and restarts, keeping the run of least cost, make it 100%. A
detector trained on patches of whole objects fires on the clutter a
sliding window shows it -- squares, and parts of things -- so its examples
come from whole scenes (the objects' windows and random windows that
barely touch them), and then from its own false alarms: one round of
mining took three false alarms a scene to one in eight scenes, a second
to none.

**The conv layer made faster.** Training the detector spent its time in
`learning.neural`'s convolution -- nested `for` and `get-in` forward, an
`update` per input element back. Loops and transients there took a
training example from 3.8 ms to 1.4 on the JVM, the same network, with
the gradients still matching numerical ones; the classifiers demo gains
too.

**The demo.** Computer Vision runs each method on a drawn scene: Canny's
stages, regions merging, k-means painting, flow arrows over a moving
disk, a disparity map filling in a row at a time, and a network learning
to find disks a few examples a frame, mining its false alarms every six
epochs. The newer modes:
- *Shading:* an ellipsoid's normals found, relit from the left beside
  the truth.
- *Texture:* the texels' ellipses and the normals they imply.
- *Contour:* overlapping shapes, their junctions marked right or wrong,
  and the depth layers they give.
- *Structure from motion:* a house's points tracked through thirty
  turning views, then the shape found turning against the truth.
- *Selective search:* regions merging and proposals accumulating.
- *Anchors and regression:* objectness learned an epoch a frame, the
  box regression an example at a time, then the best proposals moved
  onto their objects. Canny works on brightness, and two of the scene's colors are
equally bright, so no edge shows where they meet -- the caption says so.

**Shape from one image.** `allgo.vision.shape` takes the chapter's
cues in turn. The 3rd edition's text, which is what was at hand, gives
them as cues rather than procedures, so each follows the paper or the
geometry it rests on.

- *Shading:* Ikeuchi and Horn's relaxation ("Numerical Shape from
  Shading and Occluding Boundaries", Artificial Intelligence 17, 1981),
  in stereographic coordinates, which stay finite at the occluding
  contour. Brightness is one equation for two unknowns, so smoothness
  is assumed: each iteration moves every normal to the mean of its
  neighbors' and then along the brightness error. The silhouette fixes
  the normals at the contour, square to the view and pointing out, and
  so decides convex over concave. Letting the contour relax too lost
  that, and failed. On a drawn sphere the mean error falls from 41
  degrees (every normal facing the viewer) to under 5, and stays there:
  at about 4.5 the smoothness and the brightness balance, whatever the
  resolution. Weighing brightness more makes it worse -- 7.9 degrees at 4 --
  and at 8 it diverges.
- *Texture:* each round texel's ellipse is measured as a camera turned
  to look straight at it would see it. Off the axis the image plane is
  itself oblique to the line of sight, and measured in the image, the
  slant came out up to 7 degrees wrong. The minor-to-major ratio gives
  the slant to within 5 degrees of the truth, and the tilt to within 10.
  The texels' sizes, which go as 1 / distance, place them along their
  rays, and a plane through them gives the ground's normal to under 5
  degrees in the test (0.1 in the demo), needing no slant at all.
  Texels under 40 pixels are left out: those of 15 were 30 degrees off.
- *Contour:* the ground plane's depth f h / (v - cy) and height are
  checked exactly. For T-junctions, of the three boundaries meeting at a
  junction, the two that run on straightest through it are the
  occluder's outline. The first rule, taking the region covering most of
  the junction's neighborhood as the occluder, failed where the outline
  curves: a small disk covers less than half. Even the straightness rule
  misreads a junction a pixel from a box's corner, or where a small disk
  bends sharply, so a junction where no pair runs clearly straight
  (beyond about 145 degrees, and clearly straighter than the next pair)
  says nothing. On forty random rows of five overlapping shapes, 96% of
  the region pairs the junctions decide come out right, and every pair
  in 32 of the 40. Local convexity, the chapter's other figure-ground
  cue, was tried and not kept: on these scenes it was right about 80% of
  the time, fooled by the thin slivers an occluder leaves of what it hides.
  The junctions give a partial order: `in-front` is the relation, and
  `depth-order` one ordering consistent with it.

**Structure from motion.** `allgo.vision.structure` is Tomasi and
Kanade's factorization ("Shape and Motion from Image Streams under
Orthography", IJCV 9, 1992). The tracks, less their means, have rank
three, so their three largest singular directions are the camera's axes
and the points, up to a linear map. Here the singular directions are
the eigenvectors of W W^T by `allgo.numerics.eigen`'s Jacobi. The map is
fixed by the axes being unit and square to each other: six unknowns,
linear, by least squares. `lin/inverse` is for symmetric matrices only,
so Q^-1 comes from Q's own eigenvectors. Without noise, every distance
between points and the camera's turn between frames come back to 1e-9;
with 0.01 of noise, the distances are within 2% on average. A mirror
image in depth is left unknown, as under orthography it must be.

**Region proposals.** `allgo.vision.proposals` has three pieces:
- *Selective search* (Uijlings et al., IJCV 2013), merging a step from
  a graph-based segmentation, by color, size and fill; the paper's
  texture measure is left out.
- *Anchors* scored for objectness, as the region proposal network of
  Faster R-CNN (Ren et al., NeurIPS 2015) does, but by logistic
  regression on the anchor's pixels rather than shared convolutional
  features. They are sampled as that paper samples them, a fixed number
  an image at most half positive, which also spares computing features
  for every anchor.
- *R-CNN's box regression* (Girshick et al., CVPR 2014): ridge
  regressions on the proposal's pixels, with context, the normal
  equations built an example a step.

The test scenes' objects are two-toned, lighter above and darker below,
so that a segmentation splits each. The starting regions' boxes find
none of thirty objects to 0.7 IoU, and selective search finds all
thirty with about ten proposals a scene. Box regression takes jittered
boxes from 0.62 mean IoU to 0.84 and makes none worse. The anchors keep
a few dozen of some fourteen hundred and find 27 of 30 objects to 0.5
but only 15 to 0.7, their sizes being fixed; box regression makes that
29.

### The rest of chapter 25

Read now against the 4th edition's own text, the chapter left these.

**Light and color** (`allgo.vision.camera`). `shade` is section
25.2.4's model of nearly every surface, diffuse with specularities:
Lambert's term plus the ambient, and Phong's lobe for the specularity,
a point in shadow getting only the ambient. Albedos and lights may be
[r g b], Lambert's law then a channel at a time, as section 25.2.5
models color. On a lit sphere the brightest pixel is within a pixel of
where the normal bisects the light and the view. Color constancy, which
the chapter names but gives no method for, is the gray-world estimate
(Buchsbaum, 1980): the image's mean color taken for the light's. Where
the surfaces really do average to gray it recovers the light's color
and undoes it to 1e-12. With forty random albedos, which only roughly
do, every channel came within 0.035 of the scene under white light.

**Stereopsis.** The chapter's figures for human eyes, a 6 cm baseline
and 5 seconds of arc, give the depths it quotes: 0.4 mm at a meter and
0.036 mm at 30 cm (`depth-resolution`, a test).

**A moving camera** (`allgo.vision.flow`). The chapter's flow equations
for a translating camera; the focus of expansion as the point nearest
every flow line, by least squares; the time to contact as a point's
distance from the focus over its speed; and motion parallax's depth
ratio. On exact flows the focus comes back to 1e-12 and every point's
time to contact is Z / Tz. Measured, by Lucas-Kanade between two views
of a textured wall 0.5 nearer, the focus lands within half a pixel of
the center and the times average 19.1 frames against 19.

**The pose of a known object** (`camera/pose`), which section 25.6.4
calls straightforward. Scaled orthographic projection gives a linear
start, in two mirror images in depth; each is polished by Marquardt's
method (`allgo.numerics.optimize`) on the reprojection error, and the
better kept. Plain Gauss-Newton (`numerics.regression`) was tried
first: from the mirrored start some points went behind the camera and
its system turned singular. From exact images of a house's ten corners
the camera comes back to 1e-11 from three viewpoints, and with half a
pixel of noise to within 0.15, eight to ten units away.

**The horizon from pedestrians** (`shape/horizon`). People all one
height on the ground have heights in the image proportional to how far
below the horizon their feet are: a line through zero at the horizon.
The book's two ways round, the detections placing the horizon and the
horizon ruling out impossible detections, are one fit. Least squares
alone failed: one false "giant" near the horizon dragged the line more
than 50 rows off and the trimming then threw out the true people. So
the start is the consensus line, the one through the pair of detections
that most others agree with, then refit to its inliers until they
settle. With half a pixel of noise it puts the horizon within a pixel
and drops exactly the giant and the dwarf.

**Edges linked into curves** (`edges/link`): neighboring edge pixels
whose orientations agree. Within 29 degrees each object's outline is
one curve; within 17, the square's corners break it into its four sides.

**Learned boundaries** (`allgo.vision.boundaries`), section 25.3.4's
Pb(x, y, theta): at each pixel a disk is cut along a diameter at eight
orientations, and the two halves' histograms of brightness and of
gradient orientation are compared by chi-squared. Logistic regression
(`learning.models`) learns from scenes whose boundaries are known
(drawn, their labels known) how those differences mean a boundary, and
Canny's non-maximum suppression across the best orientation thins the
result. The scenes are regions told apart by texture -- stripes of two
directions, dots -- at the same mean brightness. On six new scenes,
59-81% of Canny's edges lay within two pixels of a boundary, the rest
on stripes, and it found 60-95% of the boundary. Of Pb's, 89-95% lay on
a boundary, and it found at least 97% of it. Unthinned, Pb's band was
2.6 times as many pixels as the boundary, at the same precision.

**Normalized cuts** (`segmentation/normalized-cuts`, after Shi and
Malik, PAMI 2000). The second generalized eigenvector of the graph, then
the split along it with the least Ncut, recursively. On pixels, by power
iteration, the second eigenvector needed up to 1,500 multiplications to
settle on some scenes (5 s a segmentation), and the plain background
split cheaply through its narrow parts. The chapter's own remedy for
graph methods, an over-segmentation into superpixels, fixed both: the
graph-based segmentation at a small scale leaves a few dozen
superpixels, joined by their shared boundary's length times their
similarity, and on so small a graph Jacobi gives the eigenvector
exactly. On six random scenes of two shapes, three regions each time,
each all one of them. Jacobi is cubic: 40 superpixels take 0.35 s, 200
took 54 s, so the demo's scenes are kept small.

**Detection, evaluated** (`allgo.vision.detection`): precision and
recall with detections matched best first, an object claimed once;
average precision as the area under the curve with precision taken at
its best at that recall or beyond. A hand-worked case gives 34/45 to
1e-12. The mined detector's average precision on the test scenes is
1.0. ROI pooling reduces a box of any shape to a grid of maxima.
Test-time augmentation averages the scores of a window's shifted and
mirrored copies. On these scenes it classified every test window right,
as the plain score already did, so it showed no gain here.

**The demo** gains five modes: *Boundaries* (Canny's edges and Pb's,
each pixel green if on a boundary); *Normalized cuts* (the superpixels,
then the regions a cut at a time); *Time to contact* (the flow while
flying at a wall, the focus and the estimate against the truth);
*Horizon* (detections, the horizon found, the impossible ones in red);
*Light and color* (a shiny sphere and patches under orange light, the
gray-world correction beside the scene under white light, the light
turning).

**Not implemented here:** image generation and style transfer, paired
and unpaired image translation, captioning and visual question
answering, and depth maps predicted from one image -- all networks
trained on far more data than a drawn scene; the deep networks of the
chapter's classification, for the same reason; depth from shading
integrated from the normals; shape from texture on curved surfaces; the
paper's texture measure in selective search, and a convolutional region
proposal network sharing the detector's features.

## Ethics and safety, after Russell & Norvig

Chapter 27 of *Artificial Intelligence: A Modern Approach* (4th ed.) is
mostly argument -- the limits of AI, whether machines can think, the
ethics of weapons, work and robot rights -- and that is left to the book.
These are the parts of it that are methods, in a package of their own.

| Namespace                 | What it is                                                                 |
|---------------------------|----------------------------------------------------------------------------|
| `allgo.ethics.privacy`    | k-anonymity by generalization, a step a field, with Datafly's suppression; aggregate queries and the differencing attack; the Laplace mechanism; federated averaging; secure aggregation's cancelling masks |
| `allgo.ethics.fairness`   | Selection, true and false positive rates by group; calibration; calibrated scores whose base rates differ; equal-opportunity thresholds; equal impact's utility threshold; proxy accuracy; SMOTE |
| `allgo.ethics.safety`     | Fault trees: the top event's probability and the minimal cut sets; FMEA's risk priority; a low-impact coffee fetcher on `decision.mdp`'s value iteration |

**Checked.**

- The section's differencing attack: averages of $81,199 over 13 and
  $81,234 over 12 give the one 41-year-old's salary, $80,779.
- Of 300 random records, over 250 are alone in their age, sex and zip
  code -- the point of Sweeney's 87%. Generalizing to 2-anonymity, ages
  to decades and zips to tens leave a few records still alone, and the
  generalization then deletes age outright; Sweeney's Datafly withholds
  such records instead, and with up to three withheld, the decades stay.
- The Laplace mechanism's log density moves by at most epsilon when one
  person joins, at every answer checked, and its counts are right on
  average.
- Five users each seeing a corner of the plane: federated averaging of
  their logistic regressions classifies new points over 95% right, no
  data leaving a user; one user alone does worse. Secure aggregation's
  uploads are each far from what they hide and sum exactly.
- Calibrated scores for groups with base rates 0.3 and 0.6: one
  threshold gives false positive rates 7% and 56%; thresholds per group
  that equalize the true positive rates at 0.8 leave the false positive
  rates 51% and 59% -- Kleinberg et al.'s impossibility, in the numbers.
  (The section's COMPAS figures are the real system's; these data are
  drawn.)
- Zip code predicts the deleted group 89% of the time; income alone,
  chance.
- A minority with its own rule and a tenth of the data: one logistic
  regression gets the majority 97% right and the minority 56%. SMOTE
  brings the minority to 73% and the majority down to 80% -- one linear
  model cannot fit both, and oversampling chooses how the error is
  shared rather than removing it.
- The fault tree's probability matches the product worked by hand to
  1e-15, and a cause shared by both branches of an AND is found to be a
  cut set on its own.
- The coffee fetcher walks through the vase in its way on the designer's
  reward; charged for each change it makes to the world, it goes round,
  two steps longer.

**Left out:** everything the chapter argues rather than computes; the
de-biasing of data by modeling its sources of bias (Jiang and Nachum,
2019); assistance games, which the book introduces in chapter 18 and this
library does not yet have.

## Natural language, after Russell & Norvig

The rest of chapter 23 of *Artificial Intelligence: A Modern Approach*
(4th ed.), beside the CYK parser it already had. The 3rd edition's text
was what was at hand -- its chapters 22 and 23 -- and its account of
n-grams, smoothing, perplexity, <UNK> and learning grammars is what the
code follows. The data are drawn from grammars and from
`allgo.procedural.naming`'s languages, so the right answers are known.

- **N-gram models** (`allgo.language.ngrams`): counts by order, the
  sequence padded with start markers and an end. Laplace's add-one, or
  linear interpolation of every order, with a uniform estimate below
  the unigram so that no token has probability zero. The interpolation
  weights are fitted by EM on held-out text, an iteration a step; the
  held-out likelihood never falls. On sentences from a sixty-word
  grammar: perplexity 36.5 unigram, 13.0 bigram, 13.0 trigram -- the
  grammar's words depending mostly on the one before, the trigram adds
  nothing -- while Laplace's trigram, spreading its mass over trigrams
  never seen, is 23.3, worse than its own bigram, as the book warns.
- **Unknown words:** each word's first appearance counted as <UNK>,
  the 3rd edition's trick; the 4th replaces the infrequent words instead
  (`with-unknown` with a minimum count).
- **Language identification** by character trigrams, on four languages
  drawn by the naming code: 98.5% of single words.
- **Tagging** (`allgo.language.tagging`): an HMM counted from tagged
  sentences, add-one throughout, unknown words as <UNK>, decoded by
  chapter 14's Viterbi -- `allgo.probability.temporal`'s, unchanged.
  That HMM starts from P(X0) and moves a step before the first word, so
  the sentence's start is a state of its own. On the ambiguous grammar
  every word is tagged right, where each word's most frequent tag gets
  the ambiguous ones 57%; and Viterbi, weighing the whole sentence, tags
  a word that only the next word decides, which a recurrent tagger
  reading left to right cannot.
- **A grammar from a treebank** (`learn-grammar`), by counting: from
  three thousand sampled trees every probability is within 0.015 of the
  truth, and the learned grammar parses new sentences as the true one
  does.
- **The inside-outside algorithm** (`inside-outside-seq`): EM over
  unparsed sentences, inside probabilities from the words up, outside
  from the root down, an iteration a step. With the rules given and
  their probabilities uniform, it reaches the sample's maximum in an
  iteration or two, every probability within 0.03 of the truth and the
  likelihood a little above the true grammar's on that sample. With
  every rule allowed among the categories, as the book poses it, the
  likelihood climbs every iteration from random starts but stops at
  different local maxima: on a hundred sentences one stalls about fifty
  below the true grammar, another passes it -- the chapter's warning
  about local maxima, and its cost, O(n^3 m^3): 20 s for thirty
  iterations.

`generate` now samples a tree and reads its leaves, drawing the same
random numbers in the same order, so the data the other tests draw is
unchanged.

**The demo.** Two more modes of Deep Learning for Language:
- *N-grams:* each order's samples and perplexity as EM fits the
  trigram's weights, and the HMM tagging both ducks.
- *Inside-outside:* the likelihood by iteration and each rule's learned
  probability against its true one.

**Not implemented here:** Kneser-Ney and other backoff smoothing;
text classification, which `allgo.learning.models`' naive Bayes already
covers; dependency parsing, if the 4th edition's chapter teaches it --
its text was not at hand to check.

## Deep learning for language, after Russell & Norvig

Chapter 24 of *Artificial Intelligence: A Modern Approach* (4th ed.). It
gives no figures of pseudocode; these are the methods its text teaches,
on data drawn from grammars and templates so that every right answer is
known.

| Namespace                        | What it is                                                                 |
|----------------------------------|----------------------------------------------------------------------------|
| `allgo.learning.autodiff`        | Reverse-mode differentiation over matrices; Adam and clipping              |
| `allgo.language.embeddings`      | Skip-gram with negative sampling; GloVe; similarity, nearest words, analogies, a 2-D projection |
| `allgo.language.window`          | The feed-forward tagger over a window of words, learning its embeddings as it tags |
| `allgo.language.recurrent`       | Language models with an RNN or an LSTM cell; sampling, perplexity; stacked and bidirectional taggers; sentence classification by the last state or average pooling; a bidirectional language model's contextual embeddings, probed and fine-tuned |
| `allgo.language.seq2seq`         | An encoder-decoder, basic or with attention; greedy and beam decoding      |
| `allgo.language.transformer`     | Multi-head self-attention, sinusoidal or learned positions, layer norm; tagging, masked language modeling, fine-tuning; the decoder, masked and attending to the encoder, and the whole encoder-decoder |
| `allgo.language.parsing`         | Chapter 23's CYK and treebank code, moved from `allgo.language`, now also sampling sentences from a grammar |

**One differentiation, every architecture.** Back-propagation through
time, through attention and through a transformer are each a page of
calculus by hand. `allgo.learning.autodiff` does them all: every matrix
operation records how to send a gradient back to its inputs, and
`value-and-grad` walks the graph backwards. Each model here is only its
forward pass. Every operation's gradient is checked against central
differences, through a function that uses them all. Skip-gram alone
follows its gradient directly, being two lines and many thousand updates
an epoch.

**Reused, not rewritten.** Beam decoding is `allgo.search`'s beam
search: a state a prefix, a step's cost the negative log probability of
the token added, best-first with the frontier cut to the beam. The
embeddings' projection is `allgo.numerics.eigen`'s Jacobi method on their
covariance. The training text is sampled from probabilistic grammars in
`allgo.language.parsing`, whose CYK parser then judges what the language
models write. `allgo.language` became a package for these, its parser
moving to `.parsing`.

**Steps again.** Skip-gram a chunk of pairs or an epoch; each network a
sentence, or a sentence pair, a step of Adam; beam search a prefix
expanded.

**Checked.**

- Embeddings, on a corpus where each person word is a gender, a rank and
  an age said by context words of its own: king - man + woman = queen, and
  five more analogies of the kind, on every seed tried.
- A recurrent language model of grammar sentences: perplexity 4.6 on
  unseen ones against 13.9 for word frequencies alone, and most of what it
  writes parses under the grammar.
- The last word must match the first, six words back: an LSTM learns it,
  the plain RNN stays at chance -- the vanishing gradient the chapter
  describes.
- Attention, learning to reverse digit strings: every unseen string
  reversed, greedily and by beam search, and the attention the
  anti-diagonal.
- GloVe: all six analogies on half of ten seeds, at least four on every
  one.
- Taggers: one-way at chance, 49%, on a word whose tag the next word
  decides; two-way, with one layer or two, all of them.
- The transformer tagging a grammar whose words are often two parts of
  speech -- "the duck will duck a dog": 100% of held-out words, and of the
  ambiguous ones, where tagging each word with its most frequent tag gets
  57%.

**Pretraining, honestly.** A masked language model trained on unlabeled
sentences fills each hole with a word of the right part of speech 79% of
the time, so it learns the grammar's categories without labels. But on
this small grammar it does not make a better tagger: twenty labeled
sentences already tag the held-out ones perfectly from scratch, and with
fifteen that leave four words out, fine-tuning the pretrained encoder,
or training a linear layer on its frozen output, did no better on those
words than a randomly started one (a frozen random encoder probed for
long enough reached 84%, the pretrained one 68%). The chapter's case for
pretraining is data far larger and scarcer in labels than this;
`fine-tune` and `:train-only` are here to try it, and the test claims
only what was measured.

**Faster arithmetic.** The training loop was mostly matrix products, and
`numerics.linear`'s goes through lazy sequences: a loop over the rows in
`autodiff` took a sequence-to-sequence step from 8.7 ms to 5.9.

**The demo.** Deep Learning for Language trains each model in front of
you, a little each frame within a time budget: the embeddings on their
principal components, the male-to-female arrows turning parallel; a
language model's loss and its sentences, green where they parse; the
attention grid of a reversal model becoming the anti-diagonal; and a
transformer tagging both ducks, with one head's attention. GloVe draws
the same picture as skip-gram, fitted twenty table entries a frame.
Bidirectional trains a one-way and a two-way tagger by turns, the
held-out sentences checked a couple a frame: the one-way at about 50% on
the first word, the two-way at 100%.

**GloVe.** `glove-seq` fits the log of each pair's co-occurrence count
by weighted least squares with AdaGrad, as Pennington, Socher and
Manning's paper does (EMNLP 2014). On the analogy corpus it is less sure
than skip-gram, which solves all six analogies on every seed tried. With
the paper's choice of word vector, w + w~, it solved one to six. The
wrong answers were context words: "elder" for king - man + woman.
Every word in this corpus keeps to one role, a person never a context
word about a person, so the sum mixes in what each is as a context. The
word vectors w alone, window 4, dimension 16, solve at least four of the
six on every one of ten seeds, and all six on half; that is the
default, and `:vectors :sum` gives the paper's. Centering the vectors
changed nothing.

**Stacked and bidirectional networks.** A recurrent tagger can stack
layers, each reading the outputs of the one below, and make each layer
bidirectional: one network left to right, one right to left, their
states side by side. The test grammar's first word is a noun before an
intransitive verb ("duck swim today") and a verb before an object
("duck stones today"): only the word after decides. Reading left to
right, a tagger is at chance on it, 49%, with one layer or two; reading
both ways, it gets all of them.

**Contextual embeddings.** `bilm` is a bidirectional language model
after ELMo (Peters et al., NAACL 2018), one network predicting each word
from those before and one from those after, sharing the embeddings and
the output layer. A word's contextual embedding is its embedding beside
both networks' states there. ELMo learns a weighting of its layers per
task; here they are kept side by side. The same word gets a different
vector in each sentence: on the ambiguous grammar, "duck" the noun's
vectors are nearer one another than to "duck" the verb's by more than
0.2 in cosine.

**Pretraining where it pays.** The ambiguous grammar's rules, over a
vocabulary four times its size, 61 words. The bidirectional language
model reads 1,500 unlabeled sentences. With ten labeled sentences, a
third of the test's words were never labeled:

| Tagger                                           | All words | Words never labeled |
|--------------------------------------------------|----------:|--------------------:|
| From scratch                                     |       84% |                 65% |
| A softmax probe on the frozen contextual vectors |       92% |                 84% |
| Fine-tuned from the language model               |       91% |                 81% |

With a hundred labeled sentences, every test word has been seen, and the
network trained from scratch tags all of them (the probe 97%, the
fine-tuned one 99.9%). That is the chapter's case for pretraining, here
in miniature: it pays when the labels cover little of the vocabulary and
unlabeled text covers it all. The transformer's masked language model,
above, met the other case: twenty labeled sentences that cover a small
grammar already teach it all.

### The rest of chapter 24

**The window tagger** (`allgo.language.window`) follows section 24.1's
numbered procedure: a window of w words around the one to tag, a
vocabulary of the words seen more than a minimum count with the rest
`<unk>`, an embedding matrix learned with the tagger, the window's
embeddings concatenated so that each position meets its own part of
the first layer, then two sigmoid layers and a softmax. On the grammar
whose first word only the next one decides, a window of one word is at
chance, 49%, and a window of three tags all of them. On the 61-word
grammar with 30 labeled sentences it tags the test words never labeled
95-97% right from their windows alone. Raising the minimum count, so
that rare words teach the `<unk>` embedding, did not help there: 88-89%
at one, 84-86% at two. The window already carries what an unknown word
needs.

**Classifying sentences** (`rec/classifier`, section 24.2.2) pools a
recurrent layer's states into one vector: the last state, or their
average. Reviews whose verdict is the fourth word and the ten to twelve
after it say nothing show the bias the section warns of. The last
state's classifier, a plain RNN, was at chance on seven seeds of eight
(42-50%), and averaging got every review on all eight. With two to four
words after the verdict, the last state gets all of them too.

**The basic sequence-to-sequence model** (`s2s/model ... :attention
false`), the section's Figure 24.6, knows the source only by the
encoder's last state. Reversing digits, an LSTM so built reversed 23 of
30 strings of three to five digits (26 at 1,500 steps) and none of six
to eight, at 1,000 steps or 1,500: the fixed context size limit. Its
first two outputs, the source's last digits, were right on 20 and 18 of
30, the rest near 12 -- the nearby-context bias as well. With attention
a plain RNN reversed all 30 of the long ones.

**The transformer's decoder** (`tf/encoder-decoder`), which section
24.4.2 describes in words: each layer's self-attention masked so that a
position sees only those before it (`causal-mask`, -1e9 added to the
scores of later ones), then a second attention, to the encoder's
output, then the feed-forward network, each with its residual and
normalization. Trained teacher forced, it reverses all 30 unseen strings
of six to eight digits after 2,000 pairs, greedily and by beam search.
At the taggers' rate, 0.01, it never learned; at 0.005 it did, so that
is its default. Its attention to the source is not the clean
anti-diagonal of the recurrent model's: the encoder's self-attention can
move each digit's identity to other positions first, so the test claims
none. The beam search is now `s2s/beam-search`, shared by both models.

**Learned positions.** The chapter's transformer learns an embedding for
each position, not the sinusoids of Vaswani et al.; `encoder` takes
`:max-len` to do so. Both tag every held-out word of the ambiguous
grammar.

**Left out:** coreference resolution as a softmax over
antecedents, and a character-level model, which the chapter mentions in
passing -- the recurrent models take any tokens, characters among them.

**Not implemented here:** ELMo's learned weighting of layers and its
character convolutions; anything at the scale of BERT or GPT, whose point
is the scale.

## Knowledge representation, after Russell & Norvig

Chapter 10 of *Artificial Intelligence: A Modern Approach* (4th ed.):
ontologies, categories, events, knowledge of knowledge, and reasoning
by default. The copy at hand is the 3rd edition, where the same material
is chapter 12, and the code follows its text. The chapter gives no
pseudocode, and says outright that it leaves out the truth maintenance
algorithms; what is here is the usual algorithm for each idea it
teaches.

| Namespace                    | What it is                                                                  |
|------------------------------|-----------------------------------------------------------------------------|
| `allgo.logic.categories`     | Semantic networks: MemberOf and SubsetOf, boxed and double-boxed links, inheritance with overriding, reified events; Disjoint, ExhaustiveDecomposition, Partition |
| `allgo.logic.description`    | CLASSIC without SameAs: normal forms, structural subsumption, incoherence, classification into a taxonomy, membership of individuals |
| `allgo.logic.events`         | The event calculus's T, Happens, Initiates, Terminates, Clipped, Restored; processes; the interval relations; the time scale and Date |
| `allgo.logic.modal`          | Possible worlds and accessibility, K_a with nesting and quantifiers, checked by labeling; reflexive, transitive and symmetric frames |
| `allgo.logic.defaults`       | Circumscription, plain and prioritized, as preferred models; default logic's extensions |
| `allgo.logic.tms`            | A JTMS and an ATMS                                                          |

**Steps again.** Inheritance a node a step up the links; labeling a
subsentence a step; the event calculus a time point a step; models one
DPLL call a step, with a clause blocking each one found; default logic
a process -- an order of applying the rules -- a step, depth first;
classification a concept placed a step; a JTMS a node off its queue,
an ATMS a justification.

**Reused.** Default logic and circumscription ground their first-order
sentences over the constants they mention and hand them to
`allgo.logic.propositional`'s DPLL: entailment for a rule's
prerequisite, satisfiability for its justification, and every model for
circumscription. The sentences are parsed by `allgo.logic.first-order`,
whose substitution also grounds the modal quantifiers. PartOf's
transitivity -- PartOf(Bucharest, Earth) -- is just forward chaining
there, and is tested as such rather than given code of its own.

**The book's dates are Pacific time.** Its time scale counts seconds
from midnight GMT at the start of 1900, and it gives Time(Begin(AD2001))
= 3187324800, Time(End(AD2001)) = 3218860800 and Date(0, 20, 21, 24, 1,
1995) = 3000000000. None of these is right on Greenwich time; all are
on Pacific Standard Time, eight hours behind -- midnight in California
-- and the last only with Date's arguments read seconds, minutes, hours,
not the hours, minutes, seconds the text names: 3,000,000,000 seconds
after 1900 began in Greenwich is 21:20:00 on 24 January 1995 in
California. `events/date` takes a `:utc-offset` and the seconds-first
order, and the tests reproduce every number with an offset of -8. (Only
the 3rd edition was checked.)

**Choices the text leaves open.**

- Inheritance: a value is overridden by one on any node below it, not
  just by the first one found going up, so a near general category
  cannot beat a far specific one (Opus is a penguin, and a pet, and pets
  are animals). Values that neither overrides and that disagree -- the
  Nixon diamond -- are returned as ambiguous.
- The event calculus states only when a fluent holds and when it does
  not; a fluent nothing has initiated or terminated is neither, and
  `holds` says nil. Its Clipped counts an event at the very moment of
  initiation, so an initiation and a termination at once leave the
  fluent undecided too.
- Figure 12.4's arrows were taken from its description; the figure
  itself did not survive extraction. Every claim the text makes of it
  holds: Lois knows Superman knows whether he is Clark, Superman knows
  Lois knows the forecast, and in the joined model, in w0, Superman knows
  I but not R, nor whether Lois knows R.
- CLASSIC's individuals: membership reads an individual's role fillers as
  all it has, so All(Son, Married) can be decided from the sons named.
- The truth maintenance systems are monotonic: justifications have only
  in-lists. The JTMS retracts in two passes -- everything that depended
  on the retracted node out, then whatever another justification still
  supports back in -- which is what stops two nodes that justify each
  other from keeping each other in.

**Checked**, besides the dates: Mary's two legs and John's one (figure
12.5); Fly17 (figure 12.6); the Nixon diamond as two preferred models,
religion over politics making a pacifist, and as two extensions; Tweety
flying until told it is abnormal; the bachelor, and the men with three
sons, against what follows from them and what does not; the three JTMS
cases of Q's justifications; the car that won't start; K_a P => P only
on reflexive frames, K_a P => K_a K_a P failing on a frame that is not
transitive; Bond's two readings.

**The demo.** Knowledge Representation shows each, a step at a time:
inheritance climbing figure 12.5 and the Nixon diamond; figure 12.4's
worlds lit as each part of a knowledge sentence is labeled; a wumpus
hunt's fluents through time; the Nixon diamond's processes, extensions
and preferred models; a taxonomy filling in; and the 2048 Olympic sites
under a JTMS -- the chapter's example, the consequences ours --
Athletics moving from Bucharest to Sibiu, with the car's ATMS labels
beside it.

**Not implemented here:** the ontology of measures, and of stuff and
things, which the chapter gives as axioms rather than procedures;
CLASSIC's SameAs; nonmonotonic justifications (out-lists) in the JTMS;
and the 3rd edition's internet shopping agent, which the 4th dropped.

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

Collision detection then rose to the top of what was left, for the
first time honestly rather than by default.

XPBD's velocity pass followed onto primitive doubles. Over vec3s it
built a dozen vectors per contact per substep; on the same wall settled
by `solver/settled`, at eight iterations and four substeps, it went from
13.4 ms of a 48 ms step to 0.9 ms of a 37 ms one. (That step is slower
than the 32.9 above, measured in a different session; it is the
difference that is the result.) The arithmetic is the vector version's
in the same order, and 120 steps of the wall come out identical to the
last bit.

Then the narrow phase, in `allgo.physics.contact/box-box`: the fifteen
separating axes, and the face clipping that makes a manifold, both off
vec3s and onto primitive doubles and two small buffers. On the same
wall detection went from 13.0 ms of a 36 ms step to 4.0 ms of a 29 ms
one, `box-box` itself from 11 ms to 2 over 687 pairs. One unhinted
`aset` of a long into the id array had it at 14 ms for a while -- a
reflective call per corner -- which is the reflection check earning its
place. Same answers to the bit again: the wall's 120 steps, and 200,000
random pairs of rotated boxes against the vector version, every kind of
contact among them.

`record-approach!`, the closing speed restitution is measured against,
went the same way: 1.6 ms of a 17 ms step to 0.08, measured in a warm
REPL, and the three solvers' 120 steps identical to the bit.

Then warm starting's own bookkeeping (2026-10-23). Each step copied the
three impulses back into every contact map (0.9 ms on the wall's 2,500
contacts) so that the next could build a map of them keyed by pair and
feature id (0.4 ms) and look each new contact up in it. Now the step
keeps the arrays the solve ended with, the ids beside them, and an index
from each pair of bodies to where its points start; a contact is looked
for among the three or four of its own pair. The wall went from 8.7 /
15.6 / 14.5 ms to 6.4 / 13.8 / 12.8 for sequential impulse / TGS / XPBD
-- more than the bookkeeping measured alone, which is the garbage it no
longer makes -- and 120 steps under each, from walls settled fresh,
are identical to the bit.


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
holds it there. It still cannot stack past eight courses, and tall
columns remain open work.

**XPBD's normals turn with their faces now** (2026-10-23). It measured
penetration each substep along the normal as detection found it, while
a leaning column tilted the faces under it; measured across the face as
it now stands, a column lasts about twice as long before it leans over
and is quiet enough to fall asleep first -- eleven courses with
sleeping on, seven before. Running collision detection every substep
instead did less (nine) for four times the detection.

## Every contact at once

`allgo.physics.convex`, the fourth `:solver` (2026-10-23), and the
first that is not Gauss-Seidel. The step's velocities minimize one
convex function -- kinetic energy measured from where the bodies would
go untouched, plus a term per contact whose gradient is its impulse,
projected onto its friction cone in the norm of a small regularization
-- which is Drake's SAP (Castro, Han and Permenter 2022) and the same
shape as MuJoCo's solver (Todorov 2014). Newton's method finds it: each
step one solve with `M + J' G J`, which couples every body to every
other through the contacts, and a backtracking line search. On the
bricks demo's column (0.9 x 0.45 x 0.45, friction 0.55, nothing asleep,
20 s):

    stands up to      sequential impulse  7    xpbd  9    tgs  14
                      convex             40, and 50 falls at 2.7 s

What it took, in the order it was found:

- **Stiffness is the whole of it.** SAP's "near rigid" regularization,
  `R = beta^2 w / 4 pi^2` with `beta` 1, stood twelve courses and lost
  sixteen in two seconds. `beta` 0.1, with an overlap pushed out over one
  step rather than two, stands forty. 0.03 is worse again, presumably
  the Newton matrix's conditioning. The slop was suspected first and is
  not it: none at all bought a few seconds, not courses.
- **The Newton matrix is stored by its envelope.** Dense, a 16x16 wall
  is 1,536 rows and 600 ms a step. Numbered by reverse Cuthill-McKee so
  that bodies that touch sit near each other, and kept from each row's
  first nonzero to the diagonal, it is a tenth of the dense size and a
  column is a band. Then reflection -- arrays bound without hints -- was
  most of what was left: 102 ms a step on the wall became 16, and a
  forty course column 0.6 ms, under TGS's 1.5.
- **A gap asks for nothing but itself.** A contact within the margin
  asks only that the gap close no faster than the step allows; asking
  for at least zero instead, as a missing condition on restitution did,
  let a dropped box through the floor by six centimeters and threw it
  back up at three meters a second.
- **And it found a bug in `contact/box-box` the others had hidden.** A
  box settling flat onto a much larger one, inside the margin, came back
  as a single edge contact on the larger box's far edge: the 2% handicap
  that keeps an edge pair from displacing a face multiplied a *negative*
  overlap and so favored the edge. The impulse solvers correct the next
  step and never showed it; the convex solve satisfied the bogus point
  by spinning the box rather than stopping it. Scaled by 0.98 for a gap,
  and `gap-face-test` holds it. TGS's column went from thirteen courses
  to fourteen.

On the 16x16 wall it steps in about 16 ms, against 6.4 / 13.8 / 12.8
for sequential impulse / TGS / XPBD, and it passes every test in
`solver-test` the other three do.

- **Last step's factor, kept** (2026-10-23). The factor was then most
  of a step: 4 ms each (1,536 rows averaging 75 stored), assembled in
  another 2, and two a step settled, three settling. The matrix is
  indexed by bodies, not contacts, and a sticking contact's G is the
  constant R^-1 inside the cone, so from one step to the next it barely
  changes: kept, last step's factor preconditions conjugate gradients
  on the next matrix, which reaches 1e-8 in three to eight iterations
  of 0.5 ms (a product with M + J'GJ unformed, and a substitution) on a
  factor eleven steps old. Left alone the count creeps to sixteen, past
  the dozen where a factor is cheaper, so one that took more than eight
  is replaced at the next direction, and one that fails fifteen at
  once. The wall, freshly settled: 17.1 ms a step settled becomes 11.7
  to 12.4, and 30.4 settling becomes 25.9; the ragdoll scene, whose
  islands are small, is unchanged at 2.5. Looser inner tolerances,
  1e-6 to 1e-3, only bought more Newton steps and were slower.
  `solve` returns the factors by island and takes them back; the
  callers keep them while every body has the block it had.
- **And what was around the solve** (2026-10-23), the wall settled
  after five seconds, 11.2 ms a step to 9.8. The wall is two islands --
  255 bricks, and one touching only the floor -- so the big one is
  copied into a problem of its own every step, and that copy was
  sequences and a `zipmap`: 1.5 ms, now 0.6 on arrays. In `prepare`,
  which every solver shares, `tangents` built its basis as vectors, a
  third of the 1.6 ms; written into the arrays in primitives, with the
  same arithmetic in the same order, it is bit-identical. And last
  step's impulses are looked for first at the same place in the pair's
  manifold, where they nearly always are, before the scan of vector ids.
  `prepare` is 1.0 ms now, and sequential impulse, TGS and XPBD gain
  0.4 to 0.9 ms a step on the wall with it. What is left of the convex
  step is conjugate gradients (2.7 ms, at about a nanosecond a
  multiply-add), detection (1.7) and building the problem (1.0).
- **Detection's boxes** (2026-10-23). `contact/all` built every box's
  `frame` and then `aabb` built it again, in vectors: 0.38 ms of 1.60
  on the wall. From the frame already built, into the arrays in
  primitives, it is 0.09, and the contacts are the same to the last bit.
  Two loops in `box-box` that boxed their best axis every iteration
  (a `cond` whose last branch can be nil) were unboxed as well and
  measured as nothing: the JIT had been removing it.

Its two numbers were chosen on columns, so they were then tried on
everything else the bricks demo asks (JVM, settled first, friction
0.55):

                               sequential   tgs     xpbd    convex
    jenga, 12 levels:
      all asleep at                never    4.8 s   never   0.5 s
      sleeping off, mean speed     0.031    0.000   0.008   0.000
      blocks moved over 5 cm          21        0      31       0
    keep 6 + rampart 3, 60 s:
      blocks moved                     2        2       2       2
      all asleep at                2.2 s    2.6 s   7.0 s   2.8 s
    16x16 wall, sleeping off, 20 s:
      bricks moved over 5 cm         191      211     256       1
    a ball at 26 m/s into a 5x4 wall:
      most energy at any step   all at most the shot's own 97.3 J
      fastest brick (m/s)            5.4      6.3     6.3     8.1
    rebound off the floor, asked 0 / 0.5 / 0.9:
                                0/.51/.93 0/.48/.88 0/.49/.90 0/.51/.93
    a ball rolling at 2 m/s, after 2 s:
                                   1.99     1.22    1.41    1.99
    a box sliding at 3 m/s, mu 0.5, stops at (exact 0.612 s):
                                   0.60     0.60    0.65    0.617

Nothing there asks for different numbers. The fastest brick is the one
difference worth naming, and it is the stiffness doing what it should: a
rigid wall takes the ball's momentum all at once where a springy one
takes it over a few steps.

**With articulated bodies and joints** (`allgo.physics.world`,
`:solver :convex`). The problem is now in blocks: six numbers for a
loose body, a model's whole generalized velocity for a model, each with
its dense mass matrix -- `articulated/generalized-mass-matrix`, the
composite rigid body algorithm, which agrees with the inverse of
`articulated/inverse-mass-matrix` to 1e-13 at a fiftieth of the cost
(64 us against 2.3 ms on the ragdoll). Three kinds of constraint:
contacts with their cone, joint limits as one-sided rows, pins as
two-sided ones. What came of it:

- A column of thirty cubes stands in `world`, which the sequential
  impulse holds at six. The convex solve has its own slop, 1 mm: a
  regularized contact gives a little under load every step until the
  push past the slop balances it, so every contact in a stack ends the
  whole slop deep, and at the impulse solve's 5 mm a column of thirty
  stood fifteen centimeters short.
- The ragdoll scene settles completely, figure and bricks asleep, in
  2.4 ms a step against the sequential impulse's 1.4 (6.2 before the
  three changes below).
  - The mass matrix by the composite rigid body algorithm rather than
    by inverting the articulated body algorithm's columns: 2.3 ms off.
  - Islands solved apart. One Newton direction for the whole scene
    has one step length, and a ragdoll working through its contacts
    shortened and lengthened the steps of bricks it was not touching:
    they came off the floor and back three or four times in a solve.
  - The line search on row velocities. Along a line every row's
    velocity is linear in the step, so a trial is the projections
    alone; it was a full evaluation, 81 of 185 ms over two hundred
    solves.

  The ragdoll's own island still takes ten to thirty Newton steps, and
  it is genuine: the velocities are meters a second from the answer
  until the last few. A joint limit that is exactly at its kink at the
  solution (touching, carrying nothing) flips on and off every
  iteration, each step stops at the kink, and the rest converges
  linearly. Softer limits, ten to a thousand times, cut it by a few
  percent; so does letting the line search go past the full step.
- A bike and rider is 1.3 ms a 1/240 s step against 1.0.
- It is `world`'s default (2026-10-23), and the motorcycle alone stays
  on the sequential impulse it was tuned under (`motorcycle/world-defaults`),
  for two reasons. Thrown off at 8 m/s, the rider's
  leg stays tangled with the frame under the sequential impulse, eleven
  millimeters deep, and pushing that out kicks the bike over in three
  seconds and *speeds it up* 0.2 m/s. Under the convex solve the rider
  is clear in 0.2 s, and the riderless bike does what a real one does:
  steers into the lean, recovers, and capsizes slowly, over at eight
  seconds. And a full-lever stop from 15 m/s is a knife edge -- at 1 g
  this bike has 8% of its weight on the rear, and a stoppie runs away
  once the rear is up. The sequential impulse flips it from 1.1 of full
  lever, the convex solve from 1.0, so full lever is a stoppie to
  vertical and back. Neither is obviously the wrong answer.

**Convex in a browser** (2026-10-23). The ragdoll scene, run under Node
from a `:simple` build and profiled with `--cpu-prof`, was 88 ms a step
under the convex solve against the JVM's 2.5. Three things, in order:

- **`=` on numbers.** In ClojureScript it is `cljs.core/=`, a generic
  call through `-equiv`, where Clojure on primitive longs compiles it to
  a comparison -- and every loop in `allgo.physics.convex` ended on one.
  A fifth of the step. `==` is the numeric comparison in both: 88 ms
  became 30.
- **The Newton matrix's rows found once.** `hessian!` looked up each
  entry's place in the envelope separately; a row's columns from a
  block on are one run of it. Bit-identical, and 30 became 28.
- **No response for a row the convex solve does not need.** A model's
  side of a contact was prepared as the sequential impulse wants it,
  with what an impulse does to every joint -- most of the cost of
  `prepare`. The convex solve holds the mass matrix itself, so it now
  asks only for the row, and the effective mass the regularization
  scales by is `g' M^-1 g` from its own inverse: 23 ms. Islands under
  a hundred rows are factored afresh every direction, which there is
  cheaper than conjugate gradients on a kept factor: 20.

In Chrome, from a release build, the demo is at most 13 ms a frame as
the figure lands -- the sequential impulse 8 -- and both 2 once it
sleeps, so the demo starts on convex again. The 35-45 ms that had it
start on the sequential impulse were a development build, and bricks
that started inside the slab.

The rest of physics, profiled the same way (Node, ms a step):

                                 before   after
    motorcycle and rider, 1/240     5.6     2.9
    ragdoll, sequential impulse    17.4    10-14
    10x10 wall, sequential impulse 11.4     9.6
    10x10 wall, TGS                28      14
    10x10 wall, XPBD               21      16

`=` turned out not to be the rest of it. In a browser reading a vector
is a generic `nth` and a keyword lookup on a map of more than eight
keys a hashed search, and the hot paths did both per row per
iteration: `inverse-inertia-world` built the body's rotation matrix by
rotating three vectors, `warm-start!` assembled each impulse from
vectors, XPBD's `anchor!` and `angular-mass` looked their arrays up in
the arrays map per contact per substep, the world's sweep destructured
three vectors per side per iteration, and its `prepare-rigid` ran two
rotations and three products on vectors. The articulated root's six by
six was factored, by `allgo.numerics.linear`'s seq-based Cholesky, once
an impulse rather than once a configuration, and each spherical joint's
three by three inverted the same way; and the tire's deepest point asked
`sphere-box` for a whole contact fifty times. Each is now on arrays,
with the arithmetic of what it replaced in the same order, so the
results are the same to the last bit on the JVM.

Twenty crash test dummies at once (the ragdoll demo's `dummies`) were
60-100 ms a frame in Chrome, under either solver, and the solvers were
not most of it: two hundred and twenty parts are twenty-four thousand
pairs between models, and each pair's bounding-ball test worked out
both balls' radii from their maps and vectors, the pairs filtered from
every ordered combination; and each dummy's own fifty-five pairs worked
out both parts' world boxes per pair. Each ball and box is now worked
out once a step, the pairs walked once, in the same order, and twenty
dummies are 21-45 ms a frame, convex the quicker once they land.
Under Node: 158 and 176 ms a step (sequential impulse, convex) became
80 and 93. Then a broad phase for the parts, swept with
`allgo.spatial.sweep` -- a box around each bounding ball, padded so
that it passes every pair the ball test would -- and the candidates
tried in the order every pair once was, so again the same to the bit:
74 and 62, and forty dummies 94 and 112. The ball tests, boxes and
sweep are under 2% of a step now; each dummy's own limbs against each
other are most of what detection is left with.

**Limits met ahead** (2026-10-24). The convex solve's joints in a heap
were far past their limits -- a hip 0.45 rad -- and, with `:limit-slop`,
ran away; it looked like the convex solve's. Measured beside the
sequential impulse it was not: as twenty dummies landed, the worst
limit was 0.55 rad past under the sequential impulse and 0.60 under the
convex solve, and both held to a few thousandths once settled. A limit
was a row only once passed, so a joint turning at thirty radians a
second was half a radian past before anything held it. Now a limit is a
row ahead of time, as a contact is, asking only that the gap close no
faster than it can (`limit-ahead-test`): with the rows the joint's own
turning would reach this step, 0.55; with every limit a row every step
-- twenty a ragdoll, a sixth more time -- 0.37 and 0.47, the default.
What is left is a joint turning half a radian in a step, which a row
linear in the velocity follows only so far. The runaway went with the
twist row's fix, the rows met ahead, and a Newton solve that does not
converge now giving way to projected Gauss-Seidel, which holds its
answer bounded; the slop is no longer kept from the convex solve.

Twenty dummies awake under the sequential impulse, later (2026-10-24):
96 ms a step under Node, then 57. `q/rotate` built five vectors a call
through `vec3`, and now none, in the same arithmetic; `world` built
both boxes' frames for every pair, and now each once a step;
`allgo.array`'s typed arrays came from a vector or an array by way of a
seq; and a contact side on a model asked `articulated/response-at` for
`M^-1 g`, a walk of the tree that allocates as it goes, where the
model's inverse mass matrix -- twenty-five by twenty-five, built once a
step -- answers each row in one product. That last is the same numbers
but for rounding; the rest are the same to the bit. Then 53: a model's
parts are placed every step as `rigid/collider`s, a shape and where it
is, rather than as bodies with a mass and inertia nothing read; and
`self-contacts` takes them from `world`, which had already placed them,
rather than placing them again. Left is the articulated core --
configuration, articulated inertias, forward dynamics -- which makes
many small arrays a model a step.

**Why ragdolls did not sleep** (2026-10-24). A heap of twenty never
slept under either solver, and following the fastest joints found three
things, the first two of them bugs:

- A ball joint past both its cone and its twist was held on one, the
  cone. A neck lying at its cone had its twist unheld until the swing
  dipped back inside for a step, when the twist row closed all that had
  built up at once: the head flicked round at three radians a second
  every eight steps. Now both rows.
- The twist row pushed along the bone, which is the twist's rate only
  when the joint is not swung. Its true direction, from twist =
  2 atan2(v . a, w) and q' = q (omega, 0) / 2, is
  (w^2 a + w (a x v) + (v . a) v) / ((v . a)^2 + w^2), and checked by
  differences it differs from the bone's axis by as much as its own
  size -- once in sign. Held both ways along the wrong direction, the
  twist never came back, the push never stopped, and the convex solve,
  which holds every push exactly, took the heap to 1e89 radians a second
  in sixty steps. Along its true rate, the heap stays bounded.
- And a limb with nothing touching it swings for ever: a dummy whose
  shin hung over the slab's edge never slept. `allgo.physics.world` now
  takes a passive `:damping` on a link, its own friction rather than a
  drive, so it does not keep the model awake as a `:tau` does; with
  `d dt` of armature beside it, it is stable on a forearm. At 0.2 one
  dummy sleeps in three seconds.

A heap still did not sleep, kicked at one to three radians a second.
Single-point edge contacts between box limbs looked the likely cause,
so capsules went in -- `rigid/capsule`, a segment swept by a ball, met
as balls at its ends and at its nearest point to the other shape -- and
did not help: twenty capsule dummies slept no better, a lone one took
twelve seconds to where a box one takes three, since a capsule lying
on the floor rolls, and a step cost more. Capsules stay, as a shape;
the dummy's limbs stay boxes.

Followed again, the kicks came with no contact changing at all: a limb
resting on its joint limit was pushed just inside, the row went,
gravity took it back past, and the row came back cold and kicked it.
`:limit-slop` -- a limit's own slop, held and not pushed out, as a
contact's is -- put the demo's twenty dummies to sleep in ten seconds
and ten in five. It is not a cure. In a heap on a bare floor ten slept
no better with it than without: one dummy lay with its elbow bent past
its limit and its forearm against its own chest, the limit pushing out
and the chest pushing back, and the two fought for good. It defaults to
none: the motorcycle's rider, held near its limits by muscles, was
flung forty meters with it, and under the convex solve a heap ran away.
The ragdoll demo asks for 0.05, three degrees, and that answers the
elbow too: the slop as deep as a jammed limb sits, its limit holds it
there without pushing, and nothing is left to fight. Ten dummies on a
bare floor slept in eleven seconds at 0.05 and thirty at 0.02; across a
handful of other heaps -- still, spun, fifteen of them, dropped from
three meters -- every one slept in five to eleven (`heap-sleeps-test`).

**Newton in a heap** (2026-10-24). Twenty dummies under the convex solve
take twenty Newton steps on average and up to sixty-odd, and it is not
polish: the gradient sits at one to fifty for most of them and falls to
1e-11 in the last, once the right constraints are active, and the
velocities part way are metres a second from the answer. Tried:

- The warm start was never taken. A heap's joint limits carry much of
  its weight and were not warm started, so last step's contact impulses
  alone drove the limbs into their limits, at a cost five to a hundred
  and fifty times that of starting from nothing. Keyed by model, link
  and which limit, and warm started in both solvers: taken 181 times in
  300, and a step 104 ms to 82. The iterations hardly moved, 22 to 20.
- Softer contacts, `beta` 0.3: 18 steps, 72 ms -- but `beta` 0.1 is what
  stands a tall column, and 1.0 was slower again. Not taken.
- Softer friction, the tangents' regularization the normal's: worse, 27
  steps; ten times it blew up. Stiff friction is not what churns.

What is left is the active set found one constraint at a time.

**A wedge** (2026-10-24). Two contacts nearly opposed -- a brick between
the floor and a slab tilted eighteen degrees over it, 162 degrees apart
-- each ask to be pushed apart at up to `max-push-speed`, and held near
rigid, as the convex solve holds every contact, the only velocity that
satisfies both is sideways at that over the sine of nine degrees: the
ragdoll demo's wall, started inside its slab, threw bricks at sixteen
meters a second, and deeper went to NaN. The sequential impulse never
showed it, its iterations running out first. `convex/deep-give` softens
a contact past a centimeter of overlap, by up to its own inverse
effective mass added to its regularization by two centimeters: a push
can then change its body's velocity by about half what it asks, and
pushes in any number, opposed or not, by no more than they ask between
them. Sixteen meters a second became 3.5, the NaN 2.5 (`wedge-test`).
A resting stack is millimeters deep and never reaches it.

The first cut softened every deep contact, and a ball at 120 m/s went
through a wall: an impact is deep before it is solved. The second spared
contacts closing faster than a meter a second, and threw a brick at
twenty: a quarter of a meter inside a slab, the shortest way out flips
from step to step, and last step's push reads as closing. What tells
them apart is the depth: an impact closes by more in a step than it is
deep, and a stuck body does not.

## Tried and not kept

- **Easing a joint limit's kink for Newton** (2026-10-23). The ragdoll's
  hardest solve, traced, is still falling at its thirty-step cap, and
  nearly every step exactly one constraint changes state -- a limit on
  or off, a contact sliding or separating: the exact line search lands
  on a kink and the next direction is computed as if the row there were
  off. Two remedies, neither kept. Smoothing a unilateral row's push
  over `rn y` = 1e-6 to 1e-2 either side of zero: 7.2 mean steps became
  6.0 at best with the same tail to thirty, and at 1e-3 the scene went
  to NaN. Counting a row within that band as on in the Newton matrix
  only, the gradient exact: worse at every width, 7.3 to 14.2.

- **Three items measured and not built** (2026-10-23).
  - *XPBD in `allgo.physics.world`.* It was wanted for stacking, and
    the convex solve, now the default, stands thirty cubes where the
    sequential impulse holds six. It is the weakest of the four in the
    bricks comparison above, it needs the substeps that
    `allgo.physics.world` tried and dropped (2026-10-22), and on
    articulated models it would be a new integrator, constraints on
    joint angles at the position level.
  - *Manifold reduction with per-pair stiffness.* Under the convex
    solve there is nothing to reduce: a 9x8 wall, a ten-cube column and
    a thirty-box pile settle with every pair at four points or fewer.
    Only the impulse solvers make more (TGS: 43 of 163 wall pairs at six
    or eight), and there the reduction was a wash (below), and the
    per-pair stiffness it was to come with is the load-scaled stiffness
    below, worse than uniform at every setting.
  - *SIMD and graph coloring* (Box2D v3). The demos run in a browser
    from GitHub Pages: no SIMD short of WebAssembly, and a parallel
    solve over colored batches needs workers sharing memory, which
    needs cross-origin isolation headers Pages cannot send. It could
    only ever speed up the JVM.

- **Skipping the narrow phase for pairs that have not moved** (2026-10-23,
  measured and not built). On the settled 16x16 wall it is about 1.0 ms
  of a 9.8 ms convex step or a 6 ms sequential impulse one, so that is
  the most it could save -- and only where a scene is awake and still,
  which sleeping already covers: a settled pile is not asked about at
  all, and one settling or moving moves every pair every step. It would
  also cost exactness: contacts reused within a tolerance drift by up to
  it, and 0.1 mm of depth is 6 mm/s of push in the convex solve. What
  is left of the narrow phase is spread thin -- vector lookups, the
  clipping -- with no one place for a skip to remove.

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

What it was not, at first: anything but sequential impulse, with no
sleeping and warm starting that lasted a step. It has since gained
sleeping, warm starting across steps, and the convex solve as its
default (below); it still has no substepping and no XPBD.
`allgo.physics.solver` is still the thing for a scene that is only
bricks. Two articulated models touch each other through the same sweep
and can be pinned together, which is how the motorcycle carries its
rider -- and a pin between two links of one model closes a loop in it,
a four-bar linkage as three links and a pin (`loop-in-one-model-test`).

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
missing rather than entered unverified. `phenomena/search` finds
them instead, from VSOP87, checked against every verified row from 1000
to 3000: within half an hour for Mercury and Venus and a few hours for
the outer planets and stations, about the series' own accuracy. Table 31.B, the mean elements
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
