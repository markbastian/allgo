# allgo

A grab bag of algorithms ("all algorithms") implemented in Clojure and
ClojureScript, with a browser demo for each one that has something to
show. The demos run live at <https://markbastian.github.io/allgo/>.

- **[NOTES.md](NOTES.md)**: how the finished parts work, what was
  measured along the way, and what was tried and not kept.
- **[TODOs.md](TODOs.md)**: what is still open.
- The namespace docstrings are the reference for each algorithm.

## Layout

Sources are split by platform: `src/cljc` for the algorithms themselves,
`src/cljs` for the browser demos, `src/clj` for the desktop renderers.
Everything lives under a single `allgo` root, one package per
discipline:

| Package            | Contents                                                                                                                                                                     |
|--------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `allgo.procedural` | Noise, fractals and shaping functions; caves, dungeons, diamond-square and TIN terrain, hydraulic erosion, DLA; islands with settlements, roads and names; planets and their atmospheres |
| `allgo.geometry`   | Vectors and quaternions, Delaunay triangulation, the Delaunay/Voronoi dual mesh, 3D Voronoi cells, GJK/EPA, hex grids, triangle and tetrahedral meshes                       |
| `allgo.spatial`    | Spatial hashing, sweep and prune, Morton codes, a bounding volume hierarchy                                                                                                  |
| `allgo.kinematics` | Serial chains, closed-form and numeric inverse kinematics, a redundant seven-joint human arm                                                                                 |
| `allgo.physics`    | XPBD soft bodies and cloth, rigid bodies and joints, three contact solvers, articulated bodies, mixed scenes, continuous collision, Eulerian and FLIP fluids, fire, height-field water, fluid on a sphere, skinning |
| `allgo.simulation` | Reynolds' boids, and a motorcycle with a rider and courses to ride                                                                                                           |
| `allgo.numerics`   | Runge-Kutta, Runge-Kutta-Nystrom, Adams multistep, extrapolation, FFT, tridiagonal and small dense linear solves, table interpolation, least-squares fits                     |
| `allgo.astro`      | Force models, frames, time scales, geodesy, Kepler elements, ephemerides, orbit determination, the three-body problem; calendars, coordinates, VSOP87, the Sun, Moon, planets and their satellites, eclipses, rising and setting, sundials |
| `allgo.demo`       | ClojureScript demo viewers (browser)                                                                                                                                         |
| `allgo.desktop`    | JVM-only renderers (Quil, Swing, Lanterna)                                                                                                                                   |

A few namespaces sit at the root because more than one package uses
them: `allgo.graph` (union-find, components, minimum spanning trees),
`allgo.search` (breadth-first, depth-first, Dijkstra, greedy best-first
and A* as one algorithm), and `allgo.array` (flat numeric arrays on
either platform).

Where a namespace follows a book, it says so. The main ones:

- **Montenbruck & Gill, *Satellite Orbits***: `allgo.numerics` and
  `allgo.astro`. Planetary and lunar positions follow Vallado,
  *Fundamentals of Astrodynamics and Applications*. Everything
  astrodynamical is expressed in EME2000, so the pieces can be
  visualized in one frame.
- **Meeus, *Astronomical Algorithms* (2nd ed.)**: all of its chapters
  that are algorithms, across `allgo.astro` and
  `allgo.numerics.interpolation` / `allgo.numerics.fit`, built on the
  Montenbruck & Gill pieces where the two books overlap. NOTES has the
  map from chapters to namespaces.
- **Ebert et al., *Texturing & Modeling: A Procedural Approach***:
  seven namespaces under `allgo.procedural`, built as small parts that
  compose freely rather than as finished effects.
- **Matthias Müller's [Ten Minute Physics](https://matthias-research.github.io/pages/tenMinutePhysics/)**:
  most of `allgo.physics` and `allgo.spatial`.
- **Craig, *Introduction to Robotics***: `allgo.kinematics`.
- **Featherstone**: `allgo.physics.articulated`.
- **[Red Blob Games](https://www.redblobgames.com/)**: hex grids and
  island generation.

## Demos

Live at <https://markbastian.github.io/allgo/>, or run `make serve` and
open <http://localhost:3000>. The landing page is a gallery of tiles --
searchable, filterable by topic, badged for what each demo needs or
offers -- and every tile opens its demo full window: the three with pages
of their own (`cnc/`, `moto/`, `solar/`) go there, the rest to the
player, `play/?demo=<id>`. The gallery is plain HTML and loads none of
the demos' bundle, so it is quick on a phone.

The list of demos is `resources/public/demos.json`, which both the
gallery and the player read. Adding one is an entry there and a line in
`allgo.demo.registry`. The tiles' pictures are screenshots of the demos
running, in `resources/public/thumbs/<id>.webp`; `make thumbs`, with
`make serve` running, takes them all again with a headless Chrome
(`script/thumbnails.mjs`), and `node script/thumbnails.mjs <id> ...`
just the ones named. A demo without one shows its topic's picture. The
groups:

| Group                 | Demos                                                                                          |
|-----------------------|------------------------------------------------------------------------------------------------|
| Procedural Generation | Cellular caves, dungeons, fractal terrain, hydraulic erosion, diffusion-limited aggregation, islands, planets |
| Simulation            | Fluid on a sphere                                                                              |
| Geometry              | Hex grids, Delaunay / Voronoi                                                                  |
| Spatial Queries       | Spatial hashing, broad phase (hashing, sweep and prune and a Morton BVH, raced)                |
| Robotics              | Inverse kinematics, redundancy in a human arm                                                  |
| Soft Bodies           | XPBD soft bodies, skinning, cloth                                                              |
| Rigid Bodies          | Bricks (three contact solvers), ragdoll, motorcycle, XPBD rigid bodies, joints                 |
| Fluids                | Eulerian fluid, FLIP fluid, fire and smoke, height-field water                                 |
| Flocking              | 2D, through a dungeon, with its Voronoi diagram, 3D, 3D Voronoi                                |
| Orbits & Numerics     | Lorenz, Kepler, three bodies, satellite perturbations, orbit determination, the solar system   |

The motorcycle is ridden with the keys, on a loop, an Excitebike lane or
a trials section. With no key down, the course's autopilot rides it. It
also has a full-window page of its own, `moto/`:
<https://markbastian.github.io/allgo/moto/>. So does the solar system,
with its nine thousand stars, `solar/`:
<https://markbastian.github.io/allgo/solar/>.

## Crossbows & Catapults

A game on the rigid-body engine, after Lakeside's 1983 *Crossbows and
Catapults* playset: two castles across a board, and discs flung at them
until one side's tower goes over. The catapult throws a disc through the
air; the crossbow slides one along the ground like a puck, through the
gate. Play against the computer or a second player, at
<https://markbastian.github.io/allgo/cnc/> or, under `make serve`,
<http://localhost:3000/cnc/index.html>.

The pieces are counted off a photograph of the set; the rules are a plain
version of the idea rather than Lakeside's. The rules and physics are
`allgo.simulation.crossbows`, and the page is `allgo.game.crossbows`.

## Usage

`make` on its own lists every target.

    make check     # lint + tests
    make serve     # watch and serve the demos on http://localhost:3000
    make release   # minified bundle for deployment

### Building the demo page

`make serve` rebuilds on save and serves `resources/public` at
<http://localhost:3000>. `make dev` does the same build once without
watching.

`make release` produces the deployable artifact: the whole page compiled
through Closure's advanced optimizations into a single minified
`resources/public/js/compiled/allgo.js`, plus `.gz` and `.br` copies.
Deploy `resources/public` as it stands -- six files, the page is static
and the bundle is the only script it loads. Release overwrites the same
file the development build writes, so run `make dev` afterward to get
readable names and source maps back.

The live site is GitHub Pages serving the `gh-pages` branch, which
holds nothing but the built pages: `index.html` (the gallery) and
`demos.json`, `play/index.html`, `cnc/index.html`, `moto/index.html`,
`solar/index.html`, `css/style.css` and `css/gallery.css`, `thumbs/`,
`js/compiled/allgo.js`, the three files in `data/` (the star catalogue,
the constellations and the star names, fetched by the solar system demo)
and an empty `.nojekyll`. Nothing rebuilds it; to publish, run
`make release`, copy those files onto `gh-pages`, commit and push.

    encoding      size
    identity      1.6M
    gzip          480K
    br            379K

`.br` is Brotli, which every current browser accepts and which beats gzip
by about 21% here. Both files are *hints to the server* -- nothing
requests them by name, so they are ignored unless the server is
configured to serve them in place of the `.js` with the matching
`Content-Encoding`:

    nginx   gzip_static on; brotli_static on;
    Caddy   file_server { precompressed br gzip }
    S3      upload with Content-Encoding metadata set

Netlify, Vercel, Cloudflare, GitHub Pages and the like compress on the
fly and ignore the pre-built copies; they cost nothing but are no help
there either.

Release builds delete the output directory first, which matters more than
it sounds: the development build writes ~15MB of per-namespace files into
it under `cljs-runtime/`, and those are stale the instant a release is
built but would still be deployed alongside it. Clearing first takes the
deployed tree from about 17MB to 2.6MB.

Size is dominated by three.js, at 630KB of the 1.6MB bundle -- using
`WebGLRenderer` reaches most of the library, so it does not tree-shake
usefully (`:js-provider :shadow` makes no difference). React, ReactDOM
and Reagent account for another 147KB, and are needed only by the
Delaunay demo. The two levers left, if the bundle ever
needs to be smaller, are splitting the WebGL demos into a lazily loaded
module so the first paint does not pay for three.js, and dropping Reagent
for plain DOM.

Two things to know about that build:

**It needs JDK 21 or newer.** The Closure compiler shadow-cljs shells out
to fails on 19 with a `NoClassDefFoundError` from deep inside the
compiler, which reads as a dependency conflict and is not one. The
Makefile finds a 21 with `/usr/libexec/java_home` and runs everything
under it, so `make` works regardless of what `java` on your PATH is. Pass
`JDK=/path/to/jdk` to choose a different one. If a build ever does fail
with a genuine-looking dependency conflict, `make clean` clears the
`.cpcache` a bad JVM may have left behind.

**A JS object read with `.-field` needs a `^js` hint.** `#js {:boids 30}`
compiles to a quoted key, `{"boids": 30}`, which Closure never renames;
`(.-boids controls)` compiles to a dotted access, which advanced
optimizations do rename. Without the hint the two stop agreeing, the
read comes back `undefined` in the release build alone, and the demo
renders an empty canvas with no error in the console. lil-gui is a
bystander: it reads by string, which is why the panel still shows 30
while the code sees nothing. The hint makes shadow-cljs emit externs for
the accesses so the name survives. Every `controls` object carries one:

    (def ^:private ^js controls #js {:boids 30 ...})

## Data and credits

The demos draw on other people's data, used under their terms:

- **Planet, Sun and Moon textures** (`resources/public/textures/`): the 2K
  maps from [Solar System Scope](https://www.solarsystemscope.com/textures/),
  used unmodified under the
  [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/)
  license.
- **Rotation models** (`allgo.astro.rotation`): the IAU Working Group on
  Cartographic Coordinates and Rotational Elements' 2015 report (Archinal
  et al. 2018), as distributed by NASA's NAIF in the planetary constants
  kernel `pck00011.tpc`.
- **IAU 1980 nutation** (`allgo.astro.frames`): the 106 terms of Table I
  of Seidelmann, "1980 IAU Theory of Nutation: The Final Report of the
  IAU Working Group on Nutation", *Celestial Mechanics* 27 (1982) 79,
  as the IERS also distributes them (`nut_IAU1980.dat`).
- **ERFA** ([liberfa/erfa](https://github.com/liberfa/erfa), derived from
  the IAU's SOFA library, BSD-3 license): the IAU 2006/2000A series in
  `allgo.astro.cio-data` are taken from it, with its copyright notice,
  conditions and disclaimer at the head of that file; and the FK5 and CIO
  reductions' tests compare against values computed with it.
- **Lear's and Gottlieb's gravity algorithms** (`allgo.astro.gravity`):
  transcribed from the normalized MATLAB listings in Eckman, Brown and
  Adamo, *Normalization and Implementation of Three Gravitational
  Acceleration Models*, NASA/TP-2016-218604, a work of the US government.
- **U.S. Standard Atmosphere, 1976** (`allgo.astro.us76`): its defining
  equations below 86 km and its Table I densities above, from the NOAA,
  NASA and USAF publication (NASA-TM-X-74335), a work of the US government.
- **NRLMSISE-00** (`allgo.astro.msis`, `allgo.astro.msis-data`):
  transcribed from the Fortran of Picone, Hedin and Drob (Naval Research
  Laboratory), a work of the US government, as NASA's CCMC ModelWeb
  archive distributes it; the coefficients generated from its BLOCK DATA,
  and its test driver's output in `test/data/msis/`.
- **Jacchia's models** (`allgo.astro.jacchia`, `allgo.astro.jacchia-roberts`):
  implemented from Jacchia's SAO Special Reports 313 (1970) and 332
  (1971) and Roberts's analytic form (Celestial Mechanics 4, 1971), read
  through the NASA Astrophysics Data System; the tests check them against
  those reports' tables.
- **JB2008** (`allgo.astro.jb2008`): implemented from Bowman et al.,
  AIAA 2008-6438, a work of the US government, and the JB2006 paper
  (J. Atmos. Solar-Terr. Phys. 70, 2008) for the equations it keeps; no
  code of the model's distributors was used. Its tests compare against
  densities computed with Orekit 12.2 (Apache License 2.0), in
  `test/data/jb2008/`.
- **VSOP87** (`allgo.astro.vsop87-data`): Bretagnon and Francou's
  planetary theory, version D, from CDS catalog VI/81, truncated as
  NOTES describes.
- **Meeus's tables** (the nutation, lunar, satellite, Pluto and
  phenomena series under `allgo.astro`): converted by script from Sonia
  Keys's Go port of *Astronomical Algorithms*,
  [soniakeys/meeus](https://github.com/soniakeys/meeus), MIT license.
  Checked, beyond the book's worked examples, against sources that owe
  nothing to Meeus: ERFA's eraMoon98, JPL Horizons, the US Naval
  Observatory and NASA's eclipse catalog (`test/data/meeus/`).
- **SGP4** (`allgo.astro.sgp4`): transcribed from the C++ that
  accompanies Vallado, Crawford, Hujsak and Kelso, *Revisiting Spacetrack
  Report #3* (AIAA 2006-6753), after Hoots and Roehrich's *Spacetrack
  Report No. 3*; the package's verification cases and outputs are in
  `test/data/sgp4/`. The code carries no license; its authors ask for a
  citation and a link to
  [celestrak.org/publications/AIAA/2006-6753](https://celestrak.org/publications/AIAA/2006-6753/).
- **Stars** (`resources/public/data/bsc5.tsv`): the Bright Star Catalogue,
  5th revised edition (Hoffleit & Warren 1991), VizieR catalog V/50.
- **Star names** (`resources/public/data/star-names.tsv`): the IAU Working
  Group on Star Names' catalogue.
- **Constellation figures and names**
  (`resources/public/data/constellations.tsv`): derived from
  [d3-celestial](https://github.com/ofrohn/d3-celestial), Copyright (c)
  2015 Olaf Frohn, under its BSD 3-Clause license.

## License

Copyright © 2015 Mark Bastian

Distributed under the Eclipse Public License either version 1.0 or (at
your option) any later version.
