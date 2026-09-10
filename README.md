# allgo

A grab bag of algorithms — "all algorithms" — implemented in Clojure and
ClojureScript, with a browser demo for each one that has something to show.

## Layout

Everything lives under a single `allgo` root, one package per discipline:

| Package | Contents |
| --- | --- |
| `allgo.procedural` | Cellular caves, dungeon generation, Perlin noise, diamond-square terrain, TIN, mesh shading |
| `allgo.geometry` | Delaunay triangulation, 3D Voronoi cells, GJK/EPA collision detection |
| `allgo.simulation` | Reynolds' boids |
| `allgo.numerics` | Runge-Kutta, Runge-Kutta-Nystrom, Adams multistep, Gragg-Bulirsch-Stoer extrapolation |
| `allgo.astro` | Force models, frames, time scales, geodesy, Kepler elements, ephemerides, orbit determination |
| `allgo.demo` | ClojureScript demo viewers (browser) |
| `allgo.desktop` | JVM-only renderers (Quil, Swing, Lanterna) |

`allgo.numerics` and `allgo.astro` follow Montenbruck & Gill, *Satellite
Orbits*; planetary and lunar positions follow Vallado, *Fundamentals of
Astrodynamics and Applications*. Everything astrodynamical is expressed in
EME2000 so the pieces can be visualized in one frame.

Sources are split by platform: `src/cljc` for the algorithms themselves,
`src/cljs` for the demos, `src/clj` for the desktop renderers.

## Usage

Run the tests:

    clj -M:test

Build the demo page and serve it:

    npx shadow-cljs watch app

Then open <http://localhost:3000>.

## License

Copyright © 2015 Mark Bastian

Distributed under the Eclipse Public License either version 1.0 or (at
your option) any later version.
