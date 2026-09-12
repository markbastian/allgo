# allgo

A grab bag of algorithms — "all algorithms" — implemented in Clojure and
ClojureScript, with a browser demo for each one that has something to show.

## Layout

Everything lives under a single `allgo` root, one package per discipline:

| Package            | Contents                                                                                      |
|--------------------|-----------------------------------------------------------------------------------------------|
| `allgo.procedural` | Cellular caves, dungeons, noise bases, fractals, procedural planets, terrain, TIN, atmosphere |
| `allgo.geometry`   | Delaunay triangulation, 3D Voronoi cells, GJK/EPA collision detection                         |
| `allgo.simulation` | Reynolds' boids                                                                               |
| `allgo.numerics`   | Runge-Kutta, Runge-Kutta-Nystrom, Adams multistep, Gragg-Bulirsch-Stoer extrapolation         |
| `allgo.astro`      | Force models, frames, time scales, geodesy, Kepler elements, ephemerides, orbit determination |
| `allgo.demo`       | ClojureScript demo viewers (browser)                                                          |
| `allgo.desktop`    | JVM-only renderers (Quil, Swing, Lanterna)                                                    |

`allgo.numerics` and `allgo.astro` follow Montenbruck & Gill, *Satellite
Orbits*; planetary and lunar positions follow Vallado, *Fundamentals of
Astrodynamics and Applications*. Everything astrodynamical is expressed in
EME2000 so the pieces can be visualized in one frame.

### Texturing & Modeling

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

Nothing is stored: the coastline, the colour of the ground and the cover of
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

Sources are split by platform: `src/cljc` for the algorithms themselves,
`src/cljs` for the demos, `src/clj` for the desktop renderers.

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
through Closure's advanced optimisations into a single minified
`resources/public/js/compiled/allgo.js`, plus `.gz` and `.br` copies.
Deploy `resources/public` as it stands -- eight files, the page is static
and the bundle is the only script it loads. Release overwrites the same
file the development build writes, so run `make dev` afterwards to get
readable names and source maps back.

    encoding      size
    identity      1.4M
    gzip          401K
    br            320K

`.br` is Brotli, which every current browser accepts and which beats gzip
by about 18% here. Both files are *hints to the server* -- nothing
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
deployed tree from about 17MB to 2.2MB.

Size is dominated by three.js, at 630KB of the 1.4MB bundle -- using
`WebGLRenderer` reaches most of the library, so it does not tree-shake
usefully (`:js-provider :shadow` makes no difference). React, ReactDOM
and Reagent account for another 147KB, and are needed only by the demo
picker and the Delaunay tile. The two levers left, if the bundle ever
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

**Anything a JS library reads by name needs a `^js` hint.** Advanced
optimisations rename properties, and lil-gui reaches into the demos'
control objects with string keys it cannot rename to match. Without the
hint on the `def`, shadow-cljs emits no externs for those accesses,
every read comes back `undefined` in the release build alone, and the
demo renders an empty canvas with no error in the console. Every
`controls` object therefore carries one:

    (def ^:private ^js controls #js {:boids 30 ...})

## License

Copyright © 2015 Mark Bastian

Distributed under the Eclipse Public License either version 1.0 or (at
your option) any later version.
