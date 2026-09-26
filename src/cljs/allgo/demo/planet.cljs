(ns allgo.demo.planet
  "A procedural planet, after chapter 20 of Ebert et al.

  Nothing here is stored. The terrain, the coastline, the color of the
  ground and the cover of the clouds are all evaluated from the direction
  of each vertex, so the controls below are not adjusting a model of a
  planet -- they are choosing a different function, and the planet is
  whatever that function says.

  Three meshes, and each makes a different point:

  * **The ground**, a geodesic sphere whose vertices are pushed out to
    `planet/surface-radius` and colored by `planet/surface-color`. Its
    slope comes from the mesh's own normals rather than from four more
    terrain evaluations apiece, which is both cheaper and more honest --
    it is the slope of the surface actually being drawn.
  * **The clouds**, a second shell with per-vertex opacity. Deliberately
    not a texture map: a cloud layer painted from an equirectangular
    image smears at the poles and has a seam down one meridian, and the
    whole argument of the chapter is that a function of a 3D direction has
    neither.
  * **The air**, a shell whose color at each point is what
    `allgo.procedural.atmosphere` says the sky looks like there --
    blue over the day side, red along the terminator, black at night --
    faded toward the rim by the length of the path a grazing ray takes.
    That is why the glow reddens where the sun is setting without anyone
    having said that it should."
  (:require [allgo.demo.fps :as fps]
            [allgo.geometry.tri-mesh :as tri-mesh]
            [allgo.procedural.atmosphere :as atmosphere]
            [allgo.procedural.cellular :as cellular]
            [allgo.procedural.fractal :as fractal]
            [allgo.procedural.noise :as noise]
            [allgo.procedural.planet :as planet]
            ["lil-gui" :default GUI]
            ["three" :as THREE]
            ["three/examples/jsm/controls/OrbitControls.js" :refer [OrbitControls]]))

(def basis-functions
  "The basis functions on offer, by the name the control shows.

  The two at the bottom are not lattice noises and cost roughly thirty
  times as much per sample, since each one searches a neighborhood of
  cells rather than interpolating eight corners. They are worth having
  anyway -- a planet built on `F2-F1` is made of cracked plates and looks
  like nothing summed noise can produce -- but the demo caps the octaves
  and switches the distortion off when one is chosen, or a rebuild would
  take the best part of a minute.

  Simplex is worth switching to and from on the same seed. It is the
  isotropy that shows: Perlin noise is built by fading along each axis
  and carries a faint squareness because of it, which on a planet reads
  as coastlines that prefer to run north-south and east-west. The catch
  is that it is also louder -- half again the standard deviation for the
  same nominal range -- so the terrain comes out with more relief than
  the multifractal parameters were tuned for rather than the same world
  rendered more evenly."
  {"gradient (Perlin)" #(noise/gradient-basis {:seed %})
   "simplex"           #(noise/simplex-basis {:seed %})
   "value"             #(noise/value-basis {:seed %})
   "sparse convolution" #(noise/sparse-convolution-basis {:seed %})
   "cellular F1"       #(cellular/f1 {:seed %})
   "cellular F2-F1"    #(cellular/f2-f1 {:seed %})})

(defn- expensive? [basis-name]
  (contains? #{"sparse convolution" "cellular F1" "cellular F2-F1"} basis-name))

(def constructions
  "The fractal constructions, in the order the book introduces them."
  ["fbm" "turbulence" "multifractal" "hetero-terrain"
   "hybrid-multifractal" "ridged-multifractal"])

(def config
  {:radius 1.0
   :detail 5
   ;; The clouds get a finer shell than the ground does. Their opacity is
   ;; interpolated across the triangles, and an edge you can see through
   ;; shows the triangle it crosses; the ground's color is interpolated
   ;; too but has the terrain's own shading over it to hide the seams.
   :cloud-detail 6
   :air-detail 3
   :air-radius 1.09
   :cloud-radius 1.012
   ;; Real planets are far smoother than this. Earth's relief is about a
   ;; thousandth of its radius; at that scale a planet is a billiard ball,
   ;; so every picture of one exaggerates, and this is the exaggeration.
   :relief 0.045
   :samples 2048
   :sun [0.85 0.25 0.46]
   :sun-color [1.0 0.96 0.91]
   :sun-strength 2.3
   :ambient [0.075 0.1 0.15]
   ;; How far past the terminator the light wraps. A planet with air does
   ;; not have a razor edge between day and night -- the sky over the
   ;; ground just past the terminator is still lit, and lights it.
   :light-wrap 0.12
   ;; One map, read twice: once to bake the palette and once to fill the
   ;; shader's uniforms. The snow line is the only part of the rule the
   ;; fragment shader applies itself, so it is the only part that could
   ;; drift, and this is what stops it.
   :palette-opts {:snow-line 0.55 :steep 0.6 :pole-effect 0.45}
   :palette-size [256 128]
   ;; The detail layer. `:detail-period` is the noise period the tile
   ;; spans, and `:detail-frequency` how many tiles fit in a world unit --
   ;; together they put the finest features just below what the mesh can
   ;; hold, which is exactly where the geometry stops being able to help.
   :detail-size 256
   :detail-period 8
   :detail-octaves 4.0
   ;; Picked, not guessed. The tile is 256 texels across 8 noise cells, so
   ;; its coarsest features are about 32 texels and its finest about 4. At
   ;; this frequency the coarse end lands at roughly twice the mesh
   ;; spacing -- overlapping what the geometry already does, so the two
   ;; bands meet rather than leaving a gap -- and the fine end at about
   ;; two pixels, which is as fine as anything can usefully be.
   :detail-frequency 1.8
   :bump-strength 0.015
   :detail-tint 0.07})

(defn- ^:private sun-direction []
  (let [[x y z] (:sun config)
        l (js/Math.sqrt (+ (* x x) (* y y) (* z z)))]
    [(/ x l) (/ y l) (/ z l)]))

;; ---------------------------------------------------------------------------
;; Building the meshes

(defn- vertex-directions
  "The unit directions of a geodesic sphere's vertices, as a flat array."
  [verts]
  (let [n (/ (count verts) 3)
        out (js/Float64Array. (* 3 n))]
    (dotimes [i n]
      (let [b (* 3 i)
            x (nth verts b) y (nth verts (+ b 1)) z (nth verts (+ b 2))
            l (js/Math.sqrt (+ (* x x) (* y y) (* z z)))]
        (aset out b (/ x l))
        (aset out (+ b 1) (/ y l))
        (aset out (+ b 2) (/ z l))))
    out))

(defn- index-array [tri-ids]
  (let [a (js/Uint32Array. (count tri-ids))]
    (dotimes [i (count tri-ids)] (aset a i (nth tri-ids i)))
    a))

(defn- smooth-normals!
  "Area-weighted vertex normals, straight from the triangles.

  The cross product of two edges is twice the triangle's area along its
  normal, so accumulating it unnormalized weights each face by its size
  for free -- which is what you want, since a sliver should not sway the
  average as much as a large triangle does."
  [^js positions ^js index ^js normals]
  (.fill normals 0)
  (let [tris (/ (.-length index) 3)]
    (dotimes [t tris]
      (let [ia (* 3 (aget index (* 3 t)))
            ib (* 3 (aget index (+ (* 3 t) 1)))
            ic (* 3 (aget index (+ (* 3 t) 2)))
            ax (aget positions ia) ay (aget positions (+ ia 1)) az (aget positions (+ ia 2))
            bx (- (aget positions ib) ax) by (- (aget positions (+ ib 1)) ay)
            bz (- (aget positions (+ ib 2)) az)
            cx (- (aget positions ic) ax) cy (- (aget positions (+ ic 1)) ay)
            cz (- (aget positions (+ ic 2)) az)
            nx (- (* by cz) (* bz cy))
            ny (- (* bz cx) (* bx cz))
            nz (- (* bx cy) (* by cx))]
        (doseq [i [ia ib ic]]
          (aset normals i (+ (aget normals i) nx))
          (aset normals (+ i 1) (+ (aget normals (+ i 1)) ny))
          (aset normals (+ i 2) (+ (aget normals (+ i 2)) nz)))))
    (dotimes [i (/ (.-length normals) 3)]
      (let [b (* 3 i)
            x (aget normals b) y (aget normals (+ b 1)) z (aget normals (+ b 2))
            l (js/Math.sqrt (+ (* x x) (* y y) (* z z)))
            l (if (zero? l) 1.0 l)]
        (aset normals b (/ x l))
        (aset normals (+ b 1) (/ y l))
        (aset normals (+ b 2) (/ z l))))
    normals))

(defn- srgb->linear
  "A palette entry, as a light rather than as a color.

  `allgo.procedural.planet` hands back colors in the space a person reads
  them in -- 0.5 means mid gray -- and the renderer works in linear light,
  where mid gray is about 0.21. Uploading the one as the other is the
  reason so many procedural planets come out looking like a bleached
  photograph: everything dark is lifted and all the contrast goes."
  [c] (js/Math.pow (double c) 2.2))

(def ^:private detail-texture
  "A tile of fine relief, for the detail the mesh is too coarse to carry.

  Band-limiting the terrain to the mesh removes the speckle, and it also
  removes everything finer than a triangle -- correctly, because the mesh
  could never have shown it. But the *screen* has four or five times the
  resolution of the mesh, so there is a band of detail that pixels can
  hold and vertices cannot. This is that band, baked once into a tile and
  read per pixel.

  It is one channel: a height. The shader turns it into a normal with two
  screen-space derivatives, which is cheaper than storing a normal, needs
  no tangent frame, and is automatically the right strength at any
  magnification.

  Tileable, which is the point of `:period` on the basis. The tile is
  sampled by world position, three times on the three axis planes and
  blended by the normal, so it behaves as a solid texture: no seam, no
  poles, nothing that knows where the equator is. Built once and shared,
  because it does not depend on any planet parameter."
  (delay
    (let [n (:detail-size config)
          period (:detail-period config)
          f (fractal/fbm (noise/gradient-basis {:seed 4242 :period period})
                         {:octaves (:detail-octaves config) :H 0.8 :lacunarity 2.0})
          step (/ (double period) n)
          vals (js/Float32Array. (* n n))
          data (js/Uint8Array. (* 4 n n))]
      (dotimes [j n]
        (dotimes [i n]
          (aset vals (+ (* j n) i)
                (noise/sample f (* step i) (* step j) 0.5))))
      (let [lo (areduce vals k acc js/Number.POSITIVE_INFINITY (min acc (aget vals k)))
            hi (areduce vals k acc js/Number.NEGATIVE_INFINITY (max acc (aget vals k)))
            span (max 1e-9 (- hi lo))]
        (dotimes [k (* n n)]
          (let [v (js/Math.round (* 255 (/ (- (aget vals k) lo) span)))
                b (* 4 k)]
            (aset data b v) (aset data (+ b 1) v)
            (aset data (+ b 2) v) (aset data (+ b 3) 255))))
      (let [tex (THREE/DataTexture. data n n THREE/RGBAFormat THREE/UnsignedByteType)]
        (set! (.-wrapS tex) THREE/RepeatWrapping)
        (set! (.-wrapT tex) THREE/RepeatWrapping)
        ;; Mipmaps matter more here than anywhere else on the page: the
        ;; far limb crams many texels into a pixel, and without them the
        ;; bump becomes a crawling stipple exactly where the surface is
        ;; most foreshortened.
        (set! (.-generateMipmaps tex) true)
        (set! (.-minFilter tex) THREE/LinearMipmapLinearFilter)
        (set! (.-magFilter tex) THREE/LinearFilter)
        (set! (.-anisotropy tex) 4)
        (set! (.-needsUpdate tex) true)
        tex))))

(defn- palette-texture
  "Bakes `planet/surface-color` into a lookup table the shader can read.

  The splotchiness this replaces was not the terrain, it was the color:
  a vertex-colored mesh interpolates the *palette entry* across each
  triangle, so a coastline arrives as a green-to-blue smear a triangle
  wide rather than as a line. Interpolating the elevation instead and
  looking the color up per pixel puts the coastline back where the
  terrain says it is, at no cost in vertices -- the shore is now a contour
  of the interpolated height field, and the triangles are small enough
  that its kinks land inside a pixel.

  The table is indexed by elevation across the whole relief on one axis
  and by steepness on the other, and it is filled by calling the library
  rather than by restating it, so the palette has one definition. Snow is
  the exception: it depends on latitude as well, which would need a third
  axis, so it is suppressed here -- `:snow-line 2.0` puts it above
  everything -- and applied in the fragment shader, from the same
  parameters and the same palette entry.

  Float rather than bytes because the deep ocean lives in the bottom two
  percent of the range, where eight bits are four distinct values."
  [pl]
  (let [[w h] (:palette-size config)
        relief (:relief pl)
        data (js/Float32Array. (* 4 w h))
        opts (assoc (:palette-opts config) :snow-line 2.0 :pole-effect 0.0)
        equator (planet/direction 0.0 0.0)]
    (dotimes [j h]
      (dotimes [i w]
        (let [elev (+ (- relief) (* 2.0 relief (/ (+ 0.5 i) w)))
              slope (* 0.5 js/Math.PI (/ (+ 0.5 j) h))
              [r g b] (planet/surface-color pl equator
                                            (assoc opts :elevation elev :slope slope))
              base (* 4 (+ (* j w) i))]
          ;; Straight to linear light: the shader adds this to a lit scene,
          ;; not to a screen.
          (aset data base (srgb->linear r))
          (aset data (+ base 1) (srgb->linear g))
          (aset data (+ base 2) (srgb->linear b))
          (aset data (+ base 3) 1.0))))
    (let [tex (THREE/DataTexture. data w h THREE/RGBAFormat THREE/FloatType)]
      (set! (.-minFilter tex) THREE/NearestFilter)
      (set! (.-magFilter tex) THREE/NearestFilter)
      (set! (.-wrapS tex) THREE/ClampToEdgeWrapping)
      (set! (.-wrapT tex) THREE/ClampToEdgeWrapping)
      (set! (.-needsUpdate tex) true)
      tex)))

(defn- ground-geometry
  "Displace the sphere, and hand the shader what it needs to color it.

  Two attributes beyond position and normal: the elevation at the vertex
  and how steep the mesh is there. Both are scalars, both interpolate
  meaningfully across a triangle -- unlike a color, which does not --
  and between them they are the whole input to the palette."
  [pl {:keys [verts tri-ids]}]
  (let [dirs (vertex-directions verts)
        n (/ (count verts) 3)
        positions (js/Float32Array. (* 3 n))
        normals (js/Float32Array. (* 3 n))
        elevations (js/Float32Array. n)
        steepness (js/Float32Array. n)
        index (index-array tri-ids)]
    ;; The terrain is evaluated once per vertex and no more.
    (dotimes [i n]
      (let [b (* 3 i)
            d [(aget dirs b) (aget dirs (+ b 1)) (aget dirs (+ b 2))]
            h (planet/elevation pl d)
            r (+ (:radius pl) (max (:sea-level pl) h))]
        (aset elevations i h)
        (aset positions b (* r (nth d 0)))
        (aset positions (+ b 1) (* r (nth d 1)))
        (aset positions (+ b 2) (* r (nth d 2)))))
    (smooth-normals! positions index normals)
    (dotimes [i n]
      (let [b (* 3 i)
            ;; The angle between the mesh normal and straight up is the
            ;; slope, and the mesh knows both already -- so this costs
            ;; nothing, where asking the terrain would cost four more
            ;; evaluations and would answer about a surface finer than the
            ;; one being drawn.
            c (min 1.0 (max -1.0 (+ (* (aget normals b) (aget dirs b))
                                    (* (aget normals (+ b 1)) (aget dirs (+ b 1)))
                                    (* (aget normals (+ b 2)) (aget dirs (+ b 2))))))]
        (aset steepness i (/ (js/Math.acos c) (* 0.5 js/Math.PI)))))
    (doto (THREE/BufferGeometry.)
      (.setAttribute "position" (THREE/BufferAttribute. positions 3))
      (.setAttribute "normal" (THREE/BufferAttribute. normals 3))
      (.setAttribute "aElev" (THREE/BufferAttribute. elevations 1))
      (.setAttribute "aSteep" (THREE/BufferAttribute. steepness 1))
      (.setIndex (THREE/BufferAttribute. index 1)))))

(def ^:private ground-vertex-shader
  "attribute float aElev;
   attribute float aSteep;
   varying float vElev;
   varying float vSteep;
   varying vec3 vNormal;
   varying vec3 vWorld;
   varying vec3 vObject;
   varying vec3 vObjNormal;
   void main() {
     vElev = aElev;
     vSteep = aSteep;
     // Both spaces, because they are wanted for different things. World
     // space is where the light is and where the screen-space derivatives
     // have to be resolved; object space is where the planet is, and
     // anything painted on the planet has to be anchored there or it will
     // sit still while the planet turns underneath it.
     vObject = position;
     vObjNormal = normalize(normal);
     vNormal = normalize(mat3(modelMatrix) * normal);
     vWorld = (modelMatrix * vec4(position, 1.0)).xyz;
     gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
   }")

(def ^:private ground-fragment-shader
  ;; The palette comes out of the lookup table whole; only the snow line is
  ;; computed here, because it depends on latitude and a third axis in the
  ;; table would cost more than four lines of arithmetic. The constants are
  ;; the same ones the table was baked with.
  "uniform sampler2D uPalette;
   uniform sampler2D uDetail;
   uniform float uDetailFreq;
   uniform float uBumpStrength;
   uniform float uDetailTint;
   uniform vec3 uSnow;
   uniform vec3 uSunDir;
   uniform vec3 uSunColor;
   uniform vec3 uAmbient;
   uniform float uWrap;
   uniform float uRelief;
   uniform float uSeaLevel;
   uniform float uSnowLine;
   uniform float uPoleEffect;
   uniform float uSteepLimit;
   varying float vElev;
   varying float vSteep;
   varying vec3 vNormal;
   varying vec3 vWorld;
   varying vec3 vObject;
   varying vec3 vObjNormal;

   // The tile read as a solid texture: three planar projections blended by
   // the normal, so a sphere picks up detail everywhere and nothing in it
   // knows which way is north.
   float triplanar(vec3 p, vec3 n) {
     vec3 w = abs(n);
     w = w / max(1e-5, w.x + w.y + w.z);
     return w.x * texture2D(uDetail, p.yz).r
          + w.y * texture2D(uDetail, p.zx).r
          + w.z * texture2D(uDetail, p.xy).r;
   }

   // Mikkelsen's trick: a height field becomes a normal from two
   // screen-space derivatives, with no tangent frame to build and no
   // normal map to store. It also scales itself -- the derivatives are
   // per pixel, so the bump is the same strength however close the camera
   // gets.
   vec3 perturb(vec3 n, vec3 p, float h) {
     vec3 dpdx = dFdx(p);
     vec3 dpdy = dFdy(p);
     vec3 r1 = cross(dpdy, n);
     vec3 r2 = cross(n, dpdx);
     float det = dot(dpdx, r1);
     if (abs(det) < 1e-12) return n;
     vec3 g = (r1 * dFdx(h) + r2 * dFdy(h)) / det;
     return normalize(n - g);
   }

   void main() {
     vec3 n = normalize(vNormal);
     // Latitude comes from the object, not the world, so that it survives
     // the planet being turned or its axis being tilted.
     vec3 up = normalize(vObject);

     // Water is flat and stays flat; only the ground gets the detail.
     float land = step(uSeaLevel, vElev);
     // Sampled in object space so the detail is painted *on* the planet.
     // The perturbation is still resolved in world space, because that is
     // where the surface and the light meet -- a height field is just a
     // scalar, and it does not care which space it was read from.
     float detail = triplanar(vObject * uDetailFreq, normalize(vObjNormal));
     n = mix(n, perturb(n, vWorld, detail * uBumpStrength), land);

     // Elevation across the whole relief on one axis, steepness on the
     // other. Half a texel in on both, so the clamp at the ends does not
     // return a blend of the edge with nothing.
     vec2 uv = vec2(clamp((vElev + uRelief) / (2.0 * uRelief), 0.0, 1.0),
                    clamp(vSteep, 0.0, 1.0));
     vec3 base = texture2D(uPalette, uv).rgb;

     // Snow: the line comes down toward the poles, and will not lie on
     // anything too steep to hold it.
     float lat = abs(asin(clamp(up.y, -1.0, 1.0))) / 1.5707963;
     float line = max(0.0, uSnowLine - uPoleEffect * lat);
     float t = clamp((vElev - uSeaLevel) / max(1e-9, uRelief - uSeaLevel), 0.0, 1.0);
     float snowy = smoothstep(line, line + 0.12, t)
                 * (1.0 - smoothstep(0.8 * uSteepLimit, uSteepLimit, vSteep))
                 * step(uSeaLevel, vElev);
     vec3 albedo = mix(base, uSnow, snowy);
     // And a touch of the same field in the color, so that flat ground is
     // not flat in two senses at once.
     albedo *= 1.0 + uDetailTint * land * (detail - 0.5);

     float ndl = max((dot(n, normalize(uSunDir)) + uWrap) / (1.0 + uWrap), 0.0);
     gl_FragColor = vec4(albedo * (uAmbient + ndl * uSunColor), 1.0);

     #include <tonemapping_fragment>
     #include <encodings_fragment>
   }")

(defn- ground-material [pl sun]
  (let [{:keys [snow-line steep pole-effect]} (:palette-opts config)
        [sr sg sb] (mapv srgb->linear (:snow planet/terran))
        [ar ag ab] (:ambient config)
        [cr cg cb] (:sun-color config)
        strength (:sun-strength config)]
    (THREE/ShaderMaterial.
     #js {:uniforms
          #js {:uPalette #js {:value (palette-texture pl)}
               :uDetail #js {:value @detail-texture}
               :uDetailFreq #js {:value (:detail-frequency config)}
               :uBumpStrength #js {:value (:bump-strength config)}
               :uDetailTint #js {:value (:detail-tint config)}
               :uSnow #js {:value (THREE/Vector3. sr sg sb)}
               :uSunDir #js {:value (THREE/Vector3. (nth sun 0) (nth sun 1) (nth sun 2))}
               :uSunColor #js {:value (THREE/Vector3. (* strength cr) (* strength cg)
                                                      (* strength cb))}
               :uAmbient #js {:value (THREE/Vector3. ar ag ab)}
               :uWrap #js {:value (:light-wrap config)}
               :uRelief #js {:value (:relief pl)}
               :uSeaLevel #js {:value (:sea-level pl)}
               :uSnowLine #js {:value snow-line}
               :uPoleEffect #js {:value pole-effect}
               :uSteepLimit #js {:value steep}}
          :vertexShader ground-vertex-shader
          :fragmentShader ground-fragment-shader
          ;; dFdx/dFdy are an extension under GLSL ES 1.00.
          :extensions #js {:derivatives true}})))

(defn- cloud-geometry
  "A shell whose per-vertex alpha is the cloud cover under it."
  [pl {:keys [verts tri-ids]}]
  (let [dirs (vertex-directions verts)
        n (/ (count verts) 3)
        positions (js/Float32Array. (count verts))
        normals (js/Float32Array. (count verts))
        colors (js/Float32Array. (* 4 n))]
    (dotimes [i (count verts)] (aset positions i (nth verts i)))
    (dotimes [i n]
      (let [b (* 3 i)
            d [(aget dirs b) (aget dirs (+ b 1)) (aget dirs (+ b 2))]
            a (planet/cloud-cover pl d)]
        (aset normals b (nth d 0))
        (aset normals (+ b 1) (nth d 1))
        (aset normals (+ b 2) (nth d 2))
        (aset colors (* 4 i) 1.0)
        (aset colors (+ (* 4 i) 1) 1.0)
        (aset colors (+ (* 4 i) 2) 1.0)
        ;; Never quite opaque: real cloud tops are bright but the ground
        ;; still shows faintly through the thin edges, and a layer that
        ;; reaches 1.0 reads as paint.
        (aset colors (+ (* 4 i) 3) (* 0.88 a))))
    (doto (THREE/BufferGeometry.)
      (.setAttribute "position" (THREE/BufferAttribute. positions 3))
      (.setAttribute "normal" (THREE/BufferAttribute. normals 3))
      (.setAttribute "color" (THREE/BufferAttribute. colors 4))
      (.setIndex (THREE/BufferAttribute. (index-array tri-ids) 1)))))

(def ^:private air-vertex-shader
  "attribute vec3 aSky;
   varying vec3 vSky;
   varying vec3 vWorld;
   void main() {
     vSky = aSky;
     vWorld = (modelMatrix * vec4(position, 1.0)).xyz;
     gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
   }")

(def ^:private air-fragment-shader
  ;; How far the view ray through this fragment passes from the planet's
  ;; center: at the limb it grazes the ground and has the whole atmosphere
  ;; to cross, and at the outer edge of the shell it has none. That ratio
  ;; is the glow, and it is why the halo is a ring rather than a wash.
  "uniform float uInner;
   uniform float uOuter;
   varying vec3 vSky;
   varying vec3 vWorld;
   void main() {
     vec3 d = normalize(vWorld - cameraPosition);
     float t = -dot(cameraPosition, d);
     float b = length(cameraPosition + t * d);
     float k = clamp((uOuter - b) / (uOuter - uInner), 0.0, 1.0);
     gl_FragColor = vec4(vSky * pow(k, 3.5) * 0.85, 1.0);
   }")

(defn- air-geometry
  "A shell carrying, at each vertex, the color of the sky underneath it."
  [pl {:keys [verts tri-ids]} sun]
  (let [dirs (vertex-directions verts)
        n (/ (count verts) 3)
        atm (atmosphere/atmosphere {:radius (:radius pl)})
        positions (js/Float32Array. (count verts))
        sky (js/Float32Array. (count verts))]
    (dotimes [i (count verts)] (aset positions i (nth verts i)))
    (dotimes [i n]
      (let [b (* 3 i)
            d [(aget dirs b) (aget dirs (+ b 1)) (aget dirs (+ b 2))]
            ;; Standing on the ground there and looking straight up.
            eye (mapv #(* (+ (:radius pl) 1e-4) %) d)
            [r g bl] (atmosphere/tone-map
                      1.1
                      (or (atmosphere/sky-color atm eye d sun {:steps 10 :sun-steps 4})
                          [0.0 0.0 0.0]))]
        (aset sky b r)
        (aset sky (+ b 1) g)
        (aset sky (+ b 2) bl)))
    (doto (THREE/BufferGeometry.)
      (.setAttribute "position" (THREE/BufferAttribute. positions 3))
      (.setAttribute "aSky" (THREE/BufferAttribute. sky 3))
      (.setIndex (THREE/BufferAttribute. (index-array tri-ids) 1)))))

(defn- air-material [inner outer]
  (THREE/ShaderMaterial.
   #js {:uniforms #js {:uInner #js {:value inner} :uOuter #js {:value outer}}
        :vertexShader air-vertex-shader
        :fragmentShader air-fragment-shader
        :side THREE/BackSide
        :blending THREE/AdditiveBlending
        :transparent true
        :depthWrite false}))

(defn- effective-octaves
  "What the mountain layer's octave count actually comes to once the mesh
  has had its say.

  The same call `allgo.procedural.planet/default-terrain` makes -- repeated
  here only so the number can be shown, because a control that silently
  stops mattering above three is worse than no control."
  [{:keys [octaves lacunarity distortion basis sample-spacing]}]
  (let [asked (if (expensive? basis) (min 5.0 octaves) octaves)
        distortion (if (expensive? basis) 0.0 distortion)]
    (if sample-spacing
      (min asked (fractal/octaves-for (* sample-spacing 2.2 (+ 1.0 distortion)) lacunarity))
      asked)))

(defn- build-terrain [{:keys [basis construction octaves lacunarity roughness
                              distortion seed] :as opts}]
  (let [pricey (expensive? basis)]
    (planet/default-terrain
     {:seed seed
      :basis-fn (get basis-functions basis
                     (get basis-functions "gradient (Perlin)"))
      :construction (keyword construction)
      ;; The non-lattice bases are too dear to run at nine octaves through
      ;; a distorted domain; cap them rather than let the page hang.
      :octaves (if pricey (min 5.0 octaves) octaves)
      :distortion (if pricey 0.0 distortion)
      :lacunarity lacunarity
      :H roughness
      :sample-spacing (:sample-spacing opts)})))

(defn- build-world [{:keys [sea-level clouds sample-spacing] :as opts}]
  (planet/planet {:radius (:radius config)
                  :relief (:relief config)
                  :sea-level (* sea-level (:relief config))
                  :cloud-cover clouds
                  :samples (:samples config)
                  :sample-spacing sample-spacing
                  :terrain (build-terrain opts)}))

;; ---------------------------------------------------------------------------
;; The scene

(defn- dispose! [^js mesh]
  (when mesh
    (.dispose (.-geometry mesh))
    (.dispose (.-material mesh))))

(defn init! [^js container]
  (let [scene (THREE/Scene.)
        camera (THREE/PerspectiveCamera. 45 (/ (.-clientWidth container)
                                               (max 1 (.-clientHeight container)))
                                         0.01 100)
        renderer (THREE/WebGLRenderer. #js {:antialias true})
        sun (sun-direction)
        running? (atom false)
        tick-fps! (fps/meter! container)
        group (THREE/Group.)
        state (atom {:ground nil :cloud nil :air nil :spin true :wireframe false})
        ;; lil-gui reaches into this with string keys it cannot rename, so
        ;; it needs the ^js hint or advanced optimization reads it as
        ;; undefined and every control comes back empty.
        ^js controls #js {:basis "gradient (Perlin)"
                          :construction "ridged-multifractal"
                          :octaves 9
                          :octavesUsed 0
                          :lacunarity 2.0
                          :roughness 0.9
                          :distortion 0.35
                          :sealevel 0.0
                          :clouds 0.4
                          :detail (:detail config)
                          :seed 1
                          :showClouds true
                          :showAir true
                          :wireframe false
                          :spin true}]
    (set! (.-background scene) (THREE/Color. 0x05060d))
    (.setSize renderer (.-clientWidth container) (max 1 (.-clientHeight container)))
    (.setPixelRatio renderer (min 2 (or js/window.devicePixelRatio 1)))
    ;; Filmic tone mapping and an sRGB output. Without them the sunward
    ;; side clips to white -- one hard directional light on a rough sphere
    ;; puts a wide region past 1.0 -- and the clipping takes the terrain
    ;; shading with it, so the mountains disappear exactly where they are
    ;; best lit.
    (set! (.-toneMapping renderer) THREE/ACESFilmicToneMapping)
    (set! (.-toneMappingExposure renderer) 1.0)
    (set! (.-outputEncoding renderer) THREE/sRGBEncoding)
    (.appendChild container (.-domElement renderer))
    (.add scene group)
    (.add scene (THREE/AmbientLight. 0x2a3a52 0.18))
    (let [light (THREE/DirectionalLight. 0xfff6e8 2.4)]
      (.set (.-position light) (nth sun 0) (nth sun 1) (nth sun 2))
      (.add scene light))
    (.set (.-position camera) 0 0.85 3.4)
    (let [orbit (OrbitControls. camera (.-domElement renderer))]
      (set! (.-enableDamping orbit) true)
      (set! (.-minDistance orbit) 1.15)
      (set! (.-maxDistance orbit) 12)
      (letfn [;; Variadic because lil-gui hands its handlers the new value,
              ;; and this is also called directly with none.
              (rebuild! [& _]
                (doseq [k [:ground :cloud :air]]
                  (when-let [^js m (get @state k)]
                    (.remove group m)
                    (dispose! m)))
                (let [detail (long (.-detail controls))
                      ;; Built once and measured, then handed to the
                      ;; terrain so it knows how much detail is worth
                      ;; generating, and to the geometry so it is not
                      ;; subdivided twice.
                      mesh (tri-mesh/geodesic detail)
                      opts {:basis (.-basis controls)
                            :construction (.-construction controls)
                            :octaves (double (.-octaves controls))
                            :lacunarity (double (.-lacunarity controls))
                            :roughness (double (.-roughness controls))
                            :distortion (double (.-distortion controls))
                            :seed (long (.-seed controls))
                            :sea-level (double (.-sealevel controls))
                            :clouds (double (.-clouds controls))
                            :sample-spacing (tri-mesh/mean-edge-length mesh)}
                      pl (build-world opts)
                      ground (THREE/Mesh.
                              (ground-geometry pl mesh)
                              (let [m (ground-material pl sun)]
                                (set! (.-wireframe m) (.-wireframe controls))
                                m))
                      cloud (THREE/Mesh.
                             (cloud-geometry pl (tri-mesh/geodesic
                                                 (min detail (:cloud-detail config))
                                                 (:cloud-radius config)))
                             (THREE/MeshStandardMaterial.
                              #js {:vertexColors true :transparent true
                                   :depthWrite false :roughness 1.0
                                   :side THREE/FrontSide}))
                      air (THREE/Mesh.
                           (air-geometry pl (tri-mesh/geodesic (:air-detail config)
                                                               (:air-radius config))
                                         sun)
                           ;; The glow starts at the sea surface, which
                           ;; moves with the sea level -- start it at the
                           ;; datum instead and a high sea would show a
                           ;; ring of background between water and air.
                           (air-material (+ (:radius pl) (:sea-level pl))
                                         (:air-radius config)))]
                  ;; What the octave control actually came to. Above about
                  ;; three it stops mattering at this mesh, and a number
                  ;; that visibly stops moving explains that better than a
                  ;; paragraph would.
                  (set! (.-octavesUsed controls)
                        (/ (js/Math.round (* 10 (effective-octaves opts))) 10))
                  (set! (.-visible cloud) (.-showClouds controls))
                  (set! (.-visible air) (.-showAir controls))
                  (.add group ground)
                  (.add group cloud)
                  (.add group air)
                  (swap! state assoc :ground ground :cloud cloud :air air)))
              (on-resize []
                ;; A hidden card measures 0x0, which would make the aspect NaN.
                (let [w (.-clientWidth container) h (.-clientHeight container)]
                  (when (and (pos? w) (pos? h))
                    (set! (.-aspect camera) (/ w h))
                    (.updateProjectionMatrix camera)
                    (.setSize renderer w h))))
              (animate []
                (when @running?
                  (js/requestAnimationFrame animate)
                  (let [t0 (js/performance.now)]
                    (when (:spin @state)
                      (set! (.-y (.-rotation group)) (+ (.-y (.-rotation group)) 0.0012))
                      (when-let [^js c (:cloud @state)]
                        ;; The weather moves over the ground, slowly.
                        (set! (.-y (.-rotation c)) (+ (.-y (.-rotation c)) 0.0004))))
                    (.update orbit)
                    (.render renderer scene camera)
                    (tick-fps! (- (js/performance.now) t0)))))]
        (.observe (js/ResizeObserver. (fn [& _] (on-resize))) container)
        (rebuild!)
        (let [gui (GUI. #js {:container container})
              terrain-folder (.addFolder gui "terrain")
              world-folder (.addFolder gui "world")
              view-folder (.addFolder gui "view")]
          (doto terrain-folder
            (-> (.add controls "construction" (clj->js constructions)) (.onChange rebuild!))
            (-> (.add controls "basis" (clj->js (vec (keys basis-functions)))) (.onChange rebuild!))
            (-> (.add controls "octaves" 1 12 1) (.onFinishChange rebuild!))
            (-> (.add controls "lacunarity" 1.5 3.0 0.01) (.onFinishChange rebuild!))
            (-> (.add controls "roughness" 0.1 1.5 0.01) (.onFinishChange rebuild!))
            (-> (.add controls "distortion" 0.0 1.0 0.01) (.onFinishChange rebuild!))
            (-> (.add controls "octavesUsed") (.name "octaves used") (.disable) (.listen)))
          (doto world-folder
            (-> (.add controls "sealevel" -1.0 1.0 0.01) (.onFinishChange rebuild!))
            (-> (.add controls "clouds" 0.0 1.0 0.01) (.onFinishChange rebuild!))
            (-> (.add controls "seed" 1 999 1) (.onFinishChange rebuild!))
            (-> (.add controls "detail" 3 6 1) (.onFinishChange rebuild!)))
          (doto view-folder
            (-> (.add controls "showClouds")
                (.onChange (fn [v] (when-let [^js c (:cloud @state)] (set! (.-visible c) v)))))
            (-> (.add controls "showAir")
                (.onChange (fn [v] (when-let [^js a (:air @state)] (set! (.-visible a) v)))))
            (-> (.add controls "wireframe")
                (.onChange (fn [v]
                             (when-let [^js g (:ground @state)]
                               (set! (.-wireframe (.-material g)) v)))))
            (-> (.add controls "spin")
                (.onChange (fn [v] (swap! state assoc :spin v))))))
        {:start (fn [] (when-not @running?
                         (reset! running? true)
                         (on-resize)
                         (animate)))
         :stop (fn [] (reset! running? false))}))))

(defonce ^:private controller (atom nil))

(defn start!
  "Build the scene on first selection, then resume the render loop."
  []
  (when-let [c (or @controller
                   (when-let [container (js/document.getElementById "planet")]
                     (reset! controller (init! container))))]
    ((:start c))))

(defn stop!
  "Halt the render loop, keeping the context and the framing."
  []
  (when-let [c @controller] ((:stop c))))
