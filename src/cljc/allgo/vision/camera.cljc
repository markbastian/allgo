(ns allgo.vision.camera
  "Image formation and 3D from two views: Russell and Norvig, *Artificial
  Intelligence: A Modern Approach*, 4th edition, chapter 25.

    perspective   the pinhole camera: a point at depth Z in front of it
                  lands at f X / Z, f Y / Z -- farther things smaller,
                  parallel lines meeting at vanishing points
    orthographic  scaled orthographic projection: every point divided by
                  one reference depth, which is close to perspective when
                  an object's depth varies little compared with how far
                  away it is
    shading       Lambert's law: a matte surface's brightness is its
                  albedo times the cosine of the light's angle to its
                  normal, plus what ambient light it gets
    triangulation the point two cameras both see, where the rays back
                  through its two images pass closest
    pose          where a camera must be to see a known object's points
                  where they appear
    stereopsis    depth from disparity, and how fine a difference in
                  depth fixated eyes can see
    light         surfaces diffuse with specularities, shadows getting
                  only the ambient, color a channel at a time; and color
                  constancy, the light's color estimated and divided out

  A camera is `{:position :rotation :f :center}`: where it is, the
  quaternion turning its own axes into the world's -- x right, y down, z
  out along the view, as image rows run -- its focal length in pixels,
  and the pixel the axis passes through. Rotations and vectors are
  `allgo.geometry`'s."
  (:require [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.numerics.linear :as lin]
            [allgo.numerics.optimize :as opt]
            [clojure.math :as math]))

(defn look-at
  "A camera at `position` looking at `target`, `up` (default +y) up in its
  picture."
  ([position target f center] (look-at position target [0.0 1.0 0.0] f center))
  ([position target up f center]
   (let [fwd (v/normalize (v/sub target position))
         right (v/normalize (v/cross (v/negate up) fwd))
         down (v/cross fwd right)]
     {:position position :f f :center center
      :rotation (q/from-matrix [[(right 0) (down 0) (fwd 0)]
                                [(right 1) (down 1) (fwd 1)]
                                [(right 2) (down 2) (fwd 2)]])})))

(defn to-camera
  "World point `p` in the camera's own coordinates."
  [{:keys [position rotation]} p]
  (q/rotate (q/conjugate rotation) (v/sub p position)))

(defn project
  "Where world point `p` lands in the camera's image, `[u v]`, by
  perspective; nil if it is not in front of the camera."
  [{:keys [f center] :as camera} p]
  (let [[x y z] (to-camera camera p)
        [cx cy] center]
    (when (pos? z) [(+ cx (/ (* f x) z)) (+ cy (/ (* f y) z))])))

(defn project-orthographic
  "Where `p` lands by scaled orthographic projection: every point scaled
  by f / `depth`, one reference depth for all."
  [{:keys [f center] :as camera} depth p]
  (let [[x y] (to-camera camera p)
        [cx cy] center
        s (/ f depth)]
    [(+ cx (* s x)) (+ cy (* s y))]))

(defn ray
  "The ray of world points that land at pixel `[u v]`: `[origin
  direction]`, the direction a unit vector."
  [{:keys [position rotation f center]} [u v]]
  (let [[cx cy] center]
    [position (v/normalize (q/rotate rotation [(/ (- u cx) f) (/ (- v cy) f) 1.0]))]))

(defn triangulate
  "The world point seen at `uv1` by `cam1` and `uv2` by `cam2`: the
  midpoint of the shortest segment between the two rays back through the
  images -- exact when the images are, and for noisy ones the point
  nearest both. Nil if the rays are parallel."
  [cam1 uv1 cam2 uv2]
  (let [[o1 d1] (ray cam1 uv1) [o2 d2] (ray cam2 uv2)
        w (v/sub o1 o2)
        b (v/dot d1 d2)
        den (- 1.0 (* b b))]
    (when (> den 1e-12)
      ;; the closest points o1 + s d1 and o2 + t d2: the segment between
      ;; them is perpendicular to both rays
      (let [[s t] (lin/mat-vec [[(/ -1.0 den) (/ b den)] [(/ (- b) den) (/ 1.0 den)]]
                               [(v/dot d1 w) (v/dot d2 w)])]
        (v/lerp (v/add-scaled o1 d1 s) (v/add-scaled o2 d2 t) 0.5)))))

(defn lambert
  "The brightness of a matte surface of `albedo` with unit `normal`, lit
  by a distant light of `intensity` from the unit direction `light`, with
  `ambient` light from everywhere."
  [albedo normal light intensity ambient]
  (* albedo (+ ambient (* intensity (max 0.0 (v/dot normal light))))))

(defn stereo-depth
  "The depth of a point whose images in two cameras side by side, a
  `baseline` apart with focal length `f`, are `disparity` pixels apart: f
  b / d -- the farther, the smaller the disparity."
  [f baseline disparity]
  (/ (* f baseline) disparity))

(defn angular-disparity
  "For two eyes `baseline` apart fixated at `depth`, the angular disparity
  in radians of a point `dz` farther: b dz / Z^2, for small angles."
  [baseline depth dz]
  (/ (* baseline dz) (* depth depth)))

(defn depth-resolution
  "The smallest difference in depth that eyes `baseline` apart, fixated
  at `depth`, can see when the smallest angular disparity they detect is
  `dtheta` radians: dtheta Z^2 / b -- finer the nearer, as the square of
  the distance."
  [baseline depth dtheta]
  (/ (* dtheta depth depth) baseline))

;; ---------------------------------------------------------------------------
;; Light, shading and color

(defn- channelwise
  "`f` of `xs` -- numbers, or [r g b] vectors, numbers standing for gray --
  a channel at a time if any is a color."
  [f & xs]
  (if (some vector? xs)
    (apply mapv f (map #(if (vector? %) % [% % %]) xs))
    (apply f xs)))

(defn reflect
  "The unit direction `light` mirrored about the unit `normal`: 2 (n . l)
  n - l, the way a mirror sends it."
  [normal light]
  (v/sub (v/scale normal (* 2.0 (v/dot normal light))) light))

(defn shade
  "The brightness of a surface point modeled, as the chapter says nearly
  every surface can be, as diffuse with specularities: Lambert's diffuse
  term, its `:albedo` times the light's `intensity` times the cosine of
  the light's angle to the `normal`, plus the `ambient` light its albedo
  returns from everywhere; and, where the mirror direction of the light
  comes near the `view` direction (unit, from the surface to the eye), a
  specularity -- Phong's lobe (\"Illumination for Computer Generated
  Pictures\", CACM 18, 1975), `:specular` times the intensity times the
  cosine of the angle between them to the power `:shininess`, the
  light's own color. A point `shadowed` from the light gets only the
  ambient. Albedos and intensities are numbers or [r g b]: Lambert's law
  a channel at a time, as the chapter models color."
  [{:keys [albedo specular shininess] :or {specular 0.0 shininess 1.0}} normal light view intensity ambient
   & {:keys [shadowed]}]
  (let [cos (if shadowed 0.0 (max 0.0 (v/dot normal light)))
        spec (if (pos? cos) (* specular (math/pow (max 0.0 (v/dot (reflect normal light) view)) shininess)) 0.0)]
    (channelwise (fn [rho i0 amb] (+ (* rho (+ amb (* i0 cos))) (* spec i0))) albedo intensity ambient)))

(defn sphere-image
  "A `w` by `h` image of a sphere of `radius` pixels at the center, seen
  from straight ahead by orthographic projection, of `material` (as
  `shade` takes it) lit by `intensity` from the unit direction `light`
  with `ambient` light; `background` elsewhere. Image x runs right and y
  down; normals have z toward the viewer."
  [w h radius material light intensity ambient background]
  (let [cx (/ w 2.0) cy (/ h 2.0)]
    {:w w :h h
     :px (vec (for [y (range h) x (range w)]
                (let [nx (/ (- (+ x 0.5) cx) radius) ny (/ (- (+ y 0.5) cy) radius)
                      r2 (+ (* nx nx) (* ny ny))]
                  (if (< r2 1.0)
                    (shade material [nx ny (math/sqrt (- 1.0 r2))] light [0.0 0.0 1.0] intensity ambient)
                    background))))}))

(defn gray-world
  "Color constancy by the gray-world assumption (Buchsbaum, \"A Spatial
  Processor Model for Object Colour Perception\", J. Franklin Inst. 310,
  1980): the scene's surfaces average to gray, so the image's mean color
  is the light's. `{:illuminant :image}`: the light's color, scaled to a
  mean of 1, and the image with each channel divided by it -- the scene
  as white light of the same brightness would show it."
  [{:keys [px] :as img}]
  (let [n (count px)
        means (mapv (fn [c] (/ (reduce + (map #(% c) px)) n)) (range 3))
        gray (/ (reduce + means) 3.0)
        e (mapv #(/ % gray) means)]
    {:illuminant e
     :image (assoc img :px (mapv (fn [p] (mapv / p e)) px))}))

;; ---------------------------------------------------------------------------
;; The pose of a known object

(defn- rows->camera
  "The camera whose rows `r` turn world directions into its own, at
  `position`."
  [r position f center]
  {:position position :f f :center center :rotation (q/from-matrix (lin/transpose r))})

(defn- orthographic-poses
  "The two poses scaled orthographic projection allows for `points` seen
  at `pixels`, mirror images in depth: the 2 x 3 map from the model's
  centered points to the image's by least squares, its rows made unit
  and square to each other, their scale f / Z of the model's center."
  [points pixels f [cx cy]]
  (let [n (count points)
        pm (lin/scale (reduce lin/add points) (/ 1.0 n))
        um (lin/scale (reduce lin/add pixels) (/ 1.0 n))
        P (mapv #(lin/sub % pm) points)
        U (mapv #(lin/sub % um) pixels)
        A (reduce lin/mat-add (map (fn [p] (mapv (fn [a] (mapv #(* a %) p)) p)) P))
        B (reduce lin/mat-add (map (fn [u p] (mapv (fn [a] (mapv #(* a %) p)) u)) U P))
        [m1 m2] (lin/mat-mul B (lin/spd-inverse A))
        s (/ (+ (lin/length m1) (lin/length m2)) 2.0)
        r1 (lin/normalize m1)
        r2 (lin/normalize (lin/sub m2 (lin/scale r1 (lin/dot m2 r1))))]
    (for [[a b] [[r1 r2] [(update r1 2 -) (update r2 2 -)]]]
      (let [r [a b (v/cross a b)]
            ;; the model's center in the camera's coordinates
            c [(/ (- (um 0) cx) s) (/ (- (um 1) cy) s) (/ f s)]]
        (rows->camera r (lin/sub pm (lin/mat-vec (lin/transpose r) c)) f [cx cy])))))

(defn pose
  "The camera -- position and rotation, its focal length `f` and center
  known -- that sees the model's world `points` at `pixels`: section
  25.6.4's recovery of a known object's pose from correspondences. Scaled
  orthographic projection gives a start, by linear least squares, in two
  mirror images; each is polished by Marquardt's method
  (`allgo.numerics.optimize`'s) on the squared error between the pixels
  and the points' perspective projections, over a turn and a move from
  the start, and the better kept. `{:camera :rms}`, the rms error in
  pixels."
  [points pixels f center]
  (let [polished (fn [{:keys [position rotation]}]
                   (let [at (fn [[w1 w2 w3 d1 d2 d3]]
                              (let [w [w1 w2 w3] angle (lin/length w)]
                                {:position (lin/add position [d1 d2 d3]) :f f :center center
                                 :rotation (if (< angle 1e-15) rotation (q/mul rotation (q/from-axis-angle w angle)))}))
                         ;; a point behind the camera has no image: a wall of error
                         sse (fn [p] (let [c (at p)]
                                       (reduce + (map (fn [pt px] (if-let [uv (project c pt)] (lin/distance-squared uv px) 1e12))
                                                      points pixels))))
                         {:keys [x] :as best} (opt/marquardt sse [0.0 0.0 0.0 0.0 0.0 0.0])]
                     {:camera (at x) :rms (math/sqrt (/ (:f best) (count points)))}))]
    (apply min-key :rms (map polished (orthographic-poses points pixels f center)))))
