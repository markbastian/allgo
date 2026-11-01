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

  A camera is `{:position :rotation :f :center}`: where it is, the
  quaternion turning its own axes into the world's -- x right, y down, z
  out along the view, as image rows run -- its focal length in pixels,
  and the pixel the axis passes through. Rotations and vectors are
  `allgo.geometry`'s."
  (:require [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.numerics.linear :as lin]))

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
