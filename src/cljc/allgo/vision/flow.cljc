(ns allgo.vision.flow
  "Optical flow: Russell and Norvig, *Artificial Intelligence: A Modern
  Approach*, 4th edition, chapter 25 -- how far each part of an image
  moved between two frames.

    block matching   the chapter's method: a window round the pixel in the
                     first frame compared, by the sum of squared
                     differences, with windows shifted by every whole
                     displacement within a search radius in the second;
                     the best match is the flow
    Lucas-Kanade     the brightness of a moving point stays the same, so
                     Ix u + Iy v + It = 0 at every pixel; over a window,
                     that is two unknowns from many equations, solved by
                     least squares -- subpixel, and refined a step at a
                     time by warping the second frame back by the flow
                     so far
    a moving camera  the flow radiates from the focus of expansion, the
                     point the camera is heading for; each point's
                     distance from it over its speed is the time to
                     contact, and sideways, the flows' lengths go as the
                     inverse depths

  A flow is `[u v]` in pixels, x right and y down."
  (:require [allgo.numerics.linear :as lin]
            [allgo.vision.image :as img]
            [clojure.math :as math]))

(defn- window [x y half]
  (for [dy (range (- half) (inc half)) dx (range (- half) (inc half))] [(+ x dx) (+ y dy)]))

(defn ssd
  "The sum of squared differences between `a`'s window of half-width
  `half` round `(x, y)` and `b`'s shifted by `(u, v)`."
  [a b x y half u v]
  (reduce + (map (fn [[px py]] (let [d (- (img/at a px py) (img/at b (+ px u) (+ py v)))] (* d d)))
                 (window x y half))))

(defn block-match
  "The whole-pixel flow at `(x, y)` from `a` to `b`: the shift within
  `search` pixels each way whose window best matches, smallest shift on a
  tie."
  [a b x y half search]
  (let [shifts (sort-by (fn [[u v]] (+ (* u u) (* v v)))
                        (for [v (range (- search) (inc search)) u (range (- search) (inc search))] [u v]))]
    (apply min-key #(ssd a b x y half (% 0) (% 1)) (reverse shifts))))

(defn- lucas-kanade-step
  "The flow refined by a Gauss-Newton step on `b` warped back by it."
  [a b pts Ix Iy inverse [u v]]
  (let [It (mapv (fn [[px py]] (- (img/sample b (+ px u) (+ py v)) (img/at a px py))) pts)
        bx (- (reduce + (map * Ix It))) by (- (reduce + (map * Iy It)))
        [du dv] (lin/mat-vec inverse [bx by])]
    [(+ u du) (+ v dv)]))

(defn lucas-kanade-seq
  "The flow at `(x, y)` from `a` to `b` by Lucas-Kanade over a window of
  half-width `half`, each step a least-squares correction with `b`
  warped back by the estimate so far: the estimates, from `[0 0]`. Nil
  once the window has too little texture to say -- an aperture problem.
  `grads`, `a`'s gradients, may be passed in when many points share them."
  ([a b x y half] (lucas-kanade-seq a b x y half (img/gradients a)))
  ([a b x y half [gx gy]]
   (let [pts (vec (window x y half))
         Ix (mapv (fn [[px py]] (img/at gx px py)) pts)
         Iy (mapv (fn [[px py]] (img/at gy px py)) pts)
         sxx (reduce + (map * Ix Ix)) sxy (reduce + (map * Ix Iy)) syy (reduce + (map * Iy Iy))
         det (- (* sxx syy) (* sxy sxy))]
     (when (> det 1e-12)
       (let [inverse [[(/ syy det) (/ (- sxy) det)] [(/ (- sxy) det) (/ sxx det)]]]
         (iterate (partial lucas-kanade-step a b pts Ix Iy inverse) [0.0 0.0]))))))

(defn lucas-kanade
  "Lucas-Kanade's flow at `(x, y)` after `steps` refinements (default
  10), `a`'s gradients `grads` computed if not given."
  ([a b x y half] (lucas-kanade a b x y half 10 (img/gradients a)))
  ([a b x y half steps grads]
   (some-> (lucas-kanade-seq a b x y half grads) (nth steps))))

(defn flow-field
  "The flow on a grid every `spacing` pixels, `margin` in from the edges,
  by `method` -- `:block` or `:lucas-kanade` -- `{[x y] [u v]}`."
  [a b spacing margin half search method]
  (let [grads (when (= method :lucas-kanade) (img/gradients a))]
    (into {} (for [y (range margin (- (:h a) margin) spacing)
                   x (range margin (- (:w a) margin) spacing)
                   :let [f (case method
                             :block (block-match a b x y half search)
                             :lucas-kanade (lucas-kanade a b x y half 10 grads))]
                   :when f]
               [[x y] f]))))

;; ---------------------------------------------------------------------------
;; The flow of a moving camera

(defn translational-flow
  "The optical flow at image point `(x, y)` -- focal length 1 -- of a
  scene point at depth `z`, for a camera translating with velocity `[tx
  ty tz]`: ((-Tx + x Tz) / Z, (-Ty + y Tz) / Z), the chapter's equations."
  [[tx ty tz] x y z]
  [(/ (+ (- tx) (* x tz)) z) (/ (+ (- ty) (* y tz)) z)])

(defn focus-of-expansion
  "The point the flow `field`, `{[x y] [u v]}`, of a translating camera
  radiates from -- where the flow is zero, (Tx / Tz, Ty / Tz). Every flow
  vector lies on a line through it, so it is the point nearest all those
  lines by least squares: for each, v x0 - u y0 = v x - u y."
  [field]
  (let [rows (map (fn [[[x y] [u v]]] [[v (- u)] (- (* v x) (* u y))]) field)
        a (reduce (fn [m [[p q]]] (lin/mat-add m [[(* p p) (* p q)] [(* p q) (* q q)]])) [[0.0 0.0] [0.0 0.0]] rows)
        b (reduce (fn [acc [[p q] r]] (lin/add acc [(* p r) (* q r)])) [0.0 0.0] rows)
        [[a00 a01] [_ a11]] a
        det (- (* a00 a11) (* a01 a01))]
    (when (> (abs det) 1e-12)
      [(/ (- (* a11 (b 0)) (* a01 (b 1))) det)
       (/ (- (* a00 (b 1)) (* a01 (b 0))) det)])))

(defn time-to-contact
  "Z / Tz, the time until the camera reaches the depth of the point at
  `[x y]` whose flow is `[u v]`, measured from the focus of expansion
  `foe`: the point's distance from the focus over its speed away from it.
  Neither the distance nor the speed is known alone -- the flow cannot
  tell a near slow approach from a far fast one -- but their ratio is:
  the cue a fly lands by. In the flow's time unit, frames for a flow
  between frames."
  [[x0 y0] [x y] [u v]]
  (let [dx (- x x0) dy (- y y0)
        radial (/ (+ (* u dx) (* v dy)) (math/sqrt (+ (* dx dx) (* dy dy))))]
    (/ (math/sqrt (+ (* dx dx) (* dy dy))) radial)))

(defn depth-ratio
  "Z1 / Z2 for two points seen by a camera moving sideways (Tz = 0), from
  their flows: the inverse of the ratio of the flows' lengths -- motion
  parallax, the far landscape crawling past a train window."
  [[u1 v1] [u2 v2]]
  (/ (math/sqrt (+ (* u2 u2) (* v2 v2))) (math/sqrt (+ (* u1 u1) (* v1 v1)))))
