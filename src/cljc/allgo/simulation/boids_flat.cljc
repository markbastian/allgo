(ns allgo.simulation.boids-flat
  "The flocking rules of `allgo.simulation.boids`, on flat arrays.

  Same three rules, same parameters, same results -- `boids-flat-test`
  checks the two agree step for step. What differs is the representation.
  The reference keeps a vector of `{:pos [x y] :vel [x y]}` maps and
  rebuilds a map of cells every tick; this keeps four arrays of doubles
  and rebuilds an `allgo.spatial.hash` that allocates nothing. For a few
  dozen boids the reference is the better code. For a few thousand it is
  the only thing between you and a frame budget.

  Two things beyond the representation actually matter:

  The reference walks each boid's neighbours three times, once per rule,
  because that reads clearly. Here all three sums are accumulated in one
  pass, so the neighbour list is traversed once instead of three times.

  Positions are double buffered. Every boid must see the state as it was
  at the start of the tick -- update in place and a boid steers using its
  neighbours' *new* positions, which biases the flock by iteration order.
  The reference gets this for free from persistent data; an array has to
  be explicit about it.

  Positions carry three components even for a planar flock, with z held at
  zero, so one code path serves both and the spatial hash can be used
  unchanged."
  (:require [allgo.simulation.boids :as boids]
            [allgo.simulation.flock :as flock]
            [allgo.spatial.hash :as spatial]
            [clojure.math :as math]))

;; `step` and `to-boids` are defined below but referred to by the record's
;; protocol methods, which are compiled where the record is.
(declare step to-boids)

(defn- f64
  ([n] #?(:clj (double-array n) :cljs (js/Float64Array. n)))
  ([_n coll] #?(:clj (double-array (map double coll))
                :cljs (js/Float64Array. (into-array (map double coll))))))

;; ---------------------------------------------------------------------------
;; Construction

(defrecord FlatFlock [n dims pos vel pos' vel' hash]
  ;; Implemented here rather than extended from `allgo.simulation.flock`.
  ;; Extending it there would mean naming this class from another
  ;; namespace, and a `defrecord` class only exists once its namespace has
  ;; been loaded -- which `:import` does not do, making it order-dependent
  ;; and invisible to static analysis.
  flock/Flock
  (advance [this bounds params] (step this bounds params))
  (as-boids [this] (to-boids this))
  (flock-size [this] (:n this)))

(defn flat-flock
  "A flock of `n` boids from flat `[x y z ...]` position and velocity
  arrays. `dims` is 2 or 3 and decides whether z is simulated.

  A record rather than a map so it can carry the `allgo.simulation.flock`
  protocol without claiming every map in the program."
  [n pos vel dims]
  (map->FlatFlock
   {:n    n
    :dims dims
    :pos  pos
    :vel  vel
    ;; The buffers written during a tick, swapped in at the end.
    :pos' (f64 (* 3 n))
    :vel' (f64 (* 3 n))
    :hash (spatial/spatial-hash 1.0 (max n 1))}))

(defn from-boids
  "A flat flock holding the same state as a `allgo.simulation.boids`
  flock."
  [flock]
  (let [n    (count flock)
        dims (count (:pos (first flock)))
        ^doubles pos (f64 (* 3 n))
        ^doubles vel (f64 (* 3 n))]
    (dotimes [i n]
      (let [{p :pos v :vel} (nth flock i)
            b (* 3 i)]
        (dotimes [d dims]
          (aset pos (+ b d) (double (nth p d)))
          (aset vel (+ b d) (double (nth v d))))))
    (flat-flock n pos vel dims)))

(defn to-boids
  "The same state as a vector of `{:pos :vel}`, for comparison or for the
  renderers that expect it."
  [{:keys [n dims ^doubles pos ^doubles vel]}]
  (mapv (fn [i]
          (let [b (* 3 i)]
            {:pos (mapv #(aget pos (+ b %)) (range dims))
             :vel (mapv #(aget vel (+ b %)) (range dims))}))
        (range n)))

(defn positions
  "Positions as one flat array, which is what a renderer wants."
  [{:keys [pos]}]
  pos)

;; ---------------------------------------------------------------------------
;; The step

(defn- limit!
  "Scales `[x y z]` down to `max-mag` if it is longer, in place in a
  scratch array."
  [^doubles v max-mag]
  (let [x (aget v 0) y (aget v 1) z (aget v 2)
        m2 (+ (* x x) (* y y) (* z z))]
    (when (> m2 (* max-mag max-mag))
      (let [s (/ max-mag (math/sqrt m2))]
        (aset v 0 (* x s))
        (aset v 1 (* y s))
        (aset v 2 (* z s))))))

(defn- steer!
  "Reynolds steering, written into `out`: the clamped correction that
  turns the current velocity into a full-speed run along `desired`."
  [^doubles out dx dy dz vx vy vz max-speed max-force]
  (let [m2 (+ (* dx dx) (* dy dy) (* dz dz))]
    (if (zero? m2)
      (do (aset out 0 0.0) (aset out 1 0.0) (aset out 2 0.0))
      (let [s (/ max-speed (math/sqrt m2))]
        (aset out 0 (- (* dx s) vx))
        (aset out 1 (- (* dy s) vy))
        (aset out 2 (- (* dz s) vz))
        (limit! out max-force)))))

(defn step
  "Advance the flock one tick. Same parameters as `allgo.simulation.boids/step`,
  obstacles included.

  Obstacle avoidance calls the reference implementation's own GJK code
  rather than repeating it. That means building a small map per boid, but
  only when there are obstacles at all, and a GJK query costs far more
  than the map does -- while sharing the code is what guarantees the two
  implementations steer around a barrier identically."
  ([flock bounds] (step flock bounds boids/defaults))
  ([{:keys [n dims ^doubles pos ^doubles vel ^doubles pos' ^doubles vel' hash] :as flock}
    bounds params]
   (let [{:keys [separation-radius perception-radius separation-weight
                 alignment-weight cohesion-weight max-speed max-force edges
                 obstacles avoid-weight]
          :as params}
         (merge boids/defaults params)
         avoiding? (seq obstacles)
         ;; The same index `boids/step` builds. Not only for speed: it
         ;; also fixes the order obstacles are summed in, and a different
         ;; order gives a different rounding, which is enough to make the
         ;; two implementations disagree in the sixth decimal place from
         ;; the very first tick.
         params    (cond-> params
                     avoiding?
                     (assoc :obstacle-index
                            (boids/index-obstacles
                             obstacles
                             (+ (:avoid-radius params) (:boid-radius params)))))
         sep-r2  (* separation-radius separation-radius)
         per-r2  (* perception-radius perception-radius)
         hash    (assoc hash :spacing (double perception-radius)
                        :inv-spacing (/ 1.0 (double perception-radius)))
         ;; Hinted, or every access to these is reflective -- they sit in
         ;; the innermost loop, so that alone made this ten times slower
         ;; than the implementation it was meant to beat.
         ^doubles scratch (f64 3)
         ^doubles acc     (f64 3)
         ^doubles sums    (f64 10)
         bx (double (nth bounds 0))
         by (double (nth bounds 1))
         bz (if (= dims 3) (double (nth bounds 2)) 0.0)]
     (spatial/rebuild! hash pos n)
     (dotimes [i n]
       (let [b  (* 3 i)
             px (aget pos b) py (aget pos (+ b 1)) pz (aget pos (+ b 2))
             vx (aget vel b) vy (aget vel (+ b 1)) vz (aget vel (+ b 2))
             found (spatial/query! hash pos i perception-radius)]
         ;; The three rules are summed into a scratch array rather than
         ;; carried as loop arguments. Ten accumulators through `recur`
         ;; box every one of them, and threading `[sx sy sz]` through a
         ;; conditional allocates a vector for every neighbour examined --
         ;; which is how a flat-array rewrite ends up slower than the
         ;; persistent code it was replacing.
         (dotimes [q 10] (aset sums q 0.0))
         (dotimes [k found]
           (let [j (spatial/neighbour hash k)]
             (when (not= j i)
               (let [c  (* 3 j)
                     qx (aget pos c) qy (aget pos (+ c 1)) qz (aget pos (+ c 2))
                     dx (- px qx) dy (- py qy) dz (- pz qz)
                     d2 (+ (* dx dx) (* dy dy) (* dz dz))]
                 (when (and (pos? d2) (< d2 per-r2))
                   (when (< d2 sep-r2)
                     ;; Inverse-square weighting: the unit vector away,
                     ;; over the distance squared.
                     (let [w (/ 1.0 (* (math/sqrt d2) d2))]
                       (aset sums 0 (+ (aget sums 0) (* dx w)))
                       (aset sums 1 (+ (aget sums 1) (* dy w)))
                       (aset sums 2 (+ (aget sums 2) (* dz w)))))
                   (aset sums 3 (+ (aget sums 3) (aget vel c)))
                   (aset sums 4 (+ (aget sums 4) (aget vel (+ c 1))))
                   (aset sums 5 (+ (aget sums 5) (aget vel (+ c 2))))
                   (aset sums 6 (+ (aget sums 6) qx))
                   (aset sums 7 (+ (aget sums 7) qy))
                   (aset sums 8 (+ (aget sums 8) qz))
                   (aset sums 9 (+ (aget sums 9) 1.0)))))))

         (let [near (aget sums 9)]
           (steer! scratch (aget sums 0) (aget sums 1) (aget sums 2)
                   vx vy vz max-speed max-force)
           (aset acc 0 (* (aget scratch 0) separation-weight))
           (aset acc 1 (* (aget scratch 1) separation-weight))
           (aset acc 2 (* (aget scratch 2) separation-weight))
           (when (pos? near)
             (steer! scratch (/ (aget sums 3) near) (/ (aget sums 4) near) (/ (aget sums 5) near)
                     vx vy vz max-speed max-force)
             (aset acc 0 (+ (aget acc 0) (* (aget scratch 0) alignment-weight)))
             (aset acc 1 (+ (aget acc 1) (* (aget scratch 1) alignment-weight)))
             (aset acc 2 (+ (aget acc 2) (* (aget scratch 2) alignment-weight)))
             (steer! scratch (- (/ (aget sums 6) near) px) (- (/ (aget sums 7) near) py)
                     (- (/ (aget sums 8) near) pz)
                     vx vy vz max-speed max-force)
             (aset acc 0 (+ (aget acc 0) (* (aget scratch 0) cohesion-weight)))
             (aset acc 1 (+ (aget acc 1) (* (aget scratch 1) cohesion-weight)))
             (aset acc 2 (+ (aget acc 2) (* (aget scratch 2) cohesion-weight)))))

         (when avoiding?
           (let [steer (boids/avoidance {:pos (if (= dims 3) [px py pz] [px py])
                                         :vel (if (= dims 3) [vx vy vz] [vx vy])}
                                        params)]
             (aset acc 0 (+ (aget acc 0) (* (double (nth steer 0)) avoid-weight)))
             (aset acc 1 (+ (aget acc 1) (* (double (nth steer 1)) avoid-weight)))
             (when (= dims 3)
               (aset acc 2 (+ (aget acc 2) (* (double (nth steer 2)) avoid-weight))))))

         (let [nvx (+ vx (aget acc 0)) nvy (+ vy (aget acc 1)) nvz (+ vz (aget acc 2))
               m2  (+ (* nvx nvx) (* nvy nvy) (* nvz nvz))
               s   (if (> m2 (* max-speed max-speed)) (/ max-speed (math/sqrt m2)) 1.0)
               nvx (* nvx s) nvy (* nvy s) nvz (* nvz s)
               mx  (+ px nvx) my (+ py nvy) mz (+ pz nvz)]
           (if (= :bounce edges)
             ;; Written out per axis rather than through a closure: a
             ;; closure returning [value flipped?] allocates a tuple per
             ;; axis per boid.
             (let [fx (or (neg? mx) (> mx bx))
                   fy (or (neg? my) (> my by))
                   fz (and (= dims 3) (or (neg? mz) (> mz bz)))]
               (aset pos' b (double (cond (neg? mx) (- mx) (> mx bx) (- (* 2.0 bx) mx) :else mx)))
               (aset pos' (+ b 1) (double (cond (neg? my) (- my) (> my by) (- (* 2.0 by) my) :else my)))
               (aset pos' (+ b 2) (double (if (= dims 3)
                                            (cond (neg? mz) (- mz) (> mz bz) (- (* 2.0 bz) mz) :else mz)
                                            mz)))
               (aset vel' b (if fx (- nvx) nvx))
               (aset vel' (+ b 1) (if fy (- nvy) nvy))
               (aset vel' (+ b 2) (if fz (- nvz) nvz)))
             (do
               (aset pos' b (double (cond (neg? mx) (+ mx bx) (>= mx bx) (- mx bx) :else mx)))
               (aset pos' (+ b 1) (double (cond (neg? my) (+ my by) (>= my by) (- my by) :else my)))
               (aset pos' (+ b 2) (double (if (= dims 3)
                                            (cond (neg? mz) (+ mz bz) (>= mz bz) (- mz bz) :else mz)
                                            mz)))
               (aset vel' b nvx)
               (aset vel' (+ b 1) nvy)
               (aset vel' (+ b 2) nvz)))
           ;; Steering alone cannot guarantee separation -- a boid carried
           ;; into an obstacle by its flock is pushed clear by EPA. Same
           ;; call, same place in the order, as the reference.
           (when avoiding?
             (let [three? (= dims 3)
                   [p v] (boids/resolve-contacts
                          (if three?
                            [(aget pos' b) (aget pos' (+ b 1)) (aget pos' (+ b 2))]
                            [(aget pos' b) (aget pos' (+ b 1))])
                          (if three?
                            [(aget vel' b) (aget vel' (+ b 1)) (aget vel' (+ b 2))]
                            [(aget vel' b) (aget vel' (+ b 1))])
                          params)]
               (aset pos' b (double (nth p 0)))
               (aset pos' (+ b 1) (double (nth p 1)))
               (aset vel' b (double (nth v 0)))
               (aset vel' (+ b 1) (double (nth v 1)))
               (when three?
                 (aset pos' (+ b 2) (double (nth p 2)))
                 (aset vel' (+ b 2) (double (nth v 2)))))))))
     ;; Swap the buffers: what was written becomes the state.
     (assoc flock :pos pos' :vel vel' :pos' pos :vel' vel :hash hash))))

(defn simulate
  "Run `ticks` steps, returning the flock."
  [flock bounds params ticks]
  (loop [flock flock i 0]
    (if (= i ticks) flock (recur (step flock bounds params) (inc i)))))
