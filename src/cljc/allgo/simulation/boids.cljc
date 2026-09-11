(ns allgo.simulation.boids
  "Reynolds' boids: emergent flocking from three local steering rules.

  A boid is `{:pos [..] :vel [..]}`; positions and velocities are vectors of
  any arity, so the same rules drive 2D and 3D flocks. Steering follows
  Reynolds' formulation -- each rule proposes a desired velocity, and the
  force applied is the (clamped) difference between desired and current.

  A flock may also be given `:obstacles` to steer around. Those are convex
  bodies queried through `allgo.geometry.gjk`, so anything with a support
  mapping serves."
  (:require [allgo.geometry.gjk :as gjk]
            [clojure.math :as math]))

(def defaults
  {:separation-radius 24.0
   :perception-radius 60.0
   :separation-weight 1.6
   :alignment-weight  1.0
   :cohesion-weight   0.9
   :max-speed         2.4
   :max-force         0.06
   :edges             :wrap
   :obstacles         []
   :boid-radius       1.0
   :avoid-radius      26.0
   :avoid-weight      2.2
   :restitution       0.85})

;; Vector math, arity-generic.

(defn v+ [a b] (mapv + a b))
(defn v- [a b] (mapv - a b))
(defn v* [v s] (mapv #(* % s) v))

(defn dot [a b] (reduce + (map * a b)))
(defn mag-sq [v] (reduce (fn [acc x] (+ acc (* x x))) 0 v))
(defn mag [v] (math/sqrt (mag-sq v)))

(defn normalize [v]
  (let [m (mag v)]
    (if (zero? m) v (v* v (/ m)))))

(defn limit [v max-mag]
  (if (> (mag-sq v) (* max-mag max-mag))
    (v* (normalize v) max-mag)
    v))

(defn with-magnitude [v m]
  (v* (normalize v) m))

(defn dist-sq [a b] (mag-sq (v- a b)))
(defn dist [a b] (mag (v- a b)))

;; Steering rules. Each returns a force, or a zero vector when the boid has
;; nobody to respond to.

(defn- zero-like [v] (vec (repeat (count v) 0)))

(defn- steer-toward
  "Reynolds steering: the clamped correction that turns `vel` into a
  full-speed run along `desired`."
  [desired vel {:keys [max-speed max-force]}]
  (if (zero? (mag-sq desired))
    (zero-like vel)
    (limit (v- (with-magnitude desired max-speed) vel) max-force)))

(defn separation
  "Steer away from neighbours, weighting each by inverse square distance so
  the nearest crowding dominates."
  [{:keys [pos vel]} neighbours {:keys [separation-radius] :as params}]
  (let [r2      (* separation-radius separation-radius)
        close   (filter #(< 0 (dist-sq pos (:pos %)) r2) neighbours)
        desired (reduce (fn [acc {other :pos}]
                          (let [away (v- pos other)]
                            (v+ acc (v* (normalize away) (/ (mag-sq away))))))
                        (zero-like pos)
                        close)]
    (steer-toward desired vel params)))

(defn alignment
  "Steer toward the average heading of the neighbourhood."
  [{:keys [vel]} neighbours params]
  (if (seq neighbours)
    (steer-toward (v* (reduce v+ (map :vel neighbours)) (/ (count neighbours)))
                  vel params)
    (zero-like vel)))

(defn cohesion
  "Steer toward the centre of mass of the neighbourhood."
  [{:keys [pos vel]} neighbours params]
  (if (seq neighbours)
    (steer-toward (v- (v* (reduce v+ (map :pos neighbours)) (/ (count neighbours)))
                      pos)
                  vel params)
    (zero-like vel)))

(defn acceleration
  "The weighted sum of the three rules -- the total steering force on a boid."
  [boid neighbours {:keys [separation-weight alignment-weight cohesion-weight]
                    :as params}]
  (-> (v* (separation boid neighbours params) separation-weight)
      (v+ (v* (alignment boid neighbours params) alignment-weight))
      (v+ (v* (cohesion boid neighbours params) cohesion-weight))))

;; Spatial index. A linear scan makes each tick O(n^2), which stops holding a
;; 60fps frame budget at a few hundred boids. Bucketing by a cell the size of
;; the perception radius means a boid only ever tests the 3^d cells around it.

(def ^:private cell-offsets
  (memoize (fn offsets [d]
             (if (zero? d)
               [[]]
               (vec (for [tail (offsets (dec d)) o [-1 0 1]] (conj tail o)))))))

(defn- cell-of [pos cell-size]
  (mapv #(long (math/floor (/ % cell-size))) pos))

(defn index-flock
  "Bucket the flock into cells of `cell-size` for neighbourhood lookup."
  [flock cell-size]
  {:cell-size cell-size
   :cells     (group-by #(cell-of (:pos %) cell-size) flock)})

(defn neighbours
  "Every other boid within `:perception-radius`, via a `index-flock` index."
  [{:keys [pos]} {:keys [cell-size cells]} {:keys [perception-radius]}]
  (let [r2   (* perception-radius perception-radius)
        home (cell-of pos cell-size)]
    (into []
          (comp (mapcat #(cells (mapv + home %)))
                (filter #(< 0 (dist-sq pos (:pos %)) r2)))
          (cell-offsets (count pos)))))

(defn wrap
  "Toroidal world: a boid leaving one edge re-enters at the opposite one."
  [pos bounds]
  (mapv (fn [x hi] (cond (neg? x) (+ x hi) (>= x hi) (- x hi) :else x))
        pos bounds))

(defn bounce
  "Reflect a boid off the walls of the world instead of wrapping through
  them, mirroring both the offending coordinate and its velocity."
  [pos vel bounds]
  (let [pairs (mapv (fn [x v hi]
                      (cond (neg? x)  [(- x) (- v)]
                            (> x hi)  [(- (* 2.0 hi) x) (- v)]
                            :else     [x v]))
                    pos vel bounds)]
    [(mapv first pairs) (mapv second pairs)]))

;; Obstacles. GJK is three-dimensional, so a planar flock is lifted onto
;; z = 0 for the query and the answer projected back. Obstacles must be
;; solid rather than flat -- a body with no volume drives GJK's simplex
;; degenerate -- which is what the constructors below take care of.

(defn- lift [v]
  (if (= 3 (count v)) (vec v) [(nth v 0) (nth v 1) 0.0]))

(defn- project [v n]
  (if (= n 3) (vec v) (subvec (vec v) 0 2)))

(defn sphere-obstacle
  "A ball obstacle. Its equator is the circle a planar flock meets."
  [centre radius]
  (let [c (lift centre)]
    {:support (gjk/sphere c radius) :centre c :radius radius}))

(defn box-obstacle
  "A box obstacle. Two-dimensional corners become a prism deep enough that
  the flock's plane cuts it squarely rather than grazing a flat face."
  [lo hi]
  (let [flat? (= 2 (count lo))
        depth (if flat? (max 1.0 (mag (v- (vec hi) (vec lo)))) 0.0)
        lo3   (if flat? [(nth lo 0) (nth lo 1) (- depth)] (vec lo))
        hi3   (if flat? [(nth hi 0) (nth hi 1) depth] (vec hi))]
    {:support (gjk/box lo3 hi3)
     :centre  (v* (v+ lo3 hi3) 0.5)
     ;; The bounding sphere measures the body as the flock actually meets
     ;; it. For a flat box that is the rectangle, not the prism: the depth
     ;; is invented purely to keep GJK's simplex off a degenerate face, and
     ;; the prism is centred on the flock's plane, so no boid can be near
     ;; it in z. Including that depth inflates the radius by more than a
     ;; factor of two, which defeats the cull that is supposed to keep GJK
     ;; off obstacles it will never touch.
     :radius  (* 0.5 (mag (v- (vec hi) (vec lo))))}))

(defn index-obstacles
  "Bucket obstacles into cells of `cell-size` for reach lookup, the same way
  `index-flock` buckets boids.

  Bucketing is on x and y only. A planar flock's obstacles are prisms
  spanning z, so indexing depth would put every one of them in every layer
  and buy nothing; for a spatial flock it is a coarser cull, but still a
  sound one, because the bounding-sphere test runs afterwards regardless.

  An obstacle wider than a cell lands in each cell it covers, so lookups
  have to drop duplicates."
  [obstacles cell-size]
  {:cell-size cell-size
   :cells     (reduce (fn [cells {:keys [centre radius] :as obstacle}]
                        (let [[cx cy] centre
                              lo (fn [v] (long (math/floor (/ (- v radius) cell-size))))
                              hi (fn [v] (long (math/floor (/ (+ v radius) cell-size))))]
                          (reduce (fn [cells cell] (update cells cell conj obstacle))
                                  cells
                                  (for [i (range (lo cx) (inc (hi cx)))
                                        j (range (lo cy) (inc (hi cy)))]
                                    [i j]))))
                      {}
                      obstacles)})

(defn- near-obstacles
  "Obstacles whose bounding sphere is within reach, so the GJK query is only
  run against the few that could matter.

  `step` supplies an index; without one this falls back to a linear scan,
  which keeps `avoidance` usable on its own. The scan is what makes the
  tick O(boids x obstacles) -- fine for the handful of bodies in the 2D
  flocking demo, but a flock loose in a generated dungeon meets a wall
  segment per floor tile, and at a few hundred of those it dominates the
  frame."
  [p3 {:keys [obstacles obstacle-index]} reach]
  (let [in-reach? (fn [{:keys [centre radius]}] (< (dist p3 centre) (+ radius reach)))]
    (if-let [{:keys [cell-size cells]} obstacle-index]
      (let [home (cell-of (subvec p3 0 2) cell-size)]
        (into [] (comp (mapcat #(cells (mapv + home %)))
                       (distinct)
                       (filter in-reach?))
              (cell-offsets 2)))
      (filter in-reach? obstacles))))

(defn avoidance
  "Steering away from nearby obstacles. GJK gives the gap to each convex
  body and the direction across it; the force ramps up as the gap closes,
  so a boid curves around an obstacle rather than jolting at contact."
  [{:keys [pos vel]} {:keys [obstacles avoid-radius boid-radius max-speed max-force] :as params}]
  (if (empty? obstacles)
    (zero-like pos)
    (let [n     (count pos)
          p3    (lift pos)
          me    (gjk/sphere p3 boid-radius)
          v3    (lift vel)]
      (project
       (reduce (fn [acc {:keys [support]}]
                 (let [{:keys [distance direction]} (gjk/distance me support)]
                   (if (>= distance avoid-radius)
                     acc
                     (let [away   (v* direction -1.0)
                           urgency (- 1.0 (/ distance avoid-radius))
                           steer  (limit (v- (with-magnitude away max-speed) v3) max-force)]
                       (v+ acc (v* steer urgency))))))
               [0.0 0.0 0.0]
               (near-obstacles p3 params (+ avoid-radius boid-radius)))
       n))))

(defn resolve-contacts
  "Push a boid clear of anything it has ended up inside and reflect it off
  the surface. Steering alone cannot guarantee separation -- a boid boxed in
  by its flock can be carried into an obstacle -- so EPA supplies the
  shortest way out as a backstop."
  [pos vel {:keys [obstacles boid-radius restitution] :as params}]
  (if (empty? obstacles)
    [pos vel]
    (let [n (count pos)]
      (loop [p3 (lift pos), v3 (lift vel), [o & more] (near-obstacles (lift pos) params boid-radius)]
        (if (nil? o)
          [(project p3 n) (project v3 n)]
          (let [me (gjk/sphere p3 boid-radius)]
            (if-let [{:keys [depth normal]} (and (gjk/intersects? me (:support o))
                                                 (gjk/penetration me (:support o)))]
              ;; EPA reports the shift that moves the obstacle clear, so the
              ;; boid takes the opposite; the surface normal it bounces off
              ;; points the same way.
              (let [out (v* normal -1.0)
                    vn  (dot v3 out)]
                (recur (v+ p3 (v* out depth))
                       (if (neg? vn) (v- v3 (v* out (* (+ 1.0 restitution) vn))) v3)
                       more))
              (recur p3 v3 more))))))))

(defn step-boid [boid index bounds {:keys [max-speed edges avoid-weight] :as params}]
  (let [acc       (v+ (acceleration boid (neighbours boid index params) params)
                      (v* (avoidance boid params) avoid-weight))
        vel       (limit (v+ (:vel boid) acc) max-speed)
        moved     (v+ (:pos boid) vel)
        [pos vel] (if (= :bounce edges)
                    (bounce moved vel bounds)
                    [(wrap moved bounds) vel])
        [pos vel] (resolve-contacts pos vel params)]
    (assoc boid :pos pos :vel vel)))

(defn step
  "Advance the whole flock one tick. Every boid sees the same previous state,
  so the update order cannot bias the result."
  ([flock bounds] (step flock bounds defaults))
  ([flock bounds params]
   (let [{:keys [perception-radius avoid-radius boid-radius obstacles] :as params}
         (merge defaults params)
         index  (index-flock flock perception-radius)
         ;; Built once per tick, not once per boid: the obstacles do not
         ;; move during a step.
         params (cond-> params
                  (seq obstacles)
                  (assoc :obstacle-index
                         (index-obstacles obstacles (+ avoid-radius boid-radius))))]
     (mapv #(step-boid % index bounds params) flock))))

(defn random-boid [bounds max-speed]
  {:pos (mapv rand bounds)
   :vel (with-magnitude (mapv (fn [_] (- (rand) 0.5)) bounds)
          (* max-speed (+ 0.5 (rand 0.5))))})

(defn flock
  "`n` boids at random positions and headings within `bounds`."
  ([n bounds] (flock n bounds defaults))
  ([n bounds params]
   (let [{:keys [max-speed]} (merge defaults params)]
     (vec (repeatedly n #(random-boid bounds max-speed))))))
