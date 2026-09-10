(ns procedurals.boids
  "Reynolds' boids: emergent flocking from three local steering rules.

  A boid is `{:pos [..] :vel [..]}`; positions and velocities are vectors of
  any arity, so the same rules drive 2D and 3D flocks. Steering follows
  Reynolds' formulation -- each rule proposes a desired velocity, and the
  force applied is the (clamped) difference between desired and current."
  (:require [clojure.math :as math]))

(def defaults
  {:separation-radius 24.0
   :perception-radius 60.0
   :separation-weight 1.6
   :alignment-weight  1.0
   :cohesion-weight   0.9
   :max-speed         2.4
   :max-force         0.06})

;; Vector math, arity-generic.

(defn v+ [a b] (mapv + a b))
(defn v- [a b] (mapv - a b))
(defn v* [v s] (mapv #(* % s) v))

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

(defn step-boid [boid index bounds {:keys [max-speed] :as params}]
  (let [vel (limit (v+ (:vel boid) (acceleration boid (neighbours boid index params) params))
                   max-speed)]
    (assoc boid :pos (wrap (v+ (:pos boid) vel) bounds) :vel vel)))

(defn step
  "Advance the whole flock one tick. Every boid sees the same previous state,
  so the update order cannot bias the result."
  ([flock bounds] (step flock bounds defaults))
  ([flock bounds params]
   (let [params (merge defaults params)
         index  (index-flock flock (:perception-radius params))]
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
