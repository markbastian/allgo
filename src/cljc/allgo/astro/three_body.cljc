(ns allgo.astro.three-body
  "Three bodies under their mutual gravity, in units where G = 1.

  Two bodies are solved: Kepler did it, and `allgo.astro.kepler` is the
  closed form. Add a third and there is no closed form at all -- Poincare
  showed there is not even a full set of conserved quantities to find one
  with -- so everything here is numerical integration, and the problem is
  the classic stress test for an integrator. The motion is chaotic, so an
  error is amplified rather than averaged away, and close encounters ask
  for steps a thousand times smaller than the quiet stretches between them,
  which a fixed step cannot give without wasting nearly all of its work.

  The state is flat: positions `[x0 y0 z0 x1 y1 z1 x2 y2 z2]` and
  velocities laid out the same way. `acceleration` is the second-order
  system the Nystrom methods integrate directly; `derivative` is the same
  thing as eighteen first-order equations for the Runge-Kutta family.

  The starting points in `scenarios` are the famous ones: the few special
  solutions that are known, and the problem that showed the general case is
  hopeless."
  (:require [allgo.array :as a]
            [allgo.geometry.vec3 :as v3]
            [allgo.numerics :as num]
            [allgo.random :as random]
            [clojure.math :as math]))

;; ---------------------------------------------------------------- the system

(defn body
  "Body `i`'s three components out of a flat state."
  [flat i]
  (let [k (* 3 (long i))] (subvec flat k (+ k 3))))

(defn bodies
  "A flat state as one `[x y z]` per body."
  [flat]
  (mapv #(body flat %) (range (quot (count flat) 3))))

(defn- pull
  "The acceleration of every body, flat, from positions read out of the
  first `3n` entries of `y` -- which is the whole of a position vector, or
  the front half of a stacked state, so neither has to be sliced first.

  Written as loops over a double array rather than as a sequence pipeline:
  it runs seven times a step, and in a browser the step before the JIT has
  warmed up to it costs whatever the lazy sequences cost."
  ^doubles [^doubles ms y]
  (let [n   (alength ms)
        out (a/f64 (* 3 n))]
    (dotimes [i n]
      (loop [j (inc i)]
        (when (< j n)
          (let [i3 (* 3 i) j3 (* 3 j)
                dx (- (double (nth y j3)) (double (nth y i3)))
                dy (- (double (nth y (+ j3 1))) (double (nth y (+ i3 1))))
                dz (- (double (nth y (+ j3 2))) (double (nth y (+ i3 2))))
                r2 (+ (* dx dx) (* dy dy) (* dz dz))
                k  (/ 1.0 (* r2 (math/sqrt r2)))
                ki (* k (aget ms j))
                kj (* k (aget ms i))]
            ;; equal and opposite, so each pair is visited once
            (aset out i3 (+ (aget out i3) (* ki dx)))
            (aset out (+ i3 1) (+ (aget out (+ i3 1)) (* ki dy)))
            (aset out (+ i3 2) (+ (aget out (+ i3 2)) (* ki dz)))
            (aset out j3 (- (aget out j3) (* kj dx)))
            (aset out (+ j3 1) (- (aget out (+ j3 1)) (* kj dy)))
            (aset out (+ j3 2) (- (aget out (+ j3 2)) (* kj dz))))
          (recur (inc j)))))
    out))

(defn acceleration
  "y'' = (f t y) for bodies of `masses`: each is pulled toward every other
  by m / r^2."
  [masses]
  (let [ms (a/f64 masses)]
    (fn [_ pos] (vec (pull ms pos)))))

(defn derivative
  "The same system as y' = (f t y), with positions and velocities stacked
  in one vector, for the methods that only integrate first-order systems."
  [masses]
  (let [ms (a/f64 masses)]
    (fn [_ y]
      (let [half (quot (count y) 2)
            ^doubles acc (pull ms y)]
        (persistent!
         (reduce (fn [v x] (conj! v x))
                 (reduce (fn [v k] (conj! v (nth y k))) (transient []) (range half (* 2 half)))
                 acc))))))

;; ------------------------------------------------ what should stay the same

(defn energy
  "Kinetic plus potential. Conserved exactly, so its drift is the
  integrator's error made visible."
  [masses pos vel]
  (let [ms (vec masses) ps (bodies pos) vs (bodies vel) n (count ms)]
    (- (reduce + (map (fn [m v] (* 0.5 m (v3/length-squared v))) ms vs))
       (reduce + (for [i (range n) j (range (inc i) n)]
                   (/ (* (ms i) (ms j)) (v3/distance (ps i) (ps j))))))))

(defn momentum [masses vel]
  (reduce v3/add v3/zero (map v3/scale (bodies vel) masses)))

(defn angular-momentum [masses pos vel]
  (reduce v3/add v3/zero (map (fn [m p v] (v3/scale (v3/cross p v) m))
                              masses (bodies pos) (bodies vel))))

(defn center-of-mass [masses pos]
  (v3/scale (reduce v3/add v3/zero (map v3/scale (bodies pos) masses))
            (/ 1.0 (reduce + masses))))

(defn centered
  "The same motion seen from the center of mass: its position and velocity
  subtracted from every body, so the system stays put on screen instead of
  drifting off it."
  [{:keys [masses pos vel] :as scenario}]
  (let [shift (fn [flat c] (vec (mapcat #(v3/sub % c) (bodies flat))))
        m     (reduce + masses)]
    (assoc scenario
           :pos (shift pos (center-of-mass masses pos))
           :vel (shift vel (v3/scale (momentum masses vel) (/ 1.0 m))))))

;; --------------------------------------------------------------- scenarios

(defn- flat [& vs] (vec (mapcat (fn [[x y]] [(double x) (double y) 0.0]) vs)))

(defn- lagrange
  "Three equal masses at the corners of an equilateral triangle, turning
  as a rigid body (Lagrange, 1772). Each is pulled toward the center by
  the other two with sqrt(3) m / s^2, which a circle of radius R = s /
  sqrt(3) balances when w^2 = m / (sqrt(3) R^3).

  An exact solution, and an unstable one: with equal masses the triangle
  holds only while nothing disturbs it, and rounding error is enough.
  Stability needs one mass to dominate, as the Sun does the Trojan
  asteroids at Jupiter's L4 and L5 points."
  []
  (let [w (math/pow 3.0 -0.25)
        corners (map #(math/to-radians (+ 90.0 (* 120.0 %))) (range 3))]
    {:masses [1.0 1.0 1.0]
     :pos (apply flat (map (fn [a] [(math/cos a) (math/sin a)]) corners))
     :vel (apply flat (map (fn [a] [(* w (- (math/sin a))) (* w (math/cos a))]) corners))
     :period (/ (* 2.0 math/PI) w)}))

(defn- hierarchical
  "A close binary with a light third body circling far outside it, the
  shape of most triple stars and of a planet with a moon. Treated as two
  nested two-body problems it is nearly one, which is why it is stable."
  []
  (let [m3 0.05 d 6.0
        v-bin (/ (math/sqrt 2.0) 2.0)
        v-out (math/sqrt (/ (+ 2.0 m3) d))]
    (centered {:masses [1.0 1.0 m3]
               :pos (flat [-0.5 0.0] [0.5 0.0] [d 0.0])
               :vel (flat [0.0 (- v-bin)] [0.0 v-bin] [0.0 v-out])})))

(defn free-fall
  "Three equal masses dropped from rest at random points of the unit disk.
  Almost every such start ends the same way, in a close encounter that
  flings one body out and leaves the other two bound; `seed` picks which
  start, repeatably."
  [seed]
  (let [rng (random/rng seed)
        point (fn [] (let [r (math/sqrt (rng)) a (random/uniform rng (* 2.0 math/PI))]
                       [(* r (math/cos a)) (* r (math/sin a))]))]
    (centered {:masses [1.0 1.0 1.0]
               :pos (flat (point) (point) (point))
               :vel (flat [0 0] [0 0] [0 0])})))

(def scenarios
  "The named starting points, each `{:id :name :masses :pos :vel}` and
  sometimes `:period`, with `:view` the radius worth framing."
  [(merge {:id :figure-eight :name "Figure eight" :view 1.3
           ;; Chenciner & Montgomery (2000), found numerically by Moore
           ;; (1993): three equal masses chasing each other round one
           ;; figure-eight curve, a third of a period apart.
           :masses [1.0 1.0 1.0]
           :pos (flat [0.97000436 -0.24308753] [-0.97000436 0.24308753] [0 0])
           :vel (flat [0.466203685 0.43236573] [0.466203685 0.43236573]
                      [-0.93240737 -0.86473146])
           :period 6.32591398})
   (merge {:id :lagrange :name "Lagrange triangle" :view 1.5} (lagrange))
   (merge {:id :butterfly :name "Butterfly I" :view 1.3
           ;; Suvakov & Dmitrasinovic (2013), one of the thirteen new
           ;; periodic families they found by searching numerically.
           :masses [1.0 1.0 1.0]
           :pos (flat [-1 0] [1 0] [0 0])
           :vel (flat [0.306893 0.125507] [0.306893 0.125507] [-0.613786 -0.251014])
           :period 6.235641})
   (merge {:id :pythagorean :name "Pythagorean (Burrau)" :view 4.0
           ;; Burrau (1913): masses 3, 4 and 5 at rest at the corners of a
           ;; 3-4-5 right triangle, each opposite the side of its own
           ;; length. Szebehely & Peters (1967) integrated it to the end:
           ;; some sixty time units of close encounters, then the lightest
           ;; body thrown out and the other two leaving as a binary. When
           ;; is chaotic -- tighten the tolerance and it moves -- but that
           ;; is how it ends.
           :masses [3.0 4.0 5.0]
           :pos (flat [1 3] [-2 -1] [1 -1])
           :vel (flat [0 0] [0 0] [0 0])})
   (merge {:id :hierarchical :name "Binary and outer body" :view 7.0} (hierarchical))
   (merge {:id :free-fall :name "Random free fall" :view 2.0} (free-fall 1))])

(def by-id (into {} (map (juxt :id identity)) scenarios))

;; ------------------------------------------------------------- integrating

(def integrators
  "The integrators offered, and how each is built. All but Verlet choose
  their own step, which is the only way through a close encounter; Verlet
  keeps the step it is given, and shows what that costs."
  [{:name "DOPRI5(4)" :order 1}
   {:name "RKF4(5)"   :order 1}
   {:name "RKN4"      :order 2 :adaptive? true}
   {:name "GBS8-2"    :order 2 :adaptive? true}
   {:name "Verlet"    :order 2 :h 1e-3}])

(def ^:private first-order-by-name (into {} (map (juxt :name identity)) num/first-order))
(def ^:private second-order-by-name (into {} (map (juxt :name identity)) num/second-order))

(defn integrator
  "An integrator for `scenario` by the method named `method`, holding a
  relative and absolute error of `tol` per step where it can adapt, and
  never stepping more than `h-max` so the path it leaves stays smooth."
  ([scenario method] (integrator scenario method 1e-10 0.02))
  ([{:keys [masses pos vel]} method tol h-max]
   (let [{:keys [order adaptive? h]} (some #(when (= method (:name %)) %) integrators)
         opts {:tol-abs tol :tol-rel tol :h-max h-max :adaptive? (boolean adaptive?)}]
     (if (= order 1)
       (num/integrator (first-order-by-name method) (derivative masses)
                       0.0 (into pos vel) 1e-3 (dissoc opts :adaptive?))
       (num/integrator-2 (second-order-by-name method) (acceleration masses)
                         0.0 pos vel (or h 1e-3) opts)))))

(defn positions
  "The positions an integrator is at, whichever order of system it is."
  [integ]
  (if (:dy integ) (:y integ) (subvec (:y integ) 0 (quot (count (:y integ)) 2))))

(defn velocities [integ]
  (or (:dy integ) (subvec (:y integ) (quot (count (:y integ)) 2))))

(defn advance
  "Steps `integ` to `t-end`, or until `max-steps` are spent, and returns
  `[integ path]` -- `path` the positions after each accepted step, which
  is where a trail comes from. The last step is clipped to land on
  `t-end`, without shrinking the stride the next call starts from."
  [integ t-end max-steps]
  (loop [s integ path (transient []) n 0]
    (let [remaining (- (double t-end) (:t s))]
      (if (or (<= remaining 1e-12) (>= n (long max-steps)))
        [s (persistent! path)]
        (let [natural (:h s)
              clipped (min natural remaining)
              s'      (num/step (assoc s :h clipped))
              s'      (if (and (:accepted s') (< clipped natural)) (assoc s' :h natural) s')]
          (recur s' (if (:accepted s') (conj! path (positions s')) path) (inc n)))))))

(defn nudge
  "`scenario` with its first body moved `dx` along x: the same system, a
  hair apart, to watch the two part company."
  [scenario dx]
  (update-in scenario [:pos 0] + dx))
