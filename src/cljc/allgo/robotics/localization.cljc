(ns allgo.robotics.localization
  "Monte Carlo localization: Russell and Norvig, *Artificial Intelligence:
  A Modern Approach*, 4th edition, chapter 26 (MONTE-CARLO-LOCALIZATION,
  figure 26.6), as steps.

  A robot in a known map of walls, with a range scanner and noisy
  wheels, keeps `N` guesses of where it is -- particles, each a pose `[x
  y heading]`. Each step every particle is moved by the motion model
  with its own noise, weighted by how well the scan it would see from
  there -- found by casting rays against the map -- matches the scan the
  robot did see, and the set resampled in proportion to the weights. The
  particles that survive cluster where the robot is.

  A map is a vector of wall segments `[[x0 y0] [x1 y1]]`; `scan-angles`
  are the beams, relative to the heading.

  The chapter's other localizer is the extended Kalman filter, for a map
  of landmarks the robot can tell apart: one Gaussian over the pose,
  carried through the motion model linearized about the estimate, and
  corrected by each landmark's range and bearing linearized the same way.
  `ekf-localization-seq` is that, a step a move; `allgo.robotics.slam`
  builds the map as well, with the same models."
  (:require [allgo.geometry.polygon :as polygon]
            [allgo.numerics.linear :as lin]
            [allgo.probability.temporal :as tm]
            [allgo.random :as random]
            [clojure.math :as math]))

(defn raycast
  "RAYCAST: the distance along the ray from `[x y]` at `angle` to the
  nearest wall of `walls`, or `max-range`."
  [walls [x y] angle max-range]
  (let [end [(+ x (* max-range (math/cos angle))) (+ y (* max-range (math/sin angle)))]]
    (reduce (fn [best [a b]]
              (if-let [[t u] (polygon/segment-intersection [x y] end a b)]
                (if (and (<= 0.0 t 1.0) (<= 0.0 u 1.0)) (min best (* t max-range)) best)
                best))
            max-range walls)))

(defn scan
  "The ranges the scanner sees from `pose`."
  [walls [x y heading] angles max-range]
  (mapv #(raycast walls [x y] (+ heading %) max-range) angles))

(defn drive
  "The motion model without noise: `pose` driven at speed `v` and turn
  rate `w` for `dt` -- turned first, then moved along the new heading."
  [[x y heading] v w dt]
  (let [heading' (+ heading (* w dt))]
    [(+ x (* v dt (math/cos heading'))) (+ y (* v dt (math/sin heading'))) heading']))

(defn move
  "The motion model: `pose` driven at speed `v` and turn rate `w` for
  `dt`, with Gaussian noise of `sigma-v` and `sigma-w` on each."
  [pose v w dt rng sigma-v sigma-w]
  (drive pose (random/gaussian rng v sigma-v) (random/gaussian rng w sigma-w) dt))

(defn- inside? [[[x0 y0] [x1 y1]] [x y]] (and (< x0 x x1) (< y0 y y1)))

(defn- likelihood
  "P(z | z*): each beam Gaussian about the range expected, independently."
  [z z* sigma]
  (reduce * (map (fn [a b] (let [d (- a b)] (math/exp (- (/ (* d d) (* 2 sigma sigma)))))) z z*)))

(defn resample
  "`particles` drawn again in proportion to `weights`, as many as there
  were, by low-variance resampling: one random offset, then evenly spaced
  pointers into the cumulative weights -- less noise than independent
  draws, and a particle of weight w kept about w n times. The weights
  need not sum to one, but must not all be zero."
  [particles weights rng]
  (let [total (double (reduce + weights))
        cumulative (vec (reductions + (map #(/ % total) weights)))
        n (count particles)
        r (/ (double (rng)) n)]
    (loop [i 0 j 0 out (transient [])]
      (if (= i n)
        (persistent! out)
        (let [u (+ r (/ (double i) n))
              j (loop [j j] (if (and (< j (dec n)) (> u (double (cumulative j)))) (recur (inc j)) j))]
          (recur (inc i) j (conj! out (particles j))))))))

(defn mcl-step
  "One MONTE-CARLO-LOCALIZATION update: the particles moved by `[v w]`,
  weighted by the scan `z`, and resampled."
  [{:keys [walls angles max-range dt sigma-v sigma-w sigma-z rng bounds random-fraction]
    :or {random-fraction 0.0}} particles [v w] z]
  (let [moved (mapv #(move % v w dt rng sigma-v sigma-w) particles)
        weights (mapv (fn [p] (if (inside? bounds p) (likelihood z (scan walls p angles max-range) sigma-z) 0.0)) moved)
        total (reduce + weights)]
    (if (pos? total)
      ;; A few particles replaced by random ones: the remedy for sample
      ;; impoverishment, which otherwise leaves every particle a copy of
      ;; one early guess -- right or not.
      (let [[[x0 y0] [x1 y1]] bounds]
        (mapv (fn [p] (if (< (double (rng)) random-fraction)
                        [(+ x0 (* (- x1 x0) (rng))) (+ y0 (* (- y1 y0) (rng))) (* 2 math/PI (rng))]
                        p))
              (resample moved weights rng)))
      ;; Nothing fits: start again from anywhere.
      (let [[[x0 y0] [x1 y1]] bounds]
        (vec (repeatedly (count particles) #(vector (+ x0 (* (- x1 x0) (rng))) (+ y0 (* (- y1 y0) (rng))) (* 2 math/PI (rng)))))))))

(defn mcl-seq
  "Monte Carlo localization of a robot that starts at `start` and is driven
  by `controls`, each `[v w]`: each state `{:truth :particles :scan}`, the
  robot's true pose (moved by the same noisy motion model), the
  particles, and the scan it saw. `n` particles start spread over the
  map's `bounds` at random -- the robot does not know where it is."
  [{:keys [walls angles max-range dt sigma-v sigma-w sigma-z rng bounds] :as world} start controls n]
  (let [[[x0 y0] [x1 y1]] bounds
        particles (vec (repeatedly n #(vector (+ x0 (* (- x1 x0) (rng))) (+ y0 (* (- y1 y0) (rng))) (* 2 math/PI (rng)))))]
    (reductions (fn [{:keys [truth particles]} [v w]]
                  (let [truth' (move truth v w dt rng sigma-v sigma-w)
                        ;; The scanner's own noise.
                        z (mapv #(random/gaussian rng % sigma-z) (scan walls truth' angles max-range))]
                    {:truth truth' :scan z :particles (mcl-step world particles [v w] z)}))
                {:truth start :particles particles :scan nil}
                controls)))

(defn estimate
  "Where the particles say the robot is: the mean of those within a meter
  of the densest spot -- not of them all, since while the evidence is
  ambiguous they form more than one cluster, and the mean of two clusters
  is where neither is -- and the spread about it: `[[x y] spread]`."
  [particles]
  (let [cell (fn [[x y]] [(long (math/floor (* 2 x))) (long (math/floor (* 2 y)))])
        densest (key (apply max-key val (frequencies (map cell particles))))
        [cx cy] (map #(/ (+ % 0.5) 2.0) densest)
        near (filter (fn [[x y]] (< (+ (* (- x cx) (- x cx)) (* (- y cy) (- y cy))) 1.0)) particles)
        n (count near)
        mx (/ (reduce + (map first near)) n)
        my (/ (reduce + (map second near)) n)
        particles near]
    [[mx my] (math/sqrt (/ (reduce + (map (fn [[x y]] (+ (* (- x mx) (- x mx)) (* (- y my) (- y my)))) particles)) n))]))

(defn room
  "A test map: a 10 by 6 room with an L-shaped wall and a pillar inside,
  so that few places look alike."
  []
  (let [box (fn [[x0 y0] [x1 y1]] [[[x0 y0] [x1 y0]] [[x1 y0] [x1 y1]] [[x1 y1] [x0 y1]] [[x0 y1] [x0 y0]]])]
    {:walls (vec (concat (box [0.0 0.0] [10.0 6.0])
                         [[[3.0 0.0] [3.0 2.5]] [[3.0 2.5] [5.0 2.5]]]
                         (box [7.0 3.5] [8.0 4.5])))
     :bounds [[0.0 0.0] [10.0 6.0]]
     :angles (mapv #(* % (/ math/PI 4)) (range 8))
     :max-range 12.0 :dt 0.1 :sigma-v 0.3 :sigma-w 0.2 :sigma-z 0.6 :random-fraction 0.02}))

;; ---------------------------------------------------------------------------
;; EKF localization with landmarks

(defn wrap-angle
  "`a` brought into [-pi, pi)."
  ^double [a]
  (- (mod (+ a math/PI) (* 2 math/PI)) math/PI))

(defn motion-jacobians
  "The motion model's partials at `pose` for control `[v w]`: `[G V]`, G
  with respect to the pose, V to the control."
  [[_ _ heading] [v w] dt]
  (let [h' (+ heading (* w dt)) c (math/cos h') s (math/sin h')]
    [[[1.0 0.0 (* -1.0 v dt s)]
      [0.0 1.0 (* v dt c)]
      [0.0 0.0 1.0]]
     [[(* dt c) (* -1.0 v dt s dt)]
      [(* dt s) (* v dt c dt)]
      [0.0 dt]]]))

(defn range-bearing
  "What a sensor at `pose` sees of a landmark at `[mx my]`: its range and
  its bearing relative to the heading."
  [[x y heading] [mx my]]
  (let [dx (- mx x) dy (- my y)]
    [(math/sqrt (+ (* dx dx) (* dy dy))) (wrap-angle (- (math/atan2 dy dx) heading))]))

(defn range-bearing-jacobian
  "The partials of `range-bearing` with respect to the pose (the first
  three columns) and the landmark (the last two)."
  [[x y _] [mx my]]
  (let [dx (- mx x) dy (- my y) q (+ (* dx dx) (* dy dy)) r (math/sqrt q)]
    [[(- (/ dx r)) (- (/ dy r)) 0.0 (/ dx r) (/ dy r)]
     [(/ dy q) (- (/ dx q)) -1.0 (- (/ dy q)) (/ dx q)]]))

(defn observe-landmarks
  "The noisy range and bearing of each landmark within `max-range` of
  `pose`: `[[id r bearing] ...]`."
  [{:keys [landmarks max-range sigma-r sigma-b]} pose rng]
  (vec (for [[id m] (map-indexed vector landmarks)
             :let [[r b] (range-bearing pose m)]
             :when (<= r max-range)]
         [id (random/gaussian rng r sigma-r) (wrap-angle (random/gaussian rng b sigma-b))])))

(defn ekf-correct
  "The EKF's correction of `[mean cov]` by observations `[[id r bearing]
  ...]`: each measurement row linearized about the predicted mean and
  applied a scalar at a time in Joseph form (`allgo.probability.temporal`'s
  `kalman-update`), which with independent noise is the joint update.
  `predict` gives the expected `[r bearing]` of a landmark and `jacobian`
  its rows over the whole state."
  [[mean cov] observations predict jacobian sigma-r sigma-b]
  (let [rows (for [[id r b] observations
                   :let [[er eb] (predict mean id)
                         [hr hb] (jacobian mean id)]]
               [[hr (- r er) (* sigma-r sigma-r)]
                [hb (wrap-angle (- b eb)) (* sigma-b sigma-b)]])
        {:keys [x P]} (reduce (fn [{:keys [x P]} [h innovation R]]
                                (select-keys (tm/kalman-update x P h (- innovation (lin/dot h x)) R) [:x :P]))
                              {:x (vec (repeat (count mean) 0.0)) :P cov}
                              (apply concat rows))
        mean' (lin/add mean x)]
    [(assoc mean' 2 (wrap-angle (mean' 2))) P]))

(defn ekf-predict
  "The EKF's prediction: the pose driven by `[v w]`, the covariance
  carried by the motion model's partials, the control noise added through
  them -- `temporal`'s `kalman-predict` with the linearized model."
  [{:keys [dt sigma-v sigma-w]} [mean cov] u]
  (let [[G V] (motion-jacobians mean u dt)
        M [[(* sigma-v sigma-v) 0.0] [0.0 (* sigma-w sigma-w)]]
        Q (lin/mat-mul (lin/mat-mul V M) (lin/transpose V))
        [v w] u]
    [(drive mean v w dt) (:P (tm/kalman-predict mean cov G Q))]))

(defn ekf-localization-step
  "One step of EKF localization among known, identifiable `:landmarks`:
  the belief `[mean cov]` predicted through control `[v w]` and corrected
  by the observations made after it, `[[id r bearing] ...]`."
  [{:keys [landmarks sigma-r sigma-b] :as world} belief u observations]
  (ekf-correct (ekf-predict world belief u) observations
               (fn [mean id] (range-bearing mean (landmarks id)))
               (fn [mean id] (mapv #(subvec % 0 3) (range-bearing-jacobian mean (landmarks id))))
               sigma-r sigma-b))

(defn ekf-localization-seq
  "EKF localization from the belief `[mean cov]` through `controls` and
  the observations made after each: the beliefs."
  [world belief controls observations]
  (reductions (fn [belief [u z]] (ekf-localization-step world belief u z))
              belief (map vector controls observations)))

(defn drive-and-observe
  "A simulated run: the robot driven from `start` by `controls` with the
  world's motion noise, observing the landmarks after each move --
  `{:truth [pose ...] :observations [[[id r bearing] ...] ...]}`, the
  truth including the start."
  [{:keys [dt sigma-v sigma-w] :as world} start controls rng]
  (let [truth (vec (reductions (fn [p [v w]] (move p v w dt rng sigma-v sigma-w)) start controls))]
    {:truth truth :observations (mapv #(observe-landmarks world % rng) (rest truth))}))

(defn landmark-field
  "A test world for the Kalman filters: landmarks scattered round a 20 by
  14 field, seen within 6 m."
  []
  {:landmarks [[2.0 2.0] [9.0 1.0] [17.0 2.5] [18.5 8.0] [16.0 13.0] [10.0 12.5]
               [3.0 12.0] [1.5 7.0] [7.0 6.0] [13.0 7.5] [5.0 9.5] [14.5 4.0]]
   :max-range 6.0 :dt 0.1 :sigma-v 0.1 :sigma-w 0.05 :sigma-r 0.1 :sigma-b 0.02})
