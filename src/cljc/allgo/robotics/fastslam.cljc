(ns allgo.robotics.fastslam
  "FastSLAM: Russell and Norvig, *Artificial Intelligence: A Modern
  Approach*, 4th edition, chapter 26 -- the Rao-Blackwellized particle
  filter its notes credit to Murphy and Russell, and to Montemerlo et al.
  for mapping.

  Given the robot's whole path, the landmarks' positions are independent
  of one another: each is fixed by its own sightings from known poses. So
  each particle is a path guess -- carried as its latest pose -- with, for
  every landmark it has seen, a small Kalman filter of its own: a 2-vector
  and a 2 by 2 covariance, rather than EKF-SLAM's one Gaussian over
  everything. Each step every particle draws its next pose from the motion
  model, updates the filters of the landmarks seen, and is weighted by how
  likely those sightings were under them; then the particles are
  resampled. A landmark seen for the first time is placed by inverting the
  measurement from the particle's pose.

  The motion and sensor models, and the resampling, are
  `allgo.robotics.localization`'s; known data association, as in
  `allgo.robotics.slam`."
  (:require [allgo.robotics.localization :as loc]
            [clojure.math :as math]))

(defn initial-particles
  "`n` particles at `pose`, which defines the map's frame, with no
  landmarks."
  [pose n]
  (vec (repeat n {:pose pose :landmarks {} :weight 1.0})))

(defn- mat2-mul [[[a b] [c d]] [[e f] [g h]]]
  [[(+ (* a e) (* b g)) (+ (* a f) (* b h))] [(+ (* c e) (* d g)) (+ (* c f) (* d h))]])

(defn- transpose2 [[[a b] [c d]]] [[a c] [b d]])

(defn- inverse2 [[[a b] [c d]]]
  (let [det (- (* a d) (* b c))]
    [[(/ d det) (/ (- b) det)] [(/ (- c) det) (/ a det)]]))

(defn- add2 [[[a b] [c d]] [[e f] [g h]]] [[(+ a e) (+ b f)] [(+ c g) (+ d h)]])

(defn- landmark-jacobian
  "The partials of range and bearing with respect to the landmark's
  position: the last two columns of `range-bearing-jacobian`."
  [pose m]
  (mapv #(subvec % 3 5) (loc/range-bearing-jacobian pose m)))

(defn- noise [{:keys [sigma-r sigma-b]}] [[(* sigma-r sigma-r) 0.0] [0.0 (* sigma-b sigma-b)]])

(defn- place
  "A first sighting `[r b]` from `pose`: the landmark where it says, its
  covariance the measurement noise carried back through the inverse
  of the sensor model's partials."
  [world [x y h :as pose] [r b]]
  (let [a (+ h b)
        m [(+ x (* r (math/cos a))) (+ y (* r (math/sin a)))]
        Hinv (inverse2 (landmark-jacobian pose m))]
    {:mean m :cov (mat2-mul (mat2-mul Hinv (noise world)) (transpose2 Hinv))}))

(defn- update-landmark
  "The landmark's Kalman update by a sighting `[r b]`, and the sighting's
  likelihood under it: `[landmark likelihood]`."
  [world pose {:keys [mean cov]} [r b]]
  (let [[er eb] (loc/range-bearing pose mean)
        H (landmark-jacobian pose mean)
        S (add2 (mat2-mul (mat2-mul H cov) (transpose2 H)) (noise world))
        Si (inverse2 S)
        K (mat2-mul (mat2-mul cov (transpose2 H)) Si)
        nu [(- r er) (loc/wrap-angle (- b eb))]
        [[k00 k01] [k10 k11]] K
        [n0 n1] nu
        mean' [(+ (mean 0) (* k00 n0) (* k01 n1)) (+ (mean 1) (* k10 n0) (* k11 n1))]
        KH (mat2-mul K H)
        cov' (mat2-mul [[(- 1.0 (get-in KH [0 0])) (- (get-in KH [0 1]))]
                        [(- (get-in KH [1 0])) (- 1.0 (get-in KH [1 1]))]]
                       cov)
        [[s00 s01] [s10 s11]] Si
        d2 (+ (* n0 (+ (* s00 n0) (* s01 n1))) (* n1 (+ (* s10 n0) (* s11 n1))))
        [[a b'] [c d]] S
        likelihood (/ (math/exp (* -0.5 d2)) (* 2 math/PI (math/sqrt (- (* a d) (* b' c)))))]
    [{:mean mean' :cov cov'} likelihood]))

(defn- particle-step
  [{:keys [dt sigma-v sigma-w] :as world} rng {:keys [pose landmarks]} [v w] observations]
  (let [pose (loc/move pose v w dt rng sigma-v sigma-w)]
    (reduce (fn [p [id r b]]
              (if-let [l (get-in p [:landmarks id])]
                (let [[l' lik] (update-landmark world pose l [r b])]
                  (-> p (assoc-in [:landmarks id] l') (update :weight * lik)))
                (assoc-in p [:landmarks id] (place world pose [r b]))))
            {:pose pose :landmarks landmarks :weight 1.0}
            observations)))

(defn fastslam-step
  "One FastSLAM step: every particle moved by control `u` with the motion
  noise, its landmarks updated by `observations`, `[[id r bearing] ...]`,
  and weighted; then the particles resampled by weight. The particles
  returned are the ones resampled, each with the weight it had."
  [world rng particles u observations]
  (let [moved (mapv #(particle-step world rng % u observations) particles)
        weights (mapv :weight moved)]
    (if (pos? (reduce + weights))
      (loc/resample moved weights rng)
      moved)))

(defn fastslam-seq
  "FastSLAM from `particles` through `controls` and the observations made
  after each: the particle sets."
  [world rng particles controls observations]
  (reductions (fn [ps [u z]] (fastslam-step world rng ps u z)) particles (map vector controls observations)))

(defn best
  "The particle of greatest weight: the path, and the map, most in
  keeping with the latest sightings."
  [particles]
  (apply max-key :weight particles))

(defn pose
  "The particles' mean pose, the headings averaged on the circle."
  [particles]
  (let [n (double (count particles))
        ps (map :pose particles)]
    [(/ (reduce + (map first ps)) n) (/ (reduce + (map second ps)) n)
     (math/atan2 (reduce + (map #(math/sin (% 2)) ps)) (reduce + (map #(math/cos (% 2)) ps)))]))

(defn landmarks
  "A particle's map: `{id [x y]}`."
  [particle]
  (into (sorted-map) (map (fn [[id l]] [id (:mean l)])) (:landmarks particle)))
