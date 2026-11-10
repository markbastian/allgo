(ns allgo.robotics.graph-slam
  "Graph SLAM: Russell and Norvig, *Artificial Intelligence: A Modern
  Approach*, 4th edition, chapter 26 -- mapping as one least-squares
  problem over the whole run rather than a filter over its latest pose.

  The unknowns are the robot's poses at chosen moments -- every `:every`
  steps of the run -- and the landmarks' positions. Each constraint is a
  measurement with its noise: the odometry between consecutive poses,
  the controls driven between them composed into one relative move with
  its covariance carried along (by `localization`'s `ekf-predict` from a
  pose known exactly), and each landmark's range and bearing seen from a
  pose. The map and path most probable given them all minimize the sum
  of the squared errors, each weighted by its measurement's inverse
  covariance. That is nonlinear in the headings, so it is solved by
  Gauss-Newton -- linearize every constraint about the current guess,
  solve the normal equations, move, repeat -- an iteration a step. The
  first pose is held where it is: it fixes the map's frame.

  The guess it starts from is dead reckoning, with each landmark placed
  from where it was first seen. Unlike a filter, it revisits every
  earlier pose with every later sighting, so a loop closed late corrects
  the start of the run as well."
  (:require [allgo.numerics.linear :as lin]
            [allgo.robotics.localization :as loc]
            [allgo.search :as search]
            [clojure.math :as math]))

(defn- relative
  "Pose `b` in the frame of pose `a`."
  [[ax ay ah] [bx by bh]]
  (let [dx (- bx ax) dy (- by ay) c (math/cos ah) s (math/sin ah)]
    [(+ (* c dx) (* s dy)) (+ (* -1.0 s dx) (* c dy)) (loc/wrap-angle (- bh ah))]))

(defn- relative-jacobians
  "The partials of `relative` by `a` and by `b`."
  [[ax ay ah] [bx by _]]
  (let [dx (- bx ax) dy (- by ay) c (math/cos ah) s (math/sin ah)]
    [[[(- c) (- s) (+ (* -1.0 s dx) (* c dy))]
      [s (- c) (- (* -1.0 c dx) (* s dy))]
      [0.0 0.0 -1.0]]
     [[c s 0.0] [(- s) c 0.0] [0.0 0.0 1.0]]]))

(defn constraints
  "The graph of a run: its controls and the observations after each, the
  poses kept every `every` steps. `{:poses [step ...] :odometry [{:from
  :to :z :info} ...] :sightings [{:pose :landmark :z :info} ...]}`, the
  indices into the kept poses; `:z` the measurement, `:info` its inverse
  covariance."
  [{:keys [sigma-r sigma-b] :as world} controls observations every]
  (let [n (count controls)
        kept (vec (distinct (concat (range 0 (inc n) every) [n])))
        rz [[(/ 1.0 (* sigma-r sigma-r)) 0.0] [0.0 (/ 1.0 (* sigma-b sigma-b))]]]
    {:poses kept
     :odometry (vec (for [[k [s e]] (map-indexed vector (partition 2 1 kept))
                          :let [[z cov] (reduce #(loc/ekf-predict world %1 %2)
                                                [[0.0 0.0 0.0] (vec (repeat 3 [0.0 0.0 0.0]))]
                                                (subvec (vec controls) s e))]]
                      {:from k :to (inc k) :z z :info (lin/inverse cov)}))
     ;; Observations are made after each control: those after step s are
     ;; seen from the pose at s.
     :sightings (vec (for [[k s] (map-indexed vector kept)
                           :when (pos? s)
                           [id r b] (observations (dec s))]
                       {:pose k :landmark id :z [r b] :info rz}))}))

(defn initial-guess
  "Dead reckoning from `start` through the controls, at the kept poses, and
  each landmark placed from the first pose that saw it:
  `{:poses [[x y h] ...] :landmarks {id [x y]}}`."
  [{:keys [dt]} {:keys [poses sightings]} start controls]
  (let [path (vec (reductions (fn [p [v w]] (loc/drive p v w dt)) start controls))
        ps (mapv path poses)]
    {:poses ps
     :landmarks (reduce (fn [m {:keys [pose landmark z]}]
                          (if (contains? m landmark)
                            m
                            (let [[x y h] (ps pose) [r b] z]
                              (assoc m landmark [(+ x (* r (math/cos (+ h b)))) (+ y (* r (math/sin (+ h b))))]))))
                        {} sightings)}))

(defn- layout
  "Where each unknown sits in the state vector: pose k at 3k, landmark id
  after the poses, in order of id."
  [{:keys [poses landmarks]}]
  (let [ids (sort (keys landmarks))]
    {:n (+ (* 3 (count poses)) (* 2 (count ids)))
     :landmark-at (zipmap ids (map #(+ (* 3 (count poses)) (* 2 %)) (range)))}))

(defn- quadratic
  "e^T info e."
  [info e]
  (reduce + (map (fn [ei row] (* ei (reduce + (map * row e)))) e info)))

(defn error
  "The weighted sum of squared errors of `estimate`: what Gauss-Newton
  minimizes."
  [{:keys [odometry sightings]} {:keys [poses landmarks]}]
  (+ (reduce + (for [{:keys [from to z info]} odometry
                     :let [e (mapv - (relative (poses from) (poses to)) z)
                           e (assoc e 2 (loc/wrap-angle (e 2)))]]
                 (quadratic info e)))
     (reduce + (for [{:keys [pose landmark z info]} sightings
                     :let [[er eb] (loc/range-bearing (poses pose) (landmarks landmark))
                           e [(- er (z 0)) (loc/wrap-angle (- eb (z 1)))]]]
                 (quadratic info e)))))

(defn- accumulate!
  "Adds J_a^T info J_b into H for each pair of blocks, and J_a^T info e
  into g: one constraint's share of the normal equations."
  [^objects H ^doubles g blocks info e]
  (let [m (count e)]
    (doseq [[ia Ja] blocks
            :let [JtI (mapv (fn [c] (mapv (fn [r] (reduce + (map (fn [k] (* (get-in Ja [k c]) (get-in info [k r]))) (range m))))
                                          (range m)))
                            (range (count (first Ja))))]]
      (doseq [c (range (count JtI))]
        (aset g (+ (long ia) c) (+ (aget g (+ (long ia) c)) (double (reduce + (map * (JtI c) e)))))
        (doseq [[ib Jb] blocks
                d (range (count (first Jb)))]
          (let [^doubles row (aget H (+ (long ia) c))
                v (reduce + (map (fn [r] (* (get-in JtI [c r]) (get-in Jb [r d]))) (range m)))]
            (aset row (+ (long ib) d) (+ (aget row (+ (long ib) d)) (double v)))))))))

(defn- gauss-newton-step
  [{:keys [graph estimate] :as st}]
  (let [{:keys [poses landmarks]} estimate
        {:keys [n landmark-at]} (layout estimate)
        H (object-array (repeatedly n #(double-array n)))
        g (double-array n)]
    (doseq [{:keys [from to z info]} (:odometry graph)
            :let [a (poses from) b (poses to)
                  e (mapv - (relative a b) z)
                  e (assoc e 2 (loc/wrap-angle (e 2)))
                  [Ja Jb] (relative-jacobians a b)]]
      (accumulate! H g [[(* 3 from) Ja] [(* 3 to) Jb]] info e))
    (doseq [{:keys [pose landmark z info]} (:sightings graph)
            :let [p (poses pose) m (landmarks landmark)
                  [er eb] (loc/range-bearing p m)
                  e [(- er (z 0)) (loc/wrap-angle (- eb (z 1)))]
                  J (loc/range-bearing-jacobian p m)]]
      (accumulate! H g [[(* 3 pose) (mapv #(subvec % 0 3) J)] [(landmark-at landmark) (mapv #(subvec % 3 5) J)]] info e))
    ;; The first pose held fixed: its rows and columns made the identity,
    ;; its gradient zero.
    (doseq [i (range 3)]
      (doseq [j (range n)] (aset ^doubles (aget H i) j 0.0) (aset ^doubles (aget H j) i 0.0))
      (aset ^doubles (aget H i) i 1.0)
      (aset g i 0.0))
    (when-let [dx (lin/cholesky-solve (mapv vec H) (mapv - g))]
      (let [poses' (vec (map-indexed (fn [k [x y h]]
                                       (let [i (* 3 k)]
                                         [(+ x (dx i)) (+ y (dx (+ i 1))) (loc/wrap-angle (+ h (dx (+ i 2))))]))
                                     poses))
            landmarks' (into {} (map (fn [[id [x y]]] (let [i (landmark-at id)] [id [(+ x (dx i)) (+ y (dx (inc i)))]])))
                             landmarks)
            estimate' {:poses poses' :landmarks landmarks'}]
        (assoc st :estimate estimate' :error (error graph estimate')
               :iteration (inc (long (:iteration st)))
               :change (math/sqrt (reduce + (map #(* % %) dx))))))))

(defn graph-slam-seq
  "Gauss-Newton on the run's graph, from dead reckoning, an iteration a
  step: each state `{:graph :estimate :error :iteration :change}`, the
  estimate's poses those kept every `:every` steps (default 10), `:error`
  the weighted squared error and `:change` the size of the step that led
  there. Ends when a step changes less than `:tolerance` (default 1e-6)
  or the normal equations cannot be solved."
  [world start controls observations & {:keys [every tolerance] :or {every 10 tolerance 1e-6}}]
  (let [graph (constraints world (vec controls) (vec observations) every)
        estimate (initial-guess world graph start controls)
        first-state {:graph graph :estimate estimate :error (error graph estimate) :iteration 0 :change ##Inf}]
    (concat [first-state]
            (->> (search/step-seq gauss-newton-step first-state)
                 rest
                 (reduce (fn [acc st] (if (< (double (:change st)) tolerance) (reduced (conj acc st)) (conj acc st))) [])))))

(defn graph-slam
  "The final estimate: `{:poses :landmarks}`."
  [world start controls observations & opts]
  (:estimate (last (apply graph-slam-seq world start controls observations opts))))
