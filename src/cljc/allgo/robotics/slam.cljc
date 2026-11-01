(ns allgo.robotics.slam
  "Simultaneous localization and mapping by the extended Kalman filter:
  Russell and Norvig, *Artificial Intelligence: A Modern Approach*, 4th
  edition, chapter 26, which teaches it in the text.

  The robot knows neither where it is nor where the landmarks are, only
  that there are `n` of them and which is which when it sees one. The
  state is the pose and every landmark's position, `[x y heading x1 y1
  ... xn yn]`, one Gaussian over all of it. A move carries the pose and
  its correlations with the map forward -- the landmarks do not move --
  and an observation corrects everything correlated with what was seen:
  seeing a known landmark again pulls the pose back, and with it every
  landmark mapped from the drifting pose, which is why closing a loop
  sharpens the whole map at once.

  A landmark seen for the first time is placed where its range and
  bearing say, from the current pose, and given a covariance that says
  the pose's uncertainty plus the measurement's -- its correlations with
  the pose carried over, so that it moves when the pose is corrected.

  The motion and sensor models, and the scalar Joseph-form update, are
  `allgo.robotics.localization`'s and `allgo.probability.temporal`'s."
  (:require [allgo.numerics.linear :as lin]
            [allgo.probability.temporal :as tm]
            [allgo.robotics.localization :as loc]
            [clojure.math :as math]))

(defn- landmark-index ^long [id] (+ 3 (* 2 id)))

(defn- landmark [mean id] (let [i (landmark-index id)] [(mean i) (mean (inc i))]))

(defn initial-belief
  "The belief before anything is seen: the pose known exactly at `pose`
  -- it defines the map's frame -- and `n` landmarks not yet placed."
  [pose n]
  (let [d (+ 3 (* 2 n))]
    {:mean (into (vec pose) (repeat (* 2 n) 0.0))
     :cov (vec (repeat d (vec (repeat d 0.0))))
     :seen #{}}))

(defn- predict
  "The pose moved by `[v w]`, the landmarks unchanged: the motion
  Jacobian is the identity but for its pose block."
  [{:keys [dt sigma-v sigma-w]} {:keys [mean cov] :as belief} u]
  (let [d (count mean)
        [G3 V] (loc/motion-jacobians mean u dt)
        G (mapv (fn [i] (mapv (fn [j] (if (and (< i 3) (< j 3)) (get-in G3 [i j]) (if (= i j) 1.0 0.0))) (range d))) (range d))
        M [[(* sigma-v sigma-v) 0.0] [0.0 (* sigma-w sigma-w)]]
        Q3 (lin/mat-mul (lin/mat-mul V M) (lin/transpose V))
        Q (mapv (fn [i] (mapv (fn [j] (if (and (< i 3) (< j 3)) (get-in Q3 [i j]) 0.0)) (range d))) (range d))
        [v w] u]
    (assoc belief
           :mean (into (loc/drive (subvec mean 0 3) v w dt) (subvec mean 3))
           :cov (:P (tm/kalman-predict mean cov G Q)))))

(defn- place
  "A first sighting of landmark `id` at range `r` and bearing `b`: its
  position from the pose, and its covariance -- the pose's carried
  through the placement, plus the measurement's -- with its correlation
  to everything the pose is correlated with."
  [{:keys [sigma-r sigma-b]} {:keys [mean cov] :as belief} [id r b]]
  (let [[x y h] mean
        a (+ h b) c (math/cos a) s (math/sin a)
        i (landmark-index id)
        d (count mean)
        ;; m = (x + r cos a, y + r sin a): partials by pose and by (r, b)
        Gp [[1.0 0.0 (* -1.0 r s)] [0.0 1.0 (* r c)]]
        Gz [[c (* -1.0 r s)] [s (* r c)]]
        Ppp (mapv #(subvec % 0 3) (subvec cov 0 3))
        Pmm (lin/mat-add (lin/mat-mul (lin/mat-mul Gp Ppp) (lin/transpose Gp))
                         (lin/mat-mul (lin/mat-mul Gz [[(* sigma-r sigma-r) 0.0] [0.0 (* sigma-b sigma-b)]]) (lin/transpose Gz)))
        ;; the new rows: Gp times the pose's rows of the covariance
        Pmx (lin/mat-mul Gp (subvec cov 0 3))
        cov (reduce (fn [cov k]
                      (if (or (= k i) (= k (inc i)))
                        cov
                        (-> cov
                            (assoc-in [i k] (get-in Pmx [0 k])) (assoc-in [k i] (get-in Pmx [0 k]))
                            (assoc-in [(inc i) k] (get-in Pmx [1 k])) (assoc-in [k (inc i)] (get-in Pmx [1 k])))))
                    cov (range d))
        cov (-> cov
                (assoc-in [i i] (get-in Pmm [0 0])) (assoc-in [i (inc i)] (get-in Pmm [0 1]))
                (assoc-in [(inc i) i] (get-in Pmm [1 0])) (assoc-in [(inc i) (inc i)] (get-in Pmm [1 1])))]
    (assoc belief
           :mean (-> mean (assoc i (+ x (* r c))) (assoc (inc i) (+ y (* r s))))
           :cov cov
           :seen (conj (:seen belief) id))))

(defn- correct
  "The correction by landmarks seen before: `localization`'s EKF update,
  the sensor Jacobian spread over the pose's and that landmark's columns."
  [{:keys [sigma-r sigma-b]} {:keys [mean cov] :as belief} observations]
  (if (empty? observations)
    belief
    (let [d (count mean)
          [mean' cov'] (loc/ekf-correct
                        [mean cov] observations
                        (fn [m id] (loc/range-bearing (subvec m 0 3) (landmark m id)))
                        (fn [m id]
                          (let [J (loc/range-bearing-jacobian (subvec m 0 3) (landmark m id))
                                i (landmark-index id)]
                            (mapv (fn [row]
                                    (mapv (fn [k] (cond (< k 3) (row k)
                                                        (= k i) (row 3)
                                                        (= k (inc i)) (row 4)
                                                        :else 0.0))
                                          (range d)))
                                  J)))
                        sigma-r sigma-b)]
      (assoc belief :mean mean' :cov cov'))))

(defn slam-step
  "One EKF-SLAM step: predict through control `u`, correct by the
  observations of landmarks already mapped, place those seen for the
  first time."
  [world belief u observations]
  (let [belief (predict world belief u)
        {old true new false} (group-by (fn [[id]] (contains? (:seen belief) id)) observations)
        belief (correct world belief old)]
    (reduce #(place world %1 %2) belief new)))

(defn slam-seq
  "EKF-SLAM from `belief` through `controls` and the observations made
  after each -- the beliefs, `{:mean :cov :seen}`."
  [world belief controls observations]
  (reductions (fn [b [u z]] (slam-step world b u z)) belief (map vector controls observations)))

(defn pose
  "The belief's pose."
  [{:keys [mean]}]
  (subvec mean 0 3))

(defn landmarks
  "The belief's map: `{id [x y]}` for the landmarks seen."
  [{:keys [mean seen]}]
  (into (sorted-map) (map (fn [id] [id (landmark mean id)]) seen)))

(defn landmark-cov
  "The 2 by 2 covariance of landmark `id`."
  [{:keys [cov]} id]
  (let [i (landmark-index id)]
    [[(get-in cov [i i]) (get-in cov [i (inc i)])]
     [(get-in cov [(inc i) i]) (get-in cov [(inc i) (inc i)])]]))
