(ns allgo.probability.association-test
  "Data association, after Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, section 15.3.1: two aircraft crossing,
  blips without labels. The section's figure 15.8 is drawn without
  coordinates, so the scenarios are simulated here from the model, and
  the trackers judged against the paths that made the blips -- and the
  particles against the exact posterior over histories."
  (:require [allgo.graph :as graph]
            [allgo.numerics.linear :as lin]
            [allgo.probability.association :as association]
            [allgo.random :as random]
            [clojure.test :refer [deftest is testing]]))

(def ^:private model (association/constant-velocity 0.01 0.5))

;; one starting low and climbing, one starting high and descending: their
;; paths cross near the middle
(def ^:private starts [[0.0 0.0 1.0 0.3] [0.0 6.0 1.0 -0.3]])
(def ^:private prior (mapv (fn [s] [s (lin/mat-scale (lin/eye 4) 0.25)]) starts))
(def ^:private area [[-2.0 -3.0] [22.0 9.0]])

(defn- simulate
  "Paths for `starts` under the model, and each scan's blips: each object
  seen with probability `p-detect`, through the observation noise, and
  false alarms -- on average `clutter-rate` a scan -- uniform over the
  area, all shuffled."
  [rng scans {:keys [p-detect clutter-rate] :or {p-detect 1.0 clutter-rate 0.0}}]
  (let [step (fn [x] (mapv + (lin/mat-vec (:F model) x) [0.0 0.0 (random/gaussian rng 0.0 0.05) (random/gaussian rng 0.0 0.05)]))
        paths (mapv (fn [s] (vec (take scans (rest (iterate step s))))) starts)
        r (Math/sqrt (get-in model [:R 0 0]))
        [[x0 y0] [x1 y1]] area]
    {:paths paths
     :scans (vec (for [t (range scans)]
                   (random/shuffle rng (concat (for [path paths :when (< (rng) p-detect)]
                                                 (let [[x y] (path t)] [(random/gaussian rng x r) (random/gaussian rng y r)]))
                                               (for [_ (range 10) :when (< (rng) (/ clutter-rate 10.0))]
                                                 [(+ x0 (* (rng) (- x1 x0))) (+ y0 (* (rng) (- y1 y0)))])))))}))

(defn- rms-error
  "The root mean square distance between each track's position and its
  object's, over every scan."
  [paths track-seq]
  (Math/sqrt (/ (reduce + (for [[t tracks] (map-indexed vector track-seq) [i track] (map-indexed vector tracks)]
                            (lin/distance-squared (association/position model track) (subvec ((paths i) t) 0 2))))
                (* (count track-seq) (count paths)))))

(defn- particle-tracks
  "Each object's position averaged over the particles."
  [particles]
  (vec (for [i (range (count (:tracks (first particles))))]
         [(lin/mean (map #(first (get-in % [:tracks i])) particles))])))

(defn- errors [m seed opts]
  (let [rng (random/rng seed)
        {:keys [paths scans]} (simulate rng 20 opts)]
    {:nearest (rms-error paths (association/nearest-neighbor-seq m prior scans))
     :probable (rms-error paths (association/most-probable-seq m prior scans))
     :particles (rms-error paths (map particle-tracks (association/particle-seq m prior scans rng 200)))}))

(deftest hungarian
  (testing "the Hungarian method's assignment costs what the best of every permutation does, rows fewer than columns too"
    (let [rng (random/rng 3)
          perms (fn perms [xs] (if (empty? xs) [[]] (for [x xs p (perms (remove #{x} xs))] (cons x p))))]
      (doseq [_ (range 200)]
        (let [n (inc (random/below rng 5)) m (+ n (random/below rng 3))
              cost (vec (repeatedly n #(vec (repeatedly m (fn [] (double (random/below rng 20)))))))
              cols (graph/min-cost-assignment cost)
              total (fn [cs] (reduce + (map (fn [i j] (get-in cost [i j])) (range n) cs)))]
          (is (= n (count (distinct cols))))
          (is (== (apply min (map total (perms (range m)))) (total cols))))))))

(deftest particles-against-the-exact-posterior
  (testing "on three scans the particles' posterior over association histories is the exact one, enumerated"
    (let [rng (random/rng 4)
          m (assoc model :p-detect 0.9 :clutter 0.05)
          {:keys [scans]} (simulate rng 3 {:p-detect 0.9 :clutter-rate 0.5})
          exact (association/exact-posterior m prior scans)
          estimate (association/history-posterior (last (association/particle-seq m prior scans rng 4000)))]
      (is (< (abs (- 1.0 (reduce + (vals exact)))) 1e-9))
      (doseq [[h p] (take 4 (sort-by (comp - val) exact))]
        (is (< (abs (- p (get estimate h 0.0))) 0.03) (str h))))))

(deftest crossing-paths
  (let [runs (mapv #(errors model % {}) (range 1 11))
        mean (fn [k] (/ (reduce + (map k runs)) (count runs)))]
    (testing "a single assignment a scan usually tracks well, but where the paths cross it can commit to the wrong pairing and stray"
      (is (some #(> (max (:nearest %) (:probable %)) 1.5) runs)))
    (testing "the particles keep both pairings alive until later scans settle it: never lost, and the least error on average"
      (is (every? #(< (:particles %) 1.0) runs))
      (is (< (mean :particles) (mean :nearest)))
      (is (< (mean :particles) (mean :probable))))))

(deftest clutter-and-misses
  (let [m (assoc model :p-detect 0.9 :clutter (/ 2.0 (* 24.0 12.0)))
        runs (mapv #(errors m % {:p-detect 0.9 :clutter-rate 2.0}) (range 1 9))
        mean (fn [k] (/ (reduce + (map k runs)) (count runs)))]
    (testing "with false alarms and detection failures the nearest-neighbor filter, taking whatever blip is closest, loses its tracks"
      (is (> (mean :nearest) 2.0)))
    (testing "the most probable assignment, able to put a blip down to clutter or a track down as unseen, does better"
      (is (< (mean :probable) (mean :nearest))))
    (testing "and the particles best of all"
      (is (every? #(< (:particles %) 1.1) runs))
      (is (< (mean :particles) 0.7)))))
