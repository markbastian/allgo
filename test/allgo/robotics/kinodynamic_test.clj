(ns allgo.robotics.kinodynamic-test
  (:require [allgo.random :as random]
            [allgo.robotics.kinodynamic :as kinodynamic]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

;; Planning for a double integrator, after LaValle, *Planning Algorithms*,
;; Cambridge University Press, 2006, sections 14.4.1 and 14.4.3. The exact
;; check is worked out here: the quickest rest-to-rest move along a line.

(def ^:private square [[1.0 1.0] [3.0 1.0] [3.0 3.0] [1.0 3.0]])
(def ^:private room {:bounds [[-1.0 -1.0] [5.0 5.0]] :obstacles [square] :radius 0.0})
(def ^:private opts {:amax 1.0 :vmax 1.0 :dt 0.5 :tolerance 0.1 :slow 0.3})

(defn- replays?
  "Whether the controls, propagated from the first state, give the
  states, every arc free."
  [world {:keys [states controls]}]
  (and (= (count controls) (dec (count states)))
       (every? (fn [[s u s']] (and (every? #(< (abs %) 1e-9) (map - (kinodynamic/propagate s u 0.5) s'))
                                   (kinodynamic/arc-free? world s u 0.5)))
               (map vector states controls (rest states)))))

(deftest lattice-test
  (testing "4 m from rest to rest at 1 m/s and 1 m/s^2: a second speeding up, three cruising, one braking -- 5 s"
    (let [empty-room (assoc room :obstacles [])
          p (kinodynamic/lattice-plan empty-room [0.0 2.0] [4.0 2.0] opts)]
      (is (= 5.0 (:time p)))
      (is (= [4.0 2.0 0.0 0.0] (peek (:states p))))
      (is (replays? empty-room p))))
  (testing "round the square, a little longer, and every arc free"
    (let [p (kinodynamic/lattice-plan room [0.0 2.0] [4.0 2.0] opts)]
      (is (< 5.0 (:time p) 7.0))
      (is (= [0.0 0.0] (subvec (peek (:states p)) 2)))
      (is (replays? room p))
      (is (every? (fn [[_ _ vx vy]] (and (<= (abs vx) 1.0) (<= (abs vy) 1.0))) (:states p))))))

(deftest kinodynamic-rrt-test
  (testing "the tree reaches the goal slowly enough, round the square, by controls that replay to its states"
    (doseq [seed [3 4]]
      (let [st (first (filter :goal (take 5000 (kinodynamic/rrt-seq room [0.0 2.0] [4.0 2.0] (random/rng seed) opts))))
            {:keys [states] :as p} (kinodynamic/rrt-plan st)
            [x y vx vy] (peek states)]
        (is (some? st))
        (is (<= (math/hypot (- x 4.0) (- y 2.0)) 0.1))
        (is (< (math/hypot vx vy) 0.3))
        (is (replays? room p))))))
