(ns allgo.robotics.pursuit-test
  (:require [allgo.geometry.polygon :as polygon]
            [allgo.robotics.pursuit :as pursuit]
            [clojure.set :as set]
            [clojure.test :refer [deftest is testing]]))

;; Visibility-based pursuit-evasion, after LaValle, *Planning Algorithms*,
;; Cambridge University Press, 2006, section 12.4. The environments are
;; made up here: an L, two corridors that bend twice off a hub, and three
;; -- Example 12.7's three hiding places, no two of which can watch the
;; same ground.

(def ell [[0 0] [10 0] [10 2] [2 2] [2 8] [0 8]])

(def two-hooks
  [[-3 -1] [3 -1] [3 -7] [11 -7] [11 -5] [5 -5] [5 1] [-5 1] [-5 -5] [-11 -5] [-11 -7] [-3 -7]])

(def three-hooks
  [[-3 -1] [3 -1] [3 -7] [11 -7] [11 -5] [5 -5] [5 1] [1 1] [1 3] [7 3] [7 11] [5 11] [5 5]
   [-1 5] [-1 1] [-5 1] [-5 -5] [-11 -5] [-11 -7] [-3 -7]])

(defn- dense
  "`path` with `n` points to each of its steps, so a replay moves in
  short ones."
  [path n]
  (conj (vec (mapcat (fn [[a b]] (map #(polygon/lerp a b (/ % n)) (range n))) (partition 2 1 path)))
        (peek path)))

(deftest visibility-test
  (let [env (pursuit/environment ell 0.5)]
    (testing "probes half a spacing in from the walls, a lattice step apart"
      (is (= 128 (count (:probes env))))
      (is (every? #(polygon/contains-point? ell %) (:probes env))))
    (testing "the outer corner sees the whole L; the two ends do not see each other"
      (is (= (count (:probes env)) (count (pursuit/visible env [0.25 0.25]))))
      (is (empty? (pursuit/shadow env [0.25 0.25])))
      (is (not (pursuit/sees? env [9.75 1.0] [1.0 7.75]))))
    (testing "the visibility polygon: all of a convex region, all of the L from its outer corner, a leg and a wedge from an end"
      (is (< (abs (- 16.0 (polygon/area (pursuit/visibility-polygon (pursuit/environment [[0 0] [4 0] [4 4] [0 4]] 0.5) [1 3])))) 1e-6))
      (is (< (abs (- 32.0 (polygon/area (pursuit/visibility-polygon env [0.25 0.25])))) 1e-6))
      (let [a (polygon/area (pursuit/visibility-polygon env [9.75 1.0]))]
        ;; the leg, 20, and the triangle of the other leg's mouth seen past the inner corner
        (is (< (abs (- a (+ 20.0 (* 0.5 2.0 (/ 2.0 7.75))))) 1e-4))))
    (testing "from the end of one leg the shadow is one component, the far part of the other"
      (let [[c :as s] (pursuit/shadow env [9.75 1.0])]
        (is (= 1 (count s)))
        (is (every? (fn [k] (let [[x y] ((:probes env) k)] (and (< x 2) (> y 2)))) c))))))

(deftest labeling-test
  (let [env (pursuit/environment ell 0.5)
        far #(vec (range % (+ % 2)))]
    (testing "the four critical events, read off which components overlap"
      (is (= [[:appear 1]] (pursuit/events [[1 2]] [[1 2] [9]])))
      (is (= [[:disappear 0]] (pursuit/events [[9] [1 2]] [[1 2]])))
      (is (= [[:split 0 [0 1]]] (pursuit/events [[1 2 3 4]] [[1] [4]])))
      (is (= [[:merge [0 1] 0]] (pursuit/events [[1] [4]] [[1 2 3 4]]))))
    (testing "a component that appears is cleared; a split keeps the label; a merge with a contaminated one is contaminated"
      (is (= #{0} (pursuit/relabel env [(far 40)] #{0} [(far 40) (far 100)])))
      (is (= #{0 1} (pursuit/relabel env [(vec (range 40 120))] #{0} [(far 40) (far 100)])))
      (is (= #{0} (pursuit/relabel env [(far 40) (far 100)] #{1} [(vec (range 40 120))])))
      (is (= #{} (pursuit/relabel env [(far 40) (far 100)] #{} [(vec (range 40 120))]))))))

(deftest critical-boundary-test
  (testing "figure 12.38's lesson: a loop that crosses no critical boundary changes nothing; one that does clears what it saw"
    (let [env (pursuit/environment two-hooks 0.5)
          west? (fn [k] (neg? (first ((:probes env) k))))
          small (last (pursuit/clearing-seq env (dense [[0 0] [0.5 0] [0.5 0.5] [0 0.5] [0 0]] 4)))
          large (last (pursuit/clearing-seq env (dense [[0 0] [4 0] [4 -6] [10 -6] [4 -6] [4 0] [0 0]] 8)))
          start (pursuit/start env [0 0])]
      (is (= 2 (count (:contaminated start))))
      (is (= (pursuit/contaminated-probes start) (pursuit/contaminated-probes small)))
      (is (= 2 (count (:shadow large))))
      (is (= 1 (count (:contaminated large))))
      (is (every? west? (pursuit/contaminated-probes large)))
      (is (some (complement west?) (mapcat (:shadow large) [0 1]))))))

(deftest plan-test
  (testing "two hooks: a plan from anywhere, which cleared again on a lattice twice as fine leaves nowhere to hide"
    (let [env (pursuit/environment two-hooks 0.5)
          fine (pursuit/environment two-hooks 0.25)]
      (doseq [p [[0 0] [10 -6] [-10 -6] [4 -3]]
              :let [path (pursuit/plan env p)]]
        (is (some? path))
        (is (every? (fn [[a b]] (< (abs (- (polygon/distance a b) 0.5)) 1e-9)) (partition 2 1 path)))
        (is (pursuit/cleared? (last (pursuit/clearing-seq env path))))
        (is (pursuit/cleared? (last (pursuit/clearing-seq fine (dense path 4))))))))
  (testing "the L from the end of a leg: walk to the corner, and no further"
    (let [env (pursuit/environment ell 0.5)
          path (pursuit/plan env [9.75 1.25])]
      (is (pursuit/cleared? (last (pursuit/clearing-seq env path))))
      (is (not-any? #(pursuit/cleared? (last (pursuit/clearing-seq env %))) (rest (butlast (reductions conj [] path)))))))
  (testing "where the shadow is already empty the plan is to stay"
    (is (= 1 (count (pursuit/plan (pursuit/environment ell 0.5) [0.25 0.25]))))))

(deftest three-hiding-places-test
  (let [env (pursuit/environment three-hooks 0.5)
        ends (map #(pursuit/nearest-probe env %) [[10 -6] [-10 -6] [6 10]])
        watched (fn [k] (into #{} (mapcat #(pursuit/visible env ((:probes env) %))) (pursuit/visible env ((:probes env) k))))]
    (testing "Example 12.7's condition: three places whose V(V(p)) are pairwise disjoint"
      (is (every? empty? (for [a ends b ends :when (< a b)] (set/intersection (watched a) (watched b))))))
    (testing "so one pursuer cannot clear it: the whole information space is searched without a plan, from anywhere"
      (doseq [p [[0 0] [10 -6] [6 10]]]
        (is (nil? (pursuit/plan env p)))))
    (testing "with one hook fewer it can"
      (is (some? (pursuit/plan (pursuit/environment two-hooks 0.5) [0 0]))))))
