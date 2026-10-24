(ns allgo.instancing-test
  (:require [allgo.geometry.vec3 :as v]
            [allgo.procedural.instancing :as in]
            [allgo.procedural.lsystem :as ls]
            [clojure.test :refer [deftest is testing]]))

(defn- close? [a b tol] (every? #(< (abs (double %)) tol) (map - a b)))

(deftest transform-test
  (testing "compose applies the right-hand transform first"
    (let [m (in/compose (in/translation [1 2 3]) (in/scaling 2.0))]
      (is (close? [3.0 4.0 5.0] (in/transform-point m [1 1 1]) 1e-12))))
  (testing "a rotation turns x toward y about z, and keeps lengths"
    (let [r (in/rotation [0.0 0.0 1.0] (/ Math/PI 2))]
      (is (close? [0.0 1.0 0.0] (in/transform-point r [1 0 0]) 1e-12))
      (is (< (abs (- 1.0 (in/stretch r))) 1e-12))))
  (testing "a frame's columns are where the axes go"
    (let [m (in/frame [0 1 0] [-1 0 0] [0 0 1] [5 0 0])]
      (is (close? [5.0 1.0 0.0] (in/transform-point m [1 0 0]) 1e-12))
      (is (close? [4.0 0.0 0.0] (in/transform-point m [0 1 0]) 1e-12)))))

(deftest fractal-test
  (testing "each level of a fractal multiplies its pieces by its copies"
    (doseq [d (range 1 5)]
      (is (= (long (Math/pow 4 d)) (count (:primitives (in/expand (in/sierpinski) {:max-depth d})))))
      (is (= (long (Math/pow 20 d)) (count (:primitives (in/expand (in/menger) {:max-depth d})))))))
  (testing "the pieces stay inside the bound of the whole"
    (doseq [{:keys [transform]} (:primitives (in/expand (in/menger) {:max-depth 3}))]
      ;; A box proxy's corners, at the far side of its own bound.
      (is (<= (v/length (in/transform-point transform [1 1 1])) (+ (Math/sqrt 3.0) 1e-9)))))
  (testing "Sierpinski's pieces shrink by half a level"
    (let [{:keys [primitives]} (in/expand (in/sierpinski) {:max-depth 3})]
      (is (every? #(< (abs (- 0.125 (in/stretch (:transform %)))) 1e-12) primitives)))))

(deftest level-of-detail-test
  (testing "nothing smaller than the cutoff is expanded"
    (let [{:keys [primitives]} (in/expand (in/menger) {:min-size 0.1})]
      ;; sqrt 3 / 27 = 0.064 < 0.1 < sqrt 3 / 9 = 0.19: three levels.
      (is (= 8000 (count primitives)))))
  (testing "measured from an eye, the near side gets the detail"
    (let [{:keys [primitives]} (in/expand (in/menger) {:min-size 0.02 :eye [0.0 0.0 4.0]})
          size-at (fn [z] (->> primitives
                               (filter #(let [[_ _ _ _ _ _ _ _ _ _ _ pz] (:transform %)] (< (abs (- pz z)) 0.2)))
                               (map #(in/stretch (:transform %)))
                               (apply min)))]
      (is (< (size-at 0.9) (size-at -0.9)))))
  (testing "a budget leaves the finest detail undone everywhere, not one corner"
    ;; Depth first would leave a corner at depth 6 and the rest at 1; best
    ;; first is all at 3, but for the few begun at 4 when the money ran out.
    (let [{:keys [primitives]} (in/expand (in/menger) {:limit 10000})
          depths (frequencies (map :depth primitives))]
      (is (<= (count primitives) 10100))
      (is (= #{3 4} (set (keys depths))))
      (is (< (get depths 4) (get depths 3))))))

(deftest tree-test
  (testing "a branch is the same branch however much of the tree is expanded"
    ;; The limbs of a shallow expansion are among those of a deep one.
    (let [limbs (fn [d] (->> (in/expand (in/tree {:seed 4}) {:max-depth d}) :primitives
                             (filter #(= :cylinder (:shape %)))
                             (map :transform) set))]
      (is (every? (limbs 9) (limbs 4)))))
  (testing "seeds make different trees"
    (is (not= (in/expand (in/tree {:seed 1}) {:max-depth 4})
              (in/expand (in/tree {:seed 2}) {:max-depth 4}))))
  (testing "a tree stays inside its bound"
    (let [{:keys [primitives]} (in/expand (in/tree) {})
          radius (:radius (:branch (:rules (in/tree))))]
      (is (every? #(<= (v/length (in/transform-point (:transform %) [0 0 0])) radius) primitives)))))

(deftest lsystem-test
  (testing "every symbol is rewritten at once"
    (is (= "ABAAB" (ls/derive-string {:axiom "A" :rules {\A "AB" \B "A"}} 3))))
  (testing "a stochastic system is repeatable by seed, and varies with it"
    (let [sys (:system (first (filter #(= "shrub" (:name %)) ls/plants)))]
      (is (= (ls/derive-string sys 6 {:seed 3}) (ls/derive-string sys 6 {:seed 3})))
      (is (not= (ls/derive-string sys 6 {:seed 3}) (ls/derive-string sys 6 {:seed 4})))))
  (testing "the turtle draws a limb per F, turns, and returns from a branch"
    (let [prims (ls/interpret "F[+F]F" {:angle 90 :step 1.0})
          ends (map #(in/transform-point (:transform %) [0 1 0]) prims)]
      (is (= 3 (count prims)))
      (is (close? [0.0 1.0 0.0] (nth ends 0) 1e-12))
      ;; Turned left about up: toward -x, the turtle's left.
      (is (close? [-1.0 1.0 0.0] (nth ends 1) 1e-12))
      (is (close? [0.0 2.0 0.0] (nth ends 2) 1e-12))))
  (testing "tropism bends the heading toward the pull"
    (let [end (fn [e] (in/transform-point (:transform (last (ls/interpret "+FFFF" {:angle 90 :tropism [0.0 -1.0 0.0] :elasticity e})))
                                          [0 1 0]))]
      (is (< (second (end 0.3)) (second (end 0.0))))))
  (testing "every plant grows"
    (doseq [p ls/plants]
      (is (< 100 (count (:primitives (ls/grow p {})))) (:name p)))))
