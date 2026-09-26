(ns allgo.cellular-test
  (:require [allgo.geometry.vec3 :as v3]
            [allgo.procedural.cellular :as c]
            [clojure.test :refer [deftest is testing]]))

(defn- points [n]
  (for [i (range n)] [(* 0.1379 i) (* 0.2713 i) (* 0.3137 i)]))

(deftest feature-points-test
  (testing "a cell's points lie in that cell and nowhere else"
    (doseq [i (range -2 3) j (range -2 3) k (range -2 3)]
      (doseq [[x y z _] (c/feature-points {:seed 4} i j k)]
        (is (<= i x (inc i)))
        (is (<= j y (inc j)))
        (is (<= k z (inc k))))))

  (testing "each carries a value in [-1, 1]"
    (doseq [[_ _ _ v] (c/feature-points {:seed 4} 0 0 0)]
      (is (<= -1.0 v 1.0))))

  (testing "jitter pulls them toward the center of the cell"
    ;; Jitter zero is a regular lattice, which is the thing the Poisson
    ;; scatter exists to avoid -- but it is also how you get a honeycomb
    ;; on purpose, so it has to work.
    (doseq [[x y z _] (c/feature-points {:seed 4 :jitter 0.0} 3 -1 2)]
      (is (== 3.5 x)) (is (== -0.5 y)) (is (== 2.5 z))))

  (testing "the points of a cell depend on the cell and the seed"
    (is (= (c/feature-points {:seed 4} 1 2 3) (c/feature-points {:seed 4} 1 2 3)))
    (is (not= (c/feature-points {:seed 4} 1 2 3) (c/feature-points {:seed 5} 1 2 3)))))

(deftest nearest-test
  (let [near (c/nearest {:seed 7 :n 3})]
    (testing "distances come back ascending"
      (doseq [[x y z] (points 500)]
        (let [ds (:distances (near x y z))]
          (is (= ds (sort ds))))))

    (testing "the nearest distance really is to the point reported"
      (doseq [[x y z] (points 300)]
        (let [{:keys [distances point]} (near x y z)]
          (is (< (abs (- (double (first distances))
                         (v3/distance [x y z] point)))
                 1e-9)))))

    (testing "and no feature point of the neighborhood is nearer"
      (doseq [[x y z] (points 200)]
        (let [d1 (double (first (:distances (near x y z))))
              ci (long (Math/floor x)) cj (long (Math/floor y)) ck (long (Math/floor z))
              all (for [di [-1 0 1] dj [-1 0 1] dk [-1 0 1]
                        [px py pz _] (c/feature-points {:seed 7} (+ ci di) (+ cj dj) (+ ck dk))]
                    (v3/distance [x y z] [px py pz]))]
          (is (< (abs (- d1 (apply min all))) 1e-9)))))))

(deftest metrics-test
  (testing "each metric agrees with its definition"
    (let [{:keys [euclidean manhattan chebyshev]} c/metrics]
      (is (< (abs (- 5.0 (euclidean 3.0 4.0 0.0))) 1e-12))
      (is (== 7.0 (manhattan 3.0 -4.0 0.0)))
      (is (== 4.0 (chebyshev 3.0 -4.0 1.0)))))

  (testing "and produces a different texture"
    (let [vs (fn [m] (mapv (fn [[x y z]] ((c/f1 {:seed 2 :metric m}) x y z)) (points 300)))]
      (is (not= (vs :euclidean) (vs :manhattan)))
      (is (not= (vs :euclidean) (vs :chebyshev))))))

(deftest bases-test
  (testing "F2 is never nearer than F1, so F2-F1 never goes below zero"
    (let [near (c/nearest {:seed 3 :n 2})]
      (doseq [[x y z] (points 1000)]
        (let [[d1 d2] (:distances (near x y z))]
          (is (<= (double d1) (double d2)))))))

  (testing "the scalar bases are centered and land near [-1, 1]"
    (doseq [b [(c/f1 {:seed 3}) (c/f2 {:seed 3}) (c/id-basis {:seed 3})]]
      (let [vs (mapv (fn [[x y z]] (b x y z)) (points 4000))]
        (is (> (apply min vs) -1.5))
        (is (< (apply max vs) 1.5))
        (is (< (abs (/ (reduce + vs) (count vs))) 0.15)))))

  (testing "the density changes cell size without changing the range"
    ;; Dividing by the expected spacing is what buys this: turning the
    ;; density up should make the pattern finer, not darker.
    (doseq [d [1.0 3.0 10.0]]
      (let [vs (mapv (fn [[x y z]] ((c/f1 {:seed 3 :density d}) x y z)) (points 2000))]
        (is (< 1.0 (- (apply max vs) (apply min vs)) 3.0)))))

  (testing "combination reproduces the named bases it generalizes"
    (let [opts {:seed 3 :density 3.0}
          comb (c/combination opts [-1.0 1.0])
          raw (c/distance-basis (assoc opts :n 2))]
      (doseq [[x y z] (points 200)]
        (let [[d1 d2] (raw x y z)]
          (is (< (abs (- (comb x y z) (- (double d2) (double d1)))) 1e-9))))))

  (testing "the id basis is flat inside a cell and jumps between them"
    (let [id (c/id-basis {:seed 3 :density 1.0})
          along (mapv (fn [i] (id (* 0.002 i) 0.3 0.3)) (range 2000))]
      (is (< (count (set along)) 200))
      (is (> (count (set along)) 1)))))

(deftest expected-spacing-test
  (testing "the closed form matches what the scatter actually does"
    ;; If these drifted apart, every cellular basis would quietly shift
    ;; off center as the density changed.
    (doseq [density [1.0 3.0 8.0]]
      (let [near (c/nearest {:seed 11 :density density :n 1})
            mean (/ (reduce + (map (fn [[x y z]] (first (:distances (near x y z))))
                                   (points 4000)))
                    4000.0)]
        (is (< (abs (- 1.0 (/ mean (c/expected-spacing density)))) 0.1))))))
