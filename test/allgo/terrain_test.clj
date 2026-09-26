(ns allgo.terrain-test
  (:require [allgo.procedural.terrain :as t]
            [clojure.test :refer [deftest is testing]]))

(defn- seeded
  "A repeatable stand-in for `rand`."
  [seed]
  (let [r (java.util.Random. seed)]
    (fn [] (.nextDouble r))))

(def ^:private corners [0.0 0.5 0.25 0.75])

(defn- terrain
  ([iterations] (terrain iterations 1.0 7))
  ([iterations width seed]
   (t/generate {:width width :iterations iterations
                :corners corners :rng (seeded seed)})))

(deftest size-test
  (testing "each round doubles the resolution, less the shared edge"
    (doseq [n (range 0 9)]
      (is (= (inc (long (Math/pow 2 n))) (:dim (terrain n))))))

  (testing "and the displacement halves with it"
    ;; The fractal part: each round adds a smaller bump than the last, which
    ;; is why the result reads as landscape and not as noise.
    (let [widths (map :width (take 6 (iterate t/step (t/init-grid {:width 1.0}))))]
      (is (= [1.0 0.5 0.25 0.125 0.0625 0.03125] (vec widths))))))

(deftest corners-test
  (testing "the corners given are the corners you get"
    ;; Nothing in the algorithm ever revisits a point once it has a height,
    ;; so the seeds survive to the last round.
    (doseq [n [1 2 5 7]]
      (let [g (terrain n)
            d (dec (:dim g))
            [tl tr br bl] corners]
        (is (== tl (t/height g 0 0)))
        (is (== tr (t/height g 0 d)))
        (is (== bl (t/height g d 0)))
        (is (== br (t/height g d d)))))))

(deftest determinism-test
  (testing "the same seed gives the same terrain"
    (is (= (t/cells->grid (terrain 6 1.0 42))
           (t/cells->grid (terrain 6 1.0 42)))))

  (testing "and a different one does not"
    (is (not= (t/cells->grid (terrain 6 1.0 42))
              (t/cells->grid (terrain 6 1.0 43))))))

(deftest purity-test
  ;; The grid holds a mutable array, so this is worth stating: a step must
  ;; not write into the grid it was given, or `iterate` would hand back a
  ;; row of aliases rather than a sequence of terrains.
  (testing "a step leaves the grid it was given alone"
    (let [g (terrain 4)
          before (t/cells->grid g)]
      (dotimes [_ 3] (t/step g))
      (is (= before (t/cells->grid g)))))

  (testing "so intermediate grids from iterate are independent"
    (let [gs (vec (take 5 (iterate t/step (t/init-grid {:width 1.0 :corners corners
                                                        :rng (seeded 3)}))))
          snapshots (mapv t/cells->grid gs)]
      (is (= [2 3 5 9 17] (mapv :dim gs)))
      ;; Force more work through the same arrays, then check nothing moved.
      (dotimes [_ 3] (t/step (peek gs)))
      (is (= snapshots (mapv t/cells->grid gs))))))

(deftest accessors-test
  (let [g (terrain 5)
        dim (:dim g)
        grid (t/cells->grid g)]

    (testing "cells->grid is row major and agrees with height"
      (is (= dim (count grid)))
      (is (every? #(= dim (count %)) grid))
      (is (every? true?
                  (for [i (range dim) j (range dim)]
                    (== (get-in grid [i j]) (t/height g i j))))))

    (testing "the sparse map agrees too"
      (let [m (t/cells g)]
        (is (= (* dim dim) (count m)))
        (is (every? true? (for [i (range dim) j (range dim)]
                            (== (m [i j]) (t/height g i j)))))))

    (testing "heights-seq walks every cell once, row by row"
      (is (= (apply concat grid) (t/heights-seq g))))

    (testing "bounds bracket every height"
      (let [[lo hi] (t/bounds g)
            all (t/heights-seq g)]
        (is (== lo (apply min all)))
        (is (== hi (apply max all)))))))

(deftest roughness-test
  (testing "a bigger displacement gives rougher ground"
    (let [spread (fn [width]
                   (let [[lo hi] (t/bounds (terrain 6 width 11))] (- hi lo)))]
      (is (< (spread 0.1) (spread 1.0) (spread 10.0)))))

  (testing "no displacement leaves a surface interpolated from the corners"
    ;; Every new point is then exactly the average of its neighbors, so
    ;; nothing can land outside the range the corners set.
    (let [g (t/generate {:width 0.0 :iterations 6 :corners corners :rng (seeded 5)})
          [lo hi] (t/bounds g)]
      (is (<= (- (apply min corners) 1e-12) lo))
      (is (<= hi (+ (apply max corners) 1e-12))))))

(deftest continuity-test
  (testing "neighboring cells stay close, which is what makes it terrain"
    ;; Each round's displacement is half the last, so the total any one
    ;; step can add is bounded by twice the first. A surface that failed
    ;; this would be noise with a fractal's parameters.
    (let [g (terrain 7 1.0 13)
          dim (:dim g)
          worst (reduce max 0.0
                        (for [i (range dim) j (range (dec dim))]
                          (abs (- (t/height g i j) (t/height g i (inc j))))))]
      (is (< worst 2.0) (str "neighbors differ by up to " worst)))))
