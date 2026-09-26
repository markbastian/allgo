(ns allgo.boids-flat-test
  "The flat-array flocking must agree with the reference implementation.

  That is the whole justification for keeping two of them: one is the
  readable definition of the rules, the other is fast, and the only thing
  that makes the fast one trustworthy is that it can be checked against
  the slow one."
  (:require [allgo.simulation.boids :as boids]
            [allgo.simulation.boids-flat :as flat]
            [allgo.simulation.flock :as flock]
            [clojure.test :refer [deftest is testing]]))

(defn- seeded-flock
  "A flock from a fixed seed, in place of `boids/flock`.

  `boids/flock` draws on `rand`, and a comparison to a tolerance is only
  a test if the scene it runs on is the same every time. It is not a
  detail here: on an unlucky draw two boids sit close enough to the
  perception radius that the two implementations round the same distance
  to opposite sides of it, take different neighbor sets, and part
  company by far more than double precision would explain. Seeding it
  means a failure is a real disagreement rather than that draw coming
  up."
  [seed n bounds params]
  (let [g (java.util.Random. seed)
        max-speed (:max-speed (merge boids/defaults params))]
    (vec (repeatedly n
                     #(let [dir (mapv (fn [_] (- (.nextDouble g) 0.5)) bounds)]
                        {:pos (mapv (fn [hi] (* hi (.nextDouble g))) bounds)
                         :vel (boids/with-magnitude
                                dir
                                (* max-speed (+ 0.5 (* 0.5 (.nextDouble g)))))})))))

(defn- run-both
  "Steps both implementations from the same start and returns the largest
  disagreement in any coordinate."
  [n bounds params ticks]
  (let [reference (seeded-flock 20260919 n bounds params)
        fast      (flat/from-boids reference)]
    (loop [r reference f fast i 0]
      (if (= i ticks)
        (flock/divergence r f)
        (recur (boids/step r bounds params) (flat/step f bounds params) (inc i))))))

;; Different summation order, so the two drift apart at the level of
;; double-precision rounding and no faster.
(def ^:private tolerance 1e-9)

(deftest conversion-test
  (let [reference (boids/flock 40 [300.0 300.0] {})
        fast      (flat/from-boids reference)]
    (testing "converting to flat arrays and back is lossless"
      (is (= 40 (:n fast)))
      (is (= 2 (:dims fast)))
      (is (zero? (flock/divergence reference fast))))

    (testing "three dimensions survive too"
      (let [r3 (boids/flock 25 [100.0 100.0 100.0] {})
            f3 (flat/from-boids r3)]
        (is (= 3 (:dims f3)))
        (is (zero? (flock/divergence r3 f3)))))))

(deftest agrees-with-reference-test
  (testing "wrapping edges, in two dimensions"
    (is (< (run-both 120 [400.0 400.0] {:edges :wrap} 30) tolerance)))

  (testing "bouncing edges, which is the fiddlier path"
    (is (< (run-both 120 [400.0 400.0] {:edges :bounce} 30) tolerance)))

  (testing "three dimensions"
    (is (< (run-both 80 [200.0 200.0 200.0] {:edges :wrap} 25) tolerance))
    (is (< (run-both 80 [200.0 200.0 200.0] {:edges :bounce} 25) tolerance)))

  (testing "across weightings, including ones that switch a rule off"
    (doseq [params [{:separation-weight 0.0}
                    {:alignment-weight 0.0}
                    {:cohesion-weight 0.0}
                    {:separation-weight 3.0 :cohesion-weight 0.1}
                    {:perception-radius 20.0}
                    {:perception-radius 200.0 :separation-radius 90.0}
                    {:max-speed 0.5 :max-force 0.5}]]
      (is (< (run-both 60 [300.0 300.0] (merge {:edges :wrap} params) 20) tolerance)
          (pr-str params))))

  (testing "a crowd dense enough that everything sees everything"
    (is (< (run-both 150 [60.0 60.0] {:edges :wrap :perception-radius 200.0} 15)
           tolerance)))

  (testing "and a single boid, which has no neighbors at all"
    (is (< (run-both 1 [300.0 300.0] {:edges :wrap} 20) tolerance))))

(deftest obstacles-test
  (let [obstacles [(boids/sphere-obstacle [250.0 260.0] 72.0)
                   (boids/box-obstacle [40.0 350.0] [150.0 440.0])]
        params    {:edges :bounce :obstacles obstacles :boid-radius 4.0
                   :avoid-weight 2.2}]
    (testing "steering around barriers matches the reference"
      ;; Checked over a few ticks rather than many. GJK and EPA are
      ;; iterative, so the two implementations part company in the last
      ;; bits, and flocking amplifies that: by forty ticks they have
      ;; drifted to 1e-7, which says nothing about either being wrong.
      (is (< (run-both 80 [500.0 500.0] params 5) 1e-9)))

    (testing "and the flat flock really does avoid them"
      (let [f (flat/simulate (flat/from-boids (seeded-flock 20260919 150 [500.0 500.0] params))
                             [500.0 500.0] params 80)
            inside (fn [{[x y] :pos}]
                     (or (< (Math/hypot (- x 250.0) (- y 260.0)) 68.0)
                         (and (< 44.0 x 146.0) (< 354.0 y 436.0))))]
        (is (zero? (count (filter inside (flat/to-boids f))))
            "no boid ends up inside a barrier")))

    (testing "with no obstacles the avoidance path is skipped entirely"
      (is (< (run-both 60 [400.0 400.0] {:edges :bounce} 20) tolerance)))))

(deftest protocol-test
  (let [params    {:edges :wrap}
        reference (flock/reference 50 [300.0 300.0] params)
        fast      (flat/from-boids (flock/as-boids reference))]
    (testing "both answer the same interface"
      (is (= 50 (flock/flock-size reference) (flock/flock-size fast)))
      (is (= 50 (count (flock/as-boids fast))))
      (is (every? #(and (:pos %) (:vel %)) (flock/as-boids fast))))

    (testing "and can be stepped through it interchangeably"
      (let [r (flock/simulate reference [300.0 300.0] params 10)
            f (flock/simulate fast [300.0 300.0] params 10)]
        (is (< (flock/divergence r f) tolerance))
        (is (= 50 (flock/flock-size f)))))

    (testing "divergence is zero against itself and positive against a change"
      (is (zero? (flock/divergence reference reference)))
      (is (pos? (flock/divergence reference
                                  (flock/advance reference [300.0 300.0] params)))))))

(deftest flat-mechanics-test
  (let [f (flat/from-boids (boids/flock 30 [200.0 200.0] {}))]
    (testing "positions come back as one flat array of three per boid"
      (is (= 90 (alength (flat/positions f)))))

    (testing "stepping does not lose or duplicate boids"
      (let [g (flat/simulate f [200.0 200.0] {:edges :wrap} 40)]
        (is (= 30 (:n g)))
        (is (= 30 (count (flat/to-boids g))))))

    (testing "bouncing keeps every boid inside the world"
      (let [g (flat/simulate f [200.0 200.0] {:edges :bounce} 60)]
        (is (every? (fn [{[x y] :pos}] (and (<= 0.0 x 200.0) (<= 0.0 y 200.0)))
                    (flat/to-boids g)))))

    (testing "wrapping keeps them inside too"
      (let [g (flat/simulate f [200.0 200.0] {:edges :wrap} 60)]
        (is (every? (fn [{[x y] :pos}] (and (<= 0.0 x 200.0) (<= 0.0 y 200.0)))
                    (flat/to-boids g)))))))
