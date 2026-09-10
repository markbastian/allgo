(ns procedurals.boids-test
  (:require [procedurals.boids :as b]
            [clojure.test :refer [deftest is testing]]))

(defn- seeded-flock
  "Build a flock from a fixed seed. `b/flock` draws on `rand`, which would
  make the alignment threshold below flaky from run to run."
  [seed n bounds max-speed]
  (let [g (java.util.Random. seed)]
    (vec (repeatedly n
                     #(let [dir (mapv (fn [_] (- (.nextDouble g) 0.5)) bounds)]
                        {:pos (mapv (fn [hi] (* hi (.nextDouble g))) bounds)
                         :vel (b/with-magnitude dir (* max-speed (+ 0.5 (* 0.5 (.nextDouble g)))))})))))

(defn- order
  "Polarisation: 0 when headings are scattered, 1 when the flock is aligned."
  [flock]
  (b/mag (b/v* (reduce b/v+ (map (comp b/normalize :vel) flock)) (/ 1.0 (count flock)))))

(deftest flocking-emerges
  (testing "an initially scattered flock aligns"
    (let [bounds [400 400]
          f0     (seeded-flock 20260910 60 bounds (:max-speed b/defaults))
          run    (iterate #(b/step % bounds) f0)
          start  (order f0)
          end    (order (nth run 600))]
      (is (< start 0.4) (str "starts disordered, was " start))
      (is (> end 0.6) (str "ends aligned, was " end))
      (is (> end (* 2 start)) "alignment increases substantially"))))

(deftest flock-stays-well-formed
  (let [bounds [400 400]
        final  (nth (iterate #(b/step % bounds) (seeded-flock 1 40 bounds (:max-speed b/defaults))) 200)]
    (is (= 40 (count final)))
    (testing "wrapping keeps every boid inside the world"
      (doseq [{[x y] :pos} final]
        (is (and (<= 0 x 400) (<= 0 y 400)))))
    (testing "speed never exceeds the limit"
      (doseq [{:keys [vel]} final]
        (is (<= (b/mag vel) (+ (:max-speed b/defaults) 1e-9)))))
    (testing "no boid goes non-finite"
      (doseq [{:keys [pos vel]} final]
        (is (every? #(Double/isFinite (double %)) (concat pos vel)))))))

(deftest rules-are-dimension-agnostic
  (testing "the same rules drive a 3D flock"
    (let [bounds [100 100 100]
          params {:perception-radius 22 :separation-radius 9 :max-speed 0.8 :max-force 0.02}
          f0     (seeded-flock 2 40 bounds (:max-speed params))
          final  (nth (iterate #(b/step % bounds params) f0) 400)]
      (is (= 3 (count (:pos (first f0)))))
      (is (> (order final) (order f0)) "alignment increases in 3D too")
      (doseq [{:keys [pos]} final]
        (is (every? #(<= 0 % 100) pos))))))

(deftest step-does-not-depend-on-update-order
  (testing "every boid sees the same previous state, so reversing the flock only permutes the result"
    ;; Compared within a tolerance, not bitwise. The invariant is real, but
    ;; reversing the flock also reverses the order neighbour contributions
    ;; are summed in, and floating-point addition is not associative -- the
    ;; two answers can differ in the last unit in the last place.
    (let [bounds [400 400]
          f      (seeded-flock 4242 30 bounds (:max-speed b/defaults))
          forward (sort-by :pos (b/step f bounds))
          reverse* (sort-by :pos (b/step (vec (reverse f)) bounds))]
      (is (= (count forward) (count reverse*)))
      (doseq [[x y] (map vector forward reverse*)]
        (is (< (b/mag (b/v- (:pos x) (:pos y))) 1e-9) "same positions")
        (is (< (b/mag (b/v- (:vel x) (:vel y))) 1e-9) "same velocities")))))
