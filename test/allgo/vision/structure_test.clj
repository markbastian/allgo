(ns allgo.vision.structure-test
  (:require [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.random :as random]
            [allgo.vision.structure :as structure]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

;; Structure from motion by Tomasi and Kanade's factorization, after
;; Russell and Norvig, *Artificial Intelligence: A Modern Approach*, 4th
;; edition, chapter 25: a known point cloud seen by an orthographic
;; camera turning round it, the shape recovered from the images alone.

(defn- scene [seed n]
  (let [rng (random/rng seed)]
    (vec (repeatedly n #(vector (- (* 2 (rng)) 1) (- (* 2 (rng)) 1) (- (* 2 (rng)) 1))))))

(defn- views
  "Orthographic views of `points` from a camera turning `step` radians a
  frame about a tilted axis, drifting sideways, with image noise of
  standard deviation `sigma`: frames of `[u v]`, and the true axes."
  [points frames step sigma seed]
  (let [rng (random/rng seed)
        axis (v/normalize [0.3 1.0 0.2])]
    (for [f (range frames)
          :let [r (q/from-axis-angle axis (* f step))
                i (q/rotate r [1.0 0.0 0.0]) j (q/rotate r [0.0 1.0 0.0])
                [a b] [(* 0.3 f) (* -0.2 f)]]]
      {:image (mapv (fn [p] [(+ (v/dot i p) a (* sigma (random/gaussian rng 0.0 1.0)))
                             (+ (v/dot j p) b (* sigma (random/gaussian rng 0.0 1.0)))])
                    points)
       :axes [i j]})))

(defn- distances [ps] (for [a (range (count ps)) b (range a)] (v/distance (ps a) (ps b))))

(deftest factorization
  (let [cloud (scene 1 20)
        vs (views cloud 8 0.15 0.0 2)
        states (vec (structure/factorization-seq (map :image vs)))
        {:keys [points axes residual]} (structure/factorize (map :image vs))]
    (testing "a stage a step"
      (is (= [:measured :centered :factored :metric] (map :stage states))))
    (testing "the measurements, less their means, have rank three: nothing left over
      but rounding: the root of eigenvalues zero to 1e-14"
      (is (< residual 1e-6)))
    (testing "the shape recovered exactly, up to a turn and a mirror image: every distance between points the same"
      (is (every? #(< (abs %) 1e-9) (map - (distances cloud) (distances points)))))
    (testing "each frame's camera axes unit length and square to each other, the first frame's x and y"
      (doseq [[i j] axes]
        (is (< (abs (- 1.0 (v/length i))) 1e-9))
        (is (< (abs (- 1.0 (v/length j))) 1e-9))
        (is (< (abs (v/dot i j)) 1e-9)))
      (is (< (v/distance [1.0 0.0 0.0] (first (first axes))) 1e-9)))
    (testing "and the camera's turn between frames the true 0.15 radians"
      (doseq [[[i1 j1] [i2 j2]] (partition 2 1 axes)]
        (let [k1 (v/cross i1 j1) k2 (v/cross i2 j2)
              ;; the rotation taking one frame's axes to the next's
              trace (+ (v/dot i1 i2) (v/dot j1 j2) (v/dot k1 k2))]
          (is (< (abs (- 0.15 (math/acos (/ (- trace 1.0) 2.0)))) 1e-9)))))))

(deftest with-noise
  (testing "noise in the images: a residual, and distances still within a few percent"
    (let [truth (scene 3 30)
          {:keys [points residual]} (structure/factorize (map :image (views truth 12 0.12 0.01 4)))
          errs (map #(abs (/ (- %1 %2) %2)) (distances points) (distances truth))]
      (is (< 0.01 residual))
      (is (< (/ (reduce + errs) (count errs)) 0.02))
      (is (< (apply max errs) 0.2)))))
