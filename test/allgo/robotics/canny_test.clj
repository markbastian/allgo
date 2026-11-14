(ns allgo.robotics.canny-test
  (:require [allgo.geometry.vec3 :as v]
            [allgo.numerics.multivariate :as mv]
            [allgo.robotics.canny :as canny]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

;; Canny's roadmap, after LaValle, *Planning Algorithms*, Cambridge
;; University Press, 2006, section 6.4.3: Example 6.7's sphere, Example
;; 6.8's torus, and a surface of genus two made up here -- a tube around
;; the lemniscate of Gerono, x^2 (1 - x^2) = y^2.

(def ^:private x (mv/variable 3 0))
(def ^:private y (mv/variable 3 1))
(def ^:private z (mv/variable 3 2))

(def sphere (mv/add (mv/power x 2) (mv/power y 2) (mv/power z 2) (mv/constant 3 -1)))

(defn torus [R r]
  (let [s (mv/add (mv/power x 2) (mv/power y 2) (mv/power z 2) (mv/constant 3 (- (* R R) (* r r))))]
    (mv/sub (mv/power s 2) (mv/scale (mv/add (mv/power x 2) (mv/power y 2)) (* 4 R R)))))

(def double-torus
  (mv/add (mv/power (mv/sub (mv/mul (mv/power x 2) (mv/sub (mv/constant 3 1) (mv/power x 2))) (mv/power y 2)) 2)
          (mv/power z 2)
          (mv/constant 3 -0.01)))

(defn- near? [a b] (< (v/distance a b) 1e-8))

(deftest sphere-test
  (testing "Example 6.7: with g = (x1, x2) the silhouette is the equator, and x1's critical points are (-1, 0, 0) and (1, 0, 0)"
    (let [rm (canny/roadmap sphere [[-1.5 -1.5 -1.5] [1.5 1.5 1.5]] {})
          [sil :as sils] (:silhouettes rm)]
      (is (= 1 (count sils)))
      (is (:closed? sil))
      (is (every? (fn [[_ _ c :as p]] (and (< (abs c) 1e-9) (< (abs (- 1 (v/length p))) 1e-9))) (:points sil)))
      (is (= 2 (count (:critical rm))))
      (is (some #(near? [-1.0 0.0 0.0] (:point %)) (:critical rm)))
      (is (some #(near? [1.0 0.0 0.0] (:point %)) (:critical rm)))
      (testing "the slices through them are single points, adding nothing: the equator is the roadmap"
        (is (every? (comp empty? :curves) (:slices rm)))
        (is (= 1 (count (canny/components rm))))))))

(deftest torus-test
  (let [bounds [[-1.6 -1.6 -0.6] [1.6 1.6 0.6]]
        alone (canny/roadmap (torus 1.0 0.4) bounds {:slices? false})
        rm (canny/roadmap (torus 1.0 0.4) bounds {})]
    (testing "Example 6.8: the silhouettes are the outer and inner rings, with four critical points"
      (is (= 2 (count (:silhouettes rm))))
      (is (= #{1.4 0.6} (set (map (fn [{:keys [points]}] (/ (math/round (* 1e6 (v/length (first points)))) 1e6)) (:silhouettes rm)))))
      (is (= 4 (count (:critical rm))))
      (doseq [c [-1.4 -0.6 0.6 1.4]]
        (is (some #(near? [c 0.0 0.0] (:point %)) (:critical rm)))))
    (testing "Figure 6.37: alone they leave the rings apart"
      (is (= 2 (count (canny/components alone)))))
    (testing "Figures 6.38 and 6.39: the slices through the inner critical points join them"
      (is (= 1 (count (canny/components rm))))
      (is (= 4 (count (:junctions rm))))
      (doseq [{:keys [point]} (:junctions rm)]
        ;; where the plane x = 0.6 or -0.6 meets the outer ring, radius 1.4
        (is (< (abs (- (abs (point 0)) 0.6)) 1e-9))
        (is (< (abs (- (abs (point 1)) (math/sqrt (- (* 1.4 1.4) (* 0.6 0.6))))) 1e-9)))
      (testing "the slice through an inner critical point is a figure eight: each branch ends at it"
        (doseq [{:keys [c curves]} (:slices rm) :when (< (abs c) 1)
                {:keys [points]} curves]
          (is (near? [c 0.0 0.0] (first points)))
          (is (near? [c 0.0 0.0] (peek points))))))))

(deftest route-test
  (testing "queries declared critical: a way from the top of the torus to the underside of its inner ring, on the surface throughout"
    (let [f (torus 1.0 0.4)
          s (canny/surface f)
          a (canny/onto-surface s [0.1 1.0 0.5])
          b (canny/onto-surface s [-0.6 -0.2 -0.3])
          rm (canny/roadmap f [[-1.6 -1.6 -0.6] [1.6 1.6 0.6]] {:queries [a b]})
          path (canny/route rm 0 1)]
      (is (= a (first path)))
      (is (= b (peek path)))
      (is (every? #(< (abs (canny/value s %)) 1e-9) path))
      (is (every? (fn [[p q]] (< (v/distance p q) (* 2 (:step rm)))) (partition 2 1 path))))))

(deftest general-position-test
  (testing "a turned frame is orthonormal"
    (let [{:keys [sweep across view]} (canny/turned 0.05 0.03 0.07)]
      (is (every? #(< (abs (- 1 (v/length %))) 1e-12) [sweep across view]))
      (is (every? #(< (abs %) 1e-12) [(v/dot sweep across) (v/dot sweep view) (v/dot across view)]))))
  (testing "genus two: three silhouette loops and six critical points; the slices join the loops"
    (let [bounds [[-1.3 -0.8 -0.3] [1.3 0.8 0.3]]
          frame (canny/turned 0.05 0.03 0.07)
          alone (canny/roadmap double-torus bounds {:frame frame :slices? false})
          rm (canny/roadmap double-torus bounds {:frame frame})]
      (is (= 3 (count (:silhouettes rm))))
      (is (= 6 (count (:critical rm))))
      (is (= 3 (count (canny/components alone))))
      (is (= 1 (count (canny/components rm)))))))
