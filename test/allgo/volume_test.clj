(ns allgo.volume-test
  (:require [allgo.procedural.volume :as vol]
            [clojure.test :refer [deftest is testing]]))

(defn- close? [a b tol] (every? #(< (abs (double %)) tol) (map - a b)))

(def ^:private bg [0.2 0.4 0.6])

(deftest catalog-test
  (testing "every density is finite, never negative, and gone by the edge of the unit ball"
    (doseq [{:keys [name make]} vol/catalog]
      (let [{:keys [density]} (make {})]
        (doseq [i (range 400)
                :let [r (* 0.0025 i)
                      a (* 0.7 i)
                      p [(* r (Math/cos a)) (* r (Math/sin (* 1.3 a))) (* r (Math/sin a) 0.5)]]]
          (let [d (apply density p)]
            (is (and (Double/isFinite d) (>= d 0.0)) name)))
        (doseq [i (range 200)
                :let [a (* 0.31 i) b (* 0.17 i)
                      p [(* 0.999 (Math/cos a) (Math/cos b)) (* 0.999 (Math/sin b)) (* 0.999 (Math/sin a) (Math/cos b))]]]
          (is (< (apply density p) 0.05) name))))))

(deftest march-test
  (let [ray (vol/camera-ray 0.0 0.0 3.0 0.8 0.0 0.0)]
    (testing "a ray that misses the ball, or crosses nothing, is the background"
      (is (= bg (vol/march (vol/noisy-sphere) (vol/camera-ray 0.0 0.0 3.0 0.8 0.95 0.95) {:background bg})))
      (is (close? bg (vol/march {:density (constantly 0.0) :extinction 1.0 :albedo [1 1 1]} ray
                                {:background bg})
                  1e-12)))

    (testing "through a uniform absorbing gas, the background is dimmed by Beer's law"
      ;; Straight through the middle is a diameter, 2 long.
      (let [gas {:density (constantly 1.0) :extinction 1.0 :albedo [0.0 0.0 0.0]}]
        (is (close? (mapv #(* % (Math/exp -2.0)) bg) (vol/march gas ray {:background bg :steps 200})
                    1e-9))))

    (testing "a solid's surface is found where it is, wherever the samples fall"
      ;; Lit by bisection rather than at whichever step landed inside, so
      ;; the jitter moves nothing: a hard ball looks the same at any.
      (let [ball {:density (fn [x y z] (if (< (+ (* x x) (* y y) (* z z)) 0.25) 1.0 0.0))
                  :extinction 15.0 :albedo [0.8 0.8 0.8] :surface 1.0}
            off-center (vol/camera-ray 0.0 0.0 3.0 0.8 0.2 0.1)
            at (fn [j] (vol/march ball off-center {:background bg :jitter j}))]
        (is (close? (at 0.05) (at 0.95) 1e-3))))

    (testing "fire glows hottest in its heart"
      (let [fire (vol/fireball)
            lum (fn [[r g b]] (+ r g b))]
        ;; The middle is opaque, so it is the heat of the front of it that
        ;; is seen; still more than the thin fringe gives.
        (is (> (lum (vol/march fire ray {}))
               (* 1.5 (lum (vol/march fire (vol/camera-ray 0.0 0.0 3.0 0.8 0.36 0.0) {})))))))

    (testing "a cloud lit through itself is darker underneath than on top"
      (let [cloud (vol/cumulus)
            top (vol/march cloud (vol/camera-ray 0.0 0.0 3.0 0.8 0.0 0.08) {:light [0.0 1.0 0.0] :ambient 0.1})
            bottom (vol/march cloud (vol/camera-ray 0.0 0.0 3.0 0.8 0.0 -0.06) {:light [0.0 1.0 0.0] :ambient 0.1})]
        (is (> (first top) (first bottom)))))))

(deftest camera-test
  (testing "the middle of the picture looks at the origin, from wherever the camera is"
    (doseq [[az el] [[0.0 0.0] [0.7 0.3] [-2.0 -0.5]]]
      (let [[o d] (vol/camera-ray az el 3.0 0.8 0.0 0.0)]
        (is (close? (mapv #(/ (- %) 3.0) o) d 1e-12)))))
  (testing "up on the picture is up in the world"
    (let [[_ d] (vol/camera-ray 0.4 0.0 3.0 0.8 0.0 1.0)]
      (is (pos? (second d))))))

(deftest jitter-test
  (testing "pixel jitter is in [0, 1) and differs from a pixel to the next"
    (let [js (for [x (range 40) y (range 40)] (vol/pixel-jitter x y))]
      (is (every? #(and (<= 0.0 %) (< % 1.0)) js))
      (is (> (count (distinct js)) 1500)))))
