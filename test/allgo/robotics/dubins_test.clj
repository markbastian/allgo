(ns allgo.robotics.dubins-test
  (:require [allgo.math :as am]
            [allgo.random :as random]
            [allgo.robotics.dubins :as dubins]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

;; Shortest paths for cars, after LaValle, *Planning Algorithms*,
;; Cambridge University Press, 2006, section 15.3: every path driven to
;; check that it ends where it should.

(defn- end-error [start path rho goal]
  (let [[x y th] (peek (dubins/drive start path rho 1))]
    (max (abs (- x (goal 0))) (abs (- y (goal 1))) (abs (am/wrap-angle (- th (goal 2)))))))

(def ^:private cases
  (let [rng (random/rng 2)]
    (vec (repeatedly 200 (fn [] [[(rng) (rng) (* 2 math/PI (rng))]
                                 [(- (* 8 (rng)) 4) (- (* 8 (rng)) 4) (* 2 math/PI (rng))]])))))

(deftest dubins-test
  (testing "every word the formulas give, driven, ends at the goal"
    (is (every? (fn [[s g]] (every? #(< (end-error s % 1.3 g) 1e-9) (dubins/dubins-paths s g 1.3))) cases)))
  (testing "and there is always one"
    (is (every? (fn [[s g]] (some? (dubins/dubins s g 1.3))) cases)))
  (testing "straight ahead, the path is the straight line"
    (let [{:keys [length]} (dubins/dubins [0.0 0.0 0.0] [10.0 0.0 0.0] 1.0)]
      (is (< (abs (- 10.0 length)) 1e-12))))
  (testing "turning round to face back, a radius to one side: half a circle, pi rho"
    (is (< (abs (- (* math/PI 2.0) (:length (dubins/dubins [0.0 0.0 0.0] [0.0 4.0 math/PI] 2.0)))) 1e-9)))
  (testing "never shorter than the straight-line distance"
    (is (every? (fn [[s g]] (>= (+ 1e-9 (:length (dubins/dubins s g 1.3))) (math/hypot (- (g 0) (s 0)) (- (g 1) (s 1))))) cases))))

(deftest reeds-shepp-test
  (testing "every path the formulas and their symmetries give, driven, ends at the goal"
    (is (every? (fn [[s g]] (every? #(< (end-error s % 1.0 g) 1e-9) (dubins/reeds-shepp-paths s g 1.0))) cases)))
  (testing "allowed to reverse, the car never needs longer than Dubins' car"
    (is (every? (fn [[s g]] (<= (:length (dubins/reeds-shepp s g 1.0)) (+ 1e-9 (:length (dubins/dubins s g 1.0))))) cases)))
  (testing "a goal three metres behind is three metres away, in reverse"
    (let [{:keys [length lengths]} (dubins/reeds-shepp [0.0 0.0 0.0] [-3.0 0.0 0.0] 1.0)]
      (is (< (abs (- 3.0 length)) 1e-9))
      (is (neg? (reduce + lengths)))))
  (testing "where Dubins' car must loop round, reversing is far shorter: the goal one metre behind needs 1, not 2 pi + 1"
    (is (< (abs (- 1.0 (:length (dubins/reeds-shepp [0.0 0.0 0.0] [-1.0 0.0 0.0] 1.0)))) 1e-9))
    (is (> (:length (dubins/dubins [0.0 0.0 0.0] [-1.0 0.0 0.0] 1.0)) 6.0))))
