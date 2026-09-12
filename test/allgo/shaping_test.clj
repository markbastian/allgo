(ns allgo.shaping-test
  (:require [allgo.procedural.shaping :as s]
            [clojure.test :refer [deftest is testing]]))

(defn- close? [a b] (< (abs (- (double a) (double b))) 1e-9))

(deftest selectors-test
  (testing "step is the hard edge, and it belongs to the upper side"
    (is (== 0.0 (s/step 0.5 0.4)))
    (is (== 1.0 (s/step 0.5 0.5)))
    (is (== 1.0 (s/step 0.5 0.6))))

  (testing "pulse is one step less another"
    (doseq [x [-1.0 0.0 0.3 0.5 0.9 2.0]]
      (is (== (s/pulse 0.2 0.8 x) (- (s/step 0.2 x) (s/step 0.8 x))))))

  (testing "a pulse train repeats, and only the first slice of each period"
    (is (== 1.0 (s/pulse-train 0.25 1.0 0.1)))
    (is (== 0.0 (s/pulse-train 0.25 1.0 0.5)))
    (is (== 1.0 (s/pulse-train 0.25 1.0 3.1)))))

(deftest blends-test
  (testing "every ramp is pinned at both ends and clamped outside them"
    (doseq [f [s/boxstep s/smoothstep s/smootherstep]]
      (is (== 0.0 (f 1.0 2.0 0.0)))
      (is (== 0.0 (f 1.0 2.0 1.0)))
      (is (== 1.0 (f 1.0 2.0 2.0)))
      (is (== 1.0 (f 1.0 2.0 9.0)))
      (is (close? 0.5 (f 1.0 2.0 1.5)))))

  (testing "and they are monotone in between"
    (doseq [f [s/boxstep s/smoothstep s/smootherstep]]
      (is (apply < (map #(f 0.0 1.0 (* 0.05 %)) (range 1 20))))))

  (testing "the smooth ones are flat at the ends where the linear one is not"
    ;; That flatness is the whole reason to pay for them: it is what stops
    ;; two adjacent ramps showing a crease where they meet.
    (let [slope (fn [f x] (/ (- (f 0.0 1.0 (+ x 1e-6)) (f 0.0 1.0 x)) 1e-6))]
      (is (close? 1.0 (Math/round ^double (slope s/boxstep 0.0))))
      (is (< (slope s/smoothstep 0.0) 1e-4))
      (is (< (slope s/smootherstep 0.0) 1e-8))))

  (testing "mix and remap agree with the definitions everyone remembers"
    (is (close? 7.0 (s/mix 4.0 10.0 0.5)))
    (is (close? 0.25 (s/remap 0.0 4.0 0.0 1.0 1.0)))
    (is (close? -1.0 (s/remap 0.0 1.0 -1.0 1.0 0.0)))))

(deftest warps-test
  (testing "bias and gain fix the ends of the unit interval"
    (doseq [b [0.2 0.5 0.8]]
      (is (close? 0.0 (s/bias b 0.0)))
      (is (close? 1.0 (s/bias b 1.0)))
      (is (close? 0.0 (s/gain b 0.0)))
      (is (close? 1.0 (s/gain b 1.0)))))

  (testing "bias 0.5 and gain 0.5 are the identity"
    (doseq [t [0.1 0.25 0.5 0.75 0.9]]
      (is (close? t (s/bias 0.5 t)))
      (is (close? t (s/gain 0.5 t)))))

  (testing "bias maps a half to itself"
    (doseq [b [0.2 0.35 0.5 0.7 0.9]]
      (is (close? b (s/bias b 0.5)))))

  (testing "gain also fixes the middle, which bias does not"
    (doseq [g [0.2 0.5 0.8]]
      (is (close? 0.5 (s/gain g 0.5)))))

  (testing "gain above a half adds contrast, below it takes it away"
    (is (< (s/gain 0.8 0.25) 0.25))
    (is (> (s/gain 0.8 0.75) 0.75))
    (is (> (s/gain 0.2 0.25) 0.25))
    (is (< (s/gain 0.2 0.75) 0.75))))

(deftest spline-test
  (testing "Catmull-Rom passes through its interior knots"
    ;; The first and last are phantoms that only set the end tangents, so
    ;; the curve runs from the second knot to the second to last.
    (let [knots [0.0 1.0 4.0 9.0 16.0]]
      (is (close? 1.0 (s/spline 0.0 knots)))
      (is (close? 9.0 (s/spline 1.0 knots)))
      (is (close? 4.0 (s/spline 0.5 knots)))))

  (testing "a constant run of knots is a constant curve"
    (doseq [t [0.0 0.3 0.6 1.0]]
      (is (close? 2.0 (s/spline t [2.0 2.0 2.0 2.0])))))

  (testing "and it is clamped outside [0, 1]"
    (is (close? (s/spline 0.0 [0.0 1.0 2.0 3.0]) (s/spline -5.0 [0.0 1.0 2.0 3.0])))
    (is (close? (s/spline 1.0 [0.0 1.0 2.0 3.0]) (s/spline 5.0 [0.0 1.0 2.0 3.0])))))

(deftest waves-test
  (testing "sawtooth repeats and stays in [0, 1)"
    (doseq [x [-3.7 -0.2 0.0 0.4 5.9]]
      (let [v (s/sawtooth 2.0 x)]
        (is (<= 0.0 v))
        (is (< v 1.0))
        (is (close? v (s/sawtooth 2.0 (+ x 2.0)))))))

  (testing "triangle peaks at the half period and returns"
    (is (close? 0.0 (s/triangle 1.0 0.0)))
    (is (close? 1.0 (s/triangle 1.0 0.5)))
    (is (< (s/triangle 1.0 0.999999) 1e-5))))
