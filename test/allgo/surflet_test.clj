(ns allgo.surflet-test
  (:require [allgo.procedural.noise :as noise]
            [allgo.procedural.surflet :as surflet]
            [clojure.test :refer [deftest is testing]]))

(def ^:private points
  (for [i (range 2000)] [(* 0.0731 i) (- (* 0.0419 i) 30.0) (+ 0.3 (* 0.0217 i))]))

(deftest hump-test
  (testing "1 at the point, 0 from a cell away, and flat at both ends"
    (is (= 1.0 (surflet/hump 0.0)))
    (is (= 0.0 (surflet/hump 1.0) (surflet/hump -1.5)))
    (is (< (abs (- (surflet/hump 1e-4) 1.0)) 1e-9))
    (is (< (surflet/hump 0.9999) 1e-9)))
  (testing "symmetric, and the halves sum to one: the fade's own symmetry"
    (doseq [t [0.1 0.25 0.5 0.8]]
      (is (= (surflet/hump t) (surflet/hump (- t))))
      (is (< (abs (- 1.0 (+ (surflet/hump t) (surflet/hump (- 1.0 t))))) 1e-12)))))

(deftest identity-test
  (testing "a sum of separable surflets is gradient noise"
    (doseq [seed [0 4 17]]
      (let [g (noise/gradient-basis {:seed seed})
            s (surflet/basis {:seed seed})]
        (is (every? (fn [[x y z]] (< (abs (- (g x y z) (s x y z))) 1e-12)) points)))))
  (testing "and the surflets that reach a point add up to it"
    (let [opts {:seed 4}
          s (surflet/basis opts)]
      (doseq [p (take 50 points)]
        (is (< (abs (- (apply s p) (reduce + (map second (surflet/surflets-at opts p))))) 1e-12))))))

(deftest locality-test
  (let [base {:seed 9}
        edited (assoc base :overrides {[3 2 0] [0.0 -1.4 0.0]})
        a (surflet/basis base)
        b (surflet/basis edited)]
    (testing "a changed gradient changes the noise within a cell of its point and nowhere else"
      (doseq [x (range 0.0 6.0 0.13) y (range -1.0 5.0 0.17) z [-0.6 -0.2 0.3 0.7 1.2]]
        (let [near? (and (< (abs (- x 3)) 1.0) (< (abs (- y 2)) 1.0) (< (abs z) 1.0))
              d (abs (- (a x y z) (b x y z)))]
          (if near?
            (is (< d 1.5))
            (is (zero? d))))))
    (testing "and does change it there"
      (is (> (abs (- (a 3.3 2.2 0.1) (b 3.3 2.2 0.1))) 0.05)))
    (testing "a radial surflet reaches as far as its radius, round"
      (let [r (assoc base :kernel :radial :radius 1.5)
            ra (surflet/basis r)
            rb (surflet/basis (assoc r :overrides {[3 2 0] [0.0 -1.4 0.0]}))]
        (is (zero? (- (ra 4.2 3.2 0.0) (rb 4.2 3.2 0.0))))     ; 1.70 away
        (is (not (zero? (- (ra 4.2 2.0 0.0) (rb 4.2 2.0 0.0)))))))))  ; 1.2 away

(deftest radial-test
  (let [s (surflet/basis {:seed 2})
        r (surflet/basis {:seed 2 :kernel :radial})
        rms (fn [f] (Math/sqrt (/ (reduce + (map (fn [p] (let [v (apply f p)] (* v v))) points)) (count points))))]
    (testing "the radial kind is scaled to spread as far as the separable one"
      (is (< 0.8 (/ (rms r) (rms s)) 1.25)))
    (testing "and reaches more lattice points to get there"
      (is (= 1 (surflet/reach {})))
      (is (= 2 (surflet/reach {:kernel :radial :radius 1.5})))
      (is (> (count (surflet/surflets-at {:seed 2 :kernel :radial} [0.31 0.62 0.47]))
             (count (surflet/surflets-at {:seed 2} [0.31 0.62 0.47])))))
    (testing "in two dimensions only the plane's own lattice points count"
      (let [flat (surflet/basis {:seed 2 :kernel :radial :dimensions 2})]
        (is (every? (fn [[_ ijk]] (zero? (nth ijk 2)))
                    (map (fn [[ijk v]] [v ijk]) (surflet/surflets-at {:seed 2 :kernel :radial :dimensions 2} [0.4 0.7 0.0]))))
        (is (< (abs (- (flat 0.4 0.7 0.0)
                       (reduce + (map second (surflet/surflets-at {:seed 2 :kernel :radial :dimensions 2} [0.4 0.7 0.0])))))
               1e-12))))))
