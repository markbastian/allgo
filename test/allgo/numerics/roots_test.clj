(ns allgo.numerics.roots-test
  "Bracketing root-finders and event location.

  Roots of equations and of polynomials (Chapra and Canale, *Numerical
  Methods for Engineers*, 7th ed., chapters 5-7): the book's worked
  examples, and exact roots."
  (:require [allgo.numerics.polynomial :as poly]
            [allgo.numerics.roots :as roots]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(deftest bisection
  (is (< (abs (- (roots/bisect #(- (* % %) 2.0) 0.0 2.0) (math/sqrt 2.0))) 1e-15))
  (testing "either way round"
    (is (< (abs (- (roots/bisect math/cos 3.0 0.0) (/ math/PI 2))) 1e-15))))

(deftest events
  (testing "each change of a condition, found to the tolerance"
    (let [ts (roots/transitions #(pos? (math/sin %)) 0.5 10.0 0.3 1e-9)]
      ;; each reported with the value it had before: true leaving at pi,
      ;; false entering at 2 pi, true leaving at 3 pi
      (is (= [true false true] (map second ts)))
      (is (every? true? (map #(< (abs (- (first %1) %2)) 1e-9) ts [math/PI (* 2 math/PI) (* 3 math/PI)]))))))

(defn- close? [a b tol] (< (abs (- a b)) tol))

(def ^:private f (fn [x] (- (math/exp (- x)) x)))
(def ^:private df (fn [x] (- (- (math/exp (- x))) 1.0)))
(def ^:private root 0.56714329040978387)

(defn- parachutist
  "Example 5.1: the drag coefficient c at which a 68.1 kg parachutist
  reaches 40 m/s after 10 s."
  [c]
  (- (* (/ (* 9.81 68.1) c) (- 1.0 (math/exp (- (* (/ c 68.1) 10.0))))) 40.0))

(deftest bracketing
  (testing "Example 5.5: false position from 12 and 16, the parachutist's drag coefficient"
    (is (close? (roots/false-position parachutist 12.0 16.0) 14.8011 1e-4))
    (is (close? (roots/false-position parachutist 12.0 16.0 {:modified? false}) 14.8011 1e-4)))
  (testing "the first chord at 14.9309, as the example computes it (its prose says 14.9113, a slip)"
    (is (close? (- 16.0 (/ (* (parachutist 16.0) (- 12.0 16.0)) (- (parachutist 12.0) (parachutist 16.0)))) 14.9309 1e-4)))
  (testing "incremental search finds every sign change of sin 10x + cos 3x on [3, 6]"
    (let [g (fn [x] (+ (math/sin (* 10.0 x)) (math/cos (* 3.0 x))))
          brackets (roots/incremental-search g 3.0 6.0 300)]
      (is (= 9 (count brackets)))
      (is (every? (fn [[a b]] (neg? (* (g a) (g b)))) brackets))))
  (testing "a slow case for plain false position, x^10 - 1 on [0, 1.3], which the modification fixes"
    (let [g #(- (math/pow % 10.0) 1.0)]
      (is (close? (roots/false-position g 0.0 1.3) 1.0 1e-12)))))

(deftest open-methods
  (testing "Example 6.3: Newton-Raphson for e^-x - x from 0"
    (is (close? (roots/newton f df 0.0) root 1e-15)))
  (testing "Example 6.6: the secant method from 0 and 1"
    (is (close? (roots/secant f 0.0 1.0) root 1e-15)))
  (testing "Example 6.8: the modified secant, delta 0.01, from 1"
    (is (close? (roots/modified-secant f 1.0 0.01) root 1e-15)))
  (testing "Brent's method from the bracket [0, 1], and on the parachutist"
    (is (close? (roots/brent f 0.0 1.0) root 1e-15))
    (is (close? (roots/brent parachutist 12.0 16.0) (roots/newton parachutist
                                                                  #(/ (- (parachutist (+ % 1e-6)) (parachutist (- % 1e-6))) 2e-6)
                                                                  14.0)
                1e-9)))
  (testing "and nil without a sign change"
    (is (nil? (roots/brent f 1.0 2.0)))))

(deftest multiple-roots
  (testing "Section 6.5: the double root of (x - 3)(x - 1)^2, from 0 -- quadratically"
    (let [g (fn [x] (* (- x 3.0) (- x 1.0) (- x 1.0)))
          dg (fn [x] (+ (* (- x 1.0) (- x 1.0)) (* 2.0 (- x 3.0) (- x 1.0))))
          d2g (fn [x] (+ (* 4.0 (- x 1.0)) (* 2.0 (- x 3.0))))]
      (is (close? (roots/newton-multiple g dg d2g 0.0) 1.0 1e-12)))))

(deftest systems
  (let [F (fn [[x y]] [(+ (* x x) (* x y) -10.0) (+ y (* 3.0 x y y) -57.0)])]
    (testing "Example 6.12: Newton-Raphson for x^2 + xy = 10, y + 3xy^2 = 57 from (1.5, 3.5)"
      (is (every? true? (map #(close? %1 %2 1e-12) (roots/newton-system F [1.5 3.5]) [2.0 3.0])))
      (testing "with the Jacobian given"
        (is (every? true? (map #(close? %1 %2 1e-12)
                               (roots/newton-system F (fn [[x y]] [[(+ (* 2.0 x) y) x] [(* 3.0 y y) (+ 1.0 (* 6.0 x y))]]) [1.5 3.5])
                               [2.0 3.0])))))
    (testing "Example 6.11: fixed-point iteration, the convergent rearrangement"
      (let [g (fn [[x y]] [(math/sqrt (- 10.0 (* x y))) (math/sqrt (/ (- 57.0 y) (* 3.0 x)))])]
        (is (every? true? (map #(close? %1 %2 1e-10) (roots/fixed-point-system g [1.5 3.5]) [2.0 3.0])))))
    (testing "and the divergent one gives nil"
      (let [g (fn [[x y]] [(/ (- 10.0 (* x x)) y) (- 57.0 (* 3.0 x y y))])]
        (is (nil? (roots/fixed-point-system g [1.5 3.5])))))))

(defn- agree?
  "Whether `xs` agree with the book's `table` to within `tol`, one for one."
  [table xs tol]
  (and (= (count table) (count xs)) (every? true? (map #(close? %1 %2 tol) table xs))))

(deftest iterations
  (testing "Example 5.4: bisection's estimates for the parachutist from 12 and 16"
    (is (= [14.0 15.0 14.5 14.75 14.875 14.8125] (map :x (take 6 (roots/bisect-seq parachutist 12.0 16.0))))))
  (let [g #(- (math/pow % 10.0) 1.0)]
    (testing "Example 5.6: on x^10 - 1 from 0 and 1.3, bisection is quick"
      (is (agree? [0.65 0.975 1.1375 1.05625 1.015625] (map :x (take 5 (roots/bisect-seq g 0.0 1.3))) 1e-15)))
    (testing "and plain false position crawls in from the low end"
      (is (agree? [0.09430 0.18176 0.26287 0.33811 0.40788]
                  (map :x (take 5 (roots/false-position-seq g 0.0 1.3 {:modified? false}))) 5e-6))
      (is (every? #(= 1.3 (:hi %)) (take 5 (roots/false-position-seq g 0.0 1.3 {:modified? false}))))))
  (testing "Example 5.5: false position's first two estimates for the parachutist"
    (is (agree? [14.9309 14.8151] (map :x (take 2 (roots/false-position-seq parachutist 12.0 16.0))) 5e-5)))
  (testing "Example 6.7: false position on ln x from 0.5 and 5"
    (is (agree? [1.8546 1.2163 1.0585] (map :x (take 3 (roots/false-position-seq math/log 0.5 5.0 {:modified? false}))) 5e-5)))
  (testing "Example 6.1: fixed-point iteration x <- e^-x from 0, converging linearly"
    (is (agree? [0.0 1.0 0.367879 0.692201 0.500473 0.606244 0.545396 0.579612 0.560115 0.571143 0.564879]
                (map :x (take 11 (roots/fixed-point-seq #(math/exp (- %)) 0.0))) 1e-6)))
  (testing "Example 6.3: Newton-Raphson on e^-x - x from 0, converging quadratically"
    (is (agree? [0.0 0.5 0.566311003 0.567143165 0.567143290] (map :x (take 5 (roots/newton-seq f df 0.0))) 5e-10)))
  (testing "Example 6.5: Newton-Raphson on x^10 - 1 from 0.5 overshoots and creeps back"
    (is (agree? [0.5 51.65 46.485 41.8365 37.65285 33.887565]
                (map :x (take 6 (roots/newton-seq #(- (math/pow % 10.0) 1.0) #(* 10.0 (math/pow % 9.0)) 0.5))) 5e-7)))
  (testing "Example 6.6: the secant method from 0 and 1"
    (is (agree? [1.0 0.61270 0.56384 0.56717] (map :x (take 4 (roots/secant-seq f 0.0 1.0))) 5e-6)))
  (testing "Example 6.8: the modified secant, delta 0.01, from 1"
    (is (agree? [1.0 0.537263 0.56701 0.567143] (map :x (take 4 (roots/modified-secant-seq f 1.0 0.01))) 5e-6)))
  (let [g (fn [x] (* (- x 3.0) (- x 1.0) (- x 1.0)))
        dg (fn [x] (+ (* (- x 1.0) (- x 1.0)) (* 2.0 (- x 3.0) (- x 1.0))))
        d2g (fn [x] (+ (* 4.0 (- x 1.0)) (* 2.0 (- x 3.0))))]
    (testing "Example 6.10: onto the double root of (x - 3)(x - 1)^2 from 0, Newton-Raphson linearly"
      ;; the book prints the fourth estimate 0.9133290; it is 0.91332989
      (is (agree? [0.0 0.4285714 0.6857143 0.8328654 0.9133299 0.9557833 0.9776551]
                  (map :x (take 7 (roots/newton-seq g dg 0.0))) 5e-8)))
    (testing "and the modified method quadratically"
      (is (agree? [0.0 1.105263 1.003082 1.000002] (map :x (take 4 (roots/newton-multiple-seq g dg d2g 0.0))) 5e-7)))
    (testing "and onto the simple root at 3 from 4, both"
      (is (agree? [4.0 3.4 3.1 3.008696 3.000075 3.0] (map :x (take 6 (roots/newton-seq g dg 4.0))) 5e-7))
      (is (agree? [4.0 2.636364 2.820225 2.961728 2.998479 2.999998]
                  (map :x (take 6 (roots/newton-multiple-seq g dg d2g 4.0))) 5e-7))))
  (testing "Example 6.12: the first Newton-Raphson step for x^2 + xy = 10, y + 3xy^2 = 57 from (1.5, 3.5)"
    (let [F (fn [[x y]] [(+ (* x x) (* x y) -10.0) (+ y (* 3.0 x y y) -57.0)])
          J (fn [[x y]] [[(+ (* 2.0 x) y) x] [(* 3.0 y y) (+ 1.0 (* 6.0 x y))]])]
      (is (agree? [2.03603 2.84388] (:x (second (roots/newton-system-seq F J [1.5 3.5]))) 5e-6)))))

(deftest watching-and-composing
  (testing "Brent's method on the parachutist: every bracket holds the root, and the steps say how they were taken"
    (let [states (take-while #(> (abs (- (:c %) (:b %))) 1e-12) (roots/brent-seq parachutist 12.0 16.0 1e-15))]
      (is (every? #(<= (* (:fb %) (:fc %)) 0.0) states))
      (is (every? #(<= (abs (:fb %)) (abs (:fc %))) states))
      (is (= :start (:by (first states))))
      (is (some #{:interpolation :secant} (map :by states)))))
  (testing "a hybrid by composition: bisection until the bracket is narrow, then Newton-Raphson from its midpoint"
    (let [x0 (:x (first (drop-while #(> (- (:hi %) (:lo %)) 0.1) (roots/bisect-seq f 0.0 1.0))))]
      (is (close? (:x (nth (roots/newton-seq f df x0) 3)) root 1e-15) "three Newton steps from a bracket of 0.1 suffice")))
  (testing "a step that cannot go on ends the seq: Newton-Raphson at a zero derivative"
    (is (= [0.0] (map :x (roots/newton-seq #(+ 1.0 (* % %)) #(* 2.0 %) 0.0))))
    (is (nil? (roots/newton #(+ 1.0 (* % %)) #(* 2.0 %) 0.0)))))

(deftest polynomials
  (testing "Example 7.1: (x^2 + 2x - 24)/(x - 4) = x + 6"
    (is (= {:quotient [1.0 6.0] :remainder 0.0} (poly/deflate [1.0 2.0 -24.0] 4.0))))
  (testing "long division: x^3 + 4 = (x - 1)(x^2 + x - 1) + 2x + 3"
    (is (= {:quotient [1.0 -1.0] :remainder [2.0 3.0]}
           (update (poly/long-divide [1.0 0.0 0.0 4.0] [1.0 1.0 -1.0]) :remainder #(mapv double %)))))
  (testing "value and derivatives by synthetic division"
    ;; x^3 - 2x^2 + 3x - 4 at 2: 2, 7, 8, 6
    (is (= [2.0 7.0 8.0 6.0] (poly/derivatives [1.0 -2.0 3.0 -4.0] 2.0 3))))
  (testing "Example 7.2: Muller's method on x^3 - 13x - 12 from 4.5, 5.5, 5 reaches 4"
    (let [[re im] (poly/muller [1.0 0.0 -13.0 -12.0] 4.5 5.5 5.0)]
      (is (close? re 4.0 1e-12)) (is (zero? im))))
  (testing "and finds a complex root of x^2 + 1 from real guesses"
    (let [[re im] (poly/muller [1.0 0.0 1.0] 0.5 1.0 1.5)]
      (is (close? re 0.0 1e-12)) (is (close? (abs im) 1.0 1e-12))))
  (testing "Example 7.3: Bairstow's method on x^5 - 3.5x^4 + 2.75x^3 + 2.125x^2 - 3.875x + 1.25"
    (let [rs (sort (poly/bairstow [1.0 -3.5 2.75 2.125 -3.875 1.25]))]
      (is (= 5 (count rs)))
      (is (every? true? (map (fn [[a b] [c d]] (and (close? a c 1e-10) (close? b d 1e-10)))
                             rs [[-1.0 0.0] [0.5 0.0] [1.0 -0.5] [1.0 0.5] [2.0 0.0]])))))
  (testing "Bairstow on a tenth-degree polynomial with known roots"
    (let [want [-3.0 -2.0 -1.0 0.5 1.5 2.5 4.0]
          ;; (x^2 + 1)(x^2 + 2x + 5) times the real factors
          mul (fn [a b] (reduce (fn [acc [i x]] (reduce (fn [acc [j y]] (update acc (+ i j) + (* x y))) acc (map-indexed vector b)))
                                (vec (repeat (+ (count a) (count b) -1) 0.0)) (map-indexed vector a)))
          p (reduce mul [1.0 0.0 1.0] (concat [[1.0 2.0 5.0]] (map (fn [r] [1.0 (- r)]) want)))
          got (poly/bairstow p)]
      (is (= 11 (count got)))
      (is (every? true? (map #(close? %1 %2 1e-8) (poly/real-roots got) want))))))
