(ns allgo.robotics.cylindrical-test
  (:require [allgo.numerics.multivariate :as mv]
            [allgo.robotics.cylindrical :as cad]
            [clojure.test :refer [deftest is testing]]))

;; Cylindrical algebraic decomposition, after LaValle, *Planning
;; Algorithms*, Cambridge University Press, 2006, section 6.4.2.

(def ^:private x (mv/variable 2 0))
(def ^:private y (mv/variable 2 1))
(defn- c [k] (mv/constant 2 k))

(defn- circle [cx cy r2]
  (mv/add (mv/power (mv/sub x (c cx)) 2) (mv/power (mv/sub y (c cy)) 2) (c (- r2))))

;; Example 3.1's gingerbread face, (3.11): a head, two eyes and a mouth,
;; the last three negated so that each primitive f <= 0 is the face's
;; side of it. The book gives no sizes; these put the mouth's ends
;; between the eyes' -- the head from -5 to 5 in x, the eyes from -3 to
;; -1 and 1 to 3, the mouth from -2 to 2 -- which is how Figure 6.35
;; draws it.
(def face
  [(circle 0 0 25)
   (mv/scale (circle -2 2 1) -1.0)
   (mv/scale (circle 2 2 1) -1.0)
   (mv/scale (mv/add (mv/power x 2) (mv/scale (mv/power (mv/add y (c 2)) 2) 4.0) (c -4)) -1.0)])

(defn- face? [signs] (every? #(<= % 0) signs))

(defn- signs-at [fs p] (mapv #(let [v (mv/evaluate % p)] (cond (pos? v) 1 (neg? v) -1 :else 0)) fs))

(deftest line-test
  (testing "Example 6.5: x^2 - 2x and x^2 - 4x + 3 cut the line into five 1-cells and four 0-cells, sampled at -1, 1/2, 3/2, 5/2, 4 and 0, 1, 2, 3"
    (let [cells (cad/line-cells [[1.0 -2.0 0.0] [1.0 -4.0 3.0]])]
      (is (= [-1.0 0.0 0.5 1.0 1.5 2.0 2.5 3.0 4.0] (mapv :sample cells)))
      (is (= [1 0 1 0 1 0 1 0 1] (mapv :dim cells)))
      (testing "Figure 6.31's set [1, 2] is the cells where both are at most zero, and x = 1 is the first sample found in it"
        (let [in? (fn [{:keys [signs]}] (every? #(<= % 0) signs))]
          (is (= [1.0 1.5 2.0] (mapv :sample (filter in? cells))))
          (is (= 1.0 (:sample (first (filter in? cells)))))))
      (testing "(6.21) as printed, both at least zero, is (-inf, 0] and [3, inf): true already at -1, so its inequalities are the other way round"
        (let [in? (fn [{:keys [signs]}] (every? #(>= % 0) signs))]
          (is (= [-1.0 0.0 3.0 4.0] (mapv :sample (filter in? cells)))))))))

(deftest gingerbread-test
  (let [d (cad/decompose face)]
    (testing "the critical values: where each curve folds over in y, at its ends in x; no two curves cross"
      (is (= [-5.0 -3.0 -2.0 -1.0 1.0 2.0 3.0 5.0] (mapv :x (:critical d))))
      (is (every? #(= [:fold] (map :kind (:why %))) (:critical d))))
    (testing "Figure 6.35: 37 2-cells, 64 1-cells and 28 0-cells"
      (is (= {2 37 1 64 0 28} (frequencies (map :dim (cad/cells d))))))
    (testing "Figure 6.30: nine sign-invariant regions -- five open, four curves -- so nine sign vectors among the cells"
      (is (= 9 (count (distinct (map :signs (cad/cells d))))))
      (is (= 5 (count (distinct (map :signs (filter #(= 2 (:dim %)) (cad/cells d))))))))
    (testing "sign invariance: a point anywhere in a 2-cell has the signs of its sample"
      (let [rng (java.util.Random. 7)]
        (doseq [_ (range 400)]
          (let [p [(- (* 14 (.nextDouble rng)) 7) (- (* 14 (.nextDouble rng)) 7)]
                id (cad/locate d p)]
            (is (= (:signs (cad/cell d id)) (signs-at face p)) (str p " in " id))))))
    (testing "the decision problem: there is a point in the mouth inside the head, and none in both eyes"
      (is (cad/decide d (fn [[head _ _ mouth]] (and (neg? head) (pos? mouth)))))
      (is (nil? (cad/decide d (fn [[_ left right _]] (and (pos? left) (pos? right)))))))))

(deftest roadmap-test
  (let [d (cad/decompose face)
        around (cad/roadmap d (complement face?))
        within (cad/roadmap d (fn [signs] (every? neg? signs)))]
    (testing "around the face, the free space is outside the head, in each eye and in the mouth"
      (is (= 4 (count (cad/components around)))))
    (testing "within the face, one piece with three holes"
      (is (= 1 (count (cad/components within)))))
    (testing "routes keep to the free space and end where asked"
      (doseq [[rm free? a b] [[around (complement face?) [-6.0 0.0] [6.0 0.5]]
                              [around (complement face?) [0.0 -2.0] [1.5 -2.1]]
                              [within #(every? neg? %) [-2.0 0.2] [2.5 3.8]]
                              [within #(every? neg? %) [-4.5 0.0] [0.0 -3.5]]]]
        (let [path (cad/route rm a b)]
          (is (= a (first path)))
          (is (= b (peek path)))
          (is (every? #(free? (signs-at face %)) path)))))
    (testing "none from outside the head into an eye"
      (is (nil? (cad/route around [-6.0 0.0] [-2.0 2.0]))))))

(deftest lines-test
  (testing "straight edges, vertical ones among them: a square and a triangle as half-planes"
    (let [fs [(mv/sub x (c 1)) (mv/add x (c 1)) (mv/sub y (c 1)) (mv/add y (c 1))
              (mv/sub (mv/add x y) (c 4)) (mv/sub x (c 2)) y]
          obstacle? (fn [[a b cc dd e f g]] (or (and (<= a 0) (>= b 0) (<= cc 0) (>= dd 0))
                                                (and (<= e 0) (>= f 0) (>= g 0))))
          d (cad/decompose fs)
          rm (cad/roadmap d (complement obstacle?))]
      (testing "the vertical edges are critical by their leading coefficients, the slanted one where it crosses the horizontal lines"
        (is (= [[-1.0 :leading] [1.0 :leading] [2.0 :leading] [3.0 :crossing] [4.0 :crossing] [5.0 :crossing]]
               (mapv (fn [{:keys [x why]}] [x (:kind (first why))]) (:critical d)))))
      (testing "the Euler characteristic of the plane, 1, from the counts"
        (let [{v 0 e 1 f 2} (frequencies (map :dim (cad/cells d)))]
          (is (= 1 (+ (- v e) f)))))
      (is (= 1 (count (cad/components rm))))
      (let [path (cad/route rm [1.5 -0.5] [3.5 1.5])]
        (is (= [3.5 1.5] (peek path)))
        (is (not-any? #(obstacle? (signs-at fs %)) path))))))

(deftest square-free-test
  (testing "a repeated factor makes the discriminant vanish identically, and is refused"
    (is (thrown? clojure.lang.ExceptionInfo (cad/decompose [(mv/power (circle 0 0 1) 2)])))))
