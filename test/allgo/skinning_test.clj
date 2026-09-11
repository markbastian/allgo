(ns allgo.skinning-test
  (:require [allgo.geometry.tet-mesh :as tm]
            [allgo.physics.skinning :as skin]
            [clojure.test :refer [deftest is testing]]))

(defn- ball-points
  "`n` points on a sphere, as a flat vector -- a stand-in for a detailed
  visual mesh."
  [n [cx cy cz] r seed]
  (let [rng (java.util.Random. seed)]
    (vec (for [_ (range n)
               :let [th (* 2 Math/PI (.nextDouble rng))
                     ph (Math/acos (- (* 2 (.nextDouble rng)) 1))]
               c [(+ cx (* r (Math/sin ph) (Math/cos th)))
                  (+ cy (* r (Math/cos ph)))
                  (+ cz (* r (Math/sin ph) (Math/sin th)))]]
           (double c)))))

(defn- max-error [expected actual]
  (reduce max 0.0 (map (fn [a b] (abs (- (double a) (double b)))) expected (seq actual))))

(defn- moved-cage
  "The cage's vertices put through `f`, as a flat array."
  [cage f]
  (double-array (mapcat f (partition 3 (:verts cage)))))

;; ---------------------------------------------------------------------------

(deftest barycentric-test
  (let [v [0.0 0.0 0.0  1.0 0.0 0.0  0.0 1.0 0.0  0.0 0.0 1.0]]
    (testing "a corner gets all the weight"
      (is (= [1.0 0.0 0.0 0.0] (mapv #(+ 0.0 %) (tm/barycentric v 0 1 2 3 [0.0 0.0 0.0]))))
      (is (= [0.0 1.0 0.0 0.0] (mapv #(+ 0.0 %) (tm/barycentric v 0 1 2 3 [1.0 0.0 0.0])))))

    (testing "the weights always sum to one, inside or out"
      (doseq [p [[0.25 0.25 0.25] [0.1 0.1 0.1] [5.0 -3.0 2.0] [-1.0 -1.0 -1.0]]]
        (is (< (abs (- 1.0 (reduce + (tm/barycentric v 0 1 2 3 p)))) 1e-12) (str p))))

    (testing "they reconstruct the point they came from"
      (doseq [p [[0.25 0.25 0.25] [0.9 0.05 0.02] [3.0 1.0 -2.0]]]
        (let [[b0 b1 b2 b3] (tm/barycentric v 0 1 2 3 p)
              back (mapv (fn [axis]
                           (+ (* b0 (nth (tm/vertex v 0) axis))
                              (* b1 (nth (tm/vertex v 1) axis))
                              (* b2 (nth (tm/vertex v 2) axis))
                              (* b3 (nth (tm/vertex v 3) axis))))
                         [0 1 2])]
          (is (< (max-error p back) 1e-12) (str p)))))

    (testing "inside? agrees with the sign of the weights"
      (is (tm/inside? v 0 1 2 3 [0.25 0.25 0.25]))
      (is (tm/inside? v 0 1 2 3 [0.0 0.0 0.0]) "a corner counts as inside")
      (is (not (tm/inside? v 0 1 2 3 [0.5 0.5 0.5])) "just past the slanted face")
      (is (not (tm/inside? v 0 1 2 3 [-0.01 0.1 0.1]))))

    (testing "a flat tetrahedron has no inverse and says so"
      (is (nil? (tm/barycentric [0.0 0.0 0.0 1.0 0.0 0.0 2.0 0.0 0.0 3.0 0.0 0.0]
                                0 1 2 3 [0.5 0.0 0.0]))))))

(deftest binding-test
  (let [cage (tm/lattice-box 2 2 2 0.5)
        vis  (ball-points 1500 [0.0 0.5 0.0] 0.35 42)
        sk   (skin/bind cage vis)]
    (testing "a surface well inside the cage binds entirely, and inside"
      (is (= 1500 (:n sk)))
      (is (= 1500 (:bound sk)) "every vertex found a tetrahedron containing it")
      (is (zero? (:unbound sk))))

    (testing "the rest pose is reproduced exactly"
      (is (< (max-error vis (skin/skinned-positions sk (:tet-ids cage)
                                                    (double-array (:verts cage))))
             1e-12)))

    (testing "a surface poking outside the cage still binds, by extrapolation"
      ;; Larger than the cage, so much of it is outside every tetrahedron.
      (let [big (ball-points 800 [0.0 0.5 0.0] 0.9 7)
            sk  (skin/bind cage big)]
        (is (zero? (:unbound sk)) "nothing is left without a tetrahedron")
        (is (< (:bound sk) 800) "and much of it is genuinely outside")
        (is (< (max-error big (skin/skinned-positions sk (:tet-ids cage)
                                                      (double-array (:verts cage))))
               1e-12)
            "extrapolation reproduces the rest pose too")))))

(deftest affine-invariance-test
  ;; Barycentric weights are affine invariant, so any affine motion of the
  ;; cage must carry the surface exactly -- not approximately. This is what
  ;; makes the technique correct rather than merely plausible.
  (let [cage (tm/lattice-box 2 2 2 0.5)
        vis  (ball-points 900 [0.0 0.5 0.0] 0.3 99)
        sk   (skin/bind cage vis)
        skin-with (fn [f] (skin/skinned-positions sk (:tet-ids cage) (moved-cage cage f)))]

    (testing "translation"
      (let [f (fn [[x y z]] [(+ x 10.0) (- y 3.0) (+ z 0.5)])]
        (is (< (max-error (mapcat f (partition 3 vis)) (skin-with f)) 1e-12))))

    (testing "uniform scale about the origin"
      (let [f (fn [[x y z]] [(* 3.0 x) (* 3.0 y) (* 3.0 z)])]
        (is (< (max-error (mapcat f (partition 3 vis)) (skin-with f)) 1e-12))))

    (testing "rotation"
      (let [a (/ Math/PI 5)
            f (fn [[x y z]] [(- (* x (Math/cos a)) (* z (Math/sin a)))
                             y
                             (+ (* x (Math/sin a)) (* z (Math/cos a)))])]
        (is (< (max-error (mapcat f (partition 3 vis)) (skin-with f)) 1e-12))))

    (testing "shear, which is affine but not rigid"
      (let [f (fn [[x y z]] [(+ x (* 0.4 y)) y (+ z (* 0.2 y))])]
        (is (< (max-error (mapcat f (partition 3 vis)) (skin-with f)) 1e-12))))

    (testing "a non-affine squash is followed, but only piecewise"
      ;; Bending is not affine, so the surface cannot match exactly -- it is
      ;; interpolated per tetrahedron. It must still move, and stay bounded.
      (let [f (fn [[x y z]] [x (* y (+ 1.0 (* 0.3 (Math/sin (* 3 x))))) z])
            out (skin-with f)]
        (is (pos? (max-error vis out)) "it did move")
        (is (every? #(and (not (Double/isNaN %)) (< (abs %) 10.0)) (seq out)))))))

(deftest cost-test
  (testing "the cage is far smaller than the surface it carries"
    ;; The point of the technique, as a number.
    (let [cage (tm/lattice-box 3 3 3 0.3)
          vis  (ball-points 20000 [0.0 0.45 0.0] 0.3 5)
          sk   (skin/bind cage vis)
          tets (/ (count (:tet-ids cage)) 4)]
      (is (= 162 tets))
      (is (= 20000 (:n sk)))
      (is (> (/ 20000.0 tets) 100.0) "over a hundred visual vertices per element")
      (is (zero? (:unbound sk))))))
