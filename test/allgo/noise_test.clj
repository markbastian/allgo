(ns allgo.noise-test
  (:require [allgo.procedural.noise :as n]
            [clojure.test :refer [deftest is testing]]))

(defn- points
  "A deterministic walk that avoids the lattice, so a test cannot pass by
  accidentally sampling only the points where noise is zero."
  [n]
  (for [i (range n)]
    [(* 0.1379 i) (* 0.2713 i) (* 0.3137 i)]))

(defn- values [f n] (mapv (fn [[x y z]] (f x y z)) (points n)))

(deftest poisson-test
  (testing "Poisson counts average out to the mean asked for"
    (let [mean 3.0
          counts (map #(n/poisson-count mean (/ (double %) 20000)) (range 20000))]
      (is (< (abs (- mean (/ (double (reduce + counts)) 20000))) 0.05)))))

(deftest permutation-test
  (testing "the table is a permutation of 0..255, twice"
    (let [p (n/permutation 12)]
      (is (= 512 (alength p)))
      (is (= (set (range 256)) (set (take 256 (seq p)))))
      (is (= (vec (take 256 (seq p))) (vec (drop 256 (seq p)))))))

  (testing "and a different seed shuffles it differently"
    (is (not= (vec (seq (n/permutation 1))) (vec (seq (n/permutation 2)))))))

(deftest gradient-basis-test
  (let [g (n/gradient-basis {:seed 5})]
    (testing "gradient noise is exactly zero on the lattice"
      ;; The defining property, and what keeps the grid invisible: the
      ;; value comes from the gradients, and at a lattice point the
      ;; offset the gradient is dotted with is the zero vector.
      (doseq [i (range -3 4) j (range -3 4)]
        (is (zero? (g (double i) (double j) 2.0)))))

    (testing "and lands inside [-1, 1] everywhere else"
      (let [vs (values g 20000)]
        (is (> (apply min vs) -1.0))
        (is (< (apply max vs) 1.0))
        (is (< (abs (/ (reduce + vs) (count vs))) 0.05))))

    (testing "it is continuous: a small step gives a small change"
      (doseq [[x y z] (points 200)]
        (is (< (abs (- (g x y z) (g (+ x 1e-5) y z))) 1e-3))))

    (testing "a seed is a world"
      (is (= (values g 500) (values (n/gradient-basis {:seed 5}) 500)))
      (is (not= (values g 500) (values (n/gradient-basis {:seed 6}) 500))))))

(deftest periodic-test
  (testing "a periodic basis repeats exactly, on every axis"
    ;; What makes a field safe to bake into a tile. Not bit-for-bit: a
    ;; coordinate offset by a whole period comes back with a slightly
    ;; different fractional part, so the two agree to floating point and
    ;; not to the last bit. That is the standard a seam has to meet.
    (doseq [make [n/gradient-basis n/value-basis]
            period [4 8 16]]
      (let [b (make {:seed 5 :period period})
            p (double period)
            same? (fn [a b] (< (abs (- (double a) (double b))) 1e-12))]
        (doseq [[x y z] (points 300)]
          (is (same? (b x y z) (b (+ x p) y z)))
          (is (same? (b x y z) (b x (- y (* 2 p)) z)))
          (is (same? (b x y z) (b (+ x (* 3 p)) (+ y p) (- z p))))))))

  (testing "it is still continuous across the seam"
    (let [b (n/gradient-basis {:seed 5 :period 8})]
      (doseq [t (range 0.0 1.0 0.05)]
        (is (< (abs (- (b (- 8.0 1e-6) t 0.3) (b (+ 8.0 1e-6) t 0.3))) 1e-5)))))

  (testing "a period does not flatten the field"
    (let [b (n/gradient-basis {:seed 5 :period 8})
          vs (values b 4000)]
      (is (> (- (apply max vs) (apply min vs)) 1.0))))

  (testing "and asking for none leaves the basis exactly as it was"
    (is (= (values (n/gradient-basis {:seed 5}) 500)
           (values (n/gradient-basis {:seed 5 :period nil}) 500)))))

(deftest value-basis-test
  (let [v (n/value-basis {:seed 5})]
    (testing "value noise is not zero on the lattice -- that is the point of it"
      (is (some #(not (zero? %)) (for [i (range 6)] (v (double i) 0.0 0.0)))))

    (testing "it stays inside the table's range"
      (let [vs (values v 20000)]
        (is (>= (apply min vs) -1.0))
        (is (<= (apply max vs) 1.0))))

    (testing "and repeats the lattice value exactly at the lattice"
      ;; Two visits to the same integer point have to agree, or the
      ;; hash is not a hash.
      (is (== (v 3.0 -7.0 11.0) (v 3.0 -7.0 11.0))))))

(deftest gradient-4d-test
  (let [g4 (n/gradient-basis-4d {:seed 5})]
    (testing "zero on the 4D lattice, for the same reason as in 3D"
      (doseq [i (range -2 3) l (range -2 3)]
        (is (zero? (g4 (double i) 1.0 2.0 (double l))))))

    (testing "in range"
      (let [vs (for [i (range 8000)]
                 (g4 (* 0.137 i) (* 0.271 i) (* 0.313 i) (* 0.077 i)))]
        (is (> (apply min vs) -1.0))
        (is (< (apply max vs) 1.0))))

    (testing "a slice is a 3D basis, and moving the slice moves the world"
      (let [a (n/slice g4 0.25)
            b (n/slice g4 0.25)
            c (n/slice g4 0.75)]
        (is (= (values a 200) (values b 200)))
        (is (not= (values a 200) (values c 200)))))

    (testing "and neighboring slices are near each other, which is what
              makes the fourth dimension usable for animation"
      (let [a (n/slice g4 0.500)
            b (n/slice g4 0.501)]
        (is (every? #(< % 0.02)
                    (map (fn [x y] (abs (- (double x) (double y))))
                         (values a 200) (values b 200))))))))

(deftest sparse-convolution-test
  (let [s (n/sparse-convolution-basis {:seed 5})]
    (testing "it is finite and roughly unit range"
      (let [vs (values s 8000)]
        (is (every? #(Double/isFinite (double %)) vs))
        (is (> (apply min vs) -2.0))
        (is (< (apply max vs) 2.0))
        (is (< (abs (/ (reduce + vs) (count vs))) 0.1))))

    (testing "density changes the texture, not the range"
      ;; The normalization exists so that a density knob does not double
      ;; as a brightness knob.
      (doseq [d [1.0 3.0 8.0]]
        (let [vs (values (n/sparse-convolution-basis {:seed 5 :density d}) 4000)]
          (is (< 0.5 (- (apply max vs) (apply min vs)) 3.0)))))))

(deftest vector-basis-test
  (let [vb (n/vector-basis {:seed 5})]
    (testing "three components, all in range, and none of them each other"
      (let [vs (map (fn [[x y z]] (vb x y z)) (points 2000))]
        (is (every? #(= 3 (count %)) vs))
        (is (every? (fn [v] (every? #(< (abs (double %)) 1.0) v)) vs))
        (is (not= (map first vs) (map second vs)))
        (is (not= (map second vs) (map #(nth % 2) vs)))))))

(deftest transforms-test
  (let [g (n/gradient-basis {:seed 5})]
    (testing "absolute folds the negative half up"
      (let [a (n/absolute g)]
        (doseq [[x y z] (points 500)]
          (is (== (a x y z) (abs (g x y z)))))))

    (testing "ridge is never negative and peaks where the basis crosses zero"
      (let [r (n/ridge g 1.0)]
        (doseq [[x y z] (points 500)]
          (is (>= (r x y z) 0.0)))
        ;; Exactly on the lattice the basis is zero, so the ridge is at
        ;; its offset squared.
        (is (< (abs (- 1.0 (r 2.0 3.0 4.0))) 1e-12))))

    (testing "scaled multiplies domain and range"
      (let [s (n/scaled g 2.0 3.0)]
        (doseq [[x y z] (points 200)]
          (is (< (abs (- (s x y z) (* 3.0 (g (* 2.0 x) (* 2.0 y) (* 2.0 z))))) 1e-12)))))

    (testing "translated shifts the domain"
      (let [t (n/translated g [1.5 -2.5 0.25])]
        (doseq [[x y z] (points 200)]
          (is (== (t x y z) (g (+ x 1.5) (- y 2.5) (+ z 0.25)))))))

    (testing "distortion is the identity when there is none to apply"
      (let [d (n/distorted g (n/vector-basis {:seed 9}) 0.0)]
        (doseq [[x y z] (points 200)]
          (is (== (d x y z) (g x y z))))))

    (testing "and changes things when there is"
      (let [d (n/distorted g (n/vector-basis {:seed 9}) 0.5)]
        (is (not= (values d 500) (values g 500)))))

    (testing "stepping quantizes the range"
      (let [st (n/stepped g 4 0.0)
            vs (set (values st 4000))]
        (is (<= (count vs) 5))))

    (testing "sample calls a basis and agrees with calling it directly"
      (doseq [[x y z] (points 200)]
        (is (== (n/sample g x y z) (g x y z)))))))

(deftest simplex-basis-test
  (let [s (n/simplex-basis {:seed 5})]

    (testing "it stays inside the range a basis promises"
      (let [vs (values s 20000)]
        (is (every? #(<= -1.0 % 1.0) vs))
        ;; And actually uses it -- a scaling constant set too low would
        ;; pass the bound above while quietly making everything gray.
        (is (> (reduce max vs) 0.7))
        (is (< (reduce min vs) -0.7))))

    (testing "and averages out to nothing much"
      (let [vs (values s 20000)]
        (is (< (abs (/ (reduce + vs) (count vs))) 0.02))))

    (testing "it is continuous, which is what says the simplex choice is right"
      ;; The one property a wrong corner selection cannot fake. Picking
      ;; the wrong tetrahedron gives a field that is still bounded, still
      ;; zero-mean and still looks like noise from a distance, but it
      ;; steps at the boundary between one simplex and the next. Walking
      ;; in small increments and bounding the jump catches that, and
      ;; nothing else here would.
      (let [step 1e-5
            worst (reduce max 0.0
                          (for [i (range 20000)]
                            (let [x (* 0.017 i) y (* 0.023 i) z (* 0.029 i)]
                              (abs (- (s x y z) (s (+ x step) y z))))))]
        (is (< worst (* 20.0 step)) (str "jumped by " worst))))

    (testing "it is continuous across the diagonals too"
      ;; Where the offsets tie is exactly where the ranking has to break
      ;; the tie consistently, so walk along x0 = y0 = z0 and across it.
      (let [step 1e-6
            worst (reduce max 0.0
                          (for [i (range 5000)
                                :let [t (* 0.001 i)]
                                d [[step 0.0 0.0] [0.0 step 0.0] [0.0 0.0 step]]]
                            (let [[dx dy dz] d]
                              (abs (- (s t t t) (s (+ t dx) (+ t dy) (+ t dz)))))))]
        (is (< worst (* 20.0 step)) (str "jumped by " worst))))

    (testing "the same seed gives the same field, a different one does not"
      (is (= (values s 500) (values (n/simplex-basis {:seed 5}) 500)))
      (is (not= (values s 500) (values (n/simplex-basis {:seed 6}) 500))))

    (testing "it is not the lattice noise it sits next to"
      (is (not= (values s 500) (values (n/gradient-basis {:seed 5}) 500))))

    (testing "it carries more amplitude than gradient noise, which fractals feel"
      ;; Documented on the var, and worth pinning: a multifractal tuned
      ;; against gradient noise is louder over this.
      (let [sd (fn [vs] (let [m (/ (reduce + vs) (count vs))]
                          (Math/sqrt (/ (reduce + (map #(let [d (- % m)] (* d d)) vs))
                                        (count vs)))))]
        (is (> (sd (values s 20000))
               (* 1.3 (sd (values (n/gradient-basis {:seed 5}) 20000)))))))))

(deftest simplex-basis-4d-test
  (let [s (n/simplex-basis-4d {:seed 7})
        pts (fn [n] (for [i (range n)]
                      [(* 0.0137 i) (* 0.0271 i) (* 0.0313 i) (* 0.0179 i)]))
        vals4 (fn [f n] (mapv (fn [[x y z w]] (f x y z w)) (pts n)))]

    (testing "in range, and using it"
      (let [vs (vals4 s 20000)]
        (is (every? #(<= -1.0 % 1.0) vs))
        (is (> (reduce max vs) 0.7))
        (is (< (reduce min vs) -0.7))
        (is (< (abs (/ (reduce + vs) (count vs))) 0.02))))

    (testing "continuous in every one of the four"
      (let [step 1e-5
            worst (reduce max 0.0
                          (for [i (range 6000)
                                axis (range 4)]
                            (let [p [(* 0.017 i) (* 0.023 i) (* 0.029 i) (* 0.019 i)]
                                  q (update p axis + step)]
                              (abs (- (apply s p) (apply s q))))))]
        (is (< worst (* 20.0 step)) (str "jumped by " worst))))

    (testing "holding the fourth coordinate gives an ordinary 3D basis"
      ;; What the fourth dimension is for: a slice is a field in its own
      ;; right, and neighboring slices are neighboring fields rather
      ;; than unrelated ones.
      (let [a (n/slice s 0.0)
            b (n/slice s 0.0)
            far (n/slice s 40.0)
            near (n/slice s 0.001)
            drift (fn [f g] (/ (reduce + (map #(abs (- %1 %2)) (values f 400) (values g 400)))
                               400.0))]
        (is (= (values a 400) (values b 400)))
        (is (< (drift a near) 0.02) "a slice next door is nearly the same field")
        (is (> (drift a far) 0.1) "a slice far away is a different one")))

    (testing "seeds separate it"
      (is (= (vals4 s 400) (vals4 (n/simplex-basis-4d {:seed 7}) 400)))
      (is (not= (vals4 s 400) (vals4 (n/simplex-basis-4d {:seed 8}) 400))))))
