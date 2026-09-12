(ns allgo.fractal-test
  (:require [allgo.procedural.fractal :as f]
            [allgo.procedural.noise :as n]
            [clojure.test :refer [deftest is testing]]))

(def ^:private basis (n/gradient-basis {:seed 21}))

(defn- points [n]
  (for [i (range n)] [(* 0.1379 i) (* 0.2713 i) (* 0.3137 i)]))

(defn- values [g n] (mapv (fn [[x y z]] (g x y z)) (points n)))

(deftest spectral-weights-test
  (testing "the weights are lacunarity to the minus H times the octave"
    (let [w (f/spectral-weights {:H 0.75 :lacunarity 2.0 :octaves 6})]
      (dotimes [i 7]
        (is (< (abs (- (aget w i) (Math/pow 2.0 (* -0.75 i)))) 1e-12)))))

  (testing "H = 0 gives every octave the same amplitude"
    (let [w (f/spectral-weights {:H 0.0 :lacunarity 2.0 :octaves 5})]
      (is (every? #(== 1.0 %) (seq w)))))

  (testing "and larger H makes the fall-off steeper"
    (let [gentle (f/spectral-weights {:H 0.3 :lacunarity 2.0 :octaves 5})
          steep (f/spectral-weights {:H 1.0 :lacunarity 2.0 :octaves 5})]
      (is (every? true? (map < (drop 1 (seq steep)) (drop 1 (seq gentle))))))))

(deftest octaves-for-test
  (testing "halving the sample spacing buys exactly one more octave"
    (doseq [s [0.2 0.1 0.05 0.025]]
      (is (< (abs (- 1.0 (- (f/octaves-for (/ s 2.0) 2.0) (f/octaves-for s 2.0)))) 1e-12))))

  (testing "the octave it stops at is the one the samples can still carry"
    ;; Octave i has features about lacunarity^-i across; two samples per
    ;; feature is the limit, so the answer is where those meet.
    (doseq [[spacing lacunarity] [[0.05 2.0] [0.01 2.0] [0.05 3.0] [0.02 1.93]]]
      (let [n (f/octaves-for spacing lacunarity)
            feature (Math/pow lacunarity (- n))]
        (is (< (abs (- feature (* 2.0 spacing))) 1e-9)))))

  (testing "a coarser lacunarity reaches the limit in fewer octaves"
    (is (< (f/octaves-for 0.01 3.0) (f/octaves-for 0.01 2.0))))

  (testing "it never asks for less than one octave, however coarse"
    (doseq [s [0.5 1.0 4.0 1e6]]
      (is (== 1.0 (f/octaves-for s 2.0)))))

  (testing "and an unsampled function is unlimited"
    (is (= ##Inf (f/octaves-for 0.0 2.0)))))

(deftest band-limit-test
  (testing "clamping to the sampling rate leaves the large scale alone"
    ;; The whole claim behind band-limiting: the octaves dropped were not
    ;; contributing shape, only noise. Sample both at the rate the coarse
    ;; one was built for and they should nearly agree.
    (let [spacing 0.02
          n (f/octaves-for spacing 2.0)
          full (f/fbm basis {:octaves 10.0})
          cut (f/fbm basis {:octaves n})
          at (fn [g] (mapv (fn [i] (g (* spacing i) 0.3 0.7)) (range 400)))
          a (at full) b (at cut)
          gap (reduce max (map (fn [p q] (abs (- (double p) (double q)))) a b))
          spread (- (apply max a) (apply min a))]
      (is (< gap (* 0.2 spread)))))

  (testing "and it costs less than the octaves it drops"
    (let [n (f/octaves-for 0.02 2.0)]
      (is (< n 10.0))
      (is (> n 4.0)))))

(deftest fbm-test
  (testing "one octave of fBm is the basis itself"
    (let [one (f/fbm basis {:octaves 1.0})]
      (doseq [[x y z] (points 300)]
        (is (< (abs (- (one x y z) (basis x y z))) 1e-12)))))

  (testing "the second octave arrives at the weight the spectrum says"
    (let [one (f/fbm basis {:octaves 1.0 :H 1.0 :lacunarity 2.0})
          two (f/fbm basis {:octaves 2.0 :H 1.0 :lacunarity 2.0})]
      (doseq [[x y z] (points 200)]
        (is (< (abs (- (- (two x y z) (one x y z))
                       (* 0.5 (basis (* 2.0 x) (* 2.0 y) (* 2.0 z)))))
               1e-12)))))

  (testing "a fractional octave fades the last one in rather than switching it on"
    ;; This is what lets detail appear as a camera approaches without a
    ;; visible pop, so it is worth pinning down.
    (let [at (fn [o] (values (f/fbm basis {:octaves o}) 200))
          d (fn [a b] (reduce max (map (fn [p q] (abs (- (double p) (double q)))) a b)))]
      (is (< (d (at 4.0) (at 4.05)) (d (at 4.0) (at 5.0))))
      (is (= (at 4.0) (at 4.0)))))

  (testing "more octaves add detail without running away"
    (let [range-of (fn [o] (let [vs (values (f/fbm basis {:octaves o}) 3000)]
                             (- (apply max vs) (apply min vs))))]
      (is (< (range-of 2.0) (range-of 8.0)))
      (is (< (range-of 8.0) 4.0)))))

(deftest turbulence-test
  (testing "turbulence is fBm over the folded basis, and says so in code"
    (let [t (f/turbulence basis {:octaves 5.0})
          same (f/fbm (n/absolute basis) {:octaves 5.0})]
      (is (= (values t 300) (values same 300)))))

  (testing "and so is never negative"
    (is (>= (apply min (values (f/turbulence basis {:octaves 5.0}) 3000)) 0.0))))

(deftest multifractals-test
  (testing "every construction is deterministic and finite"
    (doseq [[k g] f/constructions]
      (let [a (g basis {:octaves 6.0})
            b (g basis {:octaves 6.0})]
        (is (= (values a 200) (values b 200)) (str k " is a function"))
        (is (every? #(Double/isFinite (double %)) (values a 2000)) (str k " is finite")))))

  (testing "the ridged multifractal never goes below zero"
    ;; Each octave is a square, weighted by a clamped non-negative weight.
    (is (>= (apply min (values (f/ridged-multifractal basis {:octaves 7.0}) 4000)) 0.0)))

  (testing "hetero-terrain puts its detail where the altitude is"
    ;; The property that makes it terrain rather than noise: measure how
    ;; much the finest octave changed things, and it should be larger in
    ;; the top half of the range than in the bottom.
    (let [coarse (f/hetero-terrain basis {:octaves 3.0})
          fine (f/hetero-terrain basis {:octaves 8.0})
          samples (for [[x y z] (points 6000)]
                    [(coarse x y z) (abs (- (double (fine x y z)) (double (coarse x y z))))])
          median (nth (sort (map first samples)) 3000)
          detail (fn [pred] (let [xs (map second (filter #(pred (first %) median) samples))]
                              (/ (reduce + xs) (count xs))))]
      (is (> (detail >) (detail <)))))

  (testing "the hybrid multifractal does too"
    (let [coarse (f/hybrid-multifractal basis {:octaves 3.0})
          fine (f/hybrid-multifractal basis {:octaves 8.0})
          samples (for [[x y z] (points 6000)]
                    [(coarse x y z) (abs (- (double (fine x y z)) (double (coarse x y z))))])
          median (nth (sort (map first samples)) 3000)
          detail (fn [pred] (let [xs (map second (filter #(pred (first %) median) samples))]
                              (/ (reduce + xs) (count xs))))]
      (is (> (detail >) (detail <)))))

  (testing "fBm, being monofractal, does not -- which is the whole point"
    (let [coarse (f/fbm basis {:octaves 3.0})
          fine (f/fbm basis {:octaves 8.0})
          samples (for [[x y z] (points 6000)]
                    [(coarse x y z) (abs (- (double (fine x y z)) (double (coarse x y z))))])
          median (nth (sort (map first samples)) 3000)
          detail (fn [pred] (let [xs (map second (filter #(pred (first %) median) samples))]
                              (/ (reduce + xs) (count xs))))]
      (is (< (abs (- 1.0 (/ (detail >) (detail <)))) 0.25)))))

(deftest composition-test
  (testing "a fractal is a basis, so fractals nest"
    ;; The claim the whole namespace is arranged around. If this broke,
    ;; every composite terrain in `allgo.procedural.planet` would too.
    (let [inner (f/fbm basis {:octaves 3.0})
          outer (f/ridged-multifractal inner {:octaves 3.0})]
      (is (every? #(Double/isFinite (double %)) (values outer 500)))
      (is (not= (values outer 500) (values inner 500)))))

  (testing "build names a construction without naming its var"
    (doseq [[k g] f/constructions]
      (is (= (values (f/build k basis {:octaves 4.0}) 100)
             (values (g basis {:octaves 4.0}) 100))))))

(deftest measuring-test
  (testing "value-range brackets what the function actually does"
    (let [g (f/hybrid-multifractal basis {:octaves 6.0})
          [lo hi] (f/value-range g 4000 -4.0 4.0)]
      (is (< lo hi))
      (is (= [lo hi] (f/value-range g 4000 -4.0 4.0)))))

  (testing "normalized lands in [0, 1] and keeps the order"
    (let [g (f/hybrid-multifractal basis {:octaves 6.0})
          [lo hi] (f/value-range g 4000)
          nrm (f/normalized g lo hi)]
      (doseq [[x y z] (points 1000)]
        (is (<= 0.0 (nrm x y z) 1.0)))
      (let [a (points 300)]
        (is (= (map (fn [[x y z]] (compare (g x y z) 0.0)) a)
               (map (fn [[x y z]] (compare (g x y z) 0.0)) a))))))

  (testing "the gradient points the way the function climbs"
    (let [g (f/fbm basis {:octaves 4.0})
          e 1e-4]
      (doseq [p (take 100 (points 200))]
        (let [grad (f/gradient g p e)
              [x y z] p
              step 1e-3
              uphill (mapv (fn [c d] (+ (double c) (* step (double d)))) p grad)
              [ux uy uz] uphill]
          (when (> (reduce + (map #(* % %) grad)) 1e-6)
            (is (>= (g ux uy uz) (g x y z)))))))))
