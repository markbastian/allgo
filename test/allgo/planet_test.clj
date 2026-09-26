(ns allgo.planet-test
  (:require [allgo.geometry.vec3 :as v3]
            [allgo.procedural.fractal :as fractal]
            [allgo.procedural.noise :as noise]
            [allgo.procedural.planet :as p]
            [clojure.test :refer [deftest is testing]]))

(def ^:private world (delay (p/planet {:seed 3 :samples 2048})))

(defn- dirs [n] (p/fibonacci-directions n))

(deftest geography-test
  (testing "latitude and longitude round-trip through a direction"
    (doseq [lat [-1.4 -0.7 0.0 0.3 1.2] lon [-3.0 -1.0 0.0 2.0 3.1]]
      (let [d (p/direction lat lon)]
        (is (< (abs (- 1.0 (v3/length d))) 1e-12))
        (is (< (abs (- lat (p/latitude d))) 1e-12))
        (is (< (abs (- lon (p/longitude d))) 1e-12)))))

  (testing "the poles are where they should be"
    (is (< (abs (- (* 0.5 Math/PI) (p/latitude [0.0 1.0 0.0]))) 1e-12))
    (is (< (abs (+ (* 0.5 Math/PI) (p/latitude [0.0 -1.0 0.0]))) 1e-12)))

  (testing "Fibonacci directions are unit vectors spread over the whole sphere"
    (let [ds (dirs 2000)]
      (is (every? #(< (abs (- 1.0 (v3/length %))) 1e-9) ds))
      ;; No hemisphere is favored: the mean of an even spread is the center.
      (is (< (v3/length (v3/scale (reduce v3/add ds) (/ 1.0 2000))) 0.01))))

  (testing "the tangent basis is orthonormal everywhere, poles included"
    (doseq [d (conj (dirs 500) [0.0 1.0 0.0] [0.0 -1.0 0.0] [1.0 0.0 0.0])]
      (let [[u v] (p/tangent-basis d)]
        (is (< (abs (- 1.0 (v3/length u))) 1e-9))
        (is (< (abs (- 1.0 (v3/length v))) 1e-9))
        (is (< (abs (v3/dot u v)) 1e-9))
        (is (< (abs (v3/dot u d)) 1e-9))
        (is (< (abs (v3/dot v d)) 1e-9))))))

(deftest terrain-range-test
  (testing "the measured range really does bracket the terrain"
    (let [t (p/default-terrain {:seed 8})
          [lo hi] (p/terrain-range t 4096)]
      (is (< lo hi))
      (doseq [d (dirs 2000)]
        (let [[x y z] d]
          (is (<= (- lo 1e-9) (noise/sample t x y z) (+ hi 1e-9)))))))

  (testing "and it is the same measurement every time"
    (let [t (p/default-terrain {:seed 8})]
      (is (= (p/terrain-range t 1024) (p/terrain-range t 1024))))))

(deftest elevation-test
  (testing "elevation stays inside the relief it was given"
    (let [pl @world]
      (doseq [d (dirs 3000)]
        ;; Not an exact bound: the scaling comes from a finite sample of
        ;; directions, so one that was not sampled can overshoot a little.
        (is (<= (- (* 1.1 (:relief pl)))
                (p/elevation pl d)
                (* 1.1 (:relief pl)))))))

  (testing "the extremes are actually reached, so relief means what it says"
    (let [pl @world
          hs (map #(p/elevation pl %) (dirs 4096))]
      (is (> (apply max hs) (* 0.9 (:relief pl))))
      (is (< (apply min hs) (* -0.9 (:relief pl))))))

  (testing "a seed is a planet"
    (let [a (p/planet {:seed 1 :samples 512})
          b (p/planet {:seed 1 :samples 512})
          c (p/planet {:seed 2 :samples 512})
          hs (fn [pl] (mapv #(p/elevation pl %) (dirs 300)))]
      (is (= (hs a) (hs b)))
      (is (not= (hs a) (hs c))))))

(deftest sea-test
  (testing "raising the sea level drowns more of the planet, never less"
    (let [terrain (p/default-terrain {:seed 5})
          fraction (fn [level]
                     (let [pl (p/planet {:terrain terrain :sea-level level :samples 2048})]
                       (/ (count (filter #(p/ocean? pl %) (dirs 1500))) 1500.0)))
          levels [-0.03 -0.015 0.0 0.015 0.03]
          fs (map fraction levels)]
      (is (= fs (sort fs)))
      (is (< (first fs) 0.2))
      (is (> (last fs) 0.8))))

  (testing "the sea is flat, and the ground under it is not"
    (let [pl @world
          wet (filter #(p/ocean? pl %) (dirs 2000))]
      (is (seq wet))
      (is (every? #(== (+ (:radius pl) (:sea-level pl)) (p/surface-radius pl %)) wet))
      (is (> (count (set (map #(p/elevation pl %) wet))) 1))))

  (testing "depth is positive under water and zero on land"
    (let [pl @world]
      (doseq [d (dirs 1000)]
        (if (p/ocean? pl d)
          (is (pos? (p/depth pl d)))
          (is (zero? (p/depth pl d))))))))

(deftest surface-test
  (testing "the surface is never below the sea and never outside the relief"
    (let [{:keys [radius relief sea-level] :as pl} @world]
      (doseq [d (dirs 2000)]
        (let [r (p/surface-radius pl d)]
          (is (>= r (+ radius sea-level)))
          (is (<= r (+ radius relief 1e-9)))))))

  (testing "altitude is zero on the surface, positive above, negative below"
    (let [pl @world]
      (doseq [d (take 200 (dirs 1000))]
        (let [s (p/surface-point pl d)]
          (is (< (abs (p/altitude pl s)) 1e-9))
          (is (pos? (p/altitude pl (v3/scale s 1.05))))
          (is (neg? (p/altitude pl (v3/scale s 0.95))))))))

  (testing "normals are unit, outward, and radial over flat water"
    (let [pl @world]
      (doseq [d (take 400 (dirs 1500))]
        (let [n (p/surface-normal pl d)]
          (is (< (abs (- 1.0 (v3/length n))) 1e-6))
          (is (pos? (v3/dot n d)))
          ;; Well out to sea the surface is the flat sea, so the normal is
          ;; the radial. Not at the shore: the finite differences there
          ;; straddle the waterline and pick up the beach, which is right.
          (when (> (p/depth pl d) (* 0.05 (:relief pl)))
            (is (< (v3/distance n d) 1e-6)))))))

  (testing "slope is zero on the sea and positive on broken ground"
    (let [pl @world
          land (remove #(p/ocean? pl %) (dirs 1500))]
      (is (every? #(< (p/slope pl %) 1e-6) (take 50 (filter #(p/ocean? pl %) (dirs 1500)))))
      (is (some #(> (p/slope pl %) 0.05) land))
      (is (every? #(<= 0.0 (p/slope pl %) (* 0.5 Math/PI)) (take 200 land))))))

(deftest clouds-test
  (testing "cover is an opacity, and more cover means more cloud"
    (let [mean (fn [cover]
                 (let [pl (p/planet {:seed 4 :cloud-cover cover :samples 512})
                       vs (map #(p/cloud-cover pl %) (dirs 2000))]
                   (is (every? #(<= 0.0 % 1.0) vs))
                   (/ (reduce + vs) 2000.0)))
          ms (map mean [0.0 0.25 0.5 0.75 1.0])]
      (is (= ms (sort ms)))
      (is (< (first ms) 0.01))
      (is (> (last ms) 0.9))))

  (testing "a planet can have no clouds at all"
    (let [pl (p/planet {:seed 4 :clouds nil :samples 512})]
      (is (every? #(zero? (p/cloud-cover pl %)) (dirs 200))))))

(deftest color-test
  (testing "every color is a real color"
    (let [pl @world]
      (doseq [d (dirs 2000)]
        (let [c (p/surface-color pl d)]
          (is (= 3 (count c)))
          (is (every? #(<= 0.0 (double %) 1.0) c))))))

  (testing "water is bluer than land"
    (let [pl @world
          blueness (fn [[r _ b]] (- (double b) (double r)))
          ds (dirs 2000)
          wet (map #(blueness (p/surface-color pl %)) (filter #(p/ocean? pl %) ds))
          dry (map #(blueness (p/surface-color pl %)) (remove #(p/ocean? pl %) ds))]
      (is (> (/ (reduce + wet) (count wet)) (/ (reduce + dry) (count dry))))
      (is (every? pos? wet))))

  (testing "deeper water is darker"
    (let [pl @world
          ds (filter #(p/ocean? pl %) (dirs 2000))
          lum (fn [d] (reduce + (p/surface-color pl d)))
          by-depth (sort-by #(p/depth pl %) ds)
          shallow (take 50 by-depth)
          deep (take-last 50 by-depth)]
      (is (> (/ (reduce + (map lum shallow)) 50.0)
             (/ (reduce + (map lum deep)) 50.0)))))

  (testing "the snow line comes down toward the poles"
    ;; Testing the rule rather than a picture: hold the terrain fixed and
    ;; turn the latitude term up, and the only ground that may gain snow
    ;; is the ground away from the equator.
    (let [pl (p/planet {:seed 6 :samples 2048})
          white? (fn [[r g b]] (> (min (double r) (double g) (double b)) 0.75))
          land (remove #(p/ocean? pl %) (dirs 4000))
          polar (filter #(> (abs (p/latitude %)) 0.9) land)
          ;; Exactly on the equator the latitude term is zero, so nothing
          ;; there may change at all.
          equator (remove #(p/ocean? pl %)
                          (map #(p/direction 0.0 %) (range -3.14 3.14 0.004)))
          snowy (fn [ds effect]
                  (count (filter #(white? (p/surface-color pl % {:pole-effect effect})) ds)))]
      (is (seq polar))
      (is (seq equator))
      (is (> (snowy polar 0.5) (snowy polar 0.0)))
      (is (= (snowy equator 0.5) (snowy equator 0.0)))))

  (testing "and lowering the snow line puts snow further down the mountain"
    (let [pl @world
          white? (fn [[r g b]] (> (min (double r) (double g) (double b)) 0.75))
          land (remove #(p/ocean? pl %) (dirs 3000))
          snowy (fn [line] (count (filter #(white? (p/surface-color pl % {:snow-line line})) land)))]
      (is (> (snowy 0.3) (snowy 0.6)))
      (is (zero? (snowy 1.5))))))

(deftest equirectangular-test
  (testing "the map has the size asked for and starts at the pole"
    (let [pl @world
          m (p/equirectangular #(p/elevation pl %) 8 4)]
      (is (= 32 (count m)))
      (is (every? #(Double/isFinite (double %)) m))))

  (testing "a row near the pole barely varies, which is what the projection costs"
    (let [pl @world
          m (p/equirectangular #(p/elevation pl %) 64 32)
          row (fn [j] (subvec m (* j 64) (* (inc j) 64)))
          spread (fn [r] (- (apply max r) (apply min r)))]
      (is (< (spread (row 0)) (spread (row 16)))))))

(deftest custom-terrain-test
  (testing "a planet takes any basis at all, which is the point of the design"
    (let [pl (p/planet {:terrain (fractal/ridged-multifractal
                                  (noise/gradient-basis {:seed 2})
                                  {:octaves 4.0})
                        :samples 1024})]
      (is (every? #(Double/isFinite (p/elevation pl %)) (dirs 500)))
      (is (< (abs (- (:relief pl) (apply max (map #(p/elevation pl %) (dirs 2048)))))
             (* 0.15 (:relief pl)))))))
