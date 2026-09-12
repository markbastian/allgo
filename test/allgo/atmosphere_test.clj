(ns allgo.atmosphere-test
  (:require [allgo.procedural.atmosphere :as atm]
            [clojure.test :refer [deftest is testing]]))

(def ^:private air (atm/atmosphere))
(def ^:private ground [0.0 1.0002 0.0])
(def ^:private fine {:steps 48 :sun-steps 12})

(defn- red [c] (double (nth c 0)))
(defn- blue [c] (double (nth c 2)))

(deftest density-test
  (testing "density is 1 at the datum and falls by e every scale height"
    (let [h (:scale-height air)]
      (is (< (abs (- 1.0 (atm/density air [0.0 1.0 0.0]))) 1e-12))
      (doseq [n (range 1 6)]
        (is (< (abs (- (Math/exp (- (double n)))
                       (atm/density air [0.0 (+ 1.0 (* n h)) 0.0])))
               1e-9)))))

  (testing "and it does not go on rising below the surface"
    ;; Terrain dips under the datum, and an exponential run the other way
    ;; would give the bottom of a valley an absurd amount of air in it.
    (is (== 1.0 (atm/density air [0.0 0.98 0.0])))))

(deftest rayleigh-test
  (testing "scattering goes as one over the fourth power of the wavelength"
    (let [[r g b] (atm/rayleigh-coefficients 1.0)
          [lr lg lb] atm/wavelengths-rgb]
      (is (< (abs (- 1.0 g)) 1e-12))
      (is (< (abs (- (/ b g) (Math/pow (/ (double lg) (double lb)) 4.0))) 1e-12))
      (is (< (abs (- (/ r g) (Math/pow (/ (double lg) (double lr)) 4.0))) 1e-12))))

  (testing "which is why blue scatters several times as hard as red"
    (let [[r _ b] (atm/rayleigh-coefficients 1.0)]
      (is (< 5.0 (/ b r) 6.5))))

  (testing "and the strength argument scales all three together"
    (is (every? #(< (abs (- 3.0 %)) 1e-9)
                (map / (atm/rayleigh-coefficients 3.0) (atm/rayleigh-coefficients 1.0))))))

(deftest phase-test
  (testing "Rayleigh scatters twice as hard along the axis as across it"
    (is (< (abs (- 2.0 (/ (atm/rayleigh-phase 1.0) (atm/rayleigh-phase 0.0)))) 1e-12))
    (is (< (abs (- (atm/rayleigh-phase 1.0) (atm/rayleigh-phase -1.0))) 1e-12)))

  (testing "a phase function integrates to one over the sphere"
    ;; The conservation law that makes it a phase function and not just a
    ;; shape: scattering redistributes light, it does not create it.
    (let [n 20000
          integral (fn [f]
                     (* 2.0 Math/PI (/ 2.0 n)
                        (reduce + (for [i (range n)]
                                    (f (- (/ (* 2.0 (+ 0.5 (double i))) n) 1.0))))))]
      (is (< (abs (- 1.0 (integral atm/rayleigh-phase))) 1e-3))
      (is (< (abs (- 1.0 (integral #(atm/mie-phase 0.0 %)))) 1e-3))
      (is (< (abs (- 1.0 (integral #(atm/mie-phase 0.6 %)))) 2e-2))))

  (testing "Mie is forward-biased, and more so as g rises"
    (is (> (atm/mie-phase 0.76 1.0) (atm/mie-phase 0.76 -1.0)))
    (is (> (atm/mie-phase 0.9 1.0) (atm/mie-phase 0.5 1.0)))
    (is (< (abs (- (atm/mie-phase 0.0 1.0) (atm/mie-phase 0.0 -1.0))) 1e-12))))

(deftest optical-depth-test
  (testing "no path, no air"
    (is (zero? (atm/optical-depth air ground [0.0 1.0 0.0] 0.0 8))))

  (testing "a longer path through air is more air"
    (let [depths (map #(atm/optical-depth air ground [0.0 1.0 0.0] % 64)
                      [0.0 0.002 0.005 0.01 0.02])]
      (is (= depths (sort depths)))
      (is (apply distinct? depths))))

  (testing "and a grazing path is far more air than a vertical one"
    ;; The whole reason sunsets are a different colour from noon.
    (let [up (atm/optical-depth air ground [0.0 1.0 0.0] 0.025 128)
          along (atm/optical-depth air ground [1.0 0.0 0.0] 0.025 128)]
      (is (> along (* 5.0 up)))))

  (testing "refining the quadrature converges"
    (let [a (atm/optical-depth air ground [0.0 1.0 0.0] 0.02 64)
          b (atm/optical-depth air ground [0.0 1.0 0.0] 0.02 256)]
      (is (< (abs (- a b)) (* 0.01 b))))))

(deftest transmittance-test
  (testing "nothing is lost over no distance, and everything over a great one"
    (is (every? #(< (abs (- 1.0 (double %))) 1e-12) (atm/transmittance air 0.0)))
    (is (every? #(< (double %) 1e-6) (atm/transmittance air 10.0))))

  (testing "it only ever falls as the path grows"
    (let [ts (map #(blue (atm/transmittance air %)) [0.0 0.01 0.05 0.2 1.0])]
      (is (= ts (reverse (sort ts))))))

  (testing "and blue is taken out faster than red, which is why the sun reddens"
    (let [t (atm/transmittance air 0.05)]
      (is (< (blue t) (red t))))))

(deftest sky-test
  (testing "the daytime sky overhead is blue"
    (let [c (atm/sky-color air ground [0.0 1.0 0.0] [0.0 1.0 0.0] fine)]
      (is (> (blue c) (red c)))
      (is (every? pos? c))))

  (testing "looking away from a setting sun is redder than looking up at noon"
    ;; Both halves of chapter 18 in one assertion: the long path takes the
    ;; blue out, so what is left over is what makes the sunset.
    (let [noon (atm/sky-color air ground [0.0 1.0 0.0] [0.0 1.0 0.0] fine)
          dusk (atm/sky-color air ground [-1.0 0.0 0.0] [1.0 0.0 0.0] fine)]
      (is (> (/ (red dusk) (blue dusk)) (/ (red noon) (blue noon))))
      (is (> (red dusk) (blue dusk)))))

  (testing "the horizon is brighter than the zenith, the sun being up"
    (let [zenith (atm/sky-color air ground [0.0 1.0 0.0] [0.0 1.0 0.0] fine)
          horizon (atm/sky-color air ground [1.0 0.0 0.0] [0.0 1.0 0.0] fine)]
      (is (> (reduce + horizon) (reduce + zenith)))))

  (testing "a ray that never enters the atmosphere sees nothing"
    (is (nil? (atm/sky-color air [0.0 5.0 0.0] [0.0 1.0 0.0] [0.0 1.0 0.0]))))

  (testing "cutting the integral short gives less sky, not more"
    (let [full (atm/sky-color air ground [1.0 0.0 0.0] [0.0 1.0 0.0] fine)
          near (atm/sky-color air ground [1.0 0.0 0.0] [0.0 1.0 0.0]
                              (assoc fine :t-far 0.002))]
      (is (< (reduce + near) (reduce + full))))))

(deftest aerial-perspective-test
  (let [surface [0.3 0.4 0.2]
        sun [0.0 1.0 0.0]]
    (testing "nothing happens over no distance"
      (is (every? #(< (abs (double %)) 1e-9)
                  (map - (atm/aerial-perspective air ground ground sun surface) surface))))

    (testing "distance pulls a colour towards the haze and away from itself"
      (let [near (atm/aerial-perspective air ground [0.001 1.0002 0.0] sun surface fine)
            far (atm/aerial-perspective air ground [0.02 1.0002 0.0] sun surface fine)
            gap (fn [c] (reduce + (map (fn [a b] (abs (- (double a) (double b)))) c surface)))]
        (is (< (gap near) (gap far)))))

    (testing "and a distant surface ends up bluer than it started"
      (let [far (atm/aerial-perspective air ground [0.05 1.0002 0.0] sun surface fine)]
        (is (> (- (blue far) (blue surface)) (- (red far) (red surface))))))))

(deftest tone-map-test
  (testing "it lands in [0, 1] however bright the input"
    (doseq [c [[0.0 0.0 0.0] [0.5 1.0 2.0] [100.0 100.0 100.0]]]
      (is (every? #(<= 0.0 (double %) 1.0) (atm/tone-map c)))))

  (testing "black stays black and it never inverts the order"
    (is (every? zero? (atm/tone-map [0.0 0.0 0.0])))
    (let [vs (map #(first (atm/tone-map [% 0.0 0.0])) [0.0 0.1 0.5 1.0 4.0])]
      (is (= vs (sort vs)))))

  (testing "more exposure is brighter"
    (is (< (first (atm/tone-map 0.5 [1.0 1.0 1.0]))
           (first (atm/tone-map 2.0 [1.0 1.0 1.0]))))))
