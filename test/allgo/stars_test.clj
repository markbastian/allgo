(ns allgo.stars-test
  (:require [allgo.astro.frames :as frames]
            [allgo.astro.stars :as stars]
            [allgo.geometry.vec3 :as v]
            [clojure.java.io :as io]
            [clojure.test :refer [deftest is testing]]))

(def ^:private catalogue
  (delay (stars/parse (slurp (io/resource "public/data/bsc5.tsv")))))

(defn- hr [n] (some #(when (= n (:hr %)) %) @catalogue))

(defn- close? [a b tol] (< (abs (- (double a) (double b))) tol))

(deftest parse-test
  (testing "every star in the catalogue with a place"
    ;; 9110 entries, fourteen of which -- novae, clusters -- have none.
    (is (= 9096 (count @catalogue)))
    (is (every? #(and (:ra %) (:dec %) (:vmag %)) @catalogue)))

  (testing "and the right numbers for the ones everybody knows"
    (let [sirius (hr 2491)]
      (is (= "9Alp CMa" (:name sirius)))
      (is (close? -1.46 (:vmag sirius) 1e-9))
      (is (close? 101.287083 (Math/toDegrees (:ra sirius)) 1e-6))
      (is (close? -16.716111 (Math/toDegrees (:dec sirius)) 1e-6)))
    (is (= "M1-2Ia-Iab" (:sp (hr 2061))) "Betelgeuse")
    (is (= -0.04 (:vmag (hr 5340))) "Arcturus")))

(deftest direction-test
  (testing "a star's direction is a unit vector toward its place"
    (let [polaris (stars/direction (hr 424))]
      (is (close? 1.0 (v/length polaris) 1e-12))
      ;; Forty-four arcminutes from the pole in 2000.
      (is (close? 0.736 (Math/toDegrees (Math/acos (nth polaris 2))) 0.01))))

  (testing "and proper motion carries it along"
    ;; Arcturus moves 2.28 arcseconds a year, the fastest in the sky it
    ;; is a naked-eye star of. A century on it is 228 arcseconds away.
    (let [arcturus (hr 5340)
          then (stars/direction arcturus)
          later (stars/direction arcturus (+ 51544.5 (* 100 365.25)))
          moved (* 3600.0 (Math/toDegrees (Math/acos (v/dot then later))))]
      (is (close? 228.0 moved 1.0) (str "moved " moved " arcsec"))))

  (testing "a star with no proper motion given stays put"
    (let [s (dissoc (hr 2491) :pm-ra :pm-dec)]
      (is (= (stars/direction s) (stars/direction s 60000.0))))))

(deftest appearance-test
  (testing "the colors come out the way the stars look"
    (let [[r _ b] (stars/color (:bv (hr 2061)))]
      (is (> r b) "Betelgeuse is red"))
    (let [[r _ b] (stars/color (:bv (hr 1713)))]
      (is (> b r) "Rigel is blue-white"))
    (is (every? #(<= 0.0 % 1.0) (mapcat #(stars/color (:bv %)) @catalogue))))

  (testing "the Sun comes out near its measured temperature"
    ;; B-V 0.65 for the Sun, 5772 K measured.
    (is (close? 5772.0 (stars/temperature 0.65) 150.0)))

  (testing "five magnitudes is a factor of a hundred"
    (is (close? 100.0 (/ (stars/flux 1.0) (stars/flux 6.0)) 1e-9))))

(deftest names-and-figures-test
  (testing "the IAU's names land on the right catalogue stars"
    (let [names (stars/parse-names (slurp (io/resource "public/data/star-names.tsv")))]
      (is (= "Sirius" (names 2491)))
      (is (= "Betelgeuse" (names 2061)))
      (is (= "Polaris" (names 424)))
      (is (every? (set (map :hr @catalogue)) (keys names)) "and every one is in the catalogue")))

  (testing "every constellation, and each figure runs star to star"
    (let [cons (stars/parse-constellations (slurp (io/resource "public/data/constellations.tsv")))
          orion (some #(when (= "Ori" (:abbr %)) %) cons)
          ;; Every star, not only the bright ones: Orion's shield is an arc
          ;; of fourth- and fifth-magnitude stars.
          bright @catalogue
          nearest (fn [[ra dec]]
                    (let [u (stars/unit ra dec)]
                      (apply min (map #(Math/toDegrees (Math/acos (min 1.0 (v/dot u (stars/direction %))))) bright))))]
      ;; Serpens is two figures, head and tail, so 89.
      (is (= 89 (count cons)))
      (is (= "Orion" (:name orion)))
      ;; The figures come from a different source than the catalogue, so
      ;; this is the check that the two agree about where the stars are:
      ;; Orion's every vertex is a bright star, to a hundredth of a degree.
      (is (every? #(< (nearest %) 0.01) (apply concat (:lines orion)))
          (str "furthest from a star: " (apply max (map nearest (apply concat (:lines orion)))))))))

(deftest y-up-test
  (testing "the drawing frame is a rotation of EME2000, not a reflection"
    ;; A reflection keeps every star on the sky and every figure on its
    ;; stars, and draws the whole sky back to front: only its handedness
    ;; gives it away. A rotation keeps a x b = c between the three axes.
    (let [[x y z] (map frames/y-up [[1.0 0.0 0.0] [0.0 1.0 0.0] [0.0 0.0 1.0]])]
      (is (= z (v/cross x y)) "right-handed in, right-handed out")
      (is (= [0.0 1.0 0.0] z) "the celestial pole is up")
      (is (= [1.0 0.0 0.0] x) "the equinox is along x")))

  (testing "and Orion comes out the way round a star chart draws him"
    ;; Look toward Orion from the middle of the sky with north up, as the
    ;; demo's camera does, and project his four corners the way three.js
    ;; does: right is the view direction crossed with up. On any chart --
    ;; and in the sky -- Betelgeuse, Bellatrix, Rigel and Saiph run
    ;; clockwise: shoulders left and right, then the right knee, then the
    ;; left. A mirrored sky runs them the other way.
    (let [corner (fn [n] (frames/y-up (stars/direction (hr n))))
          [betelgeuse bellatrix rigel saiph :as cs] (map corner [2061 1790 1713 2004])
          ahead (v/normalize (reduce v/add cs))
          right (v/normalize (v/cross ahead [0.0 1.0 0.0]))
          up (v/cross right ahead)
          screen (fn [p] [(v/dot p right) (v/dot p up)])
          [p q r s] (map screen [betelgeuse bellatrix rigel saiph])
          area (reduce + (map (fn [[x1 y1] [x2 y2]] (- (* x1 y2) (* x2 y1)))
                              [p q r s] [q r s p]))]
      ;; Clockwise, with y up, is a negative signed area.
      (is (neg? area) (str "Orion's corners run " (if (neg? area) "clockwise" "counterclockwise")))
      (is (< (first (screen betelgeuse)) (first (screen bellatrix)))
          "Betelgeuse is on the left, as seen")
      (is (> (second (screen betelgeuse)) (second (screen rigel)))
          "and above Rigel"))))
