(ns allgo.vision.boundaries-test
  (:require [allgo.random :as random]
            [allgo.vision.boundaries :as bd]
            [allgo.vision.edges :as edges]
            [allgo.vision.image :as img]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

;; Boundary detection by learning, after Russell and Norvig, *Artificial
;; Intelligence: A Modern Approach*, 4th edition, section 25.3.4: scenes
;; of textured regions, where an edge detector fires on every stripe.

(deftest half-disk-test
  (testing "a disk cut in two along a diameter: two equal halves, mirror images"
    (doseq [theta [0.0 (/ math/PI 4) (/ math/PI 2) 1.1]]
      (let [[one other] (bd/half-disks 5 theta)]
        (is (= (count one) (count other)))
        (is (= (set one) (set (map (fn [[x y]] [(- x) (- y)]) other))))))))

(deftest chi-squared-test
  (testing "0 for the same histogram, 1 for two with nothing in common"
    (is (zero? (bd/chi-squared [0.5 0.5 0.0] [0.5 0.5 0.0])))
    (is (= 1.0 (bd/chi-squared [1.0 0.0] [0.0 1.0])))))

(def ^:private horizontal (fn [_ y] (if (< (mod y 4) 2) 0.25 0.75)))
(def ^:private vertical (fn [x _] (if (< (mod x 4) 2) 0.25 0.75)))
(def ^:private dots (fn [x y] (if (and (< (mod x 5) 2) (< (mod y 5) 2)) 0.9 0.4)))

(defn- scene
  "A background of horizontal stripes, a square of vertical ones, a square
  of dots and a plain disk, placed at random: `[image labels]`."
  [rng]
  (let [shapes (vec (for [t [vertical dots 0.8]]
                      (let [x (+ 4 (long (* (rng) 24))) y (+ 4 (long (* (rng) 24))) s (+ 12 (long (* (rng) 8)))]
                        (if (fn? t)
                          {:kind :rect :box [x y (+ x s) (+ y s)] :texture t}
                          {:kind :disk :center [(+ x 6) (+ y 6)] :radius 7 :value t}))))]
    [(img/add-noise (img/render 48 48 horizontal shapes) 0.03 rng)
     (img/render 48 48 0 (map-indexed (fn [i s] (-> s (dissoc :texture) (assoc :value (inc i)))) shapes))]))

(defn- within-2 [pixels [x y]] (some (fn [[a b]] (and (<= (abs (- a x)) 2) (<= (abs (- b y)) 2))) pixels))

(defn- scores [found truth]
  {:precision (/ (count (filter #(within-2 truth %) found)) (double (max 1 (count found))))
   :recall (/ (count (filter #(within-2 found %) truth)) (double (max 1 (count truth))))})

(deftest pb-test
  (let [rng (random/rng 1)
        thetas (mapv #(* % (/ math/PI 8)) (range 8))
        w (nth (bd/train-seq (bd/examples (vec (repeatedly 4 #(scene rng))) 60 3 5 thetas 8 1.0 rng) 1.0) 500)
        inner (fn [ps] (set (filter (fn [[x y]] (and (<= 5 x 42) (<= 5 y 42))) ps)))
        results (for [seed (range 10 14)]
                  (let [[image labels] (scene (random/rng seed))
                        truth (inner (bd/true-boundary labels))
                        thinned (bd/thin (bd/pb w image 5 thetas 8 1.0))]
                    {:pb (scores (inner (for [y (range 48) x (range 48) :when (> (img/at thinned x y) 0.5)] [x y])) truth)
                     :canny (scores (inner (edges/canny image)) truth)}))]
    (testing "the learned weights count a difference in texture as well as in brightness"
      (is (every? pos? (rest w))))
    (testing "on new scenes, Canny fires on the stripes inside the regions: under 85% of its edges within two pixels of a boundary"
      (is (every? #(< (get-in % [:canny :precision]) 0.85) results)))
    (testing "Pb compares the two sides' textures, alike within a region however busy: nearly 90% of what it marks is boundary, and it misses almost none"
      (is (every? #(> (get-in % [:pb :precision]) 0.88) results))
      (is (every? #(> (get-in % [:pb :recall]) 0.95) results)))))
