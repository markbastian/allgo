(ns allgo.vision.shape-test
  (:require [allgo.geometry.vec3 :as v]
            [allgo.random :as random]
            [allgo.vision.camera :as camera]
            [allgo.vision.image :as img]
            [allgo.vision.shape :as shape]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

;; Shape from shading, texture and contour, after Russell and Norvig,
;; *Artificial Intelligence: A Modern Approach*, 4th edition, chapter 25,
;; checked on drawn scenes whose shapes are known.

(defn- degrees-between [a b] (math/to-degrees (math/acos (min 1.0 (abs (v/dot a b))))))

(defn- sphere
  "A matte sphere of radius `r` pixels in an `n` by `n` image, albedo 1,
  lit from `light`: the image, its mask and its true normals."
  [n r light]
  (let [c (/ n 2.0)
        normal (fn [x y] (let [dx (/ (- (+ x 0.5) c) r) dy (/ (- (+ y 0.5) c) r) q (+ (* dx dx) (* dy dy))]
                           (when (< q 1.0) [dx dy (math/sqrt (- 1.0 q))])))]
    {:image (img/image n n (fn [x y] (if-let [nv (normal x y)] (camera/lambert 1.0 nv light 1.0 0.0) 0.0)))
     :mask (img/image n n (fn [x y] (some? (normal x y))))
     :truth (img/image n n normal)}))

(deftest shading
  (let [light (v/normalize [0.3 0.2 1.0])
        {:keys [image mask truth]} (sphere 40 16.0 light)
        boundary (shape/occluding-boundary mask)
        states (shape/shading-seq image mask 1.0 light)
        error (fn [st] (let [ns (:px (shape/normals st mask))
                             errs (for [i (range (count ns)) :when (and (ns i) (not (boundary i)))]
                                    (degrees-between (ns i) ((:px truth) i)))]
                         [(/ (reduce + errs) (count errs)) (apply max errs)]))]
    (testing "the occluding contour's normals lie in the image plane, pointing out of the silhouette"
      (doseq [[i [f g]] boundary
              :let [x (mod i 40) y (quot i 40) [nx ny nz] (shape/normal f g)]]
        (is (< (abs nz) 1e-12))
        (is (< (degrees-between [nx ny 0.0] (v/normalize [(- (+ x 0.5) 20.0) (- (+ y 0.5) 20.0) 0.0])) 15.0))))
    (testing "from every inner normal facing the viewer, 41 degrees out on average, to under 5"
      (let [[start] (error (first states))
            [mean worst] (error (nth states 1000))]
        (is (< 40.0 start 42.0))
        (is (< mean 5.0))
        (is (< worst 14.0))))
    (testing "the convex answer, not the concave one the brightness alone allows: the contour decides"
      (let [ns (:px (shape/normals (nth states 1000) mask))]
        (is (every? (fn [i] (or (nil? (ns i)) (boundary i) (pos? (v/dot (ns i) ((:px truth) i))))) (range (count ns))))))))

(defn- ground
  "A camera's `w` by `h` view of the ground plane y = 0 with round texels
  of `radius` on a square lattice `spacing` apart, bright on dark."
  [cam w h spacing radius]
  (img/image w h (fn [x y]
                   (let [[o d] (camera/ray cam [(+ x 0.5) (+ y 0.5)])]
                     (if (>= (d 1) -1e-9)
                       0.2
                       (let [[px _ pz] (v/add-scaled o d (/ (- (o 1)) (d 1)))
                             gx (* spacing (math/round (/ px spacing))) gz (* spacing (math/round (/ pz spacing)))]
                         (if (<= (math/hypot (- px gx) (- pz gz)) radius) 1.0 0.2)))))))

(deftest texture
  (let [f 120.0 center [80.0 60.0]
        cam (camera/look-at [0.0 3.0 0.0] [0.0 0.0 6.0] f center)
        ts (shape/texels (ground cam 160 120 1.0 0.3) f center 0.6 40)
        ;; the ground's normal, up, in the camera's coordinates
        up (camera/to-camera (assoc cam :position [0.0 0.0 0.0]) [0.0 1.0 0.0])]
    (is (< 12 (count ts)))
    (testing "each texel's slant, from its ellipse, within five degrees of the true angle between the ground's normal and its ray"
      (doseq [{:keys [ray slant]} ts]
        (is (< (abs (- (math/to-degrees slant) (degrees-between up ray))) 5.0))))
    (testing "and its tilt, the direction the ground falls away, within ten degrees"
      (doseq [{:keys [tilt frame]} ts
              :let [[e1 e2] frame
                    truth (mod (math/atan2 (v/dot up e2) (v/dot up e1)) math/PI)
                    off (abs (- tilt truth))]]
        (is (< (math/to-degrees (min off (- math/PI off))) 10.0))))
    (testing "the texels' sizes alone place the ground plane, its normal within five degrees"
      (is (< (degrees-between (shape/plane-from-texels ts) up) 5.0)))))

(deftest contour
  (testing "on the ground plane, the lower the nearer: distance and height from where an object stands in the image"
    (let [f 400.0 cam (camera/look-at [0.0 1.5 0.0] [0.0 1.5 10.0] f [320.0 240.0])]
      (doseq [[x z height] [[1.0 4.0 2.0] [-2.0 10.0 0.5] [0.5 25.0 3.0]]
              :let [[_ v-foot] (camera/project cam [x 0.0 z])
                    [_ v-top] (camera/project cam [x height z])]]
        (is (< (abs (- z (shape/ground-distance f 240.0 1.5 v-foot))) 1e-9))
        (is (< (abs (- height (shape/height-on-ground 240.0 1.5 v-foot v-top))) 1e-9)))
      (is (nil? (shape/ground-distance f 240.0 1.5 200.0)) "above the horizon, not on the ground")))
  (testing "T-junctions: where an outline stops against another, the other is in front"
    ;; drawn back to front: a disk, a box over it, a second disk over the box
    (let [labels (img/render 80 60 0 [{:kind :disk :center [26 30] :radius 14 :value 1}
                                      {:kind :rect :box [34 10 60 38] :value 2}
                                      {:kind :disk :center [60 44] :radius 9 :value 3}])]
      (is (= [[34 18] [60 35] [37 38] [53 38]] (map :at (shape/t-junctions labels))))
      (is (every? (fn [{:keys [front behind]}] (every? #(> (long front) (long %)) behind)) (shape/occlusions labels 4))
          "every junction puts the later-drawn region in front")
      (is (= [3 2 1 0] (shape/depth-order labels 4 0)))))
  (testing "forty random rows of five shapes, each drawn partly over the last: the junctions that decide are nearly always right"
    (let [scene (fn [seed]
                  (let [rng (random/rng seed)]
                    (img/render 96 64 0 (loop [out [] x 14.0]
                                          (if (= 5 (count out))
                                            out
                                            (let [y (+ 22 (* 22 (rng))) s (+ 9 (* 6 (rng)))]
                                              (recur (conj out (assoc (if (even? (count out))
                                                                        {:kind :disk :center [x y] :radius s}
                                                                        {:kind :rect :box [(- x s) (- y (* 0.8 s)) (+ x s) (+ y (* 0.8 s))]})
                                                                      :value (inc (count out))))
                                                     (+ x (* 1.3 s)))))))))
          pairs (map #(shape/in-front (scene %) 4) (range 1 41))
          right (fn [ps] (count (filter (fn [[a b]] (> (long a) (long b))) ps)))]
      (is (< 0.95 (/ (reduce + (map right pairs)) (reduce + (map count pairs)))))
      (is (<= 30 (count (filter #(= (right %) (count %)) pairs)))))))
