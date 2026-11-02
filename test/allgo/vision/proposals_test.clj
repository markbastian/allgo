(ns allgo.vision.proposals-test
  (:require [allgo.random :as random]
            [allgo.vision.detection :as detection]
            [allgo.vision.image :as img]
            [allgo.vision.proposals :as proposals]
            [allgo.vision.segmentation :as segmentation]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

;; Region proposals and bounding-box regression, after Russell and Norvig,
;; *Artificial Intelligence: A Modern Approach*, 4th edition, chapter 25,
;; on drawn scenes whose objects' boxes are known.

(defn- scene
  "An 80 by 60 scene of three objects apart from one another, each
  two-toned -- lighter above, darker below, as a lit object is -- so that
  a segmentation splits it, on a shaded background with noise: `{:image
  :truth}`, the objects' boxes."
  [seed]
  (let [w 80 h 60 rng (random/rng seed)
        candidate (fn [] (let [s (+ 10 (* 14 (rng))) s2 (* s (+ 0.6 (* 0.8 (rng))))
                               x (+ 2 (* (- w s 4) (rng))) y (+ 2 (* (- h s2 4) (rng)))]
                           {:box [x y (+ x s) (+ y s2)] :disk? (< (rng) 0.5)}))
        grown (fn [[x0 y0 x1 y1]] [(- x0 3) (- y0 3) (+ x1 3) (+ y1 3)])
        objects (loop [out []]
                  (if (= 3 (count out))
                    out
                    (let [c (candidate)]
                      (recur (if (every? #(zero? (detection/iou (grown (:box %)) (:box c))) out) (conj out c) out)))))
        truth (mapv (fn [{[x0 y0 x1 y1] :box disk? :disk?}]
                      (if disk?
                        (let [r (/ (min (- x1 x0) (- y1 y0)) 2) cx (/ (+ x0 x1) 2) cy (/ (+ y0 y1) 2)]
                          [(- cx r) (- cy r) (+ cx r) (+ cy r)])
                        [x0 y0 x1 y1]))
                    objects)
        shapes (mapcat (fn [{disk? :disk?} [x0 y0 x1 y1]]
                         (let [v (+ 0.6 (* 0.35 (rng))) v2 (- v 0.25) ym (+ y0 (* (- y1 y0) (+ 0.35 (* 0.3 (rng)))))]
                           (if disk?
                             [{:kind :disk :center [(/ (+ x0 x1) 2) (/ (+ y0 y1) 2)] :radius (/ (- x1 x0) 2)
                               :texture (fn [_ y] (if (< y (/ (+ y0 y1) 2)) v v2))}]
                             [{:kind :rect :box [x0 y0 x1 ym] :value v} {:kind :rect :box [x0 ym x1 y1] :value v2}])))
                       objects truth)]
    {:image (img/add-noise (img/render w h (fn [x _] (+ 0.1 (* 0.15 (/ x w)))) shapes) 0.03 rng)
     :truth truth}))

(defn- found
  "How many of the boxes `truth` some proposal overlaps by at least `t`."
  [truth proposals t]
  (count (filter (fn [g] (some #(>= (detection/iou g %) t) proposals)) truth)))

(defn- jitter
  "A box near `g` overlapping it by at least `min-iou`: shifted up to a
  quarter of its size, scaled by up to e^0.3 each way."
  [rng g min-iou]
  (let [wobble #(- (* 2 (rng)) 1)]
    (first (filter #(>= (detection/iou % g) min-iou)
                   (repeatedly (fn [] (let [[x0 y0 x1 y1] g w (- x1 x0) h (- y1 y0)
                                            cx (+ (/ (+ x0 x1) 2) (* w 0.25 (wobble))) cy (+ (/ (+ y0 y1) 2) (* h 0.25 (wobble)))
                                            sw (* w (math/exp (* 0.3 (wobble)))) sh (* h (math/exp (* 0.3 (wobble))))]
                                        [(- cx (/ sw 2)) (- cy (/ sh 2)) (+ cx (/ sw 2)) (+ cy (/ sh 2))])))))))

(defn- jittered [seeds per rng]
  (vec (for [seed seeds :let [{:keys [image truth]} (scene seed)] g truth _ (range per)]
         [image (jitter rng g 0.5) g])))

(deftest selective-search
  (let [{:keys [image truth]} (scene 3)
        labels (segmentation/graph-based image 3 8)
        states (vec (proposals/selective-search-seq image labels))]
    (testing "a merge a step, down to one region covering the image"
      (is (= (range (count (:regions (first states))) 0 -1) (map (comp count :regions) states)))
      (is (= [[0 0 80 60]] (map :box (vals (:regions (peek states)))))))
    (testing "the segmentation splits each object; grouping puts it back together"
      (is (zero? (found truth (:proposals (first states)) 0.7)))
      (is (= 3 (found truth (proposals/selective-search image labels) 0.7)))))
  (testing "over ten scenes: every object proposed, at a handful of proposals a scene"
    (let [rs (for [seed (range 10)
                   :let [{:keys [image truth]} (scene seed)
                         ps (proposals/selective-search image (segmentation/graph-based image 3 8))]]
               [(found truth ps 0.7) (count ps)])]
      (is (= 30 (reduce + (map first rs))))
      (is (every? #(<= % 12) (map second rs))))))

(deftest box-regression
  (testing "the deltas and their inverse"
    (let [p [10.0 20.0 30.0 35.0] g [12.0 18.0 34.0 40.0]]
      (is (= [0.0 0.0 0.0 0.0] (proposals/box-deltas p p)))
      (is (every? #(< (abs %) 1e-12) (map - g (proposals/apply-deltas p (proposals/box-deltas p g)))))))
  (testing "ridge regression recovers a linear function"
    (let [rng (random/rng 1)
          X (vec (repeatedly 50 (fn [] [(rng) (rng) (rng)])))
          ys (mapv (fn [[a b c]] (+ 0.5 (* 2.0 a) (* -1.0 b) (* 3.0 c))) X)]
      (is (every? #(< (abs %) 1e-6) (map - [0.5 2.0 -1.0 3.0] (proposals/ridge X ys 1e-9))))))
  (testing "learned on boxes near objects in thirty scenes, it moves boxes in ten new ones onto their objects"
    (let [rng (random/rng 7)
          regressor (proposals/box-regressor (jittered (range 100 130) 10 rng) 10 1.5 1.0)
          test (jittered (range 200 210) 5 rng)
          before (map (fn [[_ p g]] (detection/iou p g)) test)
          after (map (fn [[image p g]] (detection/iou (proposals/refine regressor image p) g)) test)
          mean #(/ (reduce + %) (count %))]
      (is (< 0.6 (mean before) 0.65))
      (is (< 0.83 (mean after)))
      (is (every? true? (map <= before after)) "none made worse"))))

(deftest anchors
  (let [k 10 context 1.5
        anchors-of (fn [image] (proposals/anchors image 4 [12 18 26] [0.6 1.0 1.6]))
        rng (random/rng 11)
        ;; 64 anchors from each of twenty scenes, at most half of them positive
        examples (for [seed (range 100 120) :let [{:keys [image truth]} (scene seed)]]
                   (proposals/objectness-examples image truth (anchors-of image) k context 0.6 0.3 64 rng))
        weights (nth (proposals/objectness-seq [(vec (mapcat first examples)) (vec (mapcat second examples))] 2.0) 50)
        regressor (proposals/box-regressor (jittered (range 100 130) 10 (random/rng 7)) k context 1.0)
        results (for [seed (range 200 210)
                      :let [{:keys [image truth]} (scene seed)
                            ps (map :box (proposals/propose weights image (anchors-of image) k context 0.5 0.5))]]
                  {:count (count ps) :anchors (count (anchors-of image))
                   :coarse (found truth ps 0.5) :tight (found truth ps 0.7)
                   :refined (found truth (map #(proposals/refine regressor image %) ps) 0.7)})
        total #(reduce + (map % results))]
    (testing "anchors lie inside the image"
      (is (every? (fn [[x0 y0 x1 y1]] (and (<= 0 x0) (<= 0 y0) (<= x1 80) (<= y1 60))) (anchors-of (:image (scene 1))))))
    (testing "learned objectness keeps a few dozen of some fourteen hundred anchors, and they find most objects roughly"
      (is (every? #(< (:count %) 50) results))
      (is (every? #(> (:anchors %) 1400) results))
      (is (<= 27 (total :coarse))))
    (testing "an anchor's fixed sizes fit few objects tightly; box regression makes nearly all fit"
      (is (<= (total :tight) 15))
      (is (<= 29 (total :refined))))))
