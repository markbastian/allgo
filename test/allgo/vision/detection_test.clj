(ns allgo.vision.detection-test
  (:require [allgo.random :as random]
            [allgo.vision.detection :as det]
            [allgo.vision.image :as img]
            [clojure.test :refer [deftest is testing]]))

;; Classification and detection, after Russell and Norvig, *Artificial
;; Intelligence: A Modern Approach*, 4th edition, chapter 25: disks to be
;; found among squares, in drawn scenes whose disks are known.

(deftest boxes-test
  (testing "intersection over union"
    (is (== 1.0 (det/iou [0 0 4 4] [0 0 4 4])))
    (is (== 0.0 (det/iou [0 0 4 4] [4 0 8 4])))
    (is (== (/ 4.0 28.0) (det/iou [0 0 4 4] [2 2 6 6]))))
  (testing "non-maximum suppression keeps the best of each overlapping cluster"
    (let [kept (det/non-maximum-suppression
                [{:box [0 0 10 10] :score 0.9} {:box [1 1 11 11] :score 0.8}
                 {:box [20 0 30 10] :score 0.7} {:box [21 0 31 10] :score 0.95}]
                0.3)]
      (is (= [[21 0 31 10] [0 0 10 10]] (map :box kept)))))
  (testing "windows tile the image at the stride"
    (is (= 9 (count (det/windows {:w 8 :h 8} 4 2))))))

(def ^:private size 12)

(defn- scene [rng]
  (let [slots (vec (random/shuffle rng (for [x [8 20 32 44] y [8 22]] [x y])))
        disks (mapv (fn [[x y]] [(+ x (- (* 2 (rng)) 1.0)) (+ y (- (* 2 (rng)) 1.0))]) (take 3 slots))
        squares (take 3 (drop 3 slots))
        shapes (concat (map (fn [c] {:kind :disk :center c :radius (+ 3.5 (rng)) :value 0.9}) disks)
                       (map (fn [[x y]] (let [r (+ 3 (rng))] {:kind :rect :box [(- x r) (- y r) (+ x r) (+ y r)] :value 0.9})) squares))]
    {:image (img/add-noise (img/render 56 32 0.1 shapes) 0.05 rng)
     :boxes (mapv (fn [[x y]] (let [x0 (long (Math/round (- x 6))) y0 (long (Math/round (- y 6)))] [x0 y0 (+ x0 size) (+ y0 size)])) disks)}))

(defn- evaluate [net scenes]
  (for [{:keys [image boxes]} scenes]
    (let [dets (det/detect net image size 1 0.5 0.3)
          hits (count (filter (fn [b] (some #(> (det/iou b (:box %)) 0.5) dets)) boxes))]
      {:hits hits :misses (- (count boxes) hits) :false (- (count dets) hits)})))

(deftest detection-test
  (let [rng (random/rng 3)
        train (vec (repeatedly 6 #(scene rng)))
        ex (map #(det/training-windows (:image %) (:boxes %) size 10 0.2 rng) train)
        pos (vec (mapcat :positives ex))
        neg (vec (mapcat :negatives ex))
        net (nth (det/train-seq (det/classifier size 4 8 (random/rng 4)) pos neg size (random/rng 5)) 10)
        hard (vec (mapcat #(det/hard-negatives net (:image %) (:boxes %) size 1 0.3 0.2) train))
        net2 (nth (det/train-seq net pos (into neg hard) size (random/rng 6)) 5)
        tests (vec (repeatedly 4 #(scene rng)))
        before (evaluate net tests)
        after (evaluate net2 tests)
        total (fn [rs k] (reduce + (map k rs)))]
    (testing "the network finds every disk"
      (is (zero? (total after :misses))))
    (testing "and mining the windows it wrongly fired on makes it far pickier"
      (is (< (total after :false) (total before :false)))
      (is (<= (total after :false) 1)))
    (testing "its average precision, every window over 0.05 ranked, is perfect: every disk scores above every false alarm"
      (is (= 1.0 (det/average-precision
                  (det/precision-recall (map (fn [{:keys [image boxes]}] [(det/detect net2 image size 1 0.05 0.3) boxes]) tests) 0.5)))))
    (testing "scoring each window's shifted and mirrored copies too, and averaging, classifies new windows as well -- here, all of them"
      (let [ex (map #(det/training-windows (:image %) (:boxes %) size 40 0.2 rng) tests)]
        (is (every? #(> (det/score-augmented net2 % size) 0.5) (mapcat :positives ex)))
        (is (every? #(<= (det/score-augmented net2 % size) 0.5) (mapcat :negatives ex)))))))

(deftest roi-pooling-test
  (let [image (img/image 4 4 (fn [x y] (+ x (* 4 y))))]
    (testing "a box cut into a grid of cells, each cell's largest value"
      (is (= [5 7 13 15] (det/roi-pool image [0 0 4 4] 2 2)))
      (testing "-- a box of any shape to the same size"
        (is (= [5 7 9 11] (det/roi-pool image [1 1 4 3] 2 2)))
        (is (= [15] (det/roi-pool image [0 0 4 4] 1 1)))))))

(deftest evaluation-test
  ;; two images, three objects. In score order: a hit, a second box on the
  ;; same object (a false alarm), a hit in the second image, a miss, and
  ;; the last object found.
  (let [results [[[{:box [0 0 10 10] :score 0.9} {:box [1 0 11 10] :score 0.8}
                   {:box [40 40 50 50] :score 0.6} {:box [20 20 30 30] :score 0.5}]
                  [[0 0 10 10] [20 20 30 30]]]
                 [[{:box [5 5 15 15] :score 0.7}] [[5 5 15 15]]]]
        curve (det/precision-recall results 0.5)]
    (testing "precision and recall after each detection, best first"
      (is (= [[1.0 (/ 1.0 3)] [0.5 (/ 1.0 3)] [(/ 2.0 3) (/ 2.0 3)] [0.5 (/ 2.0 3)] [0.6 1.0]]
             (mapv (juxt :precision :recall) curve))))
    (testing "average precision, the best precision at each recall or beyond: 1/3 + 2/9 + 1/5 = 34/45"
      (is (< (abs (- (/ 34.0 45) (det/average-precision curve))) 1e-12)))))
