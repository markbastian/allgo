(ns allgo.polygon-test
  (:require [allgo.geometry.polygon :as poly]
            [clojure.test :refer [deftest is testing]]))

(def ^:private square [[0.0 0.0] [10.0 0.0] [10.0 10.0] [0.0 10.0]])

(defn- close? [a b] (< (abs (- (double a) (double b))) 1e-9))

(deftest area-test
  (testing "the shoelace formula, signed by winding"
    (is (close? 100.0 (poly/signed-area square)))
    (is (close? -100.0 (poly/signed-area (vec (reverse square)))))
    (is (= [5.0 5.0] (poly/centroid square))))
  (testing "a right triangle's centroid is a third of the way along each leg"
    (let [[x y] (poly/centroid [[0.0 0.0] [6.0 0.0] [0.0 3.0]])]
      (is (close? 2.0 x))
      (is (close? 1.0 y)))))

(deftest predicates-test
  (is (poly/convex? square))
  (is (not (poly/convex? (vec (reverse square)))))
  (is (not (poly/convex? [[0.0 0.0] [10.0 0.0] [5.0 2.0] [10.0 10.0] [0.0 10.0]])))
  (is (poly/contains-point? square [3.0 4.0]))
  (is (not (poly/contains-point? square [11.0 4.0])))
  (testing "segments: where they cross, and where they do not"
    (is (= [0.5 0.5] (poly/segment-intersection [0.0 0.0] [2.0 2.0] [0.0 2.0] [2.0 0.0])))
    (is (nil? (poly/segment-intersection [0.0 0.0] [1.0 0.0] [0.0 1.0] [1.0 1.0])))
    (is (close? 1.0 (poly/segment-distance [5.0 1.0] [0.0 0.0] [10.0 0.0])))))

(deftest clip-test
  (testing "a square cut down the middle: two halves, the area shared, tags kept"
    (let [tags [:s :e :n :w]
          [l lt] (poly/clip square tags [5.0 0.0] [0.0 1.0] :cut)
          [r rt] (poly/clip square tags [5.0 0.0] [0.0 -1.0] :cut)]
      (is (close? 50.0 (poly/area l)))
      (is (close? 50.0 (poly/area r)))
      (is (= #{:n :w :s :cut} (set lt)))
      (is (= #{:s :e :n :cut} (set rt)))
      (is (= 1 (count (filter #{:cut} lt))))))
  (testing "a line missing the polygon leaves all or nothing"
    (is (= square (poly/clip square [20.0 0.0] [0.0 1.0])))
    (is (empty? (poly/clip square [20.0 0.0] [0.0 -1.0])))))

(deftest inset-test
  (testing "each side moved in by its own distance"
    (let [in (poly/inset square [1.0 2.0 1.0 2.0])]
      (is (close? (* 6.0 8.0) (poly/area in)))
      (is (= [[2.0 1.0] [8.0 1.0] [8.0 9.0] [2.0 9.0]] (mapv #(mapv (fn [v] (/ (Math/round (* 1e9 v)) 1e9)) %) in)))))
  (testing "and nil once it would turn inside out"
    (is (nil? (poly/inset square 6.0)))))

(deftest simplify-test
  (testing "points on a straight run and near-duplicates go, with the tags of what stays"
    (let [pts [[0.0 0.0] [5.0 0.0] [10.0 0.0] [10.0 10.0] [10.0 10.1] [0.0 10.0]]
          [p t] (poly/simplify pts [:a :b :c :d :e :f] 0.5)]
      (is (= [[0.0 0.0] [10.0 0.0] [10.0 10.0] [0.0 10.0]] p))
      (is (= [:a :c :d :f] t)))))
