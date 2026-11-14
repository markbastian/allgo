(ns allgo.robotics.clearance-test
  (:require [allgo.geometry.polygon :as polygon]
            [allgo.robotics.clearance :as clearance]
            [clojure.test :refer [deftest is testing]]))

;; The maximum-clearance roadmap, after LaValle, *Planning Algorithms*,
;; Cambridge University Press, 2006, section 6.2.3: a room with a block
;; in it, the corridor of Figure 6.8 run round an obstacle, and a room
;; with a nonconvex L and a triangle.

(def corridor {:bounds [[0 0] [10 4]] :obstacles [[[3 1.5] [7 1.5] [7 2.5] [3 2.5]]]})

(def cluttered {:bounds [[0 0] [12 8]]
                :obstacles [[[2 2] [6 2] [6 3] [3 3] [3 6] [2 6]]
                            [[8 1] [10.5 4] [7.5 5]]]})

(defn- boundary-distance
  "How far `p` is from every edge of the room and the obstacles, measured
  segment by segment -- independently of the features."
  [{:keys [bounds obstacles]} p]
  (let [[[x0 y0] [x1 y1]] bounds
        room [[x0 y0] [x1 y0] [x1 y1] [x0 y1]]]
    (reduce min (for [poly (cons room obstacles) [a b] (polygon/edges poly)] (polygon/segment-distance p a b)))))

(defn- near? [a b] (< (polygon/distance a b) 1e-6))

(deftest corridor-test
  (let [rm (clearance/roadmap corridor)]
    (testing "the pieces Figure 6.9 predicts: straight midlines either side and at the ends of the block, parabolas round its corners, and a diagonal into each corner of the room"
      (is (= 12 (count (:pieces rm))))
      (is (= {:straight 8 :parabolic 4} (frequencies (map #(if (= 2 (count (:points %))) :straight :parabolic) (:pieces rm)))))
      (is (some (fn [{:keys [points]}] (and (near? [3.0 0.75] (first points)) (near? [7.0 0.75] (peek points)))) (:pieces rm))))
    (testing "twelve nodes: the room's corners, where three walls are equally near beside the block's ends, and where the midlines turn into parabolas"
      (is (= 12 (count (:nodes rm))))
      (doseq [p [[0 0] [10 0] [10 4] [0 4] [1.5 1.5] [8.5 1.5] [1.5 2.5] [8.5 2.5] [3 0.75] [7 0.75] [3 3.25] [7 3.25]]]
        (is (some #(near? p %) (:nodes rm)) (str p))))
    (testing "retraction moves straight away from the nearest wall until another is as near"
      (is (near? [5.0 0.75] (clearance/retract rm [5.0 0.2])))
      (is (near? [9.0 3.0] (clearance/retract rm [9.0 3.5]))))
    (testing "a route keeps the most clearance any path can: half the gap beside the block"
      (let [path (clearance/route rm [1.0 1.0] [9.0 3.5])]
        (is (= [1.0 1.0] (first path)))
        (is (= [9.0 3.5] (peek path)))
        (is (< (abs (- 0.75 (reduce min (map #(boundary-distance corridor %) (butlast path))))) 1e-9))))
    (testing "a disk too wide for the gap has no way past the block; a narrower one has"
      (is (nil? (clearance/route (clearance/roadmap corridor :radius 0.8) [1.5 2.0] [8.5 2.0])))
      (is (some? (clearance/route (clearance/roadmap corridor :radius 0.7) [1.5 2.0] [8.5 2.0]))))))

(deftest equidistance-test
  (testing "every point of the roadmap, measured edge by edge, is as near its two features as to any edge at all"
    (doseq [world [corridor cluttered]
            :let [rm (clearance/roadmap world)
                  by-id (into {} (map (juxt :id identity)) (:features rm))]
            {:keys [points between]} (:pieces rm)
            p points]
      (let [[f g] (map by-id between)
            r (boundary-distance world p)]
        (is (< (abs (- r (clearance/distance-to f p))) 1e-6))
        (is (< (abs (- r (clearance/distance-to g p))) 1e-6))))))

(deftest cluttered-test
  (let [rm (clearance/roadmap cluttered)]
    (testing "the L's inner corner, on the free side, has a piece of its own between its two arms"
      (is (some (fn [{:keys [between]}] (= #{[0 :edge 2] [0 :edge 3]} (set between))) (:pieces rm))))
    (testing "the roadmap is in one piece, as the free space is"
      (let [fill (fn [start] (loop [seen #{start} todo [start]]
                               (if-let [n (peek todo)]
                                 (let [new (remove seen (keys ((:graph rm) n)))]
                                   (recur (into seen new) (into (pop todo) new)))
                                 seen)))]
        (is (= (count (:nodes rm)) (count (fill 0))))))
    (testing "routes stay free, from inside the L's crook to the far side of the triangle"
      (let [path (clearance/route rm [4.0 4.5] [11.0 2.0])]
        (is (= [11.0 2.0] (peek path)))
        (is (every? #(clearance/free-point? cluttered %) path))
        (is (every? #(> (boundary-distance cluttered %) 0.3) (butlast (rest path))))))))
