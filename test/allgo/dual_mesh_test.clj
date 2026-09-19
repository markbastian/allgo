(ns allgo.dual-mesh-test
  (:require [allgo.geometry.dual-mesh :as dm]
            [clojure.test :refer [deftest is testing]]))

(def ^:private bounds [[0.0 0.0] [1000.0 1000.0]])

(defn- points [n seed]
  (let [r (java.util.Random. seed)]
    (vec (repeatedly n #(vector (* 1000.0 (.nextDouble r)) (* 1000.0 (.nextDouble r)))))))

(defn- area [ring]
  (let [v (vec ring) n (count v)]
    (abs (* 0.5 (reduce + (for [i (range n)]
                            (let [[x0 y0] (v i)
                                  [x1 y1] (v (mod (inc i) n))]
                              (- (* x0 y1) (* x1 y0)))))))))

(defn- signed-area [ring]
  (let [v (vec ring) n (count v)]
    (* 0.5 (reduce + (for [i (range n)]
                       (let [[x0 y0] (v i)
                             [x1 y1] (v (mod (inc i) n))]
                         (- (* x0 y1) (* x1 y0))))))))

(deftest structure-test
  (let [m (dm/mesh (points 400 1) bounds)
        {:keys [centers corners edges]} m]

    (testing "every input point became a center"
      (is (= 400 (count centers)))
      (is (= (range 400) (map :id centers))))

    (testing "every cell is a polygon"
      ;; The sentinel ring exists so that this holds for every center
      ;; rather than only the interior ones: a hull cell without it is
      ;; unbounded, and there is no polygon to draw.
      (is (every? #(>= (count (:corners %)) 3) centers))
      (is (every? #(>= (count (:neighbors %)) 2) centers)))

    (testing "corners are ordered, so the cell is a simple ring"
      ;; Sorted by angle about the center, which for a convex cell is the
      ;; ring order. A wrong order gives a self-crossing star with a much
      ;; smaller area than the convex hull of the same points.
      (is (every? (fn [c] (pos? (area (dm/polygon m c)))) centers))
      (is (every? (fn [c] (pos? (signed-area (dm/polygon m c)))) centers)
          "counter-clockwise"))

    (testing "adjacency is mutual"
      (is (every? (fn [c] (every? #(some #{(:id c)} (:neighbors (centers %)))
                                  (:neighbors c)))
                  centers))
      (is (every? (fn [v] (every? #(some #{(:id v)} (:adjacent (corners %)))
                                  (:adjacent v)))
                  corners)))

    (testing "an edge joins two distinct corners and two centers"
      (is (every? (fn [e] (and (= 2 (count (:corners e)))
                               (apply not= (:corners e))
                               (= 2 (count (:centers e)))))
                  edges)))

    (testing "the dual really is dual: a center's edges are its corners' edges"
      ;; Every edge bounding a cell has both its endpoints among that
      ;; cell's corners. This is the property the whole namespace exists
      ;; to provide, and the one a naive build gets wrong.
      (is (every? (fn [c]
                    (let [cs (set (:corners c))]
                      (every? (fn [eid] (every? cs (:corners (edges eid))))
                              (:borders c))))
                  centers)))

    (testing "an edge's centers are real ids, or nil at the rim"
      ;; Sentinel centers are dropped from the mesh, so an edge that had
      ;; one must not keep an id pointing past the end of the vector --
      ;; the island passes never look, but anything drawing the map does.
      (is (every? (fn [e] (every? #(or (nil? %) (< -1 % (count centers)))
                                  (:centers e)))
                  edges))
      (is (some (fn [e] (some nil? (:centers e))) edges)
          "the rim of the map has one-sided edges")
      (is (every? (fn [e] (some some? (:centers e))) edges)
          "and no edge is sentinel on both sides"))

    (testing "ids are self-consistent"
      (is (= (range (count corners)) (map :id corners)))
      (is (= (range (count edges)) (map :id edges))))))

(deftest border-test
  (let [m (dm/mesh (points 400 2) bounds)
        {:keys [centers]} m]

    (testing "the map has an edge, and it is not everything"
      (let [border (filter :border? centers)]
        (is (pos? (count border)))
        (is (< (count border) (/ (count centers) 2)))))

    (testing "interior cells stay near the bounds they were asked for"
      ;; The sentinel ring sits one spacing outside, so nothing real
      ;; should sprawl far beyond the rectangle. A ring placed at the
      ;; distance `delaunay/bounded-cells` uses would fail this by a mile.
      (let [gap (dm/spacing bounds 400)]
        (is (every? (fn [c]
                      (every? (fn [[x y]]
                                (and (< (- (* 3 gap)) x (+ 1000.0 (* 3 gap)))
                                     (< (- (* 3 gap)) y (+ 1000.0 (* 3 gap)))))
                              (dm/polygon m c)))
                    (remove :border? centers)))))

    (testing "the border marking goes all the way round the map"
      ;; The property the ocean flood fill in `allgo.procedural.island`
      ;; depends on: whatever cell you land in at the edge of the
      ;; rectangle is marked, so the fill starts from the whole perimeter
      ;; and no stretch of coast is missed.
      ;;
      ;; Tested by asking which cell actually owns each perimeter point,
      ;; which is the nearest center. Proximity to the edge is not the
      ;; same question and gets a different answer: a cell whose site sits
      ;; 20 units in can still be screened from the edge by another, and
      ;; is then correctly not a border cell.
      (let [owner (fn [p]
                    (apply min-key
                           (fn [c] (let [[cx cy] (:point c) [x y] p]
                                     (+ (* (- cx x) (- cx x)) (* (- cy y) (- cy y)))))
                           centers))
            perimeter (concat (for [t (range 0 101)] [(* 10.0 t) 0.0])
                              (for [t (range 0 101)] [(* 10.0 t) 1000.0])
                              (for [t (range 0 101)] [0.0 (* 10.0 t)])
                              (for [t (range 0 101)] [1000.0 (* 10.0 t)]))]
        (is (every? (comp :border? owner) perimeter))))))

(deftest relax-test
  (testing "relaxation evens out the cell sizes"
    ;; Random points clump, and clumped points make cells of wildly
    ;; different sizes. This is the whole reason the pass exists, so it
    ;; is measured rather than assumed.
    (let [spread (fn [rounds]
                   (let [m (dm/relaxed (points 400 3) bounds rounds)
                         as (map #(area (dm/polygon m %)) (remove :border? (:centers m)))
                         mean (/ (reduce + as) (count as))]
                     (/ (Math/sqrt (/ (reduce + (map #(let [d (- % mean)] (* d d)) as))
                                      (count as)))
                        mean)))]
      (is (< (spread 2) (* 0.8 (spread 0))))))

  (testing "relaxed points stay inside the bounds"
    (let [m (dm/mesh (points 300 4) bounds)]
      (is (every? (fn [[x y]] (and (<= 0.0 x 1000.0) (<= 0.0 y 1000.0)))
                  (dm/relax m)))))

  (testing "zero rounds is just the mesh"
    (is (= (map :point (:centers (dm/relaxed (points 200 5) bounds 0)))
           (points 200 5)))))

(deftest determinism-test
  (testing "the same points give the same mesh"
    (is (= (dm/mesh (points 300 6) bounds)
           (dm/mesh (points 300 6) bounds)))))
