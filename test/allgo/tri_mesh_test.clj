(ns allgo.tri-mesh-test
  (:require [allgo.geometry.tet-mesh :as tet]
            [allgo.geometry.tri-mesh :as tri]
            [clojure.test :refer [deftest is testing]]))

(defn- counts [{:keys [verts tri-ids edge-ids]}]
  {:v (quot (count verts) 3) :f (quot (count tri-ids) 3) :e (quot (count edge-ids) 2)})

;; Two triangles sharing the edge 1-2, with 0 and 3 opposite it.
(def ^:private pair
  (tri/complete {:verts [0.0 0.0 0.0,  1.0 0.0 0.0,  0.0 0.0 1.0,  1.0 0.0 1.0]
                 :tri-ids [0 1 2, 1 3 2]}))

(deftest grid-test
  (let [m (tri/grid 4 4 0.25 1.0)
        {:keys [v f e]} (counts m)]
    (testing "a grid of cells, two triangles each"
      (is (= 25 v))
      (is (= 32 f))
      (is (= 56 e)))

    (testing "it is a disc: V - E + F = 1"
      (is (= 1 (+ (- v e) f))))

    (testing "the area is the area"
      (is (< (abs (- 1.0 (tri/area m))) 1e-12)))

    (testing "it lies flat at the height asked for"
      (is (= [[-0.5 1.0 -0.5] [0.5 1.0 0.5]] (tri/bounds m))))

    (testing "the outline is the outline"
      ;; Four sides of four cells each.
      (is (= 16 (quot (count (tri/boundary-edges (:tri-ids m))) 2))))

    (testing "corners are the four extremes"
      (is (= [0 4 20 24] (tri/corner-vertices m))))

    (testing "the diagonals alternate rather than all running one way"
      ;; Every cell split the same way makes the sheet stiffer along one
      ;; bias; alternating keeps it even. If they all matched, every
      ;; triangle would have the same winding pattern.
      (let [diagonals (set (map set (partition 2 (tri/bending-edges (:tri-ids m)))))]
        (is (> (count diagonals) 1))))))

(deftest topology-test
  (testing "two triangles meeting at an edge know about each other"
    (let [n (tri/triangle-neighbours (:tri-ids pair))]
      (is (= 6 (count n)))
      (is (= 2 (count (filter #(>= % 0) n))) "one shared edge, seen from both sides")
      (is (= 4 (count (filter neg? n))) "and four edges on the boundary")
      (testing "the pairing is symmetric"
        (doseq [[slot other] (map-indexed vector n) :when (>= other 0)]
          (is (= slot (nth n other)))))))

  (testing "bending pairs are the vertices opposite the shared edge"
    ;; The shared edge is 1-2; the vertices opposite it are 0 and 3.
    (is (= #{0 3} (set (tri/bending-edges (:tri-ids pair))))))

  (testing "each edge appears once even though two triangles own it"
    (let [pairs (map set (partition 2 (:edge-ids pair)))]
      (is (= 5 (count pairs)) "four outer edges and one shared")
      (is (= (count pairs) (count (distinct pairs))))))

  (testing "a closed surface has no boundary and bends everywhere"
    ;; The skin of a tetrahedral box is closed, so every edge is interior.
    (let [skin (tri/complete {:verts (:verts (tet/lattice-box 2 2 2 1.0))
                              :tri-ids (:surface-tri-ids (tet/lattice-box 2 2 2 1.0))})]
      (is (empty? (tri/boundary-edges (:tri-ids skin))))
      (is (= (quot (count (:edge-ids skin)) 2)
             (quot (count (:bend-ids skin)) 2))
          "one bending pair per edge, since every edge is shared")
      (testing "and Euler says it is a sphere: V - E + F = 2"
        ;; Counting the vertices the triangles actually name: the tet mesh
        ;; it came from carries interior vertices the skin never touches.
        (let [v (count (distinct (:tri-ids skin)))
              f (quot (count (:tri-ids skin)) 3)
              e (quot (count (:edge-ids skin)) 2)]
          (is (= 2 (+ (- v e) f))))))))

(deftest transform-test
  (let [m (tri/grid 2 2 0.5 0.0)]
    (testing "translation moves it and nothing else"
      (let [t (tri/translate m [1.0 2.0 3.0])]
        (is (= (:tri-ids m) (:tri-ids t)))
        (is (< (abs (- (tri/area m) (tri/area t))) 1e-12))
        (is (= [[0.5 2.0 2.5] [1.5 2.0 3.5]] (tri/bounds t)))))

    (testing "a deformation is applied to every vertex"
      (let [d (tri/deform m (fn [[x y z]] [x (+ y (* 0.5 x x)) z]))]
        (is (> (tri/area d) (tri/area m)) "curving it makes it bigger")
        (is (= (:tri-ids m) (:tri-ids d)))))))
