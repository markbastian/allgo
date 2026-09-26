(ns allgo.tet-mesh-test
  (:require [allgo.geometry.tet-mesh :as tm]
            [clojure.test :refer [deftest is testing]]))

(defn- tets [{:keys [tet-ids]}] (partition 4 tet-ids))

(defn- surface-edges [{:keys [surface-tri-ids]}]
  (mapcat (fn [[a b c]] [#{a b} #{b c} #{c a}]) (partition 3 surface-tri-ids)))

(deftest lattice-box-test
  (let [mesh (tm/lattice-box 3 2 4 1.0)]
    (testing "a cell becomes six tetrahedra and a corner lattice of vertices"
      (is (= (* 4 3 5) (/ (count (:verts mesh)) 3)))
      (is (= (* 3 2 4 6) (count (tets mesh)))))

    (testing "the pieces add back up to the box"
      ;; Kuhn's subdivision is exact, so this is equality up to rounding,
      ;; not an approximation.
      (is (< (abs (- 24.0 (tm/volume mesh))) 1e-9)))

    (testing "every tetrahedron is wound the same way"
      (is (every? (fn [[a b c d]] (pos? (tm/tet-volume (:verts mesh) a b c d)))
                  (tets mesh))))

    (testing "the box sits on y=0, centered in x and z"
      (is (= [[-1.5 0.0 -2.0] [1.5 2.0 2.0]] (tm/bounds mesh))))

    (testing "sizing scales the volume cubically"
      (is (< (abs (- (* 8 24.0) (tm/volume (tm/lattice-box 3 2 4 2.0)))) 1e-9)))))

(deftest edges-test
  (let [mesh (tm/lattice-box 2 2 2 1.0)]
    (testing "each edge appears once however many tetrahedra share it"
      (let [pairs (map set (partition 2 (:edge-ids mesh)))]
        (is (= (count pairs) (count (distinct pairs))))
        (is (every? #(= 2 (count %)) pairs) "no edge joins a vertex to itself")))

    (testing "every tetrahedron's six edges are all present"
      (let [have (set (map set (partition 2 (:edge-ids mesh))))]
        (is (every? (fn [tet]
                      (let [t (vec tet)]
                        (every? #(contains? have (set %))
                                (for [i (range 4) j (range (inc i) 4)] [(t i) (t j)]))))
                    (tets mesh)))))))

(deftest surface-test
  (doseq [dims [[1 1 1] [3 2 4]]]
    (let [mesh (apply tm/lattice-box dims)]
      (testing (str dims ": the skin is closed -- every edge shared by two triangles")
        ;; The test that a boundary is really a boundary: an open surface
        ;; would have edges belonging to one triangle, and a doubled one
        ;; would have edges belonging to four.
        (is (every? #(= 2 %) (vals (frequencies (surface-edges mesh))))))

      (testing (str dims ": the skin is only the outside")
        (let [[lo hi] (tm/bounds mesh)
              on-boundary? (fn [i]
                             (let [p (tm/vertex (:verts mesh) i)]
                               (some (fn [axis]
                                       (or (< (abs (- (p axis) (lo axis))) 1e-9)
                                           (< (abs (- (p axis) (hi axis))) 1e-9)))
                                     [0 1 2])))]
          (is (every? on-boundary? (:surface-tri-ids mesh))
              "no interior vertex appears in a surface triangle")))))

  (testing "a single cube's skin is twelve triangles -- two per face"
    (is (= 12 (/ (count (:surface-tri-ids (tm/lattice-box 1 1 1))) 3)))))

(deftest deform-test
  (let [mesh (tm/lattice-box 2 2 2 1.0)]
    (testing "translation moves everything and changes nothing else"
      (let [moved (tm/translate mesh [10.0 5.0 -3.0])]
        (is (< (abs (- (tm/volume mesh) (tm/volume moved))) 1e-9))
        (is (= (:tet-ids mesh) (:tet-ids moved)))
        (is (= [[9.0 5.0 -4.0] [11.0 7.0 -2.0]] (tm/bounds moved)))))

    (testing "a deformation that mirrors would invert the tetrahedra, and is fixed up"
      ;; Reflection flips the sign of every volume. Left alone the solver
      ;; would hold the body inside out, so orientation is restored.
      (let [flipped (tm/deform mesh (fn [[x y z]] [(- x) y z]))]
        (is (every? (fn [[a b c d]] (pos? (tm/tet-volume (:verts flipped) a b c d)))
                    (tets flipped)))
        (is (< (abs (- (tm/volume mesh) (tm/volume flipped))) 1e-9))))

    (testing "rounding a box toward a ball keeps every tetrahedron valid"
      (let [ball (tm/deform mesh (tm/sphere-deformation [0.0 1.0 0.0] 1.0))]
        (is (every? (fn [[a b c d]] (pos? (tm/tet-volume (:verts ball) a b c d)))
                    (tets ball)))
        (is (< (tm/volume ball) (tm/volume mesh)) "and it is smaller than its box")))))
