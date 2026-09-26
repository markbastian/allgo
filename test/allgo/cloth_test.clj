(ns allgo.cloth-test
  (:require [allgo.geometry.tri-mesh :as tri]
            [allgo.physics.xpbd :as x]
            [clojure.test :refer [deftest is testing]]))

(defn- finite? [v] (and (not (Double/isNaN v)) (not (Double/isInfinite v))))
(defn- all-finite? [body] (every? finite? (x/positions body)))
(defn- ys [body] (map second (partition 3 (x/positions body))))

(defn- hung
  "A sheet pinned by its four corners."
  ([] (hung {}))
  ([opts]
   (let [mesh (tri/grid 8 8 0.1 1.5)
         body (x/cloth mesh opts)]
     (x/pin! body (tri/corner-vertices mesh))
     [mesh body])))

(deftest construction-test
  (let [mesh (tri/grid 4 4 0.25 1.0)
        body (x/cloth mesh)]
    (testing "cloth is two distance constraints, not a new kind"
      (is (= 2 (count (:constraints body))))
      (is (= 25 (:n body))))

    (testing "mass follows area, so no particle is left massless"
      (is (every? #(pos? (x/inv-mass body %)) (range (:n body)))))

    (testing "the masses add up to the area"
      ;; The check that says the distribution is right rather than merely
      ;; plausible. Summing inverse masses -- as the tutorial does --
      ;; fails this, and inverts which particles are heavy.
      (doseq [density [1.0 7.5]]
        (let [b (x/cloth mesh {:density density})
              total (reduce + (map #(/ 1.0 (x/inv-mass b %)) (range (:n b))))]
          (is (< (abs (- (* density (tri/area mesh)) total)) 1e-9) (str "density " density)))))

    (testing "a vertex owning twice the area is twice the mass"
      (let [mass #(/ 1.0 (x/inv-mass body %))]
        ;; vertex 0 touches 2 triangles, 2 touches 4, 12 touches 8.
        (is (< (abs (- (* 2 (mass 0)) (mass 2))) 1e-9))
        (is (< (abs (- (* 4 (mass 0)) (mass 12))) 1e-9))))

    (testing "a corner belongs to fewer triangles, so it is lighter"
      (let [corner (x/inv-mass body 0)
            middle (x/inv-mass body 12)]
        (is (> corner middle) "lighter means larger inverse mass")))

    (testing "density scales mass without moving anything"
      (let [heavy (x/cloth mesh {:density 4.0})]
        (is (< (abs (- (* 4 (x/inv-mass heavy 0)) (x/inv-mass body 0))) 1e-9))
        (is (= (x/positions body) (x/positions heavy)))))

    (testing "a mesh with no interior edges gets no bending constraint"
      (let [single (tri/complete {:verts [0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 1.0]
                                  :tri-ids [0 1 2]})]
        (is (empty? (:bend-ids single)))
        (is (= 1 (count (:constraints (x/cloth single)))))))))

(deftest hanging-test
  (testing "pinned corners stay exactly where they were put"
    (let [[mesh body] (hung)
          corners (tri/corner-vertices mesh)
          before  (mapv #(x/particle body %) corners)]
      (dotimes [_ 120] (x/step! body {:floor nil}))
      (is (= before (mapv #(x/particle body %) corners)))
      (is (all-finite? body))))

  (testing "the middle sags under gravity, but does not fall away"
    (let [[_ body] (hung)
          middle (quot (:n body) 2)
          y0 (second (x/particle body middle))]
      (dotimes [_ 120] (x/step! body {:floor nil}))
      (let [y1 (second (x/particle body middle))]
        (is (< y1 y0) "it sagged")
        (is (> y1 (- y0 1.0)) "and the pins held it up")
        (is (= 4 (count (filter #(= 1.5 %) (ys body))))
            "only the four pinned particles are still at the start height"))))

  (testing "it does not stretch: edges stay near their rest length"
    (let [[mesh body] (hung {:stretch-compliance 0.0})]
      (dotimes [_ 200] (x/step! body {:floor nil}))
      (let [pos  (x/positions body)
            rest 0.1
            worst (reduce max 0.0
                          (for [[a b] (partition 2 (:edge-ids mesh))
                                :let [pa (* 3 a) pb (* 3 b)
                                      d (Math/sqrt (reduce + (for [k (range 3)]
                                                               (let [t (- (nth pos (+ pa k))
                                                                          (nth pos (+ pb k)))]
                                                                 (* t t)))))]]
                            (abs (- d rest))))]
        ;; Edges are 0.1 or 0.1*sqrt(2) long; check against the shorter.
        (is (< worst 0.06) "no edge is stretched far past its rest length")))))

(deftest bending-test
  ;; Two triangles sharing an edge, with vertices 0 and 3 opposite it.
  ;; Folding along that edge is exactly what moves 0 and 3, so the
  ;; constraint is tested where it acts rather than through the shape of a
  ;; whole sheet -- where gravity and stretching drown it out.
  (let [mesh   (tri/complete {:verts [0.0 0.0 0.0, 1.0 0.0 0.0, 0.0 0.0 1.0, 1.0 0.0 1.0]
                              :tri-ids [0 1 2, 1 3 2]})
        rest-d (Math/sqrt 2.0)
        fold   (fn [bend]
                 (let [b (x/cloth mesh {:bend-compliance bend})]
                   (x/set-particle! b 3 [1.0 0.9 1.0])
                   (dotimes [_ 200] (x/step! b {:gravity [0 0 0] :floor nil :damping 5.0}))
                   (let [[ax ay az] (x/particle b 0)
                         [dx dy dz] (x/particle b 3)]
                     (abs (- (Math/sqrt (+ (* (- ax dx) (- ax dx))
                                           (* (- ay dy) (- ay dy))
                                           (* (- az dz) (- az dz))))
                             rest-d)))))]

    (testing "the constraint spans the two vertices opposite the shared edge"
      (is (= [0 3] (:bend-ids mesh))))

    (testing "rigid bending restores the fold exactly"
      (is (< (fold 0.0) 1e-6)))

    (testing "and compliance is how far it declines to"
      (let [errors (map fold [0.0 1e-4 1e-2 1.0])]
        (is (apply <= errors) "monotone in compliance")
        (is (< (first errors) 1e-6) "rigid holds the fold out entirely")
        (is (> (last errors) 1e-3) "limp visibly does not")))

    (testing "a sheet crushed flat recovers without blowing up"
      (let [[_ body] (hung {:bend-compliance 0.0})]
        (dotimes [i (:n body)]
          (let [[px _ pz] (x/particle body i)]
            (x/set-particle! body i [px 0.0 pz])))
        (dotimes [_ 120] (x/step! body {:floor nil}))
        (is (all-finite? body))))))

(deftest floor-test
  (testing "a dropped sheet lands on the floor and stays above it"
    (let [mesh (tri/grid 6 6 0.15 1.0)
          body (x/cloth mesh {:bend-compliance 2.0})]
      (dotimes [_ 240] (x/step! body {}))
      (is (all-finite? body))
      (is (>= (apply min (ys body)) -1e-3))
      (is (< (apply max (ys body)) 1.0) "and it did fall"))))

(deftest sphere-constraint-test
  (let [center [0.0 0.55 0.0]
        r      0.33
        drop-sheet (fn [friction]
                     (let [mesh (tri/grid 20 20 0.08 1.25)
                           body (cond-> (x/cloth mesh {:bend-compliance 0.05})
                                  true (x/add-constraint (x/sphere-constraint center r friction)))]
                       (dotimes [_ 300] (x/step! body {:gravity [0.0 -9.8 0.0]
                                                       :substeps 12 :damping 0.4 :floor 0.0}))
                       body))
        distance-from-center (fn [[px py pz]]
                               (Math/sqrt (+ (* px px)
                                             (* (- py 0.55) (- py 0.55))
                                             (* pz pz))))
        on-ball (fn [body]
                  (count (filter #(< (abs (- (distance-from-center %) r)) 0.03)
                                 (partition 3 (x/positions body)))))]

    (testing "nothing ends up inside the sphere, at any friction"
      (doseq [mu [0.0 0.35 0.9]]
        (let [body (drop-sheet mu)]
          (is (every? #(>= (distance-from-center %) (- r 1e-6))
                      (partition 3 (x/positions body)))
              (str "friction " mu)))))

    (testing "frictionless, a sheet slides off and pools on the floor"
      ;; Correct, and useless as a drape -- which is what friction is for.
      (let [body (drop-sheet 0.0)]
        (is (zero? (on-ball body)))
        (is (< (apply max (ys body)) 0.5) "it all ended below the ball")))

    (testing "with friction it grips and drapes"
      (let [body (drop-sheet 0.35)]
        (is (pos? (on-ball body)))
        (is (> (apply max (ys body)) (- (+ 0.55 r) 0.02))
            "and the top of the sheet is resting on the top of the ball")))

    (testing "more friction holds more of the sheet on"
      (is (< (on-ball (drop-sheet 0.0)) (on-ball (drop-sheet 0.35)))))

    (testing "it is enforced every substep, not once a frame"
      ;; A single frame is ten-odd solver steps; a sheet dropped from high
      ;; above passes through anything checked only between frames.
      (let [mesh (tri/grid 10 10 0.1 6.0)
            body (-> (x/cloth mesh {:bend-compliance 0.05})
                     (x/add-constraint (x/sphere-constraint center r 0.3)))]
        (dotimes [_ 200] (x/step! body {:gravity [0.0 -20.0 0.0] :substeps 10 :floor 0.0}))
        (is (every? #(>= (distance-from-center %) (- r 1e-6))
                    (partition 3 (x/positions body)))
            "nothing tunneled through")))))
