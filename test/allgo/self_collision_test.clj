(ns allgo.self-collision-test
  (:require [allgo.geometry.tri-mesh :as tri]
            [allgo.physics.xpbd :as x]
            [allgo.spatial.hash :as sh]
            [clojure.test :refer [deftest is testing]]))

(defn- pair-body
  "Two free particles, resting `rest-gap` apart, currently `gap` apart."
  [rest-gap gap thickness]
  (let [body (x/body [0.0 0.0 0.0 (double rest-gap) 0.0 0.0])]
    (x/set-inv-mass! body 0 1.0)
    (x/set-inv-mass! body 1 1.0)
    (let [body (x/add-constraint body (x/self-collision-constraint body thickness 0.0))]
      (x/set-particle! body 1 [(double gap) 0.0 0.0])
      body)))

(defn- separation [body]
  (abs (- (first (x/particle body 1)) (first (x/particle body 0)))))

(deftest adjacency-test
  (testing "the neighbor list matches the pairs brute force finds"
    (let [p (double-array (for [x (range 5) y (range 5) z (range 2) c [x y z]] (double c)))
          n 50
          h (sh/spatial-hash 1.0 n)
          adj (sh/adjacency h p n 1.01)
          brute (set (for [i (range n) j (range n)
                           :when (< j i)
                           :let [a (* 3 i) b (* 3 j)
                                 d (reduce + (for [k (range 3)]
                                               (let [t (- (aget p (+ a k)) (aget p (+ b k)))] (* t t))))]
                           :when (<= d (* 1.01 1.01))]
                       [i j]))]
      (is (= brute (set (for [i (range n) j (sh/adjacent adj i)] [i j]))))
      (is (= (count brute) (:pairs adj)))))

  (testing "each pair is listed once, under the higher index"
    (let [p (double-array (repeat 30 0.0))          ; ten coincident points
          h (sh/spatial-hash 1.0 10)
          adj (sh/adjacency h p 10 1.0)]
      ;; every pair of ten points, once each
      (is (= 45 (:pairs adj)))
      (is (every? (fn [i] (every? #(< % i) (sh/adjacent adj i))) (range 10)))))

  (testing "the index is monotone and spans the whole list"
    (let [p (double-array (for [i (range 40) c [(* 0.3 i) 0.0 0.0]] (double c)))
          h (sh/spatial-hash 0.5 40)
          {:keys [^ints starts pairs]} (sh/adjacency h p 40 0.5)]
      (is (apply <= (map #(aget starts %) (range 41))))
      (is (zero? (aget starts 0)))
      (is (= pairs (aget starts 40))))))

(deftest projection-test
  (testing "an overlapping pair is separated to exactly the thickness"
    (let [body (pair-body 1.0 0.05 0.2)]
      (x/prepare! (first (:constraints body)) body)
      (x/substep! body {:gravity [0 0 0] :floor nil} (/ 1.0 600))
      (is (< (abs (- 0.2 (separation body))) 1e-9))))

  (testing "a pair already further apart than the thickness is untouched"
    (let [body (pair-body 1.0 0.5 0.2)]
      (x/prepare! (first (:constraints body)) body)
      (x/substep! body {:gravity [0 0 0] :floor nil} (/ 1.0 600))
      (is (< (abs (- 0.5 (separation body))) 1e-9))))

  (testing "particles that were always close are held at their rest gap"
    ;; The rule that keeps a fine mesh from inflating: neighbors in the
    ;; sheet are closer together than the cloth is thick, and pushing them
    ;; to the thickness would blow the sheet up.
    (let [body (pair-body 0.05 0.05 0.2)]
      (dotimes [_ 40] (x/step! body {:gravity [0 0 0] :floor nil}))
      (is (< (abs (- 0.05 (separation body))) 1e-9))))

  (testing "an immovable particle takes none of the correction"
    (let [body (pair-body 1.0 0.05 0.2)]
      (x/set-inv-mass! body 0 0.0)
      (dotimes [_ 5] (x/step! body {:gravity [0 0 0] :floor nil}))
      (is (= [0.0 0.0 0.0] (x/particle body 0))))))

(deftest speed-limit-test
  (testing "the limit is a fifth of the thickness per step"
    (is (< (abs (- 2.4 (x/speed-limit 0.2 (/ 1.0 60.0)))) 1e-9)))

  (testing "the world clamps speed to it, without turning anything"
    (let [body (x/body [0.0 0.0 0.0])]
      (x/set-inv-mass! body 0 1.0)
      (x/set-velocity! body 0 [300.0 400.0 0.0])       ; speed 500
      (x/step! body {:gravity [0 0 0] :floor nil :max-velocity 5.0 :substeps 1})
      (let [[px py _] (x/particle body 0)
            moved (Math/hypot px py)]
        ;; one step of 1/60 at the limit
        (is (< (abs (- (/ 5.0 60.0) moved)) 1e-9) "clamped to the limit")
        (is (< (abs (- 0.6 (/ px moved))) 1e-9) "and still heading the same way")))))

(deftest cloth-test
  (let [n 16 spacing 0.06 thickness (* 0.7 spacing)
        mesh (tri/grid n n spacing 0.6)
        row-a (vec (range (inc n)))
        row-b (vec (map #(+ % (* (inc n) n)) (range (inc n))))
        sheet-w (* n spacing)
        crush (fn [self?]
                (let [body (cond-> (x/cloth mesh {:bend-compliance 5.0})
                             self? (x/add-constraint
                                    (x/self-collision-constraint (x/cloth mesh {}) thickness 0.1)))
                      world {:gravity [0 -6.0 0] :substeps 10 :damping 0.8 :floor nil
                             :max-velocity (when self? (x/speed-limit thickness (/ 1.0 60.0)))}]
                  (x/pin! body (concat row-a row-b))
                  (dotimes [f 260]
                    (let [t (min 1.0 (/ (inc f) 180.0))
                          z (* 0.5 sheet-w (- 1.0 (* 0.97 t)))]
                      (doseq [i row-a] (let [[px py _] (x/particle body i)]
                                         (x/set-particle! body i [px py (- z)])))
                      (doseq [i row-b] (let [[px py _] (x/particle body i)]
                                         (x/set-particle! body i [px py z])))
                      (x/step! body world)))
                  body))
        closest-distant (fn [body]
                          (let [pts (vec (partition 3 (x/positions body)))
                                row #(quot % (inc n))]
                            (reduce min 1e9
                                    (for [i (range (count pts)) j (range (count pts))
                                          :when (and (< i j) (> (abs (- (row i) (row j))) 4))
                                          :let [[ax ay az] (pts i) [bx by bz] (pts j)]]
                                      (Math/sqrt (+ (* (- ax bx) (- ax bx))
                                                    (* (- ay by) (- ay by))
                                                    (* (- az bz) (- az bz))))))))]

    (testing "crushed without it, distant parts of the sheet pass through each other"
      (let [gap (closest-distant (crush false))]
        (is (< gap (* 0.2 thickness)) "they end up essentially coincident")))

    (testing "with it, they are held apart"
      (let [loose (closest-distant (crush false))
            tight (closest-distant (crush true))]
        (is (> tight (* 5 loose)) "an order of magnitude better")
        ;; Not the full thickness: this is an iterative projection and a
        ;; hard squeeze wins locally. It stops the sheet passing through
        ;; itself, which is what it is for.
        (is (> tight (* 0.5 thickness)))))

    (testing "and adding it does not inflate the sheet"
      ;; Without the rest-distance rule this is what goes wrong: every
      ;; neighbor is inside the thickness, so the constraint pushes the
      ;; whole mesh apart.
      (let [plain (x/cloth mesh {})
            selfy (-> (x/cloth mesh {})
                      (x/add-constraint (x/self-collision-constraint (x/cloth mesh {}) thickness 0.0)))]
        (doseq [b [plain selfy]]
          (x/pin! b row-a)
          (dotimes [_ 120] (x/step! b {:gravity [0 -9.8 0] :substeps 10 :floor nil})))
        (let [span (fn [b] (let [xs (map first (partition 3 (x/positions b)))]
                             (- (apply max xs) (apply min xs))))]
          (is (< (abs (- (span plain) (span selfy))) (* 0.05 (span plain)))
              "the sheet is the same size either way"))))))
