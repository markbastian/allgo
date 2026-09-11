(ns allgo.fluid-test
  (:require [allgo.physics.fluid :as fl]
            [clojure.test :refer [deftest is testing]]))

(defn- tunnel
  "A wind tunnel with an obstacle, ready to step."
  ([] (tunnel {}))
  ([{:keys [nx ny speed disc] :or {nx 40 ny 30 speed 2.0 disc true}}]
   (let [f (-> (fl/fluid nx ny (/ 1.0 ny) 1000.0)
               (fl/close-border! true))]
     (fl/wind-tunnel! f speed 0.1)
     (when disc (fl/disc! f 0.4 0.5 0.1))
     f)))

(defn- advance [f steps world]
  (dotimes [_ steps]
    (fl/wind-tunnel! f 2.0 0.1)
    (fl/step! f world))
  f)

(defn- finite? [^floats a n]
  (every? (fn [i] (let [x (aget a i)]
                    (and (not (Float/isNaN x)) (not (Float/isInfinite x)))))
          (range n)))

(deftest grid-test
  (let [f (fl/fluid 10 8 0.1 1000.0)]
    (testing "a border cell is added on every side"
      (is (= 12 (:nx f)))
      (is (= 10 (:ny f)))
      (is (= 120 (:n f))))

    (testing "indexing is i*ny + j"
      (is (= 0 (fl/idx f 0 0)))
      (is (= 10 (fl/idx f 1 0)))
      (is (= 13 (fl/idx f 1 3))))

    (testing "everything starts open, and the border closes"
      (is (not (fl/solid? f 0 0)))
      (fl/close-border! f)
      (is (fl/solid? f 0 0))
      (is (fl/solid? f 11 9))
      (is (not (fl/solid? f 5 5))))

    (testing "an inflow border leaves the left column open"
      (let [g (fl/close-border! (fl/fluid 10 8 0.1 1000.0) true)]
        (is (not (fl/solid? g 0 5)))
        (is (fl/solid? g 11 5) "but not the right")))))

(deftest sampling-test
  (let [f (fl/fluid 8 8 1.0 1000.0)]
    ;; Dye lives at cell centres, so sampling a centre reads that cell.
    (fl/set-smoke! f 3 4 7.0)
    (testing "sampling a cell centre returns that cell's value"
      (is (< (abs (- 7.0 (fl/sample f :smoke 3.5 4.5))) 1e-5)))

    (testing "sampling between two cells interpolates"
      (fl/set-smoke! f 4 4 3.0)
      (let [mid (fl/sample f :smoke 4.0 4.5)]
        (is (< (abs (- 5.0 mid)) 1e-5) "halfway between 7 and 3")))

    (testing "sampling outside is clamped rather than exploding"
      (is (not (Double/isNaN (fl/sample f :smoke -50.0 -50.0))))
      (is (not (Double/isNaN (fl/sample f :smoke 1e6 1e6)))))))

(deftest incompressibility-test
  (testing "projection is what drives divergence out"
    (let [f (tunnel)
          before (fl/max-divergence f)]
      (advance f 40 {})
      (is (pos? before) "the inflow starts the field divergent")
      (is (< (fl/max-divergence f) (* 0.1 before))
          "and the projection removes most of it")))

  (testing "more iterations converge further"
    ;; Gauss-Seidel: each sweep carries the correction one cell further,
    ;; so the residual falls as the sweeps accumulate.
    (let [divergence (fn [iters]
                       (let [f (tunnel)]
                         (advance f 30 {:iterations iters})
                         (fl/max-divergence f)))]
      (is (> (divergence 2) (divergence 40)))))

  (testing "over-relaxation converges faster than plain Gauss-Seidel"
    ;; No physical meaning at all; it simply overshoots each local
    ;; correction, and gets to the same place in fewer sweeps.
    (let [divergence (fn [w]
                       (let [f (tunnel)]
                         (advance f 30 {:iterations 8 :over-relaxation w})
                         (fl/max-divergence f)))]
      (is (< (divergence 1.9) (divergence 1.0)))))

  (testing "a fluid at rest with no inflow stays at rest"
    (let [f (-> (fl/fluid 20 20 0.05 1000.0) (fl/close-border!))]
      (dotimes [_ 40] (fl/step! f {:gravity 0.0}))
      (is (< (fl/max-divergence f) 1e-4))
      (is (every? (fn [[i j]] (let [[u v] (fl/velocity-at f i j)]
                                (and (< (abs u) 1e-4) (< (abs v) 1e-4))))
                  (for [i (range 2 20) j (range 2 20)] [i j]))))))

(deftest obstacle-test
  (testing "a disc marks cells solid, and clears what it left"
    (let [f (-> (fl/fluid 40 30 (/ 1.0 30)) (fl/close-border! true))]
      (fl/disc! f 0.4 0.5 0.1)
      (is (fl/solid? f (int (* 0.4 30)) (int (* 0.5 30))) "the centre is solid")
      (is (not (fl/solid? f 35 5)) "far from it is not")
      ;; Move it; the old position must open up again.
      (fl/disc! f 1.0 0.5 0.1)
      (is (not (fl/solid? f (int (* 0.4 30)) (int (* 0.5 30)))))))

  (testing "flow does not pass through a solid"
    (let [f (tunnel)]
      (advance f 60 {})
      ;; Inside the disc there is nothing moving.
      (let [ci (int (* 0.4 30)) cj (int (* 0.5 30))]
        (when (fl/solid? f ci cj)
          (let [[u v] (fl/velocity-at f ci cj)]
            (is (< (abs u) 1e-3))
            (is (< (abs v) 1e-3)))))))

  (testing "and it leaves a wake: slower behind than beside"
    (let [f (tunnel)]
      (advance f 120 {})
      (let [behind (first (fl/velocity-at f 20 15))
            beside (first (fl/velocity-at f 20 3))]
        (is (< behind beside) "the obstacle slows what is behind it")))))

(deftest stability-test
  (testing "it stays finite over a long run with an obstacle in the way"
    (let [f (tunnel)]
      (advance f 400 {})
      (is (finite? (:u f) (:n f)))
      (is (finite? (:v f) (:n f)))
      (is (finite? (:smoke f) (:n f)))))

  (testing "and under gravity with no inflow"
    (let [f (-> (fl/fluid 30 30 0.03) (fl/close-border!))]
      (dotimes [_ 200] (fl/step! f {:gravity -9.8}))
      (is (finite? (:v f) (:n f)))
      (is (< (fl/max-divergence f) 1.0))))

  (testing "dye is carried, not created"
    ;; Semi-Lagrangian advection interpolates between existing values, so
    ;; it can smear the dye but must not manufacture any.
    (let [f (-> (fl/fluid 30 20 0.05) (fl/close-border!))]
      (dotimes [i (:nx f)] (dotimes [j (:ny f)] (fl/set-smoke! f i j 0.0)))
      (doseq [i (range 8 14) j (range 8 14)] (fl/set-smoke! f i j 1.0))
      (dotimes [i (:nx f)] (fl/set-velocity! f i 10 1.0 0.0))
      (let [before (fl/total-smoke f)]
        (dotimes [_ 50] (fl/step! f {}))
        (let [after (fl/total-smoke f)]
          (is (<= after (* 1.05 before)) "no dye is created")
          (is (every? (fn [[i j]] (<= -1e-4 (fl/smoke-at f i j) 1.0001))
                      (for [i (range 1 (dec (:nx f))) j (range 1 (dec (:ny f)))] [i j]))
              "and none goes out of range"))))))
