(ns allgo.fluid-test
  (:require [allgo.physics.fluid :as fl]
            [clojure.test :refer [deftest is testing]]))

(defn- tunnel
  "A wind tunnel with an obstacle, ready to step."
  ([] (tunnel {}))
  ([{:keys [nx ny speed disk] :or {nx 40 ny 30 speed 2.0 disk true}}]
   (let [f (-> (fl/fluid nx ny (/ 1.0 ny) 1000.0)
               (fl/close-border! true))]
     (fl/wind-tunnel! f speed 0.1)
     (when disk (fl/disk! f 0.4 0.5 0.1))
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
    ;; Dye lives at cell centers, so sampling a center reads that cell.
    (fl/set-smoke! f 3 4 7.0)
    (testing "sampling a cell center returns that cell's value"
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
  (testing "a disk marks cells solid, and clears what it left"
    (let [f (-> (fl/fluid 40 30 (/ 1.0 30)) (fl/close-border! true))]
      (fl/disk! f 0.4 0.5 0.1)
      (is (fl/solid? f (int (* 0.4 30)) (int (* 0.5 30))) "the center is solid")
      (is (not (fl/solid? f 35 5)) "far from it is not")
      ;; Move it; the old position must open up again.
      (fl/disk! f 1.0 0.5 0.1)
      (is (not (fl/solid? f (int (* 0.4 30)) (int (* 0.5 30)))))))

  (testing "flow does not pass through a solid"
    (let [f (tunnel)]
      (advance f 60 {})
      ;; Inside the disk there is nothing moving.
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

(deftest staggering-test
  ;; Each velocity component lives on a different face -- u on vertical
  ;; faces, v on horizontal ones -- so advecting either needs the other
  ;; averaged from the four faces that straddle it. Take the wrong four
  ;; and the average is read from a cell diagonally away, which biases
  ;; every step along that diagonal. Nothing else here catches it: a
  ;; uniform field averages to the same value whichever four you pick.
  (let [nx 20 ny 20 h 0.1 dt 0.01
        ;; A field linear in position interpolates exactly, so whatever
        ;; comes back out names the point it was read from.
        linear-u! (fn [f fx]
                    (dotimes [i (:nx f)]
                      (dotimes [j (:ny f)]
                        (aset ^floats (:u f) (fl/idx f i j)
                              (float (fx (* i h) (* (+ j 0.5) h)))))))
        linear-v! (fn [f fx]
                    (dotimes [i (:nx f)]
                      (dotimes [j (:ny f)]
                        (aset ^floats (:v f) (fl/idx f i j)
                              (float (fx (* (+ i 0.5) h) (* j h)))))))]

    (testing "sampling a field reads it at the point asked for"
      (let [f (fl/fluid nx ny h)]
        (linear-u! f (fn [_ y] y))
        (linear-v! f (fn [x _] x))
        (is (< (abs (- (fl/sample f :u 0.97 1.04) 1.04)) 1e-5))
        (is (< (abs (- (fl/sample f :v 0.97 1.04) 0.97)) 1e-5))))

    (testing "advecting u averages v from the faces straddling it"
      ;; The u face at (i*h, (j+1/2)*h) is straddled by the v faces at
      ;; columns i-1 and i, rows j and j+1, whose centroid is the u face
      ;; itself. So with v = x, the average must be i*h -- and with the
      ;; neighboring four it would be (i+1)*h instead.
      (doseq [[label vf expected]
              [["v = x" (fn [x _] x) (fn [i _] (* i h))]
               ["v = y" (fn [_ y] y) (fn [_ j] (* (+ j 0.5) h))]]]
        (let [f (-> (fl/fluid nx ny h) (fl/close-border!))]
          (linear-u! f (fn [_ y] y))
          (linear-v! f vf)
          (fl/advect-velocity! f dt)
          (doseq [[i j] [[10 10] [7 13] [14 6]]]
            (let [u' (aget ^floats (:u f) (fl/idx f i j))
                  ;; u came back as the y it was traced to, so the
                  ;; transverse velocity used is recoverable exactly.
                  avg-v (/ (- (* (+ j 0.5) h) u') dt)]
              (is (< (abs (- avg-v (expected i j))) 1e-3)
                  (str label " at " [i j] ": averaged " avg-v
                       " instead of " (expected i j))))))))

    (testing "advecting v averages u from the faces straddling it"
      (doseq [[label uf expected]
              [["u = x" (fn [x _] x) (fn [i _] (* (+ i 0.5) h))]
               ["u = y" (fn [_ y] y) (fn [_ j] (* j h))]]]
        (let [f (-> (fl/fluid nx ny h) (fl/close-border!))]
          (linear-u! f uf)
          (linear-v! f (fn [x _] x))
          (fl/advect-velocity! f dt)
          (doseq [[i j] [[10 10] [7 13] [14 6]]]
            (let [v' (aget ^floats (:v f) (fl/idx f i j))
                  avg-u (/ (- (* (+ i 0.5) h) v') dt)]
              (is (< (abs (- avg-u (expected i j))) 1e-3)
                  (str label " at " [i j] ": averaged " avg-u
                       " instead of " (expected i j))))))))))

(defn- stirred
  "A velocity field that badly needs projecting, in a tunnel with an
  obstacle."
  [res]
  (let [f (-> (fl/fluid res res (/ 1.0 res) 1000.0) (fl/close-border! true))]
    (fl/wind-tunnel! f 2.0 0.1)
    (fl/disk! f 0.35 0.5 0.12)
    (dotimes [i res]
      (dotimes [j res]
        (let [k (fl/idx f i j)]
          (aset ^floats (:u f) k (float (Math/sin (* 0.4 i))))
          (aset ^floats (:v f) k (float (Math/cos (* 0.3 j)))))))
    f))

(defn- projected
  "Divergence left after one projection with the given solver."
  [res solver iterations]
  (let [f (stirred res)]
    (fl/project! f {:solver solver :iterations iterations
                    :dt (/ 1.0 60.0) :tolerance 0.0})
    [f (fl/max-divergence f)]))

(deftest solver-test
  (testing "every solver drives the divergence down"
    (doseq [solver [:gauss-seidel :conjugate-gradient :multigrid]]
      (let [before (fl/max-divergence (stirred 32))
            [_ after] (projected 32 solver 64)]
        (is (< after (* 0.2 before)) (str (name solver) " left " after)))))

  (testing "and they agree, because they solve the same system"
    ;; The property that makes them interchangeable. Converged, the
    ;; velocity field does not depend on how the pressure was found.
    (let [[a _] (projected 32 :gauss-seidel 4000)
          [b _] (projected 32 :multigrid 64)
          worst (reduce max 0.0
                        (for [k (range (:n a))]
                          (max (abs (- (aget ^floats (:u a) k) (aget ^floats (:u b) k)))
                               (abs (- (aget ^floats (:v a) k) (aget ^floats (:v b) k))))))]
      (is (< worst 1e-3) (str "velocity fields differ by " worst))))

  (testing "multigrid needs about the same number of cycles whatever the grid"
    ;; The whole point of it. Gauss-Seidel carries information one cell per
    ;; sweep, so a grid twice as fine wants roughly four times the sweeps;
    ;; a V-cycle carries it across the whole grid at once.
    (doseq [res [32 64 128]]
      (let [[_ d] (projected res :multigrid 32)]
        (is (< d 1e-3) (str res " grid left " d " after 32 cycles")))))

  (testing "where Gauss-Seidel cannot get there at all in sixteen times the work"
    (let [[_ d] (projected 128 :gauss-seidel 512)]
      (is (> d 1e-3) (str "512 sweeps reached " d ", so this test has gone stale"))))

  (testing "a closed tank works as well as an open one"
    ;; Closed is the harder case for a solver: with no open cell anywhere
    ;; the pressure is only determined up to a constant.
    (let [f (-> (fl/fluid 48 48 (/ 1.0 48) 1000.0) (fl/close-border!))]
      (dotimes [i 48]
        (dotimes [j 48]
          (let [k (fl/idx f i j)]
            (aset ^floats (:u f) k (float (Math/sin (* 0.4 i))))
            (aset ^floats (:v f) k (float (Math/cos (* 0.3 j)))))))
      (let [before (fl/max-divergence f)]
        (fl/project! f {:solver :multigrid :iterations 32 :dt (/ 1.0 60.0) :tolerance 0.0})
        (is (< (fl/max-divergence f) (* 0.05 before))))))

  (testing "the solvers stay finite"
    (doseq [solver [:gauss-seidel :conjugate-gradient :multigrid]]
      (let [[f _] (projected 64 solver 64)]
        (is (every? #(Float/isFinite (aget ^floats (:u f) %)) (range (:n f)))
            (str (name solver) " left something that is not a number)"))
        (is (every? #(Float/isFinite (aget ^floats (:v f) %)) (range (:n f)))))))

  (testing "the default is still the one it always was"
    (is (= :gauss-seidel (:solver fl/default-solver)))
    (let [[_ a] (projected 32 :gauss-seidel 40)
          f (stirred 32)]
      (fl/project! f {:iterations 40 :dt (/ 1.0 60.0) :tolerance 0.0})
      (is (< (abs (- a (fl/max-divergence f))) 1e-9))))

  (testing "and the positional arity still means Gauss-Seidel"
    (let [[_ a] (projected 32 :gauss-seidel 40)
          f (stirred 32)]
      (fl/project! f 40 (/ 1.0 60.0) 1.9)
      (is (< (abs (- a (fl/max-divergence f))) 1e-9))))

  (testing "an unknown solver is an error, not a silent nothing"
    (is (thrown? clojure.lang.ExceptionInfo
                 (fl/project! (stirred 16) {:solver :telepathy :dt (/ 1.0 60.0)})))))
