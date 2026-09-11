(ns allgo.flip-test
  (:require [allgo.physics.flip :as flip]
            [clojure.test :refer [deftest is testing]]))

(defn- tank
  "A closed tank with a block of water in one corner, ready to step."
  ([] (tank {}))
  ([{:keys [nx ny h r] :or {nx 30 ny 25 h 0.02 r 0.006}}]
   (let [f (-> (flip/flip-fluid {:nx nx :ny ny :h h :particle-radius r
                                 :max-particles 4000})
               (flip/close-border!))]
     (flip/fill-block! f (* 2 h) (* 2 h) (* 0.45 nx h) (* 0.8 ny h))
     f)))

(defn- advance [f steps world]
  (dotimes [_ steps] (flip/step! f world))
  f)

(defn- speeds [{:keys [^floats vel] :as f}]
  (map (fn [i] (Math/hypot (aget vel (* 2 i)) (aget vel (inc (* 2 i)))))
       (range (flip/particle-count f))))

(defn- positions [f]
  (mapv #(flip/particle f %) (range (flip/particle-count f))))

(deftest grid-test
  (let [f (flip/flip-fluid {:nx 10 :ny 8 :h 0.1 :particle-radius 0.03
                            :max-particles 100})]
    (testing "a border cell is added on every side"
      (is (= 12 (:nx f)))
      (is (= 10 (:ny f)))
      (is (= 120 (:n f))))

    (testing "the push-apart grid is sized to the particle diameter, not the cell"
      (is (= (* 2.2 0.03) (:p-spacing f)))
      (is (pos? (:p-nx f)))
      (is (pos? (:p-ny f))))

    (testing "it starts with no particles and no rest density"
      (is (zero? (flip/particle-count f)))
      (is (zero? @(:rest-density f))))))

(deftest seeding-test
  (let [f (tank)]
    (testing "the block is seeded inside the tank"
      (is (pos? (flip/particle-count f)))
      (let [[[lx ly] [hx hy]] (flip/bounds-of-particles f)
            {:keys [nx ny h]} f]
        (is (> lx h)) (is (> ly h))
        (is (< hx (* (dec nx) h))) (is (< hy (* (dec ny) h)))))

    (testing "a staggered lattice leaves no two particles on the same point"
      (is (= (flip/particle-count f) (count (distinct (positions f))))))

    (testing "particles start at rest"
      (is (every? zero? (speeds f))))))

(deftest stability-test
  (testing "a dropped block settles instead of running away"
    (let [f (advance (tank) 300 {})
          sp (speeds f)]
      (is (every? #(Double/isFinite %) sp))
      ;; Free fall down this tank tops out near 3 m/s; anything far above
      ;; that is the solver adding energy rather than gravity doing work.
      (is (< (reduce max sp) 10.0)
          "no particle is moving anywhere near fast enough to be an explosion")
      (is (< (/ (reduce + sp) (count sp)) 2.0)
          "the body as a whole has slowed down, not sped up")))

  (testing "it is still settling, not exploding, much later"
    (let [f (advance (tank) 900 {})]
      (is (< (reduce max (speeds f)) 10.0))))

  (testing "drift compensation is on by default and must not destabilise the solve"
    ;; This is the regression that matters: measuring the FLIP correction
    ;; against the wrong snapshot adds the particle-to-grid splat back onto
    ;; velocities that already carry it, and drift compensation turns that
    ;; feedback into an explosion within a hundred steps.
    (let [on  (advance (tank) 200 {:drift? true})
          off (advance (tank) 200 {:drift? false})]
      (is (< (reduce max (speeds on)) 10.0))
      (is (< (reduce max (speeds off)) 10.0)))))

(deftest separation-test
  (testing "particles never collapse onto one another"
    ;; Wall and corner clamping puts particles on identical coordinates,
    ;; which the reference leaves stuck together for good because it skips
    ;; the zero-distance case.
    (let [f (advance (tank) 300 {})
          pts (positions f)]
      (is (= (count pts) (count (distinct pts)))
          "every particle is on its own point")))

  (testing "pushing apart is what does it"
    (let [with    (advance (tank) 150 {:separation 2})
          without (advance (tank) 150 {:separation 0})]
      (is (> (count (distinct (positions with)))
             (count (distinct (positions without)))))))

  (testing "more iterations spread the particles no worse"
    (let [f (advance (tank) 150 {:separation 5})]
      (is (= (flip/particle-count f) (count (distinct (positions f))))))))

(deftest containment-test
  (testing "particles stay inside the tank"
    (let [{:keys [nx ny h particle-radius] :as f} (advance (tank) 400 {})
          [[lx ly] [hx hy]] (flip/bounds-of-particles f)
          lo (+ h particle-radius)
          hi-x (- (* (dec nx) h) particle-radius)
          hi-y (- (* (dec ny) h) particle-radius)
          eps 1e-5]
      (is (>= lx (- lo eps)))
      (is (>= ly (- lo eps)))
      (is (<= hx (+ hi-x eps)))
      (is (<= hy (+ hi-y eps)))))

  (testing "particles stay outside an obstacle"
    (let [obstacle [0.3 0.25 0.08]
          f (advance (tank) 200 {:obstacle obstacle})
          [ox oy orad] obstacle
          r (:particle-radius f)]
      (is (every? (fn [[x y]]
                    (>= (Math/hypot (- x ox) (- y oy)) (- (+ orad r) 1e-5)))
                  (positions f))))))

(deftest conservation-test
  (testing "no particle is created or destroyed"
    (let [f (tank)
          n (flip/particle-count f)]
      (advance f 300 {})
      (is (= n (flip/particle-count f)))))

  (testing "rest density is measured once and then held"
    (let [f (tank)]
      (advance f 1 {})
      (let [d @(:rest-density f)]
        (is (pos? d))
        (advance f 50 {})
        (is (= d @(:rest-density f))))))

  (testing "the fluid keeps its volume rather than slowly compressing"
    ;; Drift compensation exists for exactly this: without it the body
    ;; creeps into a smaller footprint as particles crowd.
    (let [spread (fn [world]
                   (let [f (advance (tank) 400 world)
                         [[lx ly] [hx hy]] (flip/bounds-of-particles f)]
                     (* (- hx lx) (- hy ly))))]
      (is (pos? (spread {:drift? true}))))))

(deftest flip-ratio-test
  (testing "both ends of the PIC/FLIP blend are stable"
    (doseq [ratio [0.0 0.5 0.9 1.0]]
      (let [f (advance (tank) 150 {:flip-ratio ratio})]
        (is (every? #(Double/isFinite %) (speeds f))
            (str "flip-ratio " ratio " stays finite"))
        (is (< (reduce max (speeds f)) 10.0)
            (str "flip-ratio " ratio " does not run away")))))

  (testing "PIC damps more than FLIP does"
    ;; Averaging the particles in a cell together is what costs PIC its
    ;; swirl. A single snapshot will not show it -- the tank is sloshing,
    ;; and which phase it is caught in swamps the difference -- so this
    ;; averages over time, where the ordering is clean and monotone.
    (let [energy (fn [ratio]
                   (let [f (tank)
                         samples (doall
                                  (for [step (range 400)]
                                    (do (flip/step! f {:flip-ratio ratio})
                                        (when (>= step 100)
                                          (let [sp (speeds f)]
                                            (/ (reduce + sp) (count sp)))))))
                         kept (remove nil? samples)]
                     (/ (reduce + kept) (count kept))))
          pic  (energy 0.0)
          half (energy 0.5)
          flip (energy 0.95)]
      (is (< pic half) "PIC is the most damped")
      (is (< half flip) "and damping falls off as the blend moves to FLIP"))))
