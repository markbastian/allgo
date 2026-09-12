(ns allgo.qaeb-test
  (:require [allgo.geometry.vec3 :as v3]
            [allgo.procedural.planet :as planet]
            [allgo.procedural.qaeb :as qaeb]
            [clojure.test :refer [deftest is testing]]))

(defn- sphere
  "Altitude above a sphere of radius `r` about the origin -- a surface
  whose intersections are known exactly, so the march can be graded."
  [r]
  (fn [p] (- (v3/length p) (double r))))

(deftest sphere-test
  (let [alt (sphere 1.0)
        eye [0.0 0.0 5.0]]
    (testing "a ray down the axis hits where the algebra says"
      (let [{:keys [t point]} (qaeb/intersect alt eye [0.0 0.0 -1.0]
                                              {:epsilon 1e-4 :t-min 1.0 :t-max 10.0})]
        (is (< (abs (- 4.0 t)) 1e-3))
        (is (< (abs (- 1.0 (v3/length point))) 1e-3))))

    (testing "an off-axis ray does too"
      (doseq [dx [0.05 0.1 0.15]]
        (let [d (v3/normalize [dx 0.0 -1.0])
              ;; Closest approach and the half-chord, straight out of the
              ;; quadratic.
              b (v3/dot eye d)
              c (- (v3/dot eye eye) 1.0)
              exact (- (- b) (Math/sqrt (- (* b b) c)))
              {:keys [t]} (qaeb/intersect alt eye d {:epsilon 1e-4 :t-min 1.0 :t-max 10.0})]
          (is (< (abs (- exact t)) 1e-3)))))

    (testing "a ray that passes by returns nothing"
      (is (nil? (qaeb/intersect alt eye (v3/normalize [1.0 0.0 -1.0])
                                {:epsilon 1e-3 :t-min 1.0 :t-max 20.0}))))

    (testing "a smaller epsilon costs more steps and buys more accuracy"
      (let [exact 4.0
            run (fn [e] (qaeb/intersect alt eye [0.0 0.0 -1.0] {:epsilon e :t-min 1.0 :t-max 10.0}))
            coarse (run 1e-2)
            fine (run 1e-4)]
        (is (< (:steps coarse) (:steps fine)))
        (is (<= (abs (- exact (:t fine))) (abs (- exact (:t coarse)))))))

    (testing "starting inside the surface reports the start"
      (let [{:keys [t]} (qaeb/intersect alt [0.0 0.0 0.5] [0.0 0.0 1.0] {})]
        (is (< t 1e-2))))

    (testing "a slope bound turns the march into a sphere trace, exactly"
      ;; With the bound in hand no step can overshoot, so the answer
      ;; stops depending on epsilon at all.
      (let [{:keys [t]} (qaeb/intersect alt eye [0.0 0.0 -1.0]
                                        {:epsilon 1.0 :max-slope 1.0 :t-min 1.0 :t-max 10.0
                                         :min-step 1e-7 :max-steps 200})]
        (is (< (abs (- 4.0 t)) 1e-3))))))

(deftest coherence-test
  (testing "starting the march late gives the same hit, for less"
    (let [alt (sphere 1.0)
          eye [0.0 0.0 5.0]
          cold (qaeb/intersect alt eye [0.0 0.0 -1.0] {:epsilon 1e-4 :t-min 1.0 :t-max 10.0})
          warm (qaeb/intersect alt eye [0.0 0.0 -1.0] {:epsilon 1e-4 :t-max 10.0 :t-min 3.9})]
      (is (< (abs (- (:t cold) (:t warm))) 1e-3))
      (is (< (:steps warm) (:steps cold)))))

  (testing "but starting past the surface reports the start, which is the risk"
    ;; What a too-eager `:t-min` costs: the ray is already underground, the
    ;; first sample is negative, and the hit comes back at the start of the
    ;; march rather than at the surface. In a picture that is a smear along
    ;; the silhouette, where the distance to the ground changes fastest.
    (let [alt (sphere 1.0)
          {:keys [t]} (qaeb/intersect alt [0.0 0.0 5.0] [0.0 0.0 -1.0]
                                      {:epsilon 1e-4 :t-min 4.5 :t-max 5.5})]
      (is (== 4.5 t)))))

(deftest refinement-test
  (testing "bisection tightens the interpolated crossing"
    (let [alt (sphere 1.0)
          eye [0.0 0.0 5.0]
          err (fn [r] (abs (- 4.0 (:t (qaeb/intersect alt eye [0.0 0.0 -1.0]
                                                      {:epsilon 1e-2 :t-min 1.0 :t-max 10.0 :refine r})))))]
      (is (<= (err 6) (err 0))))))

(deftest shell-test
  (testing "sphere-entry brackets the ray's passage through a shell"
    (let [[near far] (qaeb/sphere-entry [0.0 0.0 5.0] [0.0 0.0 -1.0] 2.0)]
      (is (< (abs (- 3.0 near)) 1e-12))
      (is (< (abs (- 7.0 far)) 1e-12))))

  (testing "a miss is a miss"
    (is (nil? (qaeb/sphere-entry [0.0 2.0 5.0] [0.0 0.0 -1.0] 1.0))))

  (testing "from inside, the near root is behind you"
    (let [[near far] (qaeb/sphere-entry [0.0 0.0 0.0] [0.0 0.0 1.0] 2.0)]
      (is (neg? near))
      (is (< (abs (- 2.0 far)) 1e-12)))))

(deftest planet-test
  (let [pl (planet/planet {:seed 3 :samples 1024})
        alt #(planet/altitude pl %)
        eye [0.0 0.0 4.0]]
    (testing "every hit lands on the planet's own surface"
      (doseq [dx (range -0.2 0.21 0.05)]
        (when-let [{:keys [point]} (qaeb/intersect alt eye [dx 0.0 -1.0]
                                                   {:epsilon 5e-4 :t-min 2.5 :t-max 6.0})]
          (let [d (v3/normalize point)]
            (is (< (abs (- (v3/length point) (planet/surface-radius pl d))) 1e-3))))))

    (testing "and rays that clear the planet return nothing"
      (is (nil? (qaeb/intersect alt eye [1.0 0.0 -1.0] {:epsilon 1e-3 :t-min 1.0 :t-max 20.0}))))

    (testing "a point on the day side is not in its own shadow"
      (let [{:keys [point]} (qaeb/intersect alt eye [0.0 0.0 -1.0]
                                            {:epsilon 5e-4 :t-min 2.5 :t-max 6.0})]
        (is (not (qaeb/shadowed? alt point [0.0 0.0 1.0] {:t-max 3.0})))
        ;; And the light coming from behind the planet is blocked by it.
        (is (qaeb/shadowed? alt point [0.0 0.0 -1.0] {:t-max 3.0}))))))
