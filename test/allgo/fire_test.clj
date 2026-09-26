(ns allgo.fire-test
  (:require [allgo.physics.fire :as fire]
            [allgo.physics.fluid :as fluid]
            [allgo.physics.vortex :as vortex]
            [clojure.test :refer [deftest is testing]]))

(def ^:private dt (/ 1.0 60.0))

(defn- scene
  ([] (scene {}))
  ([opts] (fire/fire (merge {:nx 60 :ny 90 :h 0.01 :vortices 200} opts))))

(def ^:private disk {:kind :disk :x 0.3 :y 0.25 :radius 0.06})

(defn- highest-hot-row
  "The topmost row holding anything warmer than `t`."
  [{:keys [nx ny ^floats smoke]} t]
  (reduce max 0 (for [i (range nx) j (range ny)
                      :when (> (aget smoke (+ (* i ny) j)) t)]
                  j)))

(defn- seeded
  "A repeatable stand-in for `rand`.

  Sources spawn vortices at random points, and vortices stir the field, so
  a scene left on `rand` gives a slightly different plume every run -- and
  a test with a threshold on it fails one run in ten."
  []
  (let [rng (java.util.Random. 20260911)]
    (fn [] (.nextDouble rng))))

(defn- advance [f steps world]
  (let [world (merge {:rng (seeded)} world)]
    (dotimes [_ steps] (fire/step! f world)))
  f)

(deftest setup-test
  (let [f (scene)]
    (testing "a fire starts cold"
      (is (zero? (fire/max-temperature f)))
      (is (zero? (fire/total-heat f))))

    (testing "the domain is left open"
      ;; A solid ceiling over a wide source makes a rising column
      ;; impossible -- there is nowhere for the fluid to go, so the
      ;; projection cancels the buoyancy and the heat pools on the floor.
      (is (not (fluid/solid? f 0 0)))
      (is (not (fluid/solid? f (dec (:nx f)) (dec (:ny f))))))

    (testing "it carries a pool of vortices"
      (is (zero? (vortex/vortex-count (:vortices f)))))))

(deftest temperature-test
  (testing "temperature is the fluid's own scalar field"
    (let [f (scene)]
      (fire/set-temperature! f 5 5 0.75)
      (is (< (abs (- (fire/temperature-at f 5 5) 0.75)) 1e-6))
      (is (< (abs (- (fluid/smoke-at f 5 5) 0.75)) 1e-6))))

  (testing "heat decays when nothing is burning"
    (let [f (scene)]
      (fire/set-temperature! f 30 20 1.0)
      (let [before (fire/total-heat f)]
        (advance f 120 {:emitters []})
        (is (< (fire/total-heat f) before))
        (is (< (fire/max-temperature f) 0.5)))))

  (testing "everything goes out eventually"
    (let [f (scene)]
      (dotimes [i 20] (fire/set-temperature! f (+ 20 i) 10 1.0))
      (advance f 400 {:emitters []})
      (is (zero? (fire/max-temperature f)))))

  (testing "flame burns out faster than the smoke it leaves"
    ;; Two cooling rates are what put a boundary between fire and smoke
    ;; instead of a plume that just dims all over.
    (let [drop-from (fn [t0]
                      (let [f (scene)]
                        (fire/set-temperature! f 30 20 t0)
                        (fire/buoyancy-and-cooling! f dt {})
                        (- t0 (fire/temperature-at f 30 20))))]
      (is (> (drop-from 0.9) (drop-from 0.1))
          "a flame-hot cell loses more per step than a smoke-cool one"))))

(deftest buoyancy-test
  (testing "hot fluid is pulled upward"
    (let [f (scene)]
      (fire/set-temperature! f 30 20 1.0)
      (fire/buoyancy-and-cooling! f dt {})
      (let [[_ v] (fluid/velocity-at f 30 20)]
        (is (pos? v)))))

  (testing "cold fluid is not"
    (let [f (scene)]
      (fire/buoyancy-and-cooling! f dt {})
      (let [[_ v] (fluid/velocity-at f 30 20)]
        (is (zero? v)))))

  (testing "hotter is faster"
    (let [lift (fn [t]
                 (let [f (scene)]
                   (fire/set-temperature! f 30 20 t)
                   (fire/buoyancy-and-cooling! f dt {})
                   (second (fluid/velocity-at f 30 20))))]
      (is (> (lift 1.0) (lift 0.5) (lift 0.1)))))

  (testing "buoyancy is the only thing lifting a fire"
    ;; There is no gravity in a fire scene, so with the lift turned off
    ;; nothing should rise at all.
    (let [f (advance (scene) 120 {:emitters [disk] :lift 0.0})
          lifted (advance (scene) 120 {:emitters [disk]})]
      (is (< (highest-hot-row f 0.05) (highest-hot-row lifted 0.05))
          "without lift the heat stays far below where buoyancy carries it"))))

(deftest plume-test
  (testing "a fire rises"
    (let [f (scene)
          source-row (do (fire/step! f {:emitters [disk]})
                         (highest-hot-row f 0.05))]
      (advance f 300 {:emitters [disk]})
      (is (> (highest-hot-row f 0.05) (+ source-row 20))
          "heat has climbed well above where it was made")))

  (testing "the plume stays inside the grid and stays finite"
    (let [f (advance (scene) 300 {:emitters [disk]})
          {:keys [^floats u ^floats v n]} f]
      (is (every? #(Float/isFinite (aget u %)) (range n)))
      (is (every? #(Float/isFinite (aget v %)) (range n)))
      (is (<= (fire/max-temperature f) 1.0))))

  (testing "a source keeps making heat and a settled fire does not run away"
    (let [f (advance (scene) 200 {:emitters [disk]})
          early (fire/total-heat f)]
      (advance f 200 {:emitters [disk]})
      (is (pos? (fire/total-heat f)))
      (is (< (fire/total-heat f) (* 3 early))
          "heat settles to a steady plume rather than accumulating")))

  (testing "sources spawn the vortices that give a plume its curl"
    (let [f (advance (scene) 120 {:emitters [disk]})]
      (is (pos? (vortex/vortex-count (:vortices f))))))

  (testing "no source, no fire"
    (let [f (advance (scene) 120 {:emitters []})]
      (is (zero? (fire/total-heat f)))
      (is (zero? (vortex/vortex-count (:vortices f)))))))

(deftest emitter-test
  (testing "a disk heats a ring and a floor heats its bottom rows"
    (let [f (scene)]
      (fire/ignite! f dt {:emitters [{:kind :floor :rows 3}]})
      (is (< (abs (- (fire/temperature-at f 30 2) 1.0)) 1e-6))
      (is (zero? (fire/temperature-at f 30 40)))))

  (testing "a floor source is not also a wind"
    (let [f (scene)]
      (fluid/set-velocity! f 30 2 5.0 5.0)
      (fire/ignite! f dt {:emitters [{:kind :floor :rows 3}]})
      (is (= [0.0 0.0] (fluid/velocity-at f 30 2)))))

  (testing "several sources burn at once"
    (let [f (scene)]
      (fire/ignite! f dt {:emitters [{:kind :floor :rows 2}
                                     {:kind :disk :x 0.3 :y 0.5 :radius 0.06}]})
      (is (pos? (fire/temperature-at f 30 1)))
      (is (pos? (fire/max-temperature f)))
      (is (> (highest-hot-row f 0.5) 20) "the disk is alight too")))

  (testing "a new kind of source is a new method, not a change to the step"
    (defmethod fire/emit! ::stripe [f _ _ _]
      (doseq [j (range 10 14)] (fire/set-temperature! f 20 j 1.0))
      [[0.2 0.11]])
    (let [f (scene)]
      (fire/ignite! f dt {:emitters [{:kind ::stripe}]})
      (is (< (abs (- (fire/temperature-at f 20 12) 1.0)) 1e-6))
      (is (= 1 (vortex/vortex-count (:vortices f))))
      (remove-method fire/emit! ::stripe)))

  (testing "an unknown source is an error, not a silent nothing"
    (is (thrown? clojure.lang.ExceptionInfo
                 (fire/ignite! (scene) dt {:emitters [{:kind :nonsense}]})))))
