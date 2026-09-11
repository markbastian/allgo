(ns allgo.vortex-test
  (:require [allgo.physics.fluid :as fluid]
            [allgo.physics.vortex :as vortex]
            [clojure.test :refer [deftest is testing]]))

(def ^:private dt (/ 1.0 60.0))

(defn- still [nx ny h]
  (-> (fluid/fluid nx ny h) (fluid/close-border!)))

(defn- uniform-flow [nx ny h speed]
  (let [f (still nx ny h)]
    (dotimes [i (:nx f)]
      (dotimes [j (:ny f)]
        (aset ^floats (:u f) (fluid/idx f i j) (float speed))))
    f))

(deftest pool-test
  (let [p (vortex/pool 3)]
    (testing "a pool starts empty"
      (is (zero? (vortex/vortex-count p))))

    (testing "vortices are added in order"
      (is (= 0 (vortex/add! p 0.1 0.2 5.0 1.0)))
      (is (= 1 (vortex/add! p 0.3 0.4 -5.0 0.5)))
      (is (= 2 (vortex/vortex-count p)))
      (let [v (vortex/vortex p 0)]
        (is (< (abs (- (:x v) 0.1)) 1e-6))
        (is (< (abs (- (:omega v) 5.0)) 1e-6))))

    (testing "a full pool refuses rather than writing past the end"
      ;; The reference leaves its capacity field undefined, so the guard
      ;; compares against nothing and always passes. The writes are then
      ;; dropped silently and the count still climbs, leaving slots that
      ;; read back as NaN.
      (is (some? (vortex/add! p 0.5 0.5 1.0 1.0)))
      (is (= 3 (vortex/vortex-count p)))
      (is (nil? (vortex/add! p 0.6 0.6 1.0 1.0)))
      (is (= 3 (vortex/vortex-count p))))))

(deftest expiry-test
  (testing "vortices age and are dropped when spent"
    ;; Lifetimes deliberately off the step boundary: one of exactly n
    ;; steps survives or not on the last bit of a float.
    (let [p (vortex/pool 8)]
      (vortex/add! p 0.1 0.1 1.0 (* 1.5 dt))
      (vortex/add! p 0.2 0.2 2.0 (* 10 dt))
      (vortex/add! p 0.3 0.3 3.0 (* 0.5 dt))
      (vortex/expire! p dt)
      (is (= 2 (vortex/vortex-count p)) "the shortest-lived one has gone")
      (vortex/expire! p dt)
      (is (= 1 (vortex/vortex-count p)) "and then the next")
      (testing "and the survivor is compacted to the front, intact"
        (let [v (vortex/vortex p 0)]
          (is (< (abs (- (:omega v) 2.0)) 1e-6))
          (is (< (abs (- (:x v) 0.2)) 1e-6))))))

  (testing "an expired vortex stops stirring"
    (let [f (still 20 20 0.1)
          p (vortex/pool 4)]
      (vortex/add! p 1.0 1.0 20.0 (* 0.5 dt))
      (vortex/step! p f dt {:radius 0.3})
      (is (zero? (vortex/vortex-count p)))
      (is (every? zero? (map #(aget ^floats (:u f) %) (range (:n f))))))))

(deftest rotation-test
  (testing "a vortex spins the velocity field around itself"
    (let [h 0.1
          f (still 20 20 h)
          p (vortex/pool 4)
          cx 1.0 cy 1.0]
      (vortex/add! p cx cy 20.0 1.0)
      (vortex/stir! p f dt {:radius 0.35 :damping 0.0})
      ;; Solid-body rotation: u is driven by how far above the centre the
      ;; face is, so it points opposite ways above and below.
      (let [u-above (aget ^floats (:u f) (fluid/idx f 10 12))
            u-below (aget ^floats (:u f) (fluid/idx f 10 7))
            v-right (aget ^floats (:v f) (fluid/idx f 12 10))
            v-left  (aget ^floats (:v f) (fluid/idx f 7 10))]
        (is (pos? u-above))
        (is (neg? u-below))
        (is (neg? v-right))
        (is (pos? v-left))
        (is (< (abs (+ u-above u-below)) 0.2) "and symmetrically so"))))

  (testing "reversing the spin reverses the flow"
    (let [spin (fn [omega]
                 (let [f (still 20 20 0.1)
                       p (vortex/pool 4)]
                   (vortex/add! p 1.0 1.0 omega 1.0)
                   (vortex/stir! p f dt {:radius 0.35 :damping 0.0})
                   (aget ^floats (:u f) (fluid/idx f 10 12))))]
      (is (pos? (spin 20.0)))
      (is (neg? (spin -20.0)))))

  (testing "nothing outside the radius is touched"
    (let [f (still 20 20 0.1)
          p (vortex/pool 4)]
      (vortex/add! p 1.0 1.0 20.0 1.0)
      (vortex/stir! p f dt {:radius 0.25 :damping 0.0})
      (is (zero? (aget ^floats (:u f) (fluid/idx f 2 2))))
      (is (zero? (aget ^floats (:u f) (fluid/idx f 18 18)))))))

(deftest blending-test
  (testing "a vortex nudges the flow it sits in rather than replacing it"
    ;; The reference assigns in the u branch and accumulates in the v one.
    ;; Assigning discards the velocity already in the cell and keeps only
    ;; the difference, scaled by the falloff -- so a spinner dropped into
    ;; a steady rightward flow reverses it at the core and erases it at
    ;; the rim, neither of which is a thing a vortex does.
    (let [speed 1.0
          f (uniform-flow 20 20 0.1 speed)
          p (vortex/pool 4)]
      ;; No spin at all: the flow should be very nearly untouched.
      (vortex/add! p 1.0 1.0 0.0 1.0)
      (vortex/stir! p f dt {:radius 0.35 :damping 0.0})
      (let [core (aget ^floats (:u f) (fluid/idx f 10 10))
            rim  (aget ^floats (:u f) (fluid/idx f 13 10))]
        (is (pos? core) "the flow still runs the way it was going")
        (is (< (abs (- core speed)) 1e-3) "and at the speed it was going")
        (is (< (abs (- rim speed)) 1e-3) "including at the rim"))))

  (testing "a spinning vortex still leaves the mean flow going the right way"
    (let [f (uniform-flow 20 20 0.1 1.0)
          p (vortex/pool 4)]
      (vortex/add! p 1.0 1.0 5.0 1.0)
      (vortex/stir! p f dt {:radius 0.35 :damping 0.0})
      (let [us (for [i (range 6 15) j (range 6 15)]
                 (aget ^floats (:u f) (fluid/idx f i j)))]
        (is (pos? (/ (reduce + us) (count us)))))))

  (testing "stirring leaves the field finite"
    (let [f (still 30 30 0.05)
          p (vortex/pool 32)]
      (dotimes [k 20] (vortex/add! p (+ 0.3 (* 0.05 k)) 0.7 (- 20.0 k) 1.0))
      (dotimes [_ 50] (vortex/step! p f dt {:radius 0.1}))
      (is (every? #(Float/isFinite (aget ^floats (:u f) %)) (range (:n f))))
      (is (every? #(Float/isFinite (aget ^floats (:v f) %)) (range (:n f)))))))
