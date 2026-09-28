(ns allgo.bplane-test
  "The B-plane against what it must be: unchanged along the hyperbola, B
  normal to S and to the reference pole's T only through R, the impact
  parameter consistent with the periapsis radius, periapsis-state its
  inverse; and targeting that lands on the B-plane point asked for."
  (:require [allgo.astro.bplane :as bp]
            [allgo.astro.universal :as u]
            [allgo.geometry.vec3 :as v3]
            [clojure.test :refer [deftest is testing]]))

(def ^:private mu-mars 42828.37)

(def ^:private arrival
  "A Mars arrival: v-infinity 3 km/s, B.T 8000 km, B.R -5000 km."
  (bp/periapsis-state mu-mars [1.2 -2.4 1.3] 8000.0 -5000.0))

(defn- inbound [dt] (let [[r v] arrival] (u/propagate mu-mars r v (- dt))))

(deftest round-trip
  (testing "periapsis-state gives the hyperbola whose B-plane is what was asked"
    (let [{:keys [bt br v-inf s]} (bp/b-plane mu-mars arrival)]
      (is (< (abs (- bt 8000.0)) 1e-8))
      (is (< (abs (- br -5000.0)) 1e-8))
      (is (< (abs (- v-inf (v3/length [1.2 -2.4 1.3]))) 1e-12))
      (is (< (v3/distance s (v3/normalize [1.2 -2.4 1.3])) 1e-12)))))

(deftest invariance
  (testing "every state on the hyperbola has the same B-plane"
    (let [ref (bp/b-plane mu-mars arrival)]
      (doseq [dt [3600.0 86400.0 (* 5 86400.0)]]
        (let [{:keys [bt br]} (bp/b-plane mu-mars (inbound dt))]
          (is (< (abs (- bt (:bt ref))) 1e-6) (str dt))
          (is (< (abs (- br (:br ref))) 1e-6) (str dt)))))))

(deftest geometry
  (let [{:keys [b s t r]} (bp/b-plane mu-mars (inbound 86400.0))]
    (testing "S, T, R orthonormal, T in the reference plane, B normal to S"
      (is (< (abs (v3/dot s t)) 1e-12))
      (is (< (abs (v3/dot s r)) 1e-12))
      (is (< (abs (nth t 2)) 1e-12))
      (is (< (abs (v3/dot b s)) 1e-8))))
  (testing "the impact parameter against the periapsis radius, b^2 = rp^2 (1 + 2 mu/(rp v^2))"
    (let [[rp] arrival
          {:keys [b v-inf]} (bp/b-plane mu-mars arrival)]
      (is (< (abs (- (bp/periapsis-radius mu-mars v-inf (v3/length b)) (v3/length rp))) 1e-8))))
  (testing "the time to periapsis is the time flown"
    (is (< (abs (- (bp/hyperbolic-time-to-periapsis mu-mars (inbound 86400.0)) 86400.0)) 1e-6))))

(deftest targeting
  (testing "a correction five days out moves the aim point by some thousands of kilometers, exactly"
    (let [s0 (inbound (* 5 86400.0))
          {:keys [dv bt br]} (bp/target mu-mars s0 (* 4 86400.0) 4000.0 2500.0)
          [r v] s0
          flown (u/propagate mu-mars r (v3/add v dv) (* 4 86400.0))
          check (bp/b-plane mu-mars flown)]
      (is (< (abs (- bt 4000.0)) 1e-6))
      (is (< (abs (- br 2500.0)) 1e-6))
      (is (< (abs (- (:bt check) 4000.0)) 1e-6))
      (is (< (abs (- (:br check) 2500.0)) 1e-6))
      ;; a few meters a second
      (is (< (v3/length dv) 0.05)))))
