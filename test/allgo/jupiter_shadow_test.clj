(ns allgo.jupiter-shadow-test
  "The Galilean satellites' shadows on Jupiter.

  There is no published table of transit times here to check against, so
  these test what the times have to satisfy: they recur at each moon's
  synodic period, their ends are where the shadow meets the limb, they
  last no longer than a shadow crossing the planet's full width would,
  and Callisto's skips the years either side of Jupiter's solstice. The
  positions underneath are checked against JPL Horizons in
  `allgo.meeus-independent-test`."
  (:require [allgo.astro.almanac :as al]
            [allgo.astro.jupiter-moons :as jm]
            [clojure.test :refer [deftest is testing]]))

(def ^:private start 61335.0)  ; 2026 October 22

(def ^:private index (zipmap jm/names (range)))

(defn- transits [moon from n]
  (take n (rest (iterate #(jm/shadow-transit moon (:mid %)) {:mid from}))))

(deftest recurrence-test
  (testing "each moon's shadow crosses once a synodic period"
    (doseq [moon [:io :europa :ganymede]]
      (let [mids (map :mid (transits moon start 12))
            p (jm/synodic-periods moon)]
        (doseq [d (map - (rest mids) mids)]
          ;; The Laplace resonance pulls each a few minutes either way.
          (is (< (abs (- d p)) 0.01) (str moon " " d " against " p)))))))

(deftest limb-test
  (testing "ingress and egress are where the shadow meets the limb"
    (let [geometry #'jm/geometry]
      (doseq [moon jm/names
              {:keys [ingress egress mid]} (transits moon start 3)]
        (is (< ingress mid egress))
        (is (< (:miss (geometry moon mid)) 1.0))
        (is (< (abs (- 1.0 (:miss (geometry moon ingress)))) 1e-3))
        (is (< (abs (- 1.0 (:miss (geometry moon egress)))) 1e-3))))))

(deftest duration-test
  (testing "no shadow takes longer to cross than the planet's full width allows"
    ;; The shadow moves across the planet as fast as the moon moves
    ;; across the line to the Sun -- measured at mid-transit, since the
    ;; orbits are eccentric enough for the mean period to be a percent or
    ;; two off. On an orbit of r radii, the moon is within a radius of that
    ;; line for an arc of 2 asin(1/r), which at that speed is the widest
    ;; crossing there is.
    (doseq [moon jm/names
            {:keys [ingress egress mid]} (transits moon start 3)]
      (let [m (nth (jm/positions-3d mid) (index moon))
            r (Math/sqrt (reduce + (map * m m)))
            across (fn [t] (let [m (nth (jm/positions-3d t) (index moon))
                                 s (jm/sun-from-jupiter t)
                                 k (reduce + (map * m s))]
                             (mapv - m (map #(* k %) s))))
            h 0.001
            speed (/ (Math/sqrt (reduce + (map #(* % %) (map - (across (+ mid h)) (across (- mid h))))))
                     (* 2.0 h))
            longest (* 1.001 (/ (* 2.0 r (Math/asin (/ 1.0 r))) speed))]
        (is (< (* 0.8 longest) (- egress ingress) longest)
            (str moon " " (- egress ingress) " of " longest))))))

(deftest callisto-test
  (testing "Callisto's shadow misses for years at a time; the others' do not"
    ;; 2005 and 2006 were Jupiter's southern solstice: the Sun stood far
    ;; enough below its equator that Callisto's shadow passed beneath.
    (let [from 53371.0 to 54101.0] ; 2005 January 1 to 2007 January 1
      (is (> (:mid (jm/shadow-transit :callisto from)) to))
      (is (< (:mid (jm/shadow-transit :io from)) (+ from 2.0))))))

(deftest almanac-test
  (testing "the almanac steps from one transit to the next and back"
    (let [a (al/next-event :europa-shadow start)
          b (al/next-event :europa-shadow a)]
      (is (< (abs (- (- b a) (jm/synodic-periods :europa))) 0.01))
      (is (< (abs (- a (al/previous-event :europa-shadow b))) 1e-4)))))
