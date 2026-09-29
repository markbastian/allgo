(ns allgo.sexagesimal-test
  (:require [allgo.astro.sexagesimal :as sx]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(defn- deg [x] (math/to-radians x))

(deftest exact-values
  (testing "values that are exact in both bases"
    (is (== (sx/dms->rad 30 30 0) (deg 30.5)))
    (is (== (sx/dms->rad 0 0 36) (deg 0.01)))
    (is (== (sx/hms->rad 1 0 0) (deg 15.0)))
    (is (== (sx/hms->rad 0 4 0) (deg 1.0)))
    (is (< (abs (- (sx/hms->rad 12 34 56) (deg 188.73333333333333))) 1e-15))
    (is (== (sx/hms->seconds 13 45 30.5) 49530.5))))

(deftest the-sign
  (testing "-0 degrees 30 minutes is half a degree south"
    (is (== (sx/dms->rad 0 -30 0) (deg -0.5)))
    (is (= {:sign -1 :d 0 :m 30 :s 0.0} (sx/rad->dms (deg -0.5)))))
  (testing "a minus anywhere makes the whole negative"
    (is (== (sx/dms->rad -12 30 0) (sx/dms->rad 12 -30 0) (deg -12.5)))))

(deftest carrying
  (testing "seconds that round to 60 carry into the minutes"
    (is (= {:sign 1 :d 11 :m 0 :s 0.0} (sx/rad->dms (deg (- 11.0 1e-12)))))
    (is (= {:h 24 :m 0 :s 0.0} (sx/seconds->hms (- 86400.0 1e-9))))))

(deftest round-trips
  (doseq [x (range -359.9 360.0 7.37)]
    (let [{:keys [sign d m s]} (sx/rad->dms (deg x) 9)
          {hs :sign :keys [h] hm :m hsec :s} (sx/rad->hms (deg x) 9)]
      (is (< (abs (- (* sign (sx/dms->rad d m s)) (deg x))) 1e-13))
      (is (< (abs (- (* hs (sx/hms->rad h hm hsec)) (deg x))) 1e-13))
      (is (<= 0 m 59)) (is (<= 0.0 s)) (is (< s 60.0))))
  (doseq [t (range 0.0 86400.0 1234.567)]
    (let [{:keys [h m s]} (sx/seconds->hms t)]
      (is (< (abs (- (sx/hms->seconds h m s) t)) 1e-6)))))
