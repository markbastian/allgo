(ns allgo.geomagnetic-test
  "The Kp-ap table against GFZ's own record of both indices since 1932
  (Kp_ap_Ap_SN_F107_since_1932.txt, GFZ Helmholtz Centre for Geosciences,
  CC BY 4.0; Matzka et al. 2021, doi:10.5880/Kp.0001): in all 276,840
  three-hour intervals with both, each Kp step comes with exactly one ap,
  and these are they."
  (:require [allgo.astro.geomagnetic :as g]
            [clojure.test :refer [deftest is testing]]))

(def ^:private gfz
  "Kp in thirds, and the ap GFZ pairs with it every time it occurs."
  [[0 0] [1 2] [2 3] [3 4] [4 5] [5 6] [6 7] [7 9] [8 12] [9 15] [10 18] [11 22]
   [12 27] [13 32] [14 39] [15 48] [16 56] [17 67] [18 80] [19 94] [20 111]
   [21 132] [22 154] [23 179] [24 207] [25 236] [26 300] [27 400]])

(deftest the-steps
  (doseq [[thirds ap] gfz]
    (is (== ap (g/kp->ap (/ thirds 3.0))))
    (is (< (abs (- (g/ap->kp ap) (/ thirds 3.0))) 1e-12))))

(deftest between-steps
  (testing "monotonic, and each the other's inverse"
    (let [kps (range 0.0 9.0 0.0137)]
      (is (apply < (map g/kp->ap kps)))
      (is (every? #(< (abs (- (g/ap->kp (g/kp->ap %)) %)) 1e-12) kps))))
  (testing "clamped to the scales' ends"
    (is (== 400.0 (g/kp->ap 9.5)))
    (is (== 0.0 (g/ap->kp -3.0)))))
