(ns allgo.coordinate-time-test
  "TCG and TCB against ERFA (BSD-3-Clause), whose eraTttcg, eraTcgtt,
  eraTdbtcb and eraTcbtdb implement the same IAU resolutions
  independently. The values, seconds, are its output at each MJD, good
  to the 6e-7 s that a double MJD resolves."
  (:require [allgo.astro.time :as time]
            [clojure.test :refer [deftest is testing]]))

(def ^:private erfa
  ;; MJD, TCG - TT from TT, TT - TCG from TCG, TCB - TDB from TDB, TDB - TCB from TCB
  [[43144.0003725 0.0 0.0 6.537884473801e-05 -6.537884473801e-05]
   [44000.0 5.154367536306e-02 -5.154367536306e-02 1.146804657765e+00 -1.146804657765e+00]
   [51544.5 5.058329785243e-01 -5.058329785243e-01 1.125378690194e+01 -1.125378690194e+01]
   [55000.25 7.139200111851e-01 -7.139200111851e-01 1.588327924255e+01 -1.588327924255e+01]
   [60000.0 1.014978275634e+00 -1.014978275634e+00 2.258119038306e+01 -2.258118975442e+01]
   [61000.75 1.075238082558e+00 -1.075238082558e+00 2.392184394412e+01 -2.392184331547e+01]
   [73000.0 1.797768613324e+00 -1.797768613324e+00 3.999662850983e+01 -3.999662725255e+01]])

(defn- secs [a b] (* 86400.0 (- a b)))

(deftest against-erfa
  (doseq [[m tcg tt tcb tdb] erfa]
    (testing (str "MJD " m)
      (is (< (abs (- (secs (time/tt->tcg m) m) tcg)) 1e-6))
      (is (< (abs (- (secs (time/tcg->tt m) m) tt)) 1e-6))
      (is (< (abs (- (secs (time/tdb->tcb m) m) tcb)) 1e-6))
      (is (< (abs (- (secs (time/tcb->tdb m) m) tdb)) 1e-6)))))

(deftest inverses
  (testing "each undoes the other, to a unit or two in the last place of
            a double MJD -- 1.3e-6 s at MJD 80000"
    (doseq [m [40000.0 51544.5 60000.0 80000.0]]
      (is (< (abs (secs (time/tcg->tt (time/tt->tcg m)) m)) 3e-6))
      (is (< (abs (secs (time/tcb->tdb (time/tdb->tcb m)) m)) 3e-6)))))
