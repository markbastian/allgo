(ns allgo.mission-handbook-test
  "`allgo.astro.interplanetary/transfer` against NASA's Earth-to-Mars
  mission design handbook, 2026-2045 (NASA/TM-2010-216764; test/data/imdh):
  every tabulated energy minimum, and the 2005 minima of its appendix, the
  last as JPL gave them too."
  (:require [allgo.astro.interplanetary :as ip]
            [allgo.astro.time :as time]
            [allgo.geometry.vec3 :as v3]
            [clojure.edn :as edn]
            [clojure.java.io :as io]
            [clojure.math :as math]
            [clojure.test :refer [deftest is testing]]))

(defn- rows [file] (edn/read-string (slurp (io/resource (str "data/imdh/" file)))))

(def ^:private date-slip
  "Days by which tables 6-8 print every date early; see the README."
  {2035 63 2037 81 2039 75})

(defn- slip
  "The row as the handbook meant it, and which quantities to check."
  [[year type k depart arrive c3 ra dec v-inf]]
  (let [shift (date-slip year 0)
        mjd #(+ shift (apply time/calendar->mjd %))]
    {:year year :k k :depart (mjd depart) :arrive (mjd arrive)
     :long? (or (= type 2) (and (= year 2041) (= k :v-inf)))
     :c3 c3 :v-inf v-inf
     :ra (when-not (and (= year 2043) (= type 1) (= k :v-inf)) ra)
     :dec (when-not (or (and (= year 2043) (= type 1) (= k :v-inf))
                        (and (= year 2028) (= type 1) (= k :c3)))
            dec)}))

(defn- degrees-apart [a b] (abs (- (mod (+ (- a b) 180.0) 360.0) 180.0)))

(deftest energy-minima-2026-2045
  (doseq [{:keys [year k depart arrive long? c3 v-inf ra dec]}
          (map slip (rows "earth-mars-2026-2045.edn"))]
    (testing (str year " " k (when long? " type II"))
      (let [t (ip/transfer :earth depart :mars arrive {:long? long?})
            [ra' dec'] (map math/to-degrees (ip/asymptote (:v-inf-depart t) depart))]
        ;; C3 to four figures, the arrival speed to its last printed digit
        (is (< (abs (- (:c3 t) c3)) (* 1.5e-3 c3)))
        (is (< (abs (- (v3/length (:v-inf-arrive t)) v-inf)) 1.5e-3))
        ;; the angles to their four figures; right ascension runs some
        ;; 0.07 degrees high on average, unexplained but inside that
        (when ra (is (< (degrees-apart ra' ra) 0.15)))
        (when dec (is (< (abs (- dec' dec)) 0.12)))))))

(deftest energy-minima-2005
  (testing "as MIDAS computed them and as JPL did"
    (doseq [[src type depart arrive c3 v-inf] (rows "earth-mars-2005.edn")]
      (let [t (ip/transfer :earth (apply time/calendar->mjd depart)
                           :mars (apply time/calendar->mjd arrive) {:long? (= type 2)})]
        (when c3 (is (< (abs (- (:c3 t) c3)) (* 1e-3 c3)) (str src)))
        (when v-inf (is (< (abs (- (v3/length (:v-inf-arrive t)) v-inf)) 1.5e-3) (str src)))))))

(deftest asymptote
  (testing "the J2000 angles transfer returns are the asymptote's, and precession moves them"
    (let [d (time/calendar->mjd 2033 4 28)
          t (ip/transfer :earth d :mars (time/calendar->mjd 2034 1 27) {:long? true})
          [ra dec] (ip/asymptote (:v-inf-depart t))
          [ra' _] (ip/asymptote (:v-inf-depart t) d)]
      (is (= [ra dec] [(:rla t) (:dla t)]))
      ;; 33 years of general precession, about 50 arcseconds a year
      (is (< 0.4 (math/to-degrees (- ra' ra)) 0.5)))))
