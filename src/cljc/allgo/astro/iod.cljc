(ns allgo.astro.iod
  "Initial orbit determination: an orbit from a few observations and no
  prior guess (Vallado, *Fundamentals of Astrodynamics and Applications*,
  chapter 7).

  A least-squares fit (`allgo.astro.estimation`) refines an orbit, but it
  must start from one. These make the first: from three position vectors
  (Gibbs's geometric construction, or Herrick-Gibbs's Taylor series when
  the three are too close together for geometry), from three sightings of
  direction alone (Gauss's method, which finds the ranges by solving an
  eighth-degree polynomial), or from two positions and the time between
  them -- Lambert's problem, which is also the targeting problem: what
  velocity takes me from here to there in this long.

  Kilometers, seconds, radians; `mu` defaults to the Earth's."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.universal :as universal]
            [allgo.geometry.vec3 :as v3]
            [allgo.math :as am]
            [clojure.math :as math]))

(def ^:private mu c/GM-earth)

(defn- angle [a b]
  (math/acos (am/clamp (/ (v3/dot a b) (* (v3/length a) (v3/length b))) -1.0 1.0)))

;; ---------------------------------------------------- three positions

(defn gibbs
  "The velocity at the middle of three positions on one orbit, by Gibbs's
  vector construction (Vallado's GIBBS): `{:v2 :theta12 :theta23 :copa}`,
  with the angles between the positions and how far the first is from the
  plane of the other two (zero for true coplanarity). Needs the positions
  well spread -- a few degrees apart at least."
  ([r1 r2 r3] (gibbs mu r1 r2 r3))
  ([mu r1 r2 r3]
   (let [m1 (v3/length r1) m2 (v3/length r2) m3 (v3/length r3)
         z12 (v3/cross r1 r2) z23 (v3/cross r2 r3) z31 (v3/cross r3 r1)
         n (v3/add (v3/add (v3/scale z23 m1) (v3/scale z31 m2)) (v3/scale z12 m3))
         d (v3/add (v3/add z12 z23) z31)
         ;; the differences of the magnitudes weight the vectors
         s (v3/add (v3/add (v3/scale r1 (- m2 m3)) (v3/scale r2 (- m3 m1)))
                   (v3/scale r3 (- m1 m2)))
         b (v3/cross d r2)
         lg (math/sqrt (/ mu (* (v3/length n) (v3/length d))))]
     {:v2 (v3/add (v3/scale b (/ lg m2)) (v3/scale s lg))
      :theta12 (angle r1 r2) :theta23 (angle r2 r3)
      :copa (math/asin (/ (v3/dot z23 r1) (* (v3/length z23) m1)))})))

(defn herrick-gibbs
  "The velocity at the middle of three closely spaced positions at times
  `t1` `t2` `t3` (seconds), by the Herrick-Gibbs Taylor expansion --
  the method for when the positions are too close for Gibbs's geometry,
  under a degree or so apart."
  ([r1 r2 r3 t1 t2 t3] (herrick-gibbs mu r1 r2 r3 t1 t2 t3))
  ([mu r1 r2 r3 t1 t2 t3]
   (let [t21 (- t2 t1) t31 (- t3 t1) t32 (- t3 t2)
         k (fn [r] (/ mu (* 12.0 (math/pow (v3/length r) 3.0))))]
     {:v2 (v3/add (v3/add (v3/scale r1 (* (- t32) (+ (/ 1.0 (* t21 t31)) (k r1))))
                          (v3/scale r2 (* (- t32 t21) (+ (/ 1.0 (* t21 t32)) (k r2)))))
                  (v3/scale r3 (* t21 (+ (/ 1.0 (* t32 t31)) (k r3)))))
      :theta12 (angle r1 r2) :theta23 (angle r2 r3)
      :copa (let [z23 (v3/cross r2 r3)]
              (math/asin (/ (v3/dot z23 r1) (* (v3/length z23) (v3/length r1)))))})))

;; ------------------------------------------------------- angles only

(defn line-of-sight
  "The unit vector toward right ascension `ra` and declination `dec`."
  [ra dec]
  [(* (math/cos dec) (math/cos ra)) (* (math/cos dec) (math/sin ra)) (math/sin dec)])

(defn- largest-positive-root
  "The largest positive root of x^8 + a x^6 + b x^3 + c, by Newton's
  method from well above it."
  [a b c]
  (let [f (fn [x] (+ (math/pow x 8) (* a (math/pow x 6)) (* b x x x) c))
        df (fn [x] (+ (* 8 (math/pow x 7)) (* 6 a (math/pow x 5)) (* 3 b x x)))]
    (loop [x (* 2.0 (math/pow (max (abs a) 1.0) 0.5)) i 0]
      (let [x' (- x (/ (f x) (df x)))]
        (if (or (< (abs (- x' x)) 1e-10) (> i 200)) x' (recur x' (inc i)))))))

(defn gauss
  "An orbit from three sightings of direction alone: right ascensions and
  declinations `[[ra dec] x3]` at times `ts` (seconds) from sites at
  `sites` (the observer's position when each was taken, in the same
  inertial frame). Gauss's method: the ranges come from the eighth-degree
  polynomial in the middle radius, the positions from the ranges, and the
  velocity from Gibbs -- or Herrick-Gibbs, when the arc is short. Returns
  `{:r2 :v2}` at the middle sighting.

  One pass, without refining the ranges; good to a percent or so over a
  short arc, which is what an initial orbit needs to be."
  ([observations ts sites] (gauss mu observations ts sites))
  ([mu observations [t1 t2 t3] [R1 R2 R3]]
   (let [[L1 L2 L3] (map (fn [[ra dec]] (line-of-sight ra dec)) observations)
         tau1 (- t1 t2) tau3 (- t3 t2) tau (- tau3 tau1)
         p1 (v3/cross L2 L3) p2 (v3/cross L1 L3) p3 (v3/cross L1 L2)
         D0 (v3/dot L1 p1)
         D (fn [R p] (v3/dot R p))
         D11 (D R1 p1) D12 (D R1 p2) D13 (D R1 p3)
         D21 (D R2 p1) D22 (D R2 p2) D23 (D R2 p3)
         D31 (D R3 p1) D32 (D R3 p2) D33 (D R3 p3)
         A (/ (+ (/ (* (- D12) tau3) tau) D22 (/ (* D32 tau1) tau)) D0)
         B (/ (+ (* D12 (- (* tau3 tau3) (* tau tau)) (/ tau3 tau))
                 (* D32 (- (* tau tau) (* tau1 tau1)) (/ tau1 tau)))
              (* 6.0 D0))
         E (v3/dot L2 R2)
         R2sq (v3/dot R2 R2)
         r2m (largest-positive-root (- (+ (* A A) (* 2.0 A E) R2sq))
                                    (* -2.0 mu B (+ A E))
                                    (- (* mu mu B B)))
         u (/ mu (* r2m r2m r2m))
         rho2 (+ A (* mu (/ B (* r2m r2m r2m))))
         rho1 (/ (- (/ (+ (* 6.0 (+ (* D31 (/ tau1 tau3)) (* D21 (/ tau tau3))) r2m r2m r2m)
                          (* mu D31 (- (* tau tau) (* tau1 tau1)) (/ tau1 tau3)))
                       (+ (* 6.0 r2m r2m r2m) (* mu (- (* tau tau) (* tau3 tau3)))))
                    D11)
                 D0)
         rho3 (/ (- (/ (+ (* 6.0 (- (* D13 (/ tau3 tau1)) (* D23 (/ tau tau1))) r2m r2m r2m)
                          (* mu D13 (- (* tau tau) (* tau3 tau3)) (/ tau3 tau1)))
                       (+ (* 6.0 r2m r2m r2m) (* mu (- (* tau tau) (* tau1 tau1)))))
                    D33)
                 D0)
         r1 (v3/add R1 (v3/scale L1 rho1))
         r2 (v3/add R2 (v3/scale L2 rho2))
         r3 (v3/add R3 (v3/scale L3 rho3))
         {:keys [theta12 theta23] :as g} (gibbs mu r1 r2 r3)
         close? (and (< theta12 (math/to-radians 1.0)) (< theta23 (math/to-radians 1.0)))]
     {:r2 r2 :v2 (:v2 (if close? (herrick-gibbs mu r1 r2 r3 t1 t2 t3) g))
      :ranges [rho1 rho2 rho3] :u u})))

;; ------------------------------------------------------------- Lambert

(defn- tof-of-psi
  "Lambert's time of flight as a function of the universal variable psi,
  with the geometry constant A, and the y(psi) it needs."
  [mu r1 r2 A psi]
  (let [[c2 c3] (universal/stumpff psi)
        y (+ r1 r2 (/ (* A (- (* psi c3) 1.0)) (math/sqrt c2)))]
    (when (and (pos? y) (pos? c2))
      (let [x (math/sqrt (/ y c2))]
        {:t (/ (+ (* x x x c3) (* A (math/sqrt y))) (math/sqrt mu)) :y y}))))

(defn- bisect-psi [f lo hi dt increasing?]
  (loop [lo lo hi hi i 0]
    (let [mid (* 0.5 (+ lo hi))
          t (:t (f mid))]
      (cond
        (or (> i 200) (< (- hi lo) 1e-12)) mid
        (nil? t) (if increasing? (recur mid hi (inc i)) (recur lo mid (inc i)))
        (= (< t dt) increasing?) (recur mid hi (inc i))
        :else (recur lo mid (inc i))))))

(defn- min-psi
  "The psi of least time of flight in the band of `revs` revolutions."
  [f lo hi]
  (let [g (/ (- (math/sqrt 5.0) 1.0) 2.0)
        t #(or (:t (f %)) ##Inf)]
    (loop [lo lo hi hi i 0]
      (if (> i 200)
        (* 0.5 (+ lo hi))
        (let [x1 (- hi (* g (- hi lo))) x2 (+ lo (* g (- hi lo)))]
          (if (< (t x1) (t x2)) (recur lo x2 (inc i)) (recur x1 hi (inc i))))))))

(defn lambert
  "The velocities `[v1 v2]` that carry a body from `r1` to `r2` in `dt`
  seconds -- Lambert's problem, by universal variables (Vallado's
  LAMBERTU). Options:

    :long?  go the long way round, through more than 180 degrees
    :revs   complete revolutions on the way (default 0)
    :high?  with revolutions there are two orbits for each time: the
            larger (high energy) or, by default, the smaller

  nil if no orbit fits -- too short a time for that many revolutions."
  ([r1 r2 dt] (lambert mu r1 r2 dt {}))
  ([r1 r2 dt opts] (lambert mu r1 r2 dt opts))
  ([mu r1 r2 dt {:keys [long? revs high?] :or {revs 0}}]
   (let [m1 (v3/length r1) m2 (v3/length r2)
         cdn (/ (v3/dot r1 r2) (* m1 m2))
         A (* (if long? -1.0 1.0) (math/sqrt (* m1 m2 (+ 1.0 cdn))))
         f (fn [psi] (tof-of-psi mu m1 m2 A psi))
         two-pi (* 2.0 math/PI)
         psi (if (zero? revs)
               ;; one branch: time rises monotonically with psi
               (bisect-psi f (* -4.0 math/PI math/PI) (* two-pi two-pi) dt true)
               ;; revs >= 1: the band between (2 pi N)^2 and (2 pi (N+1))^2,
               ;; with a minimum time inside it and an orbit either side
               (let [lo (math/pow (* two-pi revs) 2.0)
                     hi (math/pow (* two-pi (inc revs)) 2.0)
                     pm (min-psi f (+ lo 1e-9) (- hi 1e-9))]
                 (when (<= (:t (f pm)) dt)
                   (if high?
                     (bisect-psi f (+ lo 1e-9) pm dt false)
                     (bisect-psi f pm (- hi 1e-9) dt true)))))]
     (when psi
       (when-let [{:keys [t y]} (f psi)]
         (when (< (abs (- t dt)) (* 1e-6 (max 1.0 dt)))
           (let [fl (- 1.0 (/ y m1))
                 g (* A (math/sqrt (/ y mu)))
                 gd (- 1.0 (/ y m2))]
             [(v3/scale (v3/sub r2 (v3/scale r1 fl)) (/ 1.0 g))
              (v3/scale (v3/sub (v3/scale r2 gd) r1) (/ 1.0 g))])))))))

(defn lambert-min-energy
  "The transfer from `r1` to `r2` on the orbit of least energy -- the one
  whose semi-major axis is half the sum of the radii and the chord -- and
  how long it takes, the short way or the long, with `revs` extra
  revolutions: `{:v1 :a :tof :tof-parabolic}`, the last the fastest a
  parabola could make it, which no ellipse beats."
  ([r1 r2 long? revs] (lambert-min-energy mu r1 r2 long? revs))
  ([mu r1 r2 long? revs]
   (let [m1 (v3/length r1) m2 (v3/length r2)
         cdn (/ (v3/dot r1 r2) (* m1 m2))
         c (v3/length (v3/sub r2 r1))
         s (* 0.5 (+ m1 m2 c))
         a (* 0.5 s)
         p (/ (* m1 m2 (- 1.0 cdn)) c)
         beta (* 2.0 (math/asin (math/sqrt (/ (- s c) s))))
         sign (if long? 1.0 -1.0)
         sdn (* (if long? -1.0 1.0) (math/sqrt (- 1.0 (* cdn cdn))))
         v1 (v3/scale (v3/sub r2 (v3/scale r1 (- 1.0 (* (/ m2 p) (- 1.0 cdn)))))
                      (/ (math/sqrt (* mu p)) (* m1 m2 sdn)))]
     {:v1 v1 :a a
      :tof (* (math/sqrt (/ (* a a a) mu))
              (+ (* 2.0 math/PI revs) math/PI (* sign (- beta (math/sin beta)))))
      :tof-parabolic (* (/ 1.0 3.0) (math/sqrt (/ 2.0 mu))
                        (+ (math/pow s 1.5) (* sign (math/pow (- s c) 1.5))))})))
