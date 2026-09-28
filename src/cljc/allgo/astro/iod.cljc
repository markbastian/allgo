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
            [allgo.astro.kepler :as kepler]
            [allgo.astro.universal :as universal]
            [allgo.geometry.sphere :as sphere]
            [allgo.geometry.vec3 :as v3]
            [allgo.math :as am]
            [allgo.numerics.interpolation :as interp]
            [allgo.numerics.linear :as lin]
            [allgo.numerics.roots :as roots]
            [clojure.math :as math]))

(def ^:private mu c/GM-earth)

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
      :theta12 (v3/angle r1 r2) :theta23 (v3/angle r2 r3)
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
      :theta12 (v3/angle r1 r2) :theta23 (v3/angle r2 r3)
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

;; -------------------------------------------------------- Laplace

(defn laplace
  "An orbit from three sightings of direction alone, as `gauss` takes
  them, by Laplace's method: the line of sight L and the site R are
  differentiated at the middle sighting through the quadratic that passes
  their three values, and the equation of motion of r = rho L + R,

    rho'' L + 2 rho' L' + rho (L'' + mu L / r^3) + R'' + mu R / r^3 = 0,

  dotted with L x L' leaves the range, rho = -(R''.(LxL') + mu/r^3
  R.(LxL')) / L''.(LxL'), which with r^2 = rho^2 + 2 rho L.R + R^2 is an
  eighth-degree polynomial in r; dotted with L x L'' it leaves the range
  rate. Returns `{:r2 :v2}` at the middle sighting.

  The derivatives come from three points, so the answer is as good as a
  quadratic is over the arc: a percent or so, like Gauss's."
  ([observations ts sites] (laplace mu observations ts sites))
  ([mu observations [t1 t2 t3] [_ R :as sites]]
   (let [Ls (mapv (fn [[ra dec]] (line-of-sight ra dec)) observations)
         tau1 (- t1 t2) tau3 (- t3 t2)
         L (second Ls)
         [L' L''] (interp/derivatives-at-middle Ls tau1 tau3)
         [R' R''] (interp/derivatives-at-middle sites tau1 tau3)
         n1 (v3/cross L L') n2 (v3/cross L L'')
         D (v3/dot L'' n1)
         ;; rho = A + B / r^3
         A (- (/ (v3/dot R'' n1) D))
         B (- (/ (* mu (v3/dot R n1)) D))
         C (v3/dot L R)
         r (largest-positive-root (- (+ (* A A) (* 2.0 A C) (v3/dot R R)))
                                  (* -2.0 B (+ A C))
                                  (- (* B B)))
         r3 (* r r r)
         rho (+ A (/ B r3))
         rho' (- (/ (+ (v3/dot R'' n2) (* (/ mu r3) (v3/dot R n2)))
                    (* 2.0 (v3/dot L' n2))))]
     {:r2 (v3/add R (v3/scale L rho))
      :v2 (v3/add (v3/add (v3/scale L rho') (v3/scale L' rho)) R')})))

;; ------------------------------------------------------- double-r

(defn- double-r-step
  "For radii `r1m` `r2m` at the first two sightings: the three positions,
  the conic through them, and how far its timing misses the observation
  times -- `{:F [F1 F2] :r [r1 r2 r3] :conic ...}`, nil where the guess
  puts no conic through them."
  [mu [L1 L2 L3] [tau1 tau3] [R1 R2 R3] r1m r2m]
  (when-let [rho1 (sphere/far-intersection R1 L1 r1m)]
    (when-let [rho2 (sphere/far-intersection R2 L2 r2m)]
      (let [r1 (v3/add R1 (v3/scale L1 rho1))
            r2 (v3/add R2 (v3/scale L2 rho2))
            W (v3/normalize (v3/cross r1 r2))
            ;; the third position is where its line of sight meets the plane
            rho3 (- (/ (v3/dot R3 W) (v3/dot L3 W)))
            r3 (v3/add R3 (v3/scale L3 rho3))
            r3m (v3/length r3)
            ;; angles swept, signed about W
            sweep (fn [a b] (math/atan2 (v3/dot (v3/cross a b) W) (v3/dot a b)))
            dnu21 (sweep r1 r2) dnu32 (sweep r2 r3)
            ;; p, e cos nu2, e sin nu2 from p/r - 1 = e cos nu at each point
            [p X Y] (lin/solve-3 [[(/ 1.0 r1m) (- (math/cos dnu21)) (- (math/sin dnu21))]
                                  [(/ 1.0 r2m) -1.0 0.0]
                                  [(/ 1.0 r3m) (- (math/cos dnu32)) (math/sin dnu32)]]
                                 [1.0 1.0 1.0])
            e (math/hypot X Y)
            a (/ p (- 1.0 (* e e)))
            nu2 (math/atan2 Y X)
            mean (fn [nu] (second (kepler/true->anomaly-and-mean e nu)))
            anomaly (fn [nu] (first (kepler/true->anomaly-and-mean e nu)))
            n (math/sqrt (/ mu (abs (* a a a))))
            dM (fn [nu] (if (< e 1.0) (am/wrap-angle (- (mean nu) (mean nu2))) (- (mean nu) (mean nu2))))]
        (when (and (pos? p) (not= e 1.0) (pos? rho3))
          {:F [(- tau1 (/ (dM (- nu2 dnu21)) n)) (- tau3 (/ (dM (+ nu2 dnu32)) n))]
           :r [r1 r2 r3] :a a :e e
           :dx (let [d (- (anomaly (+ nu2 dnu32)) (anomaly nu2))]
                 (if (< e 1.0) (am/wrap-angle d) d))})))))

(defn double-r
  "An orbit from three sightings of direction alone, as `gauss` takes
  them, by Escobal's double-r iteration: guess the distances from the
  center at the first two sightings, which places all three positions and
  a conic through them, and adjust the two by Newton's method until that
  conic's timing matches the observations. Exact, where Gauss's method
  truncates a series, and so the way to refine it: the guesses default to
  `gauss`'s. Returns `{:r2 :v2}` at the middle sighting, nil if the
  iteration does not converge."
  ([observations ts sites] (double-r mu observations ts sites nil))
  ([mu observations [t1 t2 t3 :as ts] sites guesses]
   (let [Ls (map (fn [[ra dec]] (line-of-sight ra dec)) observations)
         taus [(- t1 t2) (- t3 t2)]
         [g1 g2] (or guesses
                     (let [{[rho1 rho2] :ranges} (gauss mu observations ts sites)]
                       [(v3/length (v3/add (first sites) (v3/scale (first Ls) rho1)))
                        (v3/length (v3/add (second sites) (v3/scale (second Ls) rho2)))]))
         step (fn [r1m r2m] (double-r-step mu Ls taus sites r1m r2m))
         scale (max (abs (first taus)) (abs (second taus)))]
     (loop [r1m g1 r2m g2 i 0]
       (when-let [{[F1 F2] :F :as here} (step r1m r2m)]
         (if (< (max (abs F1) (abs F2)) (* 1e-12 scale))
           (let [{[_ r2 r3] :r :keys [a dx]} here
                 tau3 (second taus)
                 r2m' (v3/length r2)
                 ;; f and g from the middle position to the third
                 [f g] (if (pos? a)
                         [(- 1.0 (* (/ a r2m') (- 1.0 (math/cos dx))))
                          (- tau3 (* (math/sqrt (/ (* a a a) mu)) (- dx (math/sin dx))))]
                         [(- 1.0 (* (/ a r2m') (- 1.0 (math/cosh dx))))
                          (- tau3 (* (math/sqrt (/ (- (* a a a)) mu)) (- (math/sinh dx) dx)))])]
             {:r2 r2 :v2 (v3/scale (v3/sub r3 (v3/scale r2 f)) (/ 1.0 g))})
           (when (< i 50)
             ;; Newton, with the partials by central differences
             (let [h1 (* 1e-6 r1m) h2 (* 1e-6 r2m)
                   d (fn [a b] (map #(/ (- %1 %2) (* 2.0 %3)) (:F a) (:F b) (repeat 1.0)))
                   p1 (step (+ r1m h1) r2m) m1 (step (- r1m h1) r2m)
                   p2 (step r1m (+ r2m h2)) m2 (step r1m (- r2m h2))]
               (when (and p1 m1 p2 m2)
                 (let [[a11 a21] (map #(/ % h1) (d p1 m1))
                       [a12 a22] (map #(/ % h2) (d p2 m2))
                       det (- (* a11 a22) (* a12 a21))
                       dr1 (/ (- (* a22 F1) (* a12 F2)) det)
                       dr2 (/ (- (* a11 F2) (* a21 F1)) det)]
                   (recur (- r1m dr1) (- r2m dr2) (inc i))))))))))))

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
  (roots/minimize #(or (:t (f %)) ##Inf) lo hi {:tol 0.0 :max-iter 200}))

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

;; ------------------------------------------- Lambert by Lagrange's equation

(defn- lagrange-branches
  "Lagrange's time of flight as a function of the semi-major axis, for
  each branch the orbits through the two points divide into: `[[time-fn
  p-fn] ...]` over a > s/2 for ellipses -- the upper and lower branches,
  alpha and 2 pi - alpha -- and over a < 0 for hyperbolas."
  [mu s c m1 m2 long? revs]
  (let [sgn (if long? -1.0 1.0)
        pk (/ (* 4.0 (- s m1) (- s m2)) (* c c))
        ellipse (fn [upper?]
                  (let [angles (fn [a]
                                 (let [a0 (* 2.0 (math/asin (math/sqrt (min 1.0 (/ s (* 2.0 a))))))
                                       b0 (* 2.0 (math/asin (math/sqrt (min 1.0 (/ (- s c) (* 2.0 a))))))]
                                   [(if upper? (- (* 2.0 math/PI) a0) a0) (* sgn b0)]))]
                    [(fn [a] (let [[al be] (angles a)]
                               (* (math/sqrt (/ (* a a a) mu))
                                  (- (+ (* 2.0 math/PI revs) (- al (math/sin al))) (- be (math/sin be))))))
                     (fn [a] (let [[al be] (angles a)]
                               (* pk a (math/pow (math/sin (* 0.5 (+ al be))) 2))))]))
        hyperbola (let [angles (fn [a]
                                 [(* 2.0 (am/asinh (math/sqrt (/ s (* -2.0 a)))))
                                  (* sgn 2.0 (am/asinh (math/sqrt (/ (- s c) (* -2.0 a)))))])]
                    [(fn [a] (let [[g d] (angles a)]
                               (* (math/sqrt (/ (- (* a a a)) mu))
                                  (- (- (math/sinh g) g) (- (math/sinh d) d)))))
                     (fn [a] (let [[g d] (angles a)]
                               (* (- pk) a (math/pow (math/sinh (* 0.5 (+ g d))) 2))))])]
    (concat [[:ellipse (ellipse false)] [:ellipse (ellipse true)]]
            (when (zero? revs) [[:hyperbola hyperbola]]))))

(defn lambert-lagrange
  "Every transfer from `r1` to `r2` in `dt` seconds, the short way or
  with `:long?` the long, with `:revs` complete revolutions: `[[v1 v2]
  ...]`, found from Lagrange's equation for the time of flight,

    sqrt(mu) t = a^3/2 [2 pi k + (alpha - sin alpha) -+ (beta - sin beta)],
    sin^2(alpha/2) = s / 2a,  sin^2(beta/2) = (s - c) / 2a,

  s the semiperimeter and c the chord (with sinh for hyperbolas), solved
  for the semi-major axis on each branch; the semilatus rectum, p = 4a
  (s - r1)(s - r2) / c^2 sin^2((alpha + beta) / 2), then gives the
  velocities through f and g. Independent of the universal-variable
  `lambert`, and a check on it: with revolutions it finds both orbits
  for each time, where `lambert` gives the one asked for."
  ([r1 r2 dt] (lambert-lagrange mu r1 r2 dt {}))
  ([r1 r2 dt opts] (lambert-lagrange mu r1 r2 dt opts))
  ([mu r1 r2 dt {:keys [long? revs] :or {revs 0}}]
   (let [m1 (v3/length r1) m2 (v3/length r2)
         c (v3/distance r1 r2)
         s (* 0.5 (+ m1 m2 c))
         cdn (/ (v3/dot r1 r2) (* m1 m2))
         sdn (* (if long? -1.0 1.0) (math/sqrt (max 0.0 (- 1.0 (* cdn cdn)))))
         velocities (fn [p]
                      (let [f (- 1.0 (* (/ m2 p) (- 1.0 cdn)))
                            g (/ (* m1 m2 sdn) (math/sqrt (* mu p)))
                            gd (- 1.0 (* (/ m1 p) (- 1.0 cdn)))]
                        [(v3/scale (v3/sub r2 (v3/scale r1 f)) (/ 1.0 g))
                         (v3/scale (v3/sub (v3/scale r2 gd) r1) (/ 1.0 g))]))
         ;; a on a log scale: ellipses from s/2 out; hyperbolas over every
         ;; size, their time falling to nothing as a does
         grid (fn [kind] (let [base (* 0.5 s)]
                           (if (= kind :ellipse)
                             (map #(* base (math/pow 10.0 (/ % 100.0))) (range 0 601))
                             (map #(- (* base (math/pow 10.0 (/ % 100.0)))) (range 600 -601 -1)))))]
     (vec
      (for [[kind [tf pf]] (lagrange-branches mu s c m1 m2 long? revs)
            :let [f #(- (tf %) dt) as (grid kind)]
            [lo hi] (map vector as (rest as))
            :when (not= (neg? (f lo)) (neg? (f hi)))
            :let [a (roots/bisect f lo hi)]]
        (velocities (pf a)))))))
