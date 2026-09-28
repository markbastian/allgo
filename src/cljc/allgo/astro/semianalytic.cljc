(ns allgo.astro.semianalytic
  "A semi-analytic orbit theory in the manner of the Draper Semi-analytic
  Satellite Theory (DSST; Cefola, McClain and others -- McClain, \"A
  recursively formulated first-order semianalytic artificial satellite
  theory based on the generalized method of averaging\", CSC/TR-77/6010,
  1977; Vallado, chapter 9): the mean elements carried by averaged
  equations of motion at steps of hours, and the short-periodic motion
  put back where an osculating state is wanted.

  Here the averaging, which DSST does in closed form for the zonal
  harmonics and by quadrature for drag and the third body, is done by
  quadrature for every force alike: Gauss's variational equations in the
  equinoctial elements [a af ag chi psi lambda] -- nonsingular for
  circular and equatorial orbits -- evaluated round the orbit at the mean
  elements and averaged over the mean longitude. The first-order
  short-periodic variations are the integrals over the mean longitude of
  each rate's departure from its mean, taken term by term in its Fourier
  series, the mean longitude's own including the change of mean motion
  the variation in a brings. The mean equations carry the second-order
  averages of the generalized method of averaging too -- the average of
  the rates' partials times the first-order variations, and the mean
  motion's curvature in a -- which DSST includes for J2 as its J2^2
  terms, and without which a first-order theory drifts along track by
  kilometers in days.

  `accel` is the perturbing acceleration alone, `(fn [t r v])`, km/s^2;
  times seconds."
  (:require [allgo.astro.constants :as c]
            [allgo.math :as am]
            [allgo.numerics.differentiation :as diff]
            [allgo.numerics.linear :as lin]
            [clojure.math :as math]))

(def ^:private mu c/GM-earth)

;; ------------------------------------------------ equinoctial elements

(defn- solve-kepler
  "The eccentric longitude F of mean longitude `lam`: lam = F + ag cos F - af sin F."
  [af ag lam]
  (loop [F lam k 0]
    (let [d (/ (- (+ (- F (* af (math/sin F))) (* ag (math/cos F))) lam)
               (- 1.0 (* af (math/cos F)) (* ag (math/sin F))))]
      (if (or (< (abs d) 1e-15) (> k 50)) (- F d) (recur (- F d) (inc k))))))

(defn elements->state
  "The state `[r v]` of equinoctial elements [a af ag chi psi lambda]
  (prograde: chi = tan(i/2) sin raan, psi = tan(i/2) cos raan; af, ag the
  eccentricity vector along the equinoctial axes; lambda the mean
  longitude)."
  ([x] (elements->state mu x))
  ([mu [a af ag chi psi lam]]
   (let [F (solve-kepler af ag lam)
         b (/ 1.0 (+ 1.0 (math/sqrt (- 1.0 (* af af) (* ag ag)))))
         n (math/sqrt (/ mu (* a a a)))
         cf (math/cos F) sf (math/sin F)
         r (* a (- 1.0 (* af cf) (* ag sf)))
         X (* a (+ (* (- 1.0 (* ag ag b)) cf) (* af ag b sf) (- af)))
         Y (* a (+ (* (- 1.0 (* af af b)) sf) (* af ag b cf) (- ag)))
         Xd (/ (* n a a (- (* af ag b cf) (* (- 1.0 (* ag ag b)) sf))) r)
         Yd (/ (* n a a (- (* (- 1.0 (* af af b)) cf) (* af ag b sf))) r)
         s2 (+ 1.0 (* chi chi) (* psi psi))
         f [(/ (+ (- 1.0 (* chi chi)) (* psi psi)) s2) (/ (* 2.0 chi psi) s2) (/ (* -2.0 chi) s2)]
         g [(/ (* 2.0 chi psi) s2) (/ (- (+ 1.0 (* chi chi)) (* psi psi)) s2) (/ (* 2.0 psi) s2)]
         comb (fn [u w] (mapv #(+ (* u %1) (* w %2)) f g))]
     [(comb X Y) (comb Xd Yd)])))

(defn state->elements
  "[a af ag chi psi lambda] of the state `[r v]`, the inverse of
  `elements->state`."
  ([s] (state->elements mu s))
  ([mu [r v]]
   (let [rm (lin/length r)
         h (let [[x y z] r [vx vy vz] v] [(- (* y vz) (* z vy)) (- (* z vx) (* x vz)) (- (* x vy) (* y vx))])
         hm (lin/length h)
         [wx wy wz] (lin/scale h (/ 1.0 hm))
         chi (/ wx (+ 1.0 wz)) psi (/ (- wy) (+ 1.0 wz))
         s2 (+ 1.0 (* chi chi) (* psi psi))
         f [(/ (+ (- 1.0 (* chi chi)) (* psi psi)) s2) (/ (* 2.0 chi psi) s2) (/ (* -2.0 chi) s2)]
         g [(/ (* 2.0 chi psi) s2) (/ (- (+ 1.0 (* chi chi)) (* psi psi)) s2) (/ (* 2.0 psi) s2)]
         a (/ 1.0 (- (/ 2.0 rm) (/ (lin/dot v v) mu)))
         evec (lin/sub (lin/scale r (- (/ (lin/dot v v) mu) (/ 1.0 rm))) (lin/scale v (/ (lin/dot r v) mu)))
         af (lin/dot evec f) ag (lin/dot evec g)
         X (lin/dot r f) Y (lin/dot r g)
         b (/ 1.0 (+ 1.0 (math/sqrt (- 1.0 (* af af) (* ag ag)))))
         e2 (math/sqrt (- 1.0 (* af af) (* ag ag)))
         sf (+ ag (/ (- (* (- 1.0 (* ag ag b)) Y) (* af ag b X)) (* a e2)))
         cf (+ af (/ (- (* (- 1.0 (* af af b)) X) (* af ag b Y)) (* a e2)))
         F (math/atan2 sf cf)]
     [a af ag chi psi (+ (- F (* af (math/sin F))) (* ag (math/cos F)))])))

;; ------------------------------------------------------------ averaging

(defn- rates
  "The perturbation's rates of the equinoctial elements at `x` and time
  `t`: Gauss's equations, as the elements' partials in the velocity times
  the acceleration -- the partials by central differences."
  [accel mu t x]
  (let [[r v] (elements->state mu x)
        d (accel t r v)
        J (diff/jacobian (fn [v'] (let [y (state->elements mu [r v'])]
                                    ;; the mean longitude continuous across the step
                                    (update y 5 #(+ (x 5) (am/wrap-angle (- % (x 5)))))))
                         v {:steps [1e-7 1e-7 1e-7] :richardson? false})]
    (lin/mat-vec J d)))

(defn- samples
  "The rates at `n` mean longitudes round the orbit of mean elements `x`."
  [accel mu t x n]
  (mapv (fn [k] (rates accel mu t (assoc x 5 (* 2.0 math/PI (/ k n))))) (range n)))

;; -------------------------------------------------------- short-periodic

(defn- fourier-integral
  "The zero-mean integral over lambda, divided by `nm`, of the function
  sampled as `ys` at the `n` longitudes `lams`: its Fourier series
  integrated term by term, as a function of lambda."
  [ys lams n nm]
  (let [terms (vec (for [k (range 1 (quot n 2))]
                     [k (* (/ 2.0 n) (reduce + (map #(* %1 (math/cos (* k %2))) ys lams)))
                      (* (/ 2.0 n) (reduce + (map #(* %1 (math/sin (* k %2))) ys lams)))]))]
    (fn [lam] (reduce + (for [[k ak bk] terms]
                          (/ (- (* ak (math/sin (* k lam))) (* bk (math/cos (* k lam)))) (* k nm)))))))

(defn- first-order
  "What the first order yields at mean elements `x`: the samples of the
  perturbation's rates, and the short-periodic variations as a function
  of the mean longitude -- each rate's departure from its mean integrated
  over lambda, the mean longitude's with the mean-motion change -3n/(2a)
  times a's variation."
  [accel mu t x n]
  (let [lams (mapv #(* 2.0 math/PI (/ % n)) (range n))
        rs (samples accel mu t x n)
        nm (math/sqrt (/ mu (math/pow (x 0) 3)))
        eta-a (fourier-integral (map first rs) lams n nm)
        dn (fn [l] (* -1.5 (/ nm (x 0)) (eta-a l)))
        etas (vec (concat (for [i (range 1 5)] (fourier-integral (map #(nth % i) rs) lams n nm))
                          [(fourier-integral (map (fn [r l] (+ (nth r 5) (dn l))) rs lams) lams n nm)]))]
    {:lams lams :rates rs :nm nm
     :eta (fn [lam] (into [(eta-a lam)] (map #(% lam) etas)))}))

(defn short-periodic
  "The first-order short-periodic variations of the elements, osculating
  less mean, at mean elements `x` and time `t`: each rate's departure from
  its mean over the mean longitude, its Fourier series from `n` samples
  (default 64) integrated term by term, the mean longitude's with the
  mean-motion change -3n/(2a) times a's variation integrated too."
  ([accel t x] (short-periodic accel mu t x 64))
  ([accel mu t x n] ((:eta (first-order accel mu t x n)) (x 5))))

(defn mean-rates
  "The averaged equations of motion at mean elements `x`: the Keplerian
  mean motion, the perturbation's rates averaged over the mean longitude
  (first order), and -- with `second?`, the default -- the second-order
  averages of the generalized method of averaging, <(dF/dx) eta> + 1/2
  n''(a) <eta_a^2> in the mean longitude's, F the perturbation's rates
  and eta the first-order short-periodic variations; the averages by the
  trapezoid rule on `n` points (default 64), exact to rounding for so
  smooth a periodic integrand, the partials of F by central differences."
  ([accel t x] (mean-rates accel mu t x 64 true))
  ([accel mu t x n second?]
   (let [{:keys [lams nm eta] rs :rates} (first-order accel mu t x n)
         a1 (lin/scale (reduce lin/add rs) (/ 1.0 n))
         a2 (when second?
              (let [steps [(* 1e-7 (x 0)) 1e-7 1e-7 1e-7 1e-7 1e-6]
                    terms (for [l lams
                                :let [xl (assoc x 5 l)
                                      J (diff/jacobian #(rates accel mu t %) xl {:steps steps :richardson? false})
                                      e (eta l)]]
                            (-> (lin/mat-vec J e)
                                (update 5 + (* 0.5 (/ (* 15.0 nm) (* 4.0 (x 0) (x 0))) (e 0) (e 0)))))]
                (lin/scale (reduce lin/add terms) (/ 1.0 n))))]
     (-> (if a2 (lin/add a1 a2) a1)
         (update 5 + nm)))))

(defn propagate-mean
  "The mean elements `x0` at `t0` carried to `t1` by the averaged
  equations, fourth-order Runge-Kutta at steps of about `h` seconds
  (default six hours)."
  ([accel x0 t0 t1] (propagate-mean accel mu x0 t0 t1 21600.0))
  ([accel mu x0 t0 t1 h] (propagate-mean accel mu x0 t0 t1 h true))
  ([accel mu x0 t0 t1 h second?]
   (let [steps (max 1 (long (math/ceil (/ (abs (- t1 t0)) h))))
         h (/ (- t1 t0) steps)
         f (fn [t x] (mean-rates accel mu t x 64 second?))]
     (loop [t t0 x x0 k 0]
       (if (= k steps)
         x
         (let [k1 (f t x)
               k2 (f (+ t (* 0.5 h)) (lin/add x (lin/scale k1 (* 0.5 h))))
               k3 (f (+ t (* 0.5 h)) (lin/add x (lin/scale k2 (* 0.5 h))))
               k4 (f (+ t h) (lin/add x (lin/scale k3 h)))]
           (recur (+ t h) (lin/add x (lin/scale (lin/add (lin/add k1 (lin/scale k2 2.0)) (lin/add (lin/scale k3 2.0) k4))
                                                (/ h 6.0)))
                  (inc k))))))))

(defn osculating
  "The osculating state of mean elements `x` at `t`: the mean elements
  plus their short-periodic variations."
  ([accel t x] (osculating accel mu t x))
  ([accel mu t x] (elements->state mu (lin/add x (short-periodic accel mu t x 64)))))

(defn mean-elements
  "The mean elements whose osculating state at `t` is `s`, by fixed-point
  iteration."
  ([accel t s] (mean-elements accel mu t s))
  ([accel mu t s]
   (let [target (state->elements mu s)]
     (loop [x target k 0]
       (let [x' (lin/sub target (short-periodic accel mu t x 64))
             x' (update x' 5 #(+ (x 5) (am/wrap-angle (- % (x 5)))))]
         (if (or (> k 20) (< (lin/distance x x') 1e-12))
           x'
           (recur x' (inc k))))))))

(defn propagate
  "The osculating state at `t1` of the orbit whose osculating state at
  `t0` is `s0`: to mean elements, the averaged equations at steps of `h`
  (default six hours), and back to osculating."
  ([accel s0 t0 t1] (propagate accel mu s0 t0 t1 21600.0))
  ([accel mu s0 t0 t1 h]
   (osculating accel mu t1 (propagate-mean accel mu (mean-elements accel mu t0 s0) t0 t1 h))))
