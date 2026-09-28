(ns allgo.astro.brouwer
  "J2's first-order short-period terms in full (Brouwer, \"Solution of
  the problem of artificial satellite theory without drag\", AJ 64, 1959;
  Lyddane, \"Small eccentricities or inclinations in the Brouwer theory
  of the artificial satellite\", AJ 68, 1963): every term, the ones of
  order J2 e included that Spacetrack Report No. 3 drops for SGP4
  (`allgo.astro.perturbations/j2-osculating`).

  Brouwer's theory is a canonical transformation. In Delaunay's
  variables -- L = sqrt(mu a), G = L sqrt(1 - e^2), H = G cos i and their
  angles l = M, g = argp, h = raan -- J2's Hamiltonian is

    F1 = (mu k2 / r^3) ((1 - 3 theta^2)/2 - 3/2 (1 - theta^2) cos(2f + 2g)),

  k2 = J2 R^2/2 and theta = cos i, and the generating function that
  removes its dependence on the mean anomaly is

    W1 = -(n k2/eta^3) [ (1 - 3 theta^2)/2 (f - l + e sin f)
           - 3/2 (1 - theta^2) (1/2 sin(2f + 2g) + e/2 sin(f + 2g)
                                + e/6 sin(3f + 2g)) ],

  eta = sqrt(1 - e^2) -- from n dW1/dl = <F1> - F1 with dl = r^2/(a^2
  eta) df. The osculating momenta are the mean ones plus dW1/d(angle),
  the osculating angles the mean ones less dW1/d(momentum). Here W1 is
  evaluated in closed form and its derivatives taken numerically, to
  better than ten digits of terms themselves of order J2.

  Near-circular and near-equatorial orbits make the corrections to e and
  g, and to i and h, singly large; Lyddane's remedy is to apply them
  together, to first order, in the nonsingular elements e cos(g + h),
  e sin(g + h), tan(i/2) sin h, tan(i/2) cos h and the mean longitude,
  as here. That leaves the exactly circular orbit (e = 0) and i = 180
  degrees.

  The mean elements drift at Brouwer's secular rates to second order in
  J2 (`secular-rates`); the first-order part is
  `allgo.astro.perturbations/j2-secular`'s. Second-order rates are only
  worth having if the mean L they are evaluated at is right to second
  order too -- an error of J2^2 in L is one of J2^2 n t in the mean
  anomaly -- so L alone is carried to second order, through W2's
  derivative in l, found by averaging numerically. The other elements'
  second-order terms are periodic and are left out, as are Brouwer's
  long-period terms, which over days show as a drift of order J2^2 e.
  Elements are `{:a :e :i :raan :argp :M}`, km and radians."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.geopotential :as geo]
            [allgo.astro.kepler :as kepler]
            [allgo.astro.perturbations :as pert]
            [allgo.math :as am]
            [allgo.numerics.differentiation :as diff]
            [clojure.math :as math]))

(def ^:private mu c/GM-earth)

(defn generating-function
  "Brouwer's W1, km^2/s, at the Delaunay variables [L G H l g], for
  k2 = J2 R^2/2."
  [mu k2 [L G H l g]]
  (let [eta (/ G L)
        e (/ (math/sqrt (* (- L G) (+ L G))) L)
        th2 (let [th (/ H G)] (* th th))
        n (/ (* mu mu) (* L L L))
        f (kepler/eccentric->true (kepler/kepler-equation l e) e)
        g2 (* 2.0 g)]
    (* (- (/ (* n k2) (* eta eta eta)))
       (- (* 0.5 (- 1.0 (* 3.0 th2)) (+ (am/wrap-angle (- f l)) (* e (math/sin f))))
          (* 1.5 (- 1.0 th2) (+ (* 0.5 (math/sin (+ (* 2.0 f) g2)))
                                (* 0.5 e (math/sin (+ f g2)))
                                (* (/ e 6.0) (math/sin (+ (* 3.0 f) g2)))))))))

(defn delaunay
  "[L G H l g h] of classical elements."
  ([el] (delaunay mu el))
  ([mu {:keys [a e i raan argp M]}]
   (let [L (math/sqrt (* mu a))
         G (* L (math/sqrt (- 1.0 (* e e))))]
     [L G (* G (math/cos i)) M argp raan])))

(defn- disturbing-function
  "J2's Hamiltonian F1 at the Delaunay variables [L G H l g]."
  [mu k2 [L G H l g]]
  (let [a (/ (* L L) mu)
        e (/ (math/sqrt (* (- L G) (+ L G))) L)
        th2 (let [th (/ H G)] (* th th))
        E (kepler/kepler-equation l e)
        f (kepler/eccentric->true E e)
        r (* a (- 1.0 (* e (math/cos E))))]
    (* (/ (* mu k2) (* r r r))
       (- (* 0.5 (- 1.0 (* 3.0 th2))) (* 1.5 (- 1.0 th2) (math/cos (+ (* 2.0 f) (* 2.0 g))))))))

(defn- step-LG
  "The step for differences in L and G: less than their difference, which
  keeps e real."
  [L G]
  (min (* 1e-6 L) (* 0.25 (- L G))))

(defn- w1-partials
  "[W1_L W1_G W1_H W1_l W1_g] at [L G H l g]."
  [mu k2 [L G :as x]]
  (let [h (step-LG L G)]
    (first (diff/jacobian (fn [x] [(generating-function mu k2 x)]) x {:steps [h h (* 1e-6 G) 1e-6 1e-6]}))))

(defn short-period
  "Brouwer's first-order short-period corrections, osculating less mean,
  to the Delaunay variables [L G H l g h] of mean elements `el`."
  ([el] (short-period mu c/R-earth geo/J2 el))
  ([mu R J2 el]
   (let [[WL WG WH Wl Wg] (w1-partials mu (* 0.5 J2 R R) (subvec (delaunay mu el) 0 5))]
     [Wl Wg 0.0 (- WL) (- WG) (- WH)])))

(defn- second-order-integrand
  "1/2 F0_LL W1_l^2 + F1_L W1_l + F1_G W1_g at [L G H l g]: the part of
  the Hamiltonian second order in J2 once W1 has acted."
  [mu k2 [L G H l g :as x]]
  (let [[_ _ _ Wl Wg] (w1-partials mu k2 x)
        h (step-LG L G)
        [[FL FG]] (diff/jacobian (fn [[L G]] [(disturbing-function mu k2 [L G H l g])]) [L G] {:steps [h h]})]
    (+ (* -1.5 (/ (* mu mu) (math/pow L 4)) Wl Wl) (* FL Wl) (* FG Wg))))

(defn- L-correction
  "L osculating less mean, to second order in J2, at the mean elements'
  Delaunay variables [L G H l g], the osculating angles `dl` `dg` from
  them: W1_l + W2_l, the type-2 generating function's derivative, taken
  at the osculating angles, with n W2_l the integrand's average over l
  less its value. W1_l is carried to the osculating angles linearly, by
  its derivative along (dl, dg): for a small eccentricity dl and dg are
  singly of order J2/e, and a step that large taken in full would leave
  terms of order J2^3/e^2."
  [mu k2 [L G H l g :as x] dl dg]
  (let [e (/ (math/sqrt (* (- L G) (+ L G))) L)
        n (/ (* mu mu) (* L L L))
        N (+ 32 (long (math/ceil (* 200.0 e))))
        avg (/ (reduce + (map #(second-order-integrand mu k2 [L G H (* 2.0 math/PI (/ % N)) g]) (range N))) N)
        Wl (fn [s] (nth (w1-partials mu k2 [L G H (+ l (* s dl)) (+ g (* s dg))]) 3))
        ;; a step of 0.01 rad along the shift, Richardson-refined: W1_l's
        ;; own round-off, divided by a smaller one, would be millimeters
        eps (/ 1e-2 (max (abs dl) (abs dg) 1e-300))
        d (fn [h] (/ (- (Wl h) (Wl (- h))) (* 2.0 h)))]
    (+ (Wl 0.0)
       (/ (- (* 4.0 (d (* 0.5 eps))) (d eps)) 3.0)
       (/ (- avg (second-order-integrand mu k2 x)) n))))

(defn- ->equinoctial
  "[a af ag chi psi lambda] of classical elements."
  [{:keys [a e i raan argp M]}]
  (let [lp (+ raan argp) t (math/tan (* 0.5 i))]
    [a (* e (math/cos lp)) (* e (math/sin lp)) (* t (math/sin raan)) (* t (math/cos raan)) (+ lp M)]))

(defn- equinoctial->
  [[a af ag chi psi lam]]
  (let [lp (math/atan2 ag af) node (math/atan2 chi psi)]
    {:a a :e (math/hypot af ag) :i (* 2.0 (math/atan (math/hypot chi psi)))
     :raan (am/wrap-2pi node) :argp (am/wrap-2pi (- lp node)) :M (am/wrap-2pi (- lam lp))}))

(defn osculating-elements
  "The osculating elements of mean elements `el`, `{:a :e :i :raan :argp
  :M}`, the corrections applied in Lyddane's way, first order in J2 but
  for a, which is second."
  ([el] (osculating-elements mu c/R-earth geo/J2 el))
  ([mu R J2 {:keys [a e i raan argp M] :as el}]
   (let [[L G H] (delaunay mu el)
         [dL dG _ dl dg dh] (short-period mu R J2 el)
         ;; L to second order, which sets the mean motion: the other
         ;; elements' second-order errors are periodic, L's would be a drift
         L' (+ L (L-correction mu (* 0.5 J2 R R) [L G H M argp] dl dg))
         ;; e^2 = 1 - G^2/L^2; cos i = H/G with H unchanged
         da (- (/ (* L' L') mu) a)
         de (/ (- (/ (* G G dL) (* L L L)) (/ (* G dG) (* L L))) e)
         si (math/sin i)
         di (if (zero? si) 0.0 (/ (* (math/cos i) dG) (* G si)))
         ;; to first order in the nonsingular elements
         lp (+ raan argp) dlp (+ dg dh)
         t (math/tan (* 0.5 i)) dt (/ di (* 2.0 (math/pow (math/cos (* 0.5 i)) 2)))
         [_ af ag chi psi lam] (->equinoctial el)]
     (equinoctial-> [(+ a da)
                     (+ af (* de (math/cos lp)) (- (* e (math/sin lp) dlp)))
                     (+ ag (* de (math/sin lp)) (* e (math/cos lp) dlp))
                     (+ chi (* dt (math/sin raan)) (* t (math/cos raan) dh))
                     (+ psi (* dt (math/cos raan)) (- (* t (math/sin raan) dh)))
                     (+ lam dl dlp)]))))

(defn osculating-state
  "The osculating state `[r v]` of mean elements `el`."
  ([el] (osculating-state mu c/R-earth geo/J2 el))
  ([mu R J2 el] (kepler/elements->state mu (osculating-elements mu R J2 el))))

(defn mean-elements
  "The mean elements whose osculating state is `s`, by fixed-point
  iteration in the nonsingular elements: each step moves the guess by how
  far its own osculating elements miss the target's, and J2 being small,
  a few settle it."
  ([s] (mean-elements mu c/R-earth geo/J2 s))
  ([mu R J2 [r v]]
   (let [target (->equinoctial (kepler/state->elements mu r v))
         miss (fn [x y] (let [d (mapv - x y)] (update d 5 am/wrap-angle)))]
     (loop [guess target k 0]
       (let [step (miss target (->equinoctial (osculating-elements mu R J2 (equinoctial-> guess))))
             guess' (mapv + guess step)]
         (if (or (>= k 50) (< (apply max (abs (/ (step 0) (guess 0))) (map abs (rest step))) 1e-13))
           (equinoctial-> guess')
           (recur guess' (inc k))))))))

(defn secular-rates
  "Brouwer's secular rates to second order in J2, `{:raan :argp :M}`,
  rad/s, the last the whole mean-anomaly rate, of mean elements `a` `e`
  `i`: J2's first-order rates and, with gamma' = k2/(a^2 eta^4), eta =
  sqrt(1 - e^2) and theta = cos i,

    dl/dt += 3/32 n gamma'^2 eta [-15 + 16 eta + 25 eta^2
               + (30 - 96 eta - 90 eta^2) theta^2 + (105 + 144 eta + 25 eta^2) theta^4]
    dg/dt += 3/32 n gamma'^2 [-35 + 24 eta + 25 eta^2
               + (90 - 192 eta - 126 eta^2) theta^2 + (385 + 360 eta + 45 eta^2) theta^4]
    dh/dt += 3/8 n gamma'^2 [(-5 + 12 eta + 9 eta^2) theta + (-35 - 36 eta - 5 eta^2) theta^3]

  -- the derivatives of the mean Hamiltonian's second-order secular part,
  the average over l and g of 1/2 F0_LL W1_l^2 + F1_L W1_l + F1_G W1_g.
  At eta = 1 they are SGP4's J2^2 terms."
  ([a e i] (secular-rates mu c/R-earth geo/J2 a e i))
  ([mu R J2 a e i]
   (let [{:keys [raan argp M]} (pert/j2-secular mu R J2 a e i)
         n (math/sqrt (/ mu (* a a a)))
         eta (math/sqrt (- 1.0 (* e e))) eta2 (* eta eta)
         gp (/ (* 0.5 J2 R R) (* a a eta2 eta2))
         k (* n gp gp)
         th (math/cos i) th2 (* th th) th4 (* th2 th2)]
     {:raan (+ raan (* 0.375 k (+ (* (+ -5.0 (* 12.0 eta) (* 9.0 eta2)) th)
                                  (* (+ -35.0 (* -36.0 eta) (* -5.0 eta2)) th th2))))
      :argp (+ argp (* (/ 3.0 32.0) k (+ -35.0 (* 24.0 eta) (* 25.0 eta2)
                                         (* (+ 90.0 (* -192.0 eta) (* -126.0 eta2)) th2)
                                         (* (+ 385.0 (* 360.0 eta) (* 45.0 eta2)) th4))))
      :M (+ M (* (/ 3.0 32.0) k eta (+ -15.0 (* 16.0 eta) (* 25.0 eta2)
                                       (* (+ 30.0 (* -96.0 eta) (* -90.0 eta2)) th2)
                                       (* (+ 105.0 (* 144.0 eta) (* 25.0 eta2)) th4))))})))

(defn propagate
  "The osculating state `dt` seconds after an epoch where the mean
  elements are `el`: the mean elements carried by Brouwer's secular
  rates, then the short-period terms put back."
  ([el dt] (propagate mu c/R-earth geo/J2 el dt))
  ([mu R J2 {:keys [a e i raan argp M] :as el} dt]
   (let [rates (secular-rates mu R J2 a e i)]
     (osculating-state mu R J2 (assoc el
                                      :raan (+ raan (* (:raan rates) dt))
                                      :argp (+ argp (* (:argp rates) dt))
                                      :M (+ M (* (:M rates) dt)))))))
