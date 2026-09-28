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
  second-order short-period terms are periodic and are left out.

  What the short-period averaging leaves in the second-order Hamiltonian
  still turns with the perigee, and Brouwer's second transformation takes
  that out: the long-period generator W*, of first order in J2 --
  second-order terms over the first-order perigee rate -- here from a
  Fourier series in g of the Hamiltonian averaged numerically over l.
  Brouwer's mean elements are doubly averaged, the long-period terms
  applied first and the short-period after. The perigee rate vanishes at
  the critical inclination, 63.4 degrees, and there the long-period terms
  fail, as Brouwer's do: where they grow too large to trust, the
  conversions and the propagator throw (`long-period-limit`).

  With all of it, what is left against J2's motion integrated
  numerically is third order in J2. Elements are `{:a :e :i :raan :argp
  :M}`, km and radians."
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

(defn- eccentricity [L G] (/ (math/sqrt (* (- L G) (+ L G))) L))

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
  the Hamiltonian second order in J2 once W1 has acted, every partial in
  closed form -- the long-period part of its average is of order e^2,
  which differences in G, their steps bounded by L - G, would bury in
  round-off. With A = (1 - 3 theta^2)/2, B = 3/2 (1 - theta^2) and u2 =
  2f + 2g, at fixed l: df/dl = (a/r)^2 eta, df/de = sin f (2 + e cos f)/
  eta^2, dr/de = -a cos f; F1 goes as L^-6 at fixed e and theta; and
  de/dL = eta^2/(L e), de/dG = -eta/(L e), dtheta/dG = -theta/G."
  [mu k2 [L G H l g]]
  (let [a (/ (* L L) mu)
        e (eccentricity L G)
        eta (/ G L) eta2 (* eta eta)
        th (/ H G) th2 (* th th)
        n (/ (* mu mu) (* L L L))
        E (kepler/kepler-equation l e)
        r (* a (- 1.0 (* e (math/cos E))))
        f (kepler/eccentric->true E e)
        cf (math/cos f) sf (math/sin f)
        A (* 0.5 (- 1.0 (* 3.0 th2)))
        B (* 1.5 (- 1.0 th2))
        u2 (+ (* 2.0 f) (* 2.0 g))
        c1 (math/cos (+ f (* 2.0 g))) c2 (math/cos u2) c3 (math/cos (+ (* 3.0 f) (* 2.0 g)))
        ;; F1 and its partials
        P (/ (* mu k2) (* r r r))
        S (- A (* B c2))
        fe (/ (* sf (+ 2.0 (* e cf))) eta2)
        Fe (* P (+ (/ (* 3.0 a cf S) r) (* 2.0 B (math/sin u2) fe)))
        Fth (* P 3.0 th (- c2 1.0))
        FL (+ (/ (* -6.0 P S) L) (* Fe (/ eta2 (* L e))))
        FG (+ (* Fe (/ (- eta) (* L e))) (* Fth (/ (- th) G)))
        ;; W1's partials in the angles
        pre (- (/ (* n k2) (* eta2 eta)))
        fl (* (/ (* a a) (* r r)) eta)
        Wl (* pre (- (* A (+ (- fl 1.0) (* e cf fl)))
                     (* B fl (+ c2 (* 0.5 e c1) (* 0.5 e c3)))))
        Wg (* pre (- (* B (+ c2 (* e c1) (* (/ e 3.0) c3)))))]
    (+ (* -1.5 (/ (* mu mu) (math/pow L 4)) Wl Wl) (* FL Wl) (* FG Wg))))

(defn- l-points
  "How many points the average over l takes, for eccentricity `e`."
  [e]
  (+ 32 (long (math/ceil (* 200.0 e)))))

(defn- l-averaged-integrand
  "The second-order integrand averaged over l, by the trapezoid rule on
  `N` points, at momenta [L G H] and perigee `g`: the Hamiltonian's
  second-order part, secular and long-period."
  [mu k2 [L G H] g N]
  (/ (reduce + (map #(second-order-integrand mu k2 [L G H (* 2.0 math/PI (/ % N)) g]) (range N))) N))

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
  (let [n (/ (* mu mu) (* L L L))
        avg (l-averaged-integrand mu k2 [L G H] g (l-points (eccentricity L G)))
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

(defn- corrected
  "Elements `el` with the Delaunay corrections applied: L to `L'` in full,
  the others to first order in the nonsingular elements, Lyddane's way;
  `dL` is L's first-order part, which is what e's correction takes."
  [mu {:keys [a e i raan argp] :as el} L' dL dG dl dg dh]
  (let [[L G] (delaunay mu el)
        ;; e^2 = 1 - G^2/L^2; cos i = H/G with H unchanged
        da (- (/ (* L' L') mu) a)
        de (/ (- (/ (* G G dL) (* L L L)) (/ (* G dG) (* L L))) e)
        si (math/sin i)
        di (if (zero? si) 0.0 (/ (* (math/cos i) dG) (* G si)))
        lp (+ raan argp) dlp (+ dg dh)
        t (math/tan (* 0.5 i)) dt (/ di (* 2.0 (math/pow (math/cos (* 0.5 i)) 2)))
        [_ af ag chi psi lam] (->equinoctial el)]
    (equinoctial-> [(+ a da)
                    (+ af (* de (math/cos lp)) (- (* e (math/sin lp) dlp)))
                    (+ ag (* de (math/sin lp)) (* e (math/cos lp) dlp))
                    (+ chi (* dt (math/sin raan)) (* t (math/cos raan) dh))
                    (+ psi (* dt (math/cos raan)) (- (* t (math/sin raan) dh)))
                    (+ lam dl dlp)])))

(defn short-period-elements
  "The osculating elements of the singly averaged elements `el` -- the
  long-period terms still in them: the short-period corrections, first
  order in J2 but for L, which is second. L sets the mean motion: the
  other elements' second-order errors are periodic, L's would be a drift."
  ([el] (short-period-elements mu c/R-earth geo/J2 el))
  ([mu R J2 {:keys [argp M] :as el}]
   (let [[L G H] (delaunay mu el)
         [dL dG _ dl dg dh] (short-period mu R J2 el)
         L' (+ L (L-correction mu (* 0.5 J2 R R) [L G H M argp] dl dg))]
     (corrected mu el L' dL dG dl dg dh))))

(defn- long-period-generator
  "Brouwer's long-period generator W* and its derivative in g, at momenta
  [L G H] and perigee `g`, the l average taken on `N` points: the
  Hamiltonian's second-order part K2, averaged over l, is sampled at eight
  perigees and expanded in a Fourier series in g, and

    W* = -(1/g1') integral of (K2 - <K2>_g) dg,    W*_g = -(K2 - <K2>_g)/g1'

  g1' the first-order perigee rate -- which vanishes at the critical
  inclination, where the long-period terms, like Brouwer's, fail."
  [mu R J2 [L e th] g N]
  (let [k2 (* 0.5 J2 R R)
        eta (math/sqrt (- 1.0 (* e e)))
        G (* L eta) H (* G th)
        a (/ (* L L) mu)
        p (* a eta eta)
        ;; J2's first-order perigee rate, in theta rather than i: steps in
        ;; theta may pass 1
        gdot (* 0.75 (math/sqrt (/ mu (* a a a))) J2 (/ (* R R) (* p p)) (- (* 5.0 th th) 1.0))
        gs (map #(* 2.0 math/PI (/ % 8)) (range 8))
        K (mapv #(l-averaged-integrand mu k2 [L G H] % N) gs)
        coef (fn [f k w] (* w (reduce + (map (fn [Kj gj] (* Kj (f (* k gj)))) K gs))))
        ab (for [k [1 2 3]] [k (coef math/cos k 0.25) (coef math/sin k 0.25)])
        a4 (coef math/cos 4 0.125)]
    {:W (/ (- (+ (reduce + (map (fn [[k a b]] (/ (- (* a (math/sin (* k g))) (* b (math/cos (* k g)))) k)) ab))
                 (* 0.25 a4 (math/sin (* 4.0 g)))))
           gdot)
     :Wg (/ (- (+ (reduce + (map (fn [[k a b]] (+ (* a (math/cos (* k g))) (* b (math/sin (* k g))))) ab))
                  (* a4 (math/cos (* 4.0 g)))))
            gdot)}))

(defn long-period
  "Brouwer's first-order long-period corrections, singly averaged less
  doubly, to the Delaunay variables [L G H l g h] of mean elements `el`:
  W*_g, -W*_L, -W*_G and -W*_H, L and H unchanged. W* is differenced in
  L, e and theta -- steps in e in proportion to e, where steps in G would
  be bounded by L - G, of order e^2 -- and carried to L, G and H by
  de/dL = eta^2/(L e), de/dG = -eta/(L e) and dtheta/dG = -theta/G."
  ([el] (long-period mu c/R-earth geo/J2 el))
  ([mu R J2 {:keys [e i argp] :as el}]
   (let [[L G] (delaunay mu el)
         eta (/ G L)
         th (math/cos i)
         N (l-points e)
         W (fn [L e th] (:W (long-period-generator mu R J2 [L e th] argp N)))
         d (fn [f h] (/ (- (f h) (f (- h))) (* 2.0 h)))
         WL (d #(W (+ L %) e th) (* 1e-4 L))
         We (d #(W L (+ e %) th) (* 1e-3 e))
         Wth (d #(W L e (+ th %)) 1e-4)]
     [0.0 (:Wg (long-period-generator mu R J2 [L e th] argp N)) 0.0
      (- (+ WL (* We (/ (* eta eta) (* L e)))))
      (- (+ (* We (/ (- eta) (* L e))) (* Wth (/ (- th) G))))
      (- (/ Wth G))])))

(def long-period-limit
  "The largest long-period correction, in any of the nonsingular elements
  e cos(g + h), e sin(g + h), tan(i/2) sin h, tan(i/2) cos h and the mean
  longitude, the theory accepts. The terms are first order in it, so what
  they leave out goes as its square: at 1e-3 thirty orbits are still
  within a few meters of J2's motion integrated numerically, at 2e-3 they
  are some twenty out and at 4e-3 a hundred. Ordinary orbits keep below
  2e-4; near the critical inclination the correction grows without bound,
  the sooner the more eccentric the orbit -- within a degree of it at e =
  0.3, a twentieth of one at e = 0.01."
  1e-3)

(defn long-period-elements
  "The singly averaged elements of the doubly averaged `el`: the
  long-period corrections, first order in J2. Throws, with `:type`
  `::critical-inclination` in its data, where the correction passes
  `long-period-limit` -- near the critical inclination, 63.4 degrees or
  116.6, where the first-order perigee rate the terms divide by vanishes."
  ([el] (long-period-elements mu c/R-earth geo/J2 el))
  ([mu R J2 el]
   (let [[L] (delaunay mu el)
         [_ dG _ dl dg dh] (long-period mu R J2 el)
         el' (corrected mu el L 0.0 dG dl dg dh)
         ;; the mean longitude's change wrapped: `el`'s angles may run on
         ;; past 2 pi, as `propagate` carries them, and `el'`'s are reduced
         change (update (mapv - (->equinoctial el') (->equinoctial el)) 5 am/wrap-angle)
         size (apply max (map abs (rest change)))]
     (when-not (<= size long-period-limit)
       (throw (ex-info (str "Too near the critical inclination for Brouwer's long-period terms: a correction of "
                            size " passes the limit of " long-period-limit)
                       {:type ::critical-inclination :i (:i el) :e (:e el) :correction size
                        :limit long-period-limit})))
     el')))

(defn osculating-elements
  "The osculating elements of Brouwer's mean elements `el` -- doubly
  averaged, the long- and short-period terms both taken out -- `{:a :e :i
  :raan :argp :M}`."
  ([el] (osculating-elements mu c/R-earth geo/J2 el))
  ([mu R J2 el]
   (short-period-elements mu R J2 (long-period-elements mu R J2 el))))

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
