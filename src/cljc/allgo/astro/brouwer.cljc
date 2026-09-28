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

  The mean elements drift at J2's first-order secular rates
  (`allgo.astro.perturbations/j2-secular`, which are Brouwer's).
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

(defn short-period
  "Brouwer's first-order short-period corrections, osculating less mean,
  to the Delaunay variables [L G H l g h] of mean elements `el`."
  ([el] (short-period mu c/R-earth geo/J2 el))
  ([mu R J2 el]
   (let [k2 (* 0.5 J2 R R)
         [L G H l g] (delaunay mu el)
         ;; L and G stepped by less than their difference, which keeps e real
         hLG (min (* 1e-6 L) (* 0.25 (- L G)))
         [[WL WG WH Wl Wg]] (diff/jacobian (fn [x] [(generating-function mu k2 x)]) [L G H l g]
                                           {:steps [hLG hLG (* 1e-6 G) 1e-6 1e-6]})]
     [Wl Wg 0.0 (- WL) (- WG) (- WH)])))

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
  :M}`, the corrections applied in Lyddane's way."
  ([el] (osculating-elements mu c/R-earth geo/J2 el))
  ([mu R J2 {:keys [a e i raan argp] :as el}]
   (let [[L G] (delaunay mu el)
         [dL dG _ dl dg dh] (short-period mu R J2 el)
         ;; a = L^2/mu; e^2 = 1 - G^2/L^2; cos i = H/G with H unchanged
         da (/ (* 2.0 L dL) mu)
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

(defn propagate
  "The osculating state `dt` seconds after an epoch where the mean
  elements are `el`: the mean elements carried by J2's secular rates, then
  the short-period terms put back."
  ([el dt] (propagate mu c/R-earth geo/J2 el dt))
  ([mu R J2 {:keys [a e i raan argp M] :as el} dt]
   (let [rates (pert/j2-secular mu R J2 a e i)]
     (osculating-state mu R J2 (assoc el
                                      :raan (+ raan (* (:raan rates) dt))
                                      :argp (+ argp (* (:argp rates) dt))
                                      :M (+ M (* (:M rates) dt)))))))
