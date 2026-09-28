(ns allgo.astro.covariance
  "An orbit's uncertainty carried between the ways of writing it
  (Vallado, *Fundamentals of Astrodynamics and Applications*, chapter
  10's covariance transformations): Cartesian, classical and equinoctial
  elements, flight elements, and the satellite's own radial/along-track/
  cross-track and velocity-aligned frames.

  To first order a covariance P in one set becomes J P J^T in another, J
  the Jacobian of the conversion. For the local frames J is a rotation and
  exact; for the classical elements it is derived analytically
  (`classical-partials`); for the equinoctial and flight elements it is
  taken numerically (`allgo.numerics.differentiation`), to ten digits, the
  angles among the outputs wrapped.

  States are `[r v]`, km and km/s; element vectors are ordered as their
  functions below say, radians for the angles."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.kepler :as kepler]
            [allgo.astro.states :as states]
            [allgo.geometry.rotation :as rot]
            [allgo.geometry.vec3 :as v3]
            [allgo.numerics.differentiation :as diff]
            [allgo.numerics.linear :as lin]
            [clojure.math :as math]))

(def ^:private mu c/GM-earth)

;; ------------------------------------------------------------ machinery

(defn- flat [[r v]] (vec (concat r v)))
(defn- state [x] [(subvec x 0 3) (subvec x 3 6)])

;; ------------------------------------------------------ element sets

(defn classical-vector
  "[a e i raan argp M] of a flat state [x y z vx vy vz]."
  ([x] (classical-vector mu x))
  ([mu x] (let [[r v] (state x) {:keys [a e i raan argp M]} (kepler/state->elements mu r v)]
            [a e i raan argp M])))

(defn classical->flat
  "The flat state [x y z vx vy vz] of classical elements [a e i raan argp M]."
  [mu [a e i raan argp M]]
  (flat (kepler/elements->state mu {:a a :e e :i i :raan raan :argp argp :M M})))

(defn classical-partials
  "The Jacobian of the state [x y z vx vy vz] with respect to the classical
  elements [a e i raan argp M], analytically, one column per element:

    a      r/a and -v/2a -- at fixed mean anomaly the orbit scales, and the
           speed with a n, as a^-1/2
    e      through the perifocal coordinates, with dE/de = sin E/(1 - e cos E)
    i      N x r and N x v, N the line of nodes -- a turn about it
    raan   z x r and z x v, a turn about the pole
    argp   W x r and W x v, a turn about the orbit's normal
    M      v/n and -mu r/(n r^3), the motion itself"
  ([el] (classical-partials mu el))
  ([mu [a e i raan argp M]]
   (let [n (kepler/mean-motion mu a)
         E (kepler/kepler-equation M e)
         sE (math/sin E) cE (math/cos E)
         b (math/sqrt (- 1.0 (* e e)))
         D (- 1.0 (* e cE))
         Ee (/ sE D)
         De (+ (- cE) (* e sE Ee))
         ;; the perifocal state and its derivative in e
         rpf [(* a (- cE e)) (* a b sE) 0.0]
         vpf [(/ (* (- a) n sE) D) (/ (* a n b cE) D) 0.0]
         drpf [(* a (- (- (* sE Ee)) 1.0)) (* a (+ (* (/ (- e) b) sE) (* b cE Ee))) 0.0]
         dvpf [(/ (* (- a) n (- (* cE Ee D) (* sE De))) (* D D))
               (/ (* a n (- (* (/ (- e) b) cE D) (* b sE Ee D) (* b cE De))) (* D D))
               0.0]
         Q (rot/chain (rot/rotate-z raan) (rot/rotate-x i) (rot/rotate-z argp))
         to (fn [w] (lin/mat-vec Q w))
         r (to rpf) v (to vpf)
         N [(math/cos raan) (math/sin raan) 0.0]
         W (to [0.0 0.0 1.0])
         z [0.0 0.0 1.0]
         col (fn [dr dv] (vec (concat dr dv)))
         cols [(col (v3/scale r (/ 1.0 a)) (v3/scale v (/ -0.5 a)))
               (col (to drpf) (to dvpf))
               (col (v3/cross N r) (v3/cross N v))
               (col (v3/cross z r) (v3/cross z v))
               (col (v3/cross W r) (v3/cross W v))
               (col (v3/scale v (/ 1.0 n)) (v3/scale r (/ (- mu) (* n (math/pow (v3/length r) 3)))))]]
     (lin/transpose cols))))

(defn cartesian->classical
  "Carry the Cartesian covariance `P` of state `s` into classical
  elements [a e i raan argp M], through the inverse of
  `classical-partials`. Singular where the elements are: circular or
  equatorial orbits."
  ([P s] (cartesian->classical mu P s))
  ([mu P s] (lin/congruence (lin/inverse-general (classical-partials mu (classical-vector mu (flat s)))) P)))

(defn classical->cartesian
  "Carry the covariance `P` of classical elements [a e i raan argp M]
  `el` into Cartesian, through `classical-partials`."
  ([P el] (classical->cartesian mu P el))
  ([mu P el] (lin/congruence (classical-partials mu (vec el)) P)))

(defn equinoctial-vector
  "[a af ag chi psi meanlon] of a flat state -- `allgo.astro.states`'s
  equinoctial set, nonsingular for circular and equatorial orbits."
  ([x] (equinoctial-vector mu x))
  ([mu x] (let [{:keys [a af ag chi psi meanlon]} (states/state->equinoctial mu (state x))]
            [a af ag chi psi meanlon])))

(defn cartesian->equinoctial
  "Carry the Cartesian covariance `P` of state `s` into equinoctial
  elements [a af ag chi psi meanlon], retrograde orbits in their own
  set as `states/state->equinoctial` chooses."
  ([P s] (cartesian->equinoctial mu P s))
  ([mu P s] (lin/congruence (diff/jacobian #(equinoctial-vector mu %) (flat s) {:angles #{5}}) P)))

(defn equinoctial->cartesian
  "Carry the covariance `P` of equinoctial elements `eq`, a map as
  `states/state->equinoctial` gives, into Cartesian."
  ([P eq] (equinoctial->cartesian mu P eq))
  ([mu P {:keys [fr] :as eq}]
   (let [x (mapv eq [:a :af :ag :chi :psi :meanlon])
         f (fn [[a af ag chi psi meanlon]]
             (flat (states/equinoctial->state mu {:a a :af af :ag ag :chi chi :psi psi :meanlon meanlon :fr fr})))]
     (lin/congruence (diff/jacobian f x {:steps (mapv #(* 1e-6 (max 1e-3 (abs %))) x)}) P))))

(defn cartesian->flight
  "Carry the Cartesian (inertial) covariance `P` of state `s` into flight
  elements [rm vm latgc lon fpa az] at TT `mjd-tt` and UT1 `mjd-ut1`."
  [P s mjd-tt mjd-ut1 eop]
  (lin/congruence (diff/jacobian #(vec (take 6 (states/state->flight (state %) mjd-tt mjd-ut1 eop))) (flat s)
                                 {:angles #{2 3 4 5}})
                  P))

;; ------------------------------------------------------- local frames

(defn- block [m] (let [z [0.0 0.0 0.0]]
                   (vec (concat (map #(vec (concat % z)) m) (map #(vec (concat z %)) m)))))

(defn cartesian->rsw
  "Carry the Cartesian covariance `P` of state `s` into the satellite's
  radial, along-track and cross-track axes: a rotation of position and
  velocity alike (the frame's own turning left out, as is usual for
  comparing uncertainties)."
  [P s]
  (lin/congruence (block (states/rsw s)) P))

(defn rsw->cartesian [P s] (lin/congruence (lin/transpose (block (states/rsw s))) P))

(defn cartesian->ntw
  "The same into the velocity-aligned N, T, W axes."
  [P s]
  (lin/congruence (block (states/ntw s)) P))

(defn ntw->cartesian [P s] (lin/congruence (lin/transpose (block (states/ntw s))) P))
