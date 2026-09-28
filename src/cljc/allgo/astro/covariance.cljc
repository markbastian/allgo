(ns allgo.astro.covariance
  "An orbit's uncertainty carried between the ways of writing it
  (Vallado, *Fundamentals of Astrodynamics and Applications*, chapter
  10's covariance transformations): Cartesian, classical and equinoctial
  elements, flight elements, and the satellite's own radial/along-track/
  cross-track and velocity-aligned frames.

  To first order a covariance P in one set becomes J P J^T in another, J
  the Jacobian of the conversion. For the local frames J is a rotation and
  exact. For the element sets it is taken numerically
  (`allgo.numerics.differentiation`), to ten digits for any conversion the
  library has, the angles among the outputs wrapped.

  States are `[r v]`, km and km/s; element vectors are ordered as their
  functions below say, radians for the angles."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.kepler :as kepler]
            [allgo.astro.states :as states]
            [allgo.numerics.differentiation :as diff]
            [allgo.numerics.linear :as lin]))

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

(defn- from-classical [mu [a e i raan argp M]]
  (flat (kepler/elements->state mu {:a a :e e :i i :raan raan :argp argp :M M})))

(defn cartesian->classical
  "Carry the Cartesian covariance `P` of state `s` into classical
  elements [a e i raan argp M]. Singular where the elements are: circular
  or equatorial orbits."
  ([P s] (cartesian->classical mu P s))
  ([mu P s] (lin/congruence (diff/jacobian #(classical-vector mu %) (flat s) {:angles #{2 3 4 5}}) P)))

(defn classical->cartesian
  "Carry the covariance `P` of classical elements [a e i raan argp M]
  `el` into Cartesian."
  ([P el] (classical->cartesian mu P el))
  ([mu P el] (lin/congruence (diff/jacobian #(from-classical mu %) (vec el) {:steps (mapv #(* 1e-6 (max 1e-3 (abs %))) el)}) P)))

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
