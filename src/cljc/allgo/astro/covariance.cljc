(ns allgo.astro.covariance
  "An orbit's uncertainty carried between the ways of writing it
  (Vallado, *Fundamentals of Astrodynamics and Applications*, chapter
  10's covariance transformations): Cartesian, classical and equinoctial
  elements, flight elements, and the satellite's own radial/along-track/
  cross-track and velocity-aligned frames.

  To first order a covariance P in one set becomes J P J^T in another, J
  the Jacobian of the conversion. For the local frames J is a rotation and
  exact. For the element sets it is taken here numerically, by central
  differences refined by Richardson extrapolation, which gets it to ten
  digits for any conversion the library has -- the angles among the
  outputs wrapped, so a difference straddling 2 pi is not a jump.

  States are `[r v]`, km and km/s; element vectors are ordered as their
  functions below say, radians for the angles."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.kepler :as kepler]
            [allgo.astro.states :as states]
            [allgo.math :as am]
            [allgo.numerics.linear :as lin]))

(def ^:private mu c/GM-earth)

;; ------------------------------------------------------------ machinery

(defn jacobian
  "The Jacobian of `f`, a function from a vector to a vector, at `x`: each
  column a central difference at step `h_j` and half of it, combined by
  Richardson extrapolation to cancel the h^2 error. `:steps` gives the
  h_j (default 1e-5 of each component, or 1e-5); `:angles` the output
  indices to difference modulo 2 pi."
  ([f x] (jacobian f x {}))
  ([f x {:keys [steps angles] :or {angles #{}}}]
   (let [n (count x)
         steps (or steps (mapv #(* 1e-5 (max 1.0 (abs %))) x))
         diff (fn [a b] (vec (map-indexed (fn [k [p q]] (let [d (- p q)] (if (angles k) (am/wrap-angle d) d)))
                                          (map vector a b))))
         column (fn [j]
                  (let [d (fn [h] (lin/scale (diff (f (update x j + h)) (f (update x j - h))) (/ 1.0 (* 2.0 h))))
                        h (steps j)]
                    (lin/scale (lin/sub (lin/scale (d (* 0.5 h)) 4.0) (d h)) (/ 1.0 3.0))))]
     (lin/transpose (mapv column (range n))))))

(defn transform
  "J P J^T: the covariance `P` carried through the Jacobian `J`."
  [J P]
  (lin/mat-mul (lin/mat-mul J P) (lin/transpose J)))

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
  ([mu P s] (transform (jacobian #(classical-vector mu %) (flat s) {:angles #{2 3 4 5}}) P)))

(defn classical->cartesian
  "Carry the covariance `P` of classical elements [a e i raan argp M]
  `el` into Cartesian."
  ([P el] (classical->cartesian mu P el))
  ([mu P el] (transform (jacobian #(from-classical mu %) (vec el) {:steps (mapv #(* 1e-6 (max 1e-3 (abs %))) el)}) P)))

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
  ([mu P s] (transform (jacobian #(equinoctial-vector mu %) (flat s) {:angles #{5}}) P)))

(defn equinoctial->cartesian
  "Carry the covariance `P` of equinoctial elements `eq`, a map as
  `states/state->equinoctial` gives, into Cartesian."
  ([P eq] (equinoctial->cartesian mu P eq))
  ([mu P {:keys [fr] :as eq}]
   (let [x (mapv eq [:a :af :ag :chi :psi :meanlon])
         f (fn [[a af ag chi psi meanlon]]
             (flat (states/equinoctial->state mu {:a a :af af :ag ag :chi chi :psi psi :meanlon meanlon :fr fr})))]
     (transform (jacobian f x {:steps (mapv #(* 1e-6 (max 1e-3 (abs %))) x)}) P))))

(defn cartesian->flight
  "Carry the Cartesian (inertial) covariance `P` of state `s` into flight
  elements [rm vm latgc lon fpa az] at TT `mjd-tt` and UT1 `mjd-ut1`."
  [P s mjd-tt mjd-ut1 eop]
  (transform (jacobian #(vec (take 6 (states/state->flight (state %) mjd-tt mjd-ut1 eop))) (flat s)
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
  (transform (block (states/rsw s)) P))

(defn rsw->cartesian [P s] (transform (lin/transpose (block (states/rsw s))) P))

(defn cartesian->ntw
  "The same into the velocity-aligned N, T, W axes."
  [P s]
  (transform (block (states/ntw s)) P))

(defn ntw->cartesian [P s] (transform (lin/transpose (block (states/ntw s))) P))
