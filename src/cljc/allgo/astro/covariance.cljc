(ns allgo.astro.covariance
  "An orbit's uncertainty carried between the ways of writing it
  (Vallado, *Fundamentals of Astrodynamics and Applications*, chapter
  10's covariance transformations): Cartesian, classical and equinoctial
  elements, flight elements, and the satellite's own radial/along-track/
  cross-track and velocity-aligned frames.

  To first order a covariance P in one set becomes J P J^T in another, J
  the Jacobian of the conversion. For the local frames J is a rotation and
  exact; for the classical and equinoctial elements it is derived
  analytically (`classical-partials`, `equinoctial-partials`); for the
  flight elements it is taken numerically
  (`allgo.numerics.differentiation`), to ten digits, the angles among the
  outputs wrapped.

  States are `[r v]`, km and km/s; element vectors are ordered as their
  functions below say, radians for the angles."
  (:require [allgo.astro.constants :as c]
            [allgo.astro.kepler :as kepler]
            [allgo.astro.states :as states]
            [allgo.geometry.rotation :as rot]
            [allgo.geometry.vec3 :as v3]
            [allgo.numerics.differentiation :as diff]
            [allgo.numerics.linear :as lin]
            [allgo.numerics.linear-systems :as ls]
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
  ([mu P s] (lin/congruence (ls/inverse (classical-partials mu (classical-vector mu (flat s)))) P)))

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

(defn equinoctial-partials
  "The Jacobian of the state [x y z vx vy vz] with respect to the
  equinoctial elements [a af ag chi psi meanlon] of `eq`, a map as
  `states/state->equinoctial` gives, analytically and directly (after
  Broucke and Cefola, \"On the equinoctial orbit elements\", Celestial
  Mechanics 5, 1972), one column per element. With k = af, h = ag,
  p = chi, q = psi, I = fr and F the eccentric longitude,

    lambda = F - k sin F + h cos F
    r = X f + Y g,  v = X' f + Y' g
    X = a ((1 - h^2 b) cos F + h k b sin F - k)
    Y = a ((1 - k^2 b) sin F + h k b cos F - h),   b = 1/(1 + sqrt(1 - h^2 - k^2))
    f = (1 - p^2 + q^2, 2pq, -2Ip)/(1 + p^2 + q^2)
    g = (2Ipq, (1 + p^2 - q^2) I, 2q)/(1 + p^2 + q^2)

  and so

    a       r/a and -v/2a, as for the classical elements
    k, h    the in-plane coordinates' partials at fixed lambda, F moving
            with them: dF/dk = (a/r) sin F, dF/dh = -(a/r) cos F
    p, q    the frame's own partials
    lambda  v/n and -mu r/(n r^3), the motion itself

  none singular for circular or equatorial orbits."
  ([eq] (equinoctial-partials mu eq))
  ([mu {:keys [a af ag chi psi meanlon fr] :or {fr 1.0}}]
   (let [k af h ag p chi q psi I fr
         n (math/sqrt (/ mu (* a a a)))
         F (loop [F meanlon j 0]
             (let [dF (/ (- (+ (- F (* k (math/sin F))) (* h (math/cos F))) meanlon)
                         (- 1.0 (* k (math/cos F)) (* h (math/sin F))))]
               (if (or (< (abs dF) 1e-15) (> j 50)) (- F dF) (recur (- F dF) (inc j)))))
         c (math/cos F) s (math/sin F)
         eta (math/sqrt (- 1.0 (* h h) (* k k)))
         b (/ 1.0 (+ 1.0 eta))
         bh (/ (* b b h) eta) bk (/ (* b b k) eta)
         D (- 1.0 (* k c) (* h s))                      ; r/a
         ;; X/a, Y/a and their partials in F, h and k at fixed F
         x (+ (* (- 1.0 (* h h b)) c) (* h k b s) (- k))
         y (+ (* (- 1.0 (* k k b)) s) (* h k b c) (- h))
         Px (- (* h k b c) (* (- 1.0 (* h h b)) s))       ; dx/dF
         Py (- (* (- 1.0 (* k k b)) c) (* h k b s))       ; dy/dF
         PxF (- (+ (* (- 1.0 (* h h b)) c) (* h k b s)))
         PyF (- (+ (* (- 1.0 (* k k b)) s) (* h k b c)))
         xh (+ (* -1.0 (+ (* 2.0 h b) (* h h bh)) c) (* (+ (* k b) (* h k bh)) s))
         xk (+ (* -1.0 h h bk c) (* (+ (* h b) (* h k bk)) s) -1.0)
         yh (+ (* -1.0 k k bh s) (* (+ (* k b) (* h k bh)) c) -1.0)
         yk (+ (* -1.0 (+ (* 2.0 k b) (* k k bk)) s) (* (+ (* h b) (* h k bk)) c))
         Pxh (+ (* (+ (* 2.0 h b) (* h h bh)) s) (* (+ (* k b) (* h k bh)) c))
         Pxk (+ (* h h bk s) (* (+ (* h b) (* h k bk)) c))
         Pyh (- (* -1.0 k k bh c) (* (+ (* k b) (* h k bh)) s))
         Pyk (- (* -1.0 (+ (* 2.0 k b) (* k k bk)) c) (* (+ (* h b) (* h k bk)) s))
         DF (- (* k s) (* h c)) Dh (- s) Dk (- c)
         Fh (/ (- c) D) Fk (/ s D)
         ;; X, Y, X', Y' in an element e at fixed lambda
         in-plane (fn [xe ye Pxe Pye De Fe]
                    (let [Dt (+ De (* DF Fe))]
                      [(* a (+ xe (* Px Fe))) (* a (+ ye (* Py Fe)))
                       (* n a (- (/ (+ Pxe (* PxF Fe)) D) (/ (* Px Dt) (* D D))))
                       (* n a (- (/ (+ Pye (* PyF Fe)) D) (/ (* Py Dt) (* D D))))]))
         s2 (+ 1.0 (* p p) (* q q))
         fhat (v3/scale [(+ (- 1.0 (* p p)) (* q q)) (* 2.0 p q) (* -2.0 I p)] (/ 1.0 s2))
         ghat (v3/scale [(* 2.0 I p q) (* I (- (+ 1.0 (* p p)) (* q q))) (* 2.0 q)] (/ 1.0 s2))
         d-frame (fn [dN w frame] (v3/scale (v3/sub dN (v3/scale frame (* 2.0 w))) (/ 1.0 s2)))
         fp (d-frame [(* -2.0 p) (* 2.0 q) (* -2.0 I)] p fhat)
         fq (d-frame [(* 2.0 q) (* 2.0 p) 0.0] q fhat)
         gp (d-frame [(* 2.0 I q) (* 2.0 I p) 0.0] p ghat)
         gq (d-frame [(* 2.0 I p) (* -2.0 I q) 2.0] q ghat)
         X (* a x) Y (* a y) Xd (/ (* n a Px) D) Yd (/ (* n a Py) D)
         in (fn [u w] (v3/add (v3/scale fhat u) (v3/scale ghat w)))
         r (in X Y) v (in Xd Yd)
         rm (v3/length r)
         col (fn [dr dv] (vec (concat dr dv)))
         plane-col (fn [[dX dY dXd dYd]] (col (in dX dY) (in dXd dYd)))
         frame-col (fn [df dg] (col (v3/add (v3/scale df X) (v3/scale dg Y))
                                    (v3/add (v3/scale df Xd) (v3/scale dg Yd))))]
     (lin/transpose
      [(col (v3/scale r (/ 1.0 a)) (v3/scale v (/ -0.5 a)))
       (plane-col (in-plane xk yk Pxk Pyk Dk Fk))
       (plane-col (in-plane xh yh Pxh Pyh Dh Fh))
       (frame-col fp gp)
       (frame-col fq gq)
       (col (v3/scale v (/ 1.0 n)) (v3/scale r (/ (- mu) (* n rm rm rm))))]))))

(defn cartesian->equinoctial
  "Carry the Cartesian covariance `P` of state `s` into equinoctial
  elements [a af ag chi psi meanlon], retrograde orbits in their own
  set as `states/state->equinoctial` chooses, through the inverse of
  `equinoctial-partials`."
  ([P s] (cartesian->equinoctial mu P s))
  ([mu P s] (lin/congruence (ls/inverse (equinoctial-partials mu (states/state->equinoctial mu s))) P)))

(defn equinoctial->cartesian
  "Carry the covariance `P` of equinoctial elements `eq`, a map as
  `states/state->equinoctial` gives, into Cartesian, through
  `equinoctial-partials`."
  ([P eq] (equinoctial->cartesian mu P eq))
  ([mu P eq] (lin/congruence (equinoctial-partials mu eq) P)))

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
