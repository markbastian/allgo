(ns allgo.astro.cr3bp
  "The circular restricted three-body problem (Vallado, *Fundamentals of
  Astrodynamics and Applications*, chapter 2; Szebehely, *Theory of
  Orbits*, 1967): a body too small to pull on two others that circle their
  common center of mass -- a spacecraft between the Earth and the Moon.

  In units where the primaries' separation, their total mass and their
  angular rate are all 1, and in the frame that turns with them, the
  primaries sit still on the x axis: the larger, of mass 1 - mu, at -mu,
  and the smaller, of mass mu, at 1 - mu. There the motion obeys

    x'' - 2y' = dU/dx,   y'' + 2x' = dU/dy,   z'' = dU/dz,
    U = (x^2 + y^2)/2 + (1 - mu)/r1 + mu/r2,

  which conserves the Jacobi constant C = 2U - v^2 -- the one integral the
  problem has, and the fence it builds: where 2U < C the body cannot go.
  U has five stationary points, the Lagrange points: three on the axis
  and two at the tips of equilateral triangles."
  (:require [allgo.geometry.vec3 :as v3]
            [clojure.math :as math]))

(defn mass-ratio
  "mu for primaries of parameters `mu1` (the larger) and `mu2`."
  [mu1 mu2]
  (/ mu2 (+ mu1 mu2)))

(defn- distances [mu [x y z]]
  [(math/sqrt (+ (* (+ x mu) (+ x mu)) (* y y) (* z z)))
   (math/sqrt (+ (* (- x (- 1.0 mu)) (- x (- 1.0 mu))) (* y y) (* z z)))])

(defn potential
  "U, the effective potential of the rotating frame."
  [mu [x y :as r]]
  (let [[r1 r2] (distances mu r)]
    (+ (* 0.5 (+ (* x x) (* y y))) (/ (- 1.0 mu) r1) (/ mu r2))))

(defn gradient
  "dU/dx, dU/dy, dU/dz."
  [mu [x y z :as r]]
  (let [[r1 r2] (distances mu r)
        k1 (/ (- 1.0 mu) (* r1 r1 r1))
        k2 (/ mu (* r2 r2 r2))]
    [(- x (* k1 (+ x mu)) (* k2 (- x (- 1.0 mu))))
     (- y (* k1 y) (* k2 y))
     (- (+ (* k1 z) (* k2 z)))]))

(defn jacobi
  "The Jacobi constant 2U - v^2 of the state `[r v]`."
  [mu [r v]]
  (- (* 2.0 (potential mu r)) (v3/dot v v)))

(defn derivative
  "d/dt of the state [x y z vx vy vz], for the integrators."
  [mu]
  (fn [_ [x y z vx vy vz]]
    (let [[ux uy uz] (gradient mu [x y z])]
      [vx vy vz (+ ux (* 2.0 vy)) (- uy (* 2.0 vx)) uz])))

;; ---------------------------------------------------- the Lagrange points

(defn- axis-root
  "The x between `lo` and `hi` where dU/dx vanishes on the axis, by
  bisection: one sign change in each of the three intervals the primaries
  cut the axis into."
  [mu lo hi]
  (let [f #(first (gradient mu [% 0.0 0.0]))]
    (loop [lo lo hi hi k 0]
      (let [mid (* 0.5 (+ lo hi))]
        (if (or (> k 200) (<= (- hi lo) 1e-15))
          mid
          (if (= (neg? (f lo)) (neg? (f mid))) (recur mid hi (inc k)) (recur lo mid (inc k))))))))

(defn lagrange-points
  "The five Lagrange points: `{:L1 :L2 :L3 :L4 :L5}`. L1 between the
  primaries, L2 beyond the smaller, L3 beyond the larger, and L4 and L5
  leading and trailing the smaller by 60 degrees."
  [mu]
  (let [eps 1e-9
        p2 (- 1.0 mu)]
    {:L1 [(axis-root mu (+ (- mu) eps) (- p2 eps)) 0.0 0.0]
     :L2 [(axis-root mu (+ p2 eps) 2.0) 0.0 0.0]
     :L3 [(axis-root mu -2.0 (- (- mu) eps)) 0.0 0.0]
     :L4 [(- 0.5 mu) (/ (math/sqrt 3.0) 2.0) 0.0]
     :L5 [(- 0.5 mu) (- (/ (math/sqrt 3.0) 2.0)) 0.0]}))

(defn- hessian-in-plane
  "Uxx, Uyy, Uxy at `r`, in the plane."
  [mu [x y :as r]]
  (let [[r1 r2] (distances mu r)
        a (+ x mu) b (- x (- 1.0 mu))
        k1 (/ (- 1.0 mu) (math/pow r1 3)) k2 (/ mu (math/pow r2 3))
        m1 (/ (* 3.0 (- 1.0 mu)) (math/pow r1 5)) m2 (/ (* 3.0 mu) (math/pow r2 5))]
    [(+ (- 1.0 k1 k2) (* m1 a a) (* m2 b b))
     (+ (- 1.0 k1 k2) (* m1 y y) (* m2 y y))
     (+ (* m1 a y) (* m2 b y))]))

(defn planar-eigenvalues
  "The eigenvalues, as `[re im]` pairs, of the motion linearized about the
  equilibrium `r` in the plane: the roots of
  lambda^4 + (4 - Uxx - Uyy) lambda^2 + Uxx Uyy - Uxy^2 = 0. An
  equilibrium is linearly stable when all four are imaginary."
  [mu r]
  (let [[uxx uyy uxy] (hessian-in-plane mu r)
        b (- 4.0 uxx uyy)
        cc (- (* uxx uyy) (* uxy uxy))
        disc (- (* b b) (* 4.0 cc))
        ;; lambda^2 = s, each s gives +-sqrt(s)
        roots (if (>= disc 0.0)
                (let [d (math/sqrt disc)] [[(* 0.5 (+ (- b) d)) 0.0] [(* 0.5 (- (- b) d)) 0.0]])
                (let [d (math/sqrt (- disc))] [[(* -0.5 b) (* 0.5 d)] [(* -0.5 b) (* -0.5 d)]]))
        csqrt (fn [[re im]]
                (let [m (math/hypot re im)
                      a (math/sqrt (* 0.5 (+ m re)))
                      b (* (if (neg? im) -1.0 1.0) (math/sqrt (* 0.5 (- m re))))]
                  [a b]))]
    (vec (mapcat (fn [s] (let [[a b] (csqrt s)] [[a b] [(- a) (- b)]])) roots))))

(defn stable?
  "Whether the equilibrium `r` is linearly stable in the plane."
  [mu r]
  (every? (fn [[re _]] (< (abs re) 1e-12)) (planar-eigenvalues mu r)))

(def routh-critical
  "The mass ratio above which L4 and L5 are unstable: (1 - sqrt(23/27))/2,
  Routh's, 0.0385 -- the Sun and Jupiter, the Earth and the Moon, are far
  below it, which is why Trojans gather there."
  (* 0.5 (- 1.0 (math/sqrt (/ 23.0 27.0)))))

(defn forbidden?
  "Whether a body of Jacobi constant `C` is barred from `r`: 2U(r) < C,
  where it would need an imaginary speed."
  [mu C r]
  (< (* 2.0 (potential mu r)) C))
