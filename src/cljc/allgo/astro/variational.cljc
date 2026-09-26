(ns allgo.astro.variational
  "Linearized dynamics: how a trajectory responds to a change in where it
  started (Montenbruck & Gill chapter 7).

  Orbit determination never solves the real problem. It guesses an initial
  state, propagates it, compares against observations and corrects -- and
  the correction requires knowing how the trajectory at time t would change
  had the initial state been slightly different. That derivative is the
  state transition matrix, and it is the bridge between chapter 4's
  propagation and chapter 8's estimation.

    dPhi/dt = A(t) Phi,   Phi(t0, t0) = I

  where A is the Jacobian of the equations of motion. It is integrated
  alongside the trajectory, 36 extra states for the six real ones, because
  A depends on where the satellite actually is."
  (:require [clojure.math :as math]))

;; ------------------------------------------------------------ small matrices

(defn identity-matrix [n]
  (mapv (fn [i] (mapv (fn [j] (if (= i j) 1.0 0.0)) (range n))) (range n)))

(defn mat-mul [a b]
  (let [m (count (first b))]
    (mapv (fn [row] (mapv (fn [j] (reduce + (map-indexed (fn [k v] (* v (nth (nth b k) j))) row)))
                          (range m)))
          a)))

(defn mat-vec [m v] (mapv (fn [row] (reduce + (map * row v))) m))

(defn det-3 [[[a b cc] [d e f] [g h i]]]
  (- (+ (* a e i) (* b f g) (* cc d h)) (+ (* cc e g) (* b d i) (* a f h))))

;; --------------------------------------------------------------- gradients

(defn two-body-gradient
  "The exact Jacobian of two-body acceleration with respect to position:

    da_i/dr_j = mu (3 r_i r_j / r^5 - delta_ij / r^3)

  Worth having in closed form as a check on the numerical version, since
  everything else in the force model has to be differentiated numerically
  and there is otherwise nothing to test that machinery against."
  [mu r]
  (let [r2 (reduce + (map * r r))
        rm (math/sqrt r2)
        r3 (* r2 rm)
        r5 (* r3 r2)]
    (mapv (fn [i]
            (mapv (fn [j]
                    (- (/ (* 3.0 mu (nth r i) (nth r j)) r5)
                       (if (= i j) (/ mu r3) 0.0)))
                  (range 3)))
          (range 3))))

(defn- perturb [v i h] (update (vec v) i + h))

(defn acceleration-gradients
  "Numerical Jacobians of an acceleration with respect to position and
  velocity, by central differences.

  Numerical rather than analytic because the alternative is differentiating
  a spherical harmonic recursion, a density table and a shadow function by
  hand -- and a mistake in any of them would be invisible, since a wrong
  Jacobian degrades convergence rather than breaking it. The step is scaled
  to the magnitude of what is being perturbed, so it works as well at
  geostationary radius as in low orbit."
  [accel t r v]
  (let [step (fn [x] (max 1e-4 (* 1e-7 (abs x))))
        col  (fn [f x i] (let [h (step (nth x i))]
                           (mapv #(/ % (* 2.0 h))
                                 (mapv - (f (perturb x i h)) (f (perturb x i (- h)))))))
        dr   (mapv (fn [i] (col (fn [r'] (accel t r' v)) r i)) (range 3))
        dv   (mapv (fn [i] (col (fn [v'] (accel t r v')) v i)) (range 3))]
    ;; the columns above are gradients; transpose into Jacobians
    {:d-dr (mapv (fn [i] (mapv (fn [j] (nth (nth dr j) i)) (range 3))) (range 3))
     :d-dv (mapv (fn [i] (mapv (fn [j] (nth (nth dv j) i)) (range 3))) (range 3))}))

(defn jacobian
  "The 6x6 matrix A = d(state-rate)/d(state).

    A = [  0    I  ]
        [ da/dr  da/dv ]

  The top half says position rate is velocity, exactly and always. All the
  physics is in the bottom half."
  [accel t r v]
  (let [{:keys [d-dr d-dv]} (acceleration-gradients accel t r v)]
    (vec (concat
          (mapv (fn [i] (vec (concat (repeat 3 0.0)
                                     (map #(if (= i %) 1.0 0.0) (range 3)))))
                (range 3))
          (mapv (fn [i] (vec (concat (nth d-dr i) (nth d-dv i)))) (range 3))))))

;; ------------------------------------------------ the variational equations

(defn pack
  "State and transition matrix as one 42-vector, which is what a general
  integrator can carry."
  [r v phi]
  (vec (concat r v (apply concat phi))))

(defn unpack [y]
  {:r   (subvec (vec y) 0 3)
   :v   (subvec (vec y) 3 6)
   :phi (mapv (fn [i] (subvec (vec y) (+ 6 (* i 6)) (+ 12 (* i 6)))) (range 6))})

(defn rhs
  "Right-hand side for the combined system: the equations of motion, and the
  variational equations riding alongside them.

  `accel` is `(fn [t r v])` returning acceleration."
  [accel]
  (fn [t y]
    (let [{:keys [r v phi]} (unpack y)
          a   (accel t r v)
          A   (jacobian accel t r v)
          dphi (mat-mul A phi)]
      (pack v a dphi))))

(defn initial
  "Start the combined system: the given state, and the identity, since at
  the initial epoch the trajectory is exactly as sensitive to its own start
  as it could possibly be."
  [r v]
  (pack r v (identity-matrix 6)))
