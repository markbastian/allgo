(ns allgo.astro.od
  "Orbit determination from a stream of observations (Vallado, chapter
  10; Tapley, Schutz and Born, *Statistical Orbit Determination*, 2004):
  the sequential-batch least squares, and the extended Kalman filter.

  Both linearize about a reference trajectory and carry the state
  transition matrix along it (`allgo.astro.variational`). The sequential
  batch takes the observations in batches, each batch's solution the
  a-priori estimate for the next: the information of every batch
  accumulates, and the last solution is the one the whole batch at once
  would give. The extended filter takes them one epoch at a time and
  moves its reference to each new estimate, so the linearization stays
  good however far the first guess was out; process noise, as state-noise
  compensation, keeps it from growing so sure of itself that it stops
  listening.

  States are `[r v]`, km and km/s; `accel` is `(fn [t r v])`, the
  acceleration, two-body included."
  (:require [allgo.astro.variational :as var]
            [allgo.numerics.core :as core]
            [allgo.numerics.differentiation :as diff]
            [allgo.numerics.linear :as lin]
            [allgo.numerics.rk :as rk]))

(defn- flat [r v] (vec (concat r v)))

(defn propagate-with-stm
  "The state at `t1` of the trajectory through `[r v]` at `t0`, and the
  state transition matrix from `t0` to it: `{:r :v :phi}`, the variational
  equations integrated alongside the motion by Dormand-Prince at tight
  tolerance."
  [accel t0 [r v] t1]
  (if (== t0 t1)
    {:r r :v v :phi (lin/eye 6)}
    (let [integ (rk/integrator rk/dopri54 (var/rhs accel) t0 (var/initial r v) (* 0.01 (- t1 t0))
                               {:tol-abs 1e-11 :tol-rel 1e-11})]
      (var/unpack (:y (core/step-until integ t1))))))

;; ------------------------------------------------------- sequential batch

(defn sequential-batch
  "The sequential-batch least squares: `batches`, each a collection of
  rows as `allgo.astro.estimation/normal-equations` takes them -- `{:H
  :residual :weight}`, the partial mapped to the common epoch -- folded
  one after another into the information matrix and vector, starting from
  the a-priori `prior`, `{:x :P}` (the deviation and its covariance) or
  nil for none. Each batch's estimate is

    x = (P^-1 + H^T W H)^-1 (P^-1 x-bar + H^T W y)

  with the previous batch's the a-priori. Returns the estimate after each
  batch, `[{:x :P} ...]`."
  [prior batches]
  (let [n (count (:H (ffirst (filter seq batches))))
        [N0 b0] (if prior
                  (let [Pinv (lin/inverse-general (:P prior))]
                    [Pinv (lin/mat-vec Pinv (:x prior))])
                  [(lin/mat-scale (lin/eye n) 0.0) (vec (repeat n 0.0))])]
    (->> batches
         (reductions (fn [[N b] rows]
                       [(reduce (fn [N {:keys [H weight]}]
                                  (lin/mat-add N (lin/mat-scale (mapv (fn [hi] (mapv #(* hi %) H)) H) (or weight 1.0))))
                                N rows)
                        (reduce (fn [b {:keys [H residual weight]}]
                                  (lin/add b (lin/scale H (* (or weight 1.0) residual))))
                                b rows)])
                     [N0 b0])
         rest
         (mapv (fn [[N b]] (let [P (lin/inverse-general N)] {:x (lin/mat-vec P b) :P P}))))))

;; ------------------------------------------------------------- batch

(defn- measurement-partials
  "The partials of `(model t r v)` in the state at `x`, by central
  differences: one row per measured value."
  [model t x]
  (let [predict (fn [x] (vec (model t (subvec x 0 3) (subvec x 3 6))))]
    (diff/jacobian predict x {:steps (vec (concat (repeat 3 1e-4) (repeat 3 1e-7))) :richardson? false})))

(defn batch-rows
  "The rows of the batch least squares about the epoch state `[r v]` at
  `t0`: for each of `observations` (as `ekf` takes them, in time order),
  each measured value's residual, weight and partial mapped back to the
  epoch through the transition matrix -- chained from one observation to
  the next, so the trajectory is flown once."
  [accel t0 [r v] observations]
  (loop [obs observations t t0 r r v v phi (lin/eye 6) out []]
    (if-let [[{tk :t :keys [z model sigma]} & more] (seq obs)]
      (let [{r1 :r v1 :v step :phi} (propagate-with-stm accel t [r v] tk)
            phi1 (lin/mat-mul step phi)
            x1 (flat r1 v1)
            H (lin/mat-mul (measurement-partials model tk x1) phi1)
            res (mapv - z (model tk r1 v1))]
        (recur more tk r1 v1 phi1
               (into out (map (fn [k] {:H (H k) :residual (res k) :weight (/ 1.0 (* (sigma k) (sigma k)))})
                              (range (count z))))))
      out)))

(defn batch
  "The batch least squares: the epoch state at `t0` that best fits
  `observations`, iterated from `[r v]` until the correction is below
  `:tol` in each component (default 1e-9): `{:r :v :P :iterations}`, P
  the covariance of the epoch state. With `:batches` n, the rows are
  taken n at a time by `sequential-batch`, to the same answer."
  [accel t0 [r v] observations {:keys [tol max-iter batches] :or {tol 1e-9 max-iter 15}}]
  (loop [r r v v k 0]
    (let [rows (batch-rows accel t0 [r v] observations)
          {:keys [x P]} (peek (sequential-batch nil (if batches (partition-all batches rows) [rows])))
          r' (lin/add r (subvec x 0 3)) v' (lin/add v (subvec x 3 6))]
      (if (or (< (apply max (map abs x)) tol) (>= k max-iter))
        {:r r' :v v' :P P :iterations (inc k)}
        (recur r' v' (inc k))))))

;; ----------------------------------------------------------- extended

(defn process-noise
  "State-noise compensation: the covariance a white acceleration of
  spectral density `q` (km^2/s^3, the same on each axis) adds over `dt`
  seconds -- [[q dt^3/3 I, q dt^2/2 I], [q dt^2/2 I, q dt I]]."
  [q dt]
  (let [a (* q (/ (* dt dt dt) 3.0)) b (* q (/ (* dt dt) 2.0)) c (* q dt)]
    (vec (for [i (range 6)]
           (vec (for [j (range 6)]
                  (cond (not= (mod i 3) (mod j 3)) 0.0
                        (and (< i 3) (< j 3)) a
                        (and (>= i 3) (>= j 3)) c
                        :else b)))))))

(defn ekf
  "The extended Kalman filter from state `[r v]` with covariance `P` at
  time `t` through `observations`, each `{:t :z :model :sigma}`: at time
  t, measured values `z`, the function `(model t r v)` that predicts them,
  and their standard deviations. Between observations the estimate and
  covariance are carried by `accel` and the transition matrix, the process
  noise of spectral density `:q` (default 0) added; at each, the
  measurements are taken one component at a time, the partials of the
  model by central differences, the covariance updated in Joseph form.
  Returns the estimate after each observation, `[{:t :r :v :P
  :residuals} ...]`, the residuals before the update."
  [accel [r v] P t observations {:keys [q] :or {q 0.0}}]
  (->> observations
       (reductions
        (fn [{:keys [r v P t]} {tk :t :keys [z model sigma]}]
          (let [{r1 :r v1 :v :keys [phi]} (propagate-with-stm accel t [r v] tk)
                P1 (lin/mat-add (lin/mat-mul (lin/mat-mul phi P) (lin/transpose phi)) (process-noise q (abs (- tk t))))
                x1 (flat r1 v1)
                H (lin/transpose (measurement-partials model tk x1))
                residuals (mapv - z (model tk r1 v1))
                ;; one scalar update per component, the reference fixed at x1
                [dx P'] (reduce (fn [[dx P] k]
                                  (let [h (mapv #(nth % k) H)
                                        pht (lin/mat-vec P h)
                                        s (+ (lin/dot h pht) (* (sigma k) (sigma k)))
                                        gain (lin/scale pht (/ 1.0 s))
                                        innov (- (residuals k) (lin/dot h dx))
                                        ikh (lin/mat-sub (lin/eye 6) (mapv (fn [gi] (mapv #(* gi %) h)) gain))]
                                    [(lin/add dx (lin/scale gain innov))
                                     (lin/mat-add (lin/mat-mul (lin/mat-mul ikh P) (lin/transpose ikh))
                                                  (mapv (fn [gi] (mapv #(* gi % (sigma k) (sigma k)) gain)) gain))]))
                                [(vec (repeat 6 0.0)) P1]
                                (range (count z)))
                x (lin/add x1 dx)]
            {:t tk :r (subvec x 0 3) :v (subvec x 3 6) :P P' :residuals residuals}))
        {:t t :r r :v v :P P})
       rest
       vec))
