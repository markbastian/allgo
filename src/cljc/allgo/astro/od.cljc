(ns allgo.astro.od
  "Orbit determination from a stream of observations (Vallado, chapter
  10; Tapley, Schutz and Born, *Statistical Orbit Determination*, 2004):
  the sequential-batch least squares, and the extended and unscented
  Kalman filters.

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
  (:require [allgo.astro.estimation :as est]
            [allgo.astro.variational :as var]
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

(defn lkf
  "The linearized Kalman filter (Vallado's algorithm 68): as `ekf`, but
  about a reference trajectory flown once from `[r v]` and never
  corrected, the filter estimating only the deviation from it -- carried
  by the reference's transition matrix and updated by the residuals from
  the reference, their partials taken there. Cheaper and, while the
  deviation stays small, as good; it is the batch least squares' first
  iteration made sequential. Returns `[{:t :r :v :P :dx :residuals} ...]`,
  the estimate the reference plus the deviation."
  [accel [r v] P t observations {:keys [q] :or {q 0.0}}]
  (->> observations
       (reductions
        (fn [{:keys [ref-r ref-v dx P t]} {tk :t :keys [z model sigma]}]
          (let [{r1 :r v1 :v :keys [phi]} (propagate-with-stm accel t [ref-r ref-v] tk)
                dx1 (lin/mat-vec phi dx)
                P1 (lin/mat-add (lin/mat-mul (lin/mat-mul phi P) (lin/transpose phi)) (process-noise q (abs (- tk t))))
                x1 (flat r1 v1)
                H (measurement-partials model tk x1)
                residuals (mapv - z (model tk r1 v1))
                {dx' :x P' :P} (reduce (fn [{:keys [x P]} k]
                                         (select-keys (est/kalman-update x P (H k) (- (residuals k) (lin/dot (H k) x))
                                                                         (* (sigma k) (sigma k)))
                                                      [:x :P]))
                                       {:x dx1 :P P1}
                                       (range (count z)))
                x (lin/add x1 dx')]
            {:t tk :ref-r r1 :ref-v v1 :dx dx' :P P' :r (subvec x 0 3) :v (subvec x 3 6) :residuals residuals}))
        {:t t :ref-r r :ref-v v :dx (vec (repeat 6 0.0)) :P P})
       rest
       (mapv #(dissoc % :ref-r :ref-v))))

;; ----------------------------------------------------------- unscented

(defn- propagate-state
  "The state `dt` from `[r v]` under `accel`, without the transition
  matrix."
  [accel t0 x t1]
  (if (== t0 t1)
    x
    (let [first-order (fn [t y] (into (subvec y 3 6) (accel t (subvec y 0 3) (subvec y 3 6))))
          integ (rk/integrator rk/dopri54 first-order t0 x (* 0.01 (- t1 t0)) {:tol-abs 1e-11 :tol-rel 1e-11})]
      (:y (core/step-until integ t1)))))

(defn ukf
  "The unscented Kalman filter (Julier and Uhlmann, \"A new extension of
  the Kalman filter to nonlinear systems\", 1997, in the scaled form of
  Wan and van der Merwe, 2000), taking the arguments and returning what
  `ekf` does. No transition matrix and no partials: the covariance is
  represented by 2n + 1 sigma points, x and x +/- the columns of
  sqrt((n + lambda) P), each flown through the dynamics and the
  measurement model in full, and the mean and covariance rebuilt from
  them with the weights W0m = lambda/(n + lambda), W0c = W0m + 1 -
  alpha^2 + beta, the rest 1/(2(n + lambda)); lambda = alpha^2 (n +
  kappa) - n. Options `:q` as for `ekf`, and `:alpha` (default 1e-3),
  `:beta` (2, right for Gaussian errors) and `:kappa` (0)."
  [accel [r v] P t observations {:keys [q alpha beta kappa] :or {q 0.0 alpha 1e-3 beta 2.0 kappa 0.0}}]
  (let [n 6
        lambda (- (* alpha alpha (+ n kappa)) n)
        wm (into [(/ lambda (+ n lambda))] (repeat (* 2 n) (/ 1.0 (* 2.0 (+ n lambda)))))
        wc (assoc wm 0 (+ (wm 0) (- 1.0 (* alpha alpha)) beta))
        mean (fn [xs] (reduce lin/add (map lin/scale xs wm)))
        cov (fn [xs mx ys my] (reduce lin/mat-add (map (fn [x y w] (let [dx (lin/sub x mx) dy (lin/sub y my)]
                                                                     (mapv (fn [a] (mapv #(* w a %) dy)) dx)))
                                                       xs ys wc)))]
    (->> observations
         (reductions
          (fn [{:keys [r v P t]} {tk :t :keys [z model sigma]}]
            (let [x (flat r v)
                  L (lin/cholesky (lin/mat-scale P (+ n lambda)))
                  cols (lin/transpose L)
                  points (into [x] (concat (map #(lin/add x %) cols) (map #(lin/sub x %) cols)))
                  ;; the time update: each point flown to tk
                  flown (mapv #(propagate-state accel t % tk) points)
                  x1 (mean flown)
                  P1 (lin/mat-add (cov flown x1 flown x1) (process-noise q (abs (- tk t))))
                  ;; the measurement update, from points redrawn about x1
                  L1 (lin/cholesky (lin/mat-scale P1 (+ n lambda)))
                  cols1 (lin/transpose L1)
                  pts (into [x1] (concat (map #(lin/add x1 %) cols1) (map #(lin/sub x1 %) cols1)))
                  zs (mapv #(vec (model tk (subvec % 0 3) (subvec % 3 6))) pts)
                  zbar (mean zs)
                  noise (vec (for [i (range (count z))]
                               (vec (for [j (range (count z))] (if (= i j) (* (sigma i) (sigma i)) 0.0)))))
                  Pzz (lin/mat-add (cov zs zbar zs zbar) noise)
                  Pxz (cov pts x1 zs zbar)
                  K (lin/mat-mul Pxz (lin/inverse-general Pzz))
                  residuals (mapv - z zbar)
                  xk (lin/add x1 (lin/mat-vec K residuals))
                  Pk (lin/mat-sub P1 (lin/mat-mul (lin/mat-mul K Pzz) (lin/transpose K)))
                  ;; kept symmetric against rounding
                  Pk (lin/mat-scale (lin/mat-add Pk (lin/transpose Pk)) 0.5)]
              {:t tk :r (subvec xk 0 3) :v (subvec xk 3 6) :P Pk :residuals residuals}))
          {:t t :r r :v v :P P})
         rest
         vec)))
