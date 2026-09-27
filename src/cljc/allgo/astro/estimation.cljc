(ns allgo.astro.estimation
  "Orbit determination: recovering a trajectory from tracking data
  (Montenbruck & Gill chapter 8).

  Everything else in this package builds to here. A truth orbit is never
  known -- what exists is a pile of noisy ranges and angles, and the task is
  the state that best explains them. Two approaches, with different
  characters:

  Batch least squares takes all the observations at once and iterates on the
  initial state until the residuals stop improving. It is what you use when
  the data has already been collected.

  A Kalman filter processes observations one at a time, carrying a covariance
  that says how well it currently knows the state. It is what you use when
  the data is still arriving, and it never revisits a measurement.

  Both need the state transition matrix of chapter 7 to relate a correction
  at the epoch to its effect at the time of an observation."
  (:require [allgo.numerics.linear :as lin]
            [clojure.math :as math]))

;; --------------------------------------------------------- batch estimation

(defn normal-equations
  "Accumulate the weighted normal equations from a set of observations.

  Each observation contributes H^T W H to the information matrix and
  H^T W r to the right-hand side, where H is the partial of the measurement
  with respect to the epoch state -- already including the state transition
  matrix that carries it back from the observation time."
  [rows n]
  (reduce (fn [{:keys [N b]} {:keys [H residual weight]}]
            (let [w (or weight 1.0)]
              {:N (lin/mat-add N (lin/mat-scale (mapv (fn [hi] (mapv #(* hi %) H)) H) w))
               :b (lin/add b (lin/scale H (* w residual)))}))
          {:N (lin/mat-scale (lin/eye n) 0.0) :b (vec (repeat n 0.0))}
          rows))

(defn solve-batch
  "One Gauss-Newton step: the correction that minimizes the weighted sum of
  squared residuals, given their partials.

  Returns nil when the normal matrix is not positive definite, which means
  the observations leave some direction of the state undetermined -- range
  alone from a single station cannot fix an orbit, however much of it there
  is."
  [rows n]
  (let [{:keys [N b]} (normal-equations rows n)]
    (when-let [dx (lin/cholesky-solve N b)]
      {:correction dx :information N :covariance (lin/inverse N)})))

(defn rms
  "Root-mean-square residual, the number that says whether an iteration
  helped."
  [rows]
  (math/sqrt (/ (reduce + (map #(let [r (:residual %)] (* r r)) rows))
                (max 1 (count rows)))))

;; -------------------------------------------------------- sequential filter

(defn kalman-predict
  "Time update: carry the state and its covariance forward.

    P = Phi P Phi^T + Q

  The process noise Q is what stops a filter becoming so confident that it
  ignores new data. Without it the covariance shrinks forever and the filter
  quietly stops learning, which is the classic way to make one diverge."
  [x P phi Q]
  {:x (lin/mat-vec phi x)
   :P (lin/mat-add (lin/mat-mul (lin/mat-mul phi P) (lin/transpose phi)) Q)})

(defn kalman-update
  "Measurement update for a scalar observation.

  The covariance is updated in Joseph form, which stays symmetric and
  positive definite under roundoff where the shorter (I - KH)P does not.
  That matters over a long arc: a covariance that loses positive
  definiteness produces a negative variance, and the filter is finished."
  [x P H z-residual R]
  (let [PHt  (lin/mat-vec P H)                      ; P is symmetric
        S    (+ (lin/dot H PHt) R)
        K     (lin/scale PHt (/ 1.0 S))
        x'    (lin/add x (lin/scale K z-residual))
        n     (count x)
        IKH   (lin/mat-sub (lin/eye n) (mapv (fn [ki] (mapv #(* ki %) H)) K))
        joseph (lin/mat-add (lin/mat-mul (lin/mat-mul IKH P) (lin/transpose IKH))
                            (lin/mat-scale (mapv (fn [ki] (mapv #(* ki % R) K)) K) 1.0))]
    {:x x' :P joseph :gain K :innovation-variance S}))

;; ------------------------------------------------- orthogonal least squares
;;
;; The normal equations are convenient and lossy. Forming A^T A squares the
;; condition number, so a problem solvable in double precision at cond(A) =
;; 1e8 becomes singular at cond(A^T A) = 1e16 -- and the information is
;; destroyed by the forming, not by the problem. Reducing A directly by
;; orthogonal transformations never squares anything, which is why M&G
;; prefers it and why a square-root information filter carries R rather
;; than P.

(defn- householder-column
  "Apply one Householder reflection, zeroing below the diagonal of column
  `k`. Returns the updated rows, each row being its A-part and b-part."
  [rows k n]
  (let [m     (count rows)
        col   (mapv #(nth (first %) k) rows)
        tail  (subvec col k)
        norm  (lin/length tail)]
    (if (< norm 1e-300)
      rows
      (let [x0    (nth col k)
            alpha (if (neg? x0) norm (- norm))
            v     (assoc (vec (repeat m 0.0)) k (- x0 alpha))
            v     (reduce (fn [v i] (assoc v i (nth col i))) v (range (inc k) m))
            vv    (lin/length-squared v)]
        (if (< vv 1e-300)
          rows
          (let [;; H = I - 2 v v^T / (v^T v), applied to every column at once
                dot-col (fn [get-el]
                          (reduce + (map-indexed (fn [i r] (* (nth v i) (get-el r))) rows)))
                factors (mapv (fn [j] (/ (* 2.0 (dot-col #(nth (first %) j))) vv)) (range n))
                bfactor (/ (* 2.0 (dot-col second)) vv)]
            (mapv (fn [i [a b]]
                    [(mapv (fn [j aj] (- aj (* (nth v i) (nth factors j)))) (range n) a)
                     (- b (* (nth v i) bfactor))])
                  (range m) rows)))))))

(defn qr-least-squares
  "Solve min ||A x - b|| by Householder reflections, without ever forming
  A^T A.

  `rows` are the same maps `solve-batch` takes. Weights enter as a scaling
  of each row by the square root of the weight, which is the same weighted
  problem written so that the solver never sees a weight at all."
  [rows n]
  (let [scaled (mapv (fn [{:keys [H residual weight]}]
                       (let [s (math/sqrt (or weight 1.0))]
                         [(mapv #(* % s) H) (* residual s)]))
                     rows)]
    (when (>= (count scaled) n)
      (let [reduced (reduce (fn [rws k] (householder-column rws k n)) scaled (range n))
            R (mapv #(first (nth reduced %)) (range n))
            y (mapv #(second (nth reduced %)) (range n))
            ;; Rank test against the size of the factor, not an absolute
            ;; floor: a diagonal of 1e-15 is negligible beside a norm of 1
            ;; and perfectly significant beside a norm of 1e-20.
            scale (apply max 1e-300 (map (fn [i] (abs (nth (nth R i) i))) (range n)))
            tol   (* scale n 1e-14)
            back  (fn [rhs]
                    (reduce (fn [x i]
                              (if (nil? x)
                                nil
                                (let [sum (reduce + (map (fn [j] (* (nth (nth R i) j) (nth x j)))
                                                         (range (inc i) n)))
                                      d   (nth (nth R i) i)]
                                  (if (< (abs d) tol)
                                    nil
                                    (assoc x i (/ (- (nth rhs i) sum) d))))))
                            (vec (repeat n 0.0))
                            (reverse (range n))))]
        (when-let [x (back y)]
          ;; R^-1 column by column, each column solving R c = e_i
          (let [cols  (mapv (fn [i] (back (mapv #(if (= i %) 1.0 0.0) (range n)))) (range n))]
            (when (every? some? cols)
              (let [Rinv (lin/transpose cols)]
                {:correction x
                 :r-factor R
                 ;; (A^T A)^-1 = (R^T R)^-1 = R^-1 R^-T. The order matters:
                 ;; the transpose on the right, not the left.
                 :covariance (lin/mat-mul Rinv (lin/transpose Rinv))}))))))))

(defn solve
  "Solve the batch problem, by orthogonal reduction unless asked otherwise.

  `:method :normal` uses the normal equations, which are cheaper and quite
  adequate for a well-conditioned arc. `:method :qr` is the default because
  the cost of being wrong is silent: an ill-conditioned problem does not
  announce itself, it returns fewer correct digits than it appears to."
  ([rows n] (solve rows n {}))
  ([rows n {:keys [method] :or {method :qr}}]
   (case method
     :qr     (qr-least-squares rows n)
     :normal (solve-batch rows n))))
