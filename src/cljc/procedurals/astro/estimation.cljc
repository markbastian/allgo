(ns procedurals.astro.estimation
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
  (:require [clojure.math :as math]))

;; -------------------------------------------------------- linear algebra

(defn transpose [m] (apply mapv vector m))

(defn mat-mul [a b]
  (let [bt (transpose b)]
    (mapv (fn [row] (mapv (fn [col] (reduce + (map * row col))) bt)) a)))

(defn mat-vec [m v] (mapv (fn [row] (reduce + (map * row v))) m))

(defn mat-add [a b] (mapv (fn [r s] (mapv + r s)) a b))
(defn mat-sub [a b] (mapv (fn [r s] (mapv - r s)) a b))
(defn mat-scale [m s] (mapv (fn [r] (mapv #(* % s) r)) m))
(defn eye [n] (mapv (fn [i] (mapv #(if (= i %) 1.0 0.0) (range n))) (range n)))

(defn cholesky
  "Lower-triangular L with L L^T = A, for symmetric positive definite A.

  Returns nil when A is not positive definite. For a normal matrix that is
  not a numerical mishap but a statement about the data: some direction of
  the state is not determined by the observations, and no amount of solving
  will invent it."
  [A]
  (let [n (count A)]
    (loop [i 0 L (vec (repeat n (vec (repeat n 0.0))))]
      (cond
        (nil? L) nil
        (= i n)  L
        :else
        (recur (inc i)
               (loop [j 0 L L]
                 (cond
                   (nil? L) nil
                   (> j i)  L
                   :else
                   (let [s (reduce + (map * (take j (nth L i)) (take j (nth L j))))
                         v (- (nth (nth A i) j) s)]
                     (if (= i j)
                       (if (<= v 0.0)
                         nil
                         (recur (inc j) (assoc-in L [i j] (math/sqrt v))))
                       (recur (inc j) (assoc-in L [i j] (/ v (nth (nth L j) j)))))))))))))

(defn cholesky-solve
  "Solve A x = b for symmetric positive definite A, by forward then back
  substitution through the Cholesky factor."
  [A b]
  (when-let [L (cholesky A)]
    (let [n (count A)
          y (reduce (fn [y i]
                      (conj y (/ (- (nth b i) (reduce + (map * (take i (nth L i)) y)))
                                 (nth (nth L i) i))))
                    [] (range n))
          x (reduce (fn [x i]
                      (let [s (reduce + (map (fn [k] (* (nth (nth L k) i) (get x k 0.0)))
                                             (range (inc i) n)))]
                        (assoc x i (/ (- (nth y i) s) (nth (nth L i) i)))))
                    (vec (repeat n 0.0)) (reverse (range n)))]
      x)))

(defn inverse
  "Inverse of a symmetric positive definite matrix, column by column."
  [A]
  (let [n (count A)]
    (when (cholesky A)
      (transpose (mapv (fn [i] (cholesky-solve A (mapv #(if (= i %) 1.0 0.0) (range n))))
                       (range n))))))

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
              {:N (mat-add N (mat-scale (mapv (fn [hi] (mapv #(* hi %) H)) H) w))
               :b (mapv + b (mapv #(* % w residual) H))}))
          {:N (mat-scale (eye n) 0.0) :b (vec (repeat n 0.0))}
          rows))

(defn solve-batch
  "One Gauss-Newton step: the correction that minimises the weighted sum of
  squared residuals, given their partials.

  Returns nil when the normal matrix is not positive definite, which means
  the observations leave some direction of the state undetermined -- range
  alone from a single station cannot fix an orbit, however much of it there
  is."
  [rows n]
  (let [{:keys [N b]} (normal-equations rows n)]
    (when-let [dx (cholesky-solve N b)]
      {:correction dx :information N :covariance (inverse N)})))

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
  {:x (mat-vec phi x)
   :P (mat-add (mat-mul (mat-mul phi P) (transpose phi)) Q)})

(defn kalman-update
  "Measurement update for a scalar observation.

  The covariance is updated in Joseph form, which stays symmetric and
  positive definite under roundoff where the shorter (I - KH)P does not.
  That matters over a long arc: a covariance that loses positive
  definiteness produces a negative variance, and the filter is finished."
  [x P H z-residual R]
  (let [PHt  (mat-vec P H)                      ; P is symmetric
        S    (+ (reduce + (map * H PHt)) R)
        K     (mapv #(/ % S) PHt)
        x'    (mapv + x (mapv #(* % z-residual) K))
        n     (count x)
        IKH   (mat-sub (eye n) (mapv (fn [ki] (mapv #(* ki %) H)) K))
        joseph (mat-add (mat-mul (mat-mul IKH P) (transpose IKH))
                        (mat-scale (mapv (fn [ki] (mapv #(* ki % R) K)) K) 1.0))]
    {:x x' :P joseph :gain K :innovation-variance S}))
