(ns allgo.probability.temporal
  "Reasoning over time: Russell and Norvig, *Artificial Intelligence: A
  Modern Approach*, 4th edition, chapter 14, as steps.

    filtering       the forward message, a step an observation
    smoothing       FORWARD-BACKWARD (figure 14.4)
    fixed-lag       FIXED-LAG-SMOOTHING (figure 14.6): online smoothing
    smoothing       `d` steps back, a step an observation, the backward
                    message kept as a matrix and updated in constant time
    Viterbi         the most likely sequence of states, which the chapter
                    gives without a figure
    particle        PARTICLE-FILTERING (figure 14.17): N samples moved by
    filtering       the transition model, weighted by the evidence, and
                    resampled, a step an observation
    Kalman filter   the chapter's linear-Gaussian filter, also without a
                    figure: predict through the transition model, correct
                    by the observation

  A hidden Markov model is `{:prior [...] :transition [[...]] :sensor
  {e [...]}}`: P(X0), T[i][j] = P(X_t = j | X_t-1 = i), and for each
  observation the vector P(e | X = i). `umbrella` is the chapter's
  umbrella world as aima-python (MIT) states it."
  (:require [allgo.numerics.linear :as lin]
            [allgo.numerics.linear-systems :as ls]
            [allgo.probability.bayes :as bayes]))

(def umbrella
  "The umbrella world: rain persists with probability 0.7, and the
  director carries an umbrella on 90% of rainy days and 20% of dry ones.
  State 0 is rain, 1 no rain."
  {:prior [0.5 0.5]
   :transition [[0.7 0.3] [0.3 0.7]]
   :sensor {true [0.9 0.2] false [0.1 0.8]}})

(defn- predict
  "The distribution one step on: sum over i of f(i) T[i][j]."
  [{:keys [transition]} f]
  (let [n (count f)] (mapv (fn [j] (reduce + (map (fn [i] (* (f i) (get-in transition [i j]))) (range n)))) (range n))))

(defn forward
  "FORWARD: the filtered distribution after observing `e`, from `f`."
  [hmm f e]
  (bayes/normalize (mapv * ((:sensor hmm) e) (predict hmm f))))

(defn backward
  "BACKWARD: the backward message one step earlier, b'(i) = sum over j of
  P(e | j) b(j) T[i][j]."
  [{:keys [transition sensor]} b e]
  (let [n (count b) o (sensor e)]
    (mapv (fn [i] (reduce + (map (fn [j] (* (o j) (b j) (get-in transition [i j]))) (range n)))) (range n))))

(defn filtering-seq
  "The filtered distribution after each observation in `evidence`, prior
  first."
  [hmm evidence]
  (reductions #(forward hmm %1 %2) (:prior hmm) evidence))

(defn forward-backward
  "FORWARD-BACKWARD: the smoothed distribution of each step 1 to t given
  all of `evidence`."
  [hmm evidence]
  (let [fv (vec (filtering-seq hmm evidence))
        t (count evidence)]
    (loop [i t b (vec (repeat (count (:prior hmm)) 1.0)) sv ()]
      (if (zero? i)
        (vec sv)
        (recur (dec i) (backward hmm b (nth evidence (dec i)))
               (cons (bayes/normalize (mapv * (fv i) b)) sv))))))

(defn- diag [v] (vec (map-indexed (fn [i x] (assoc (vec (repeat (count v) 0.0)) i x)) v)))

(defn fixed-lag-smoothing-seq
  "FIXED-LAG-SMOOTHING over `evidence`: for each observation, P(X_t-d |
  e_1:t) once `t > d`, nil before. The backward message is kept as the
  matrix B = product over the last d steps of T O, updated each step by
  B <- O_t-d^-1 T^-1 B T O_t -- so each step costs the same however long
  the lag; the forward message runs d steps behind."
  [{:keys [prior transition sensor] :as hmm} evidence d]
  (let [T transition
        T-inv (ls/inverse T)
        n (count prior)]
    (map :smoothed
         (rest (reductions
                (fn [{:keys [t f B window]} e]
                  (let [O (diag (sensor e))
                        window (conj window e)
                        [f B window] (if (> t d)
                                       (let [old (first window)
                                             f (forward hmm f old)
                                             window (vec (rest window))
                                             O-old-inv (lin/spd-inverse (diag (sensor old)))]
                                         [f (lin/mat-mul (lin/mat-mul (lin/mat-mul O-old-inv T-inv) B) (lin/mat-mul T O)) window])
                                       [f (lin/mat-mul B (lin/mat-mul T O)) window])
                        t' (inc t)]
                    {:t t' :f f :B B :window window
                     :smoothed (when (> t' (inc d))
                                 (bayes/normalize (mapv * f (lin/mat-vec B (vec (repeat n 1.0))))))}))
                {:t 1 :f prior :B (lin/eye n) :window []}
                evidence)))))

(defn viterbi
  "The most likely sequence of states given `evidence`, and at each step
  the probability of the most likely path ending in the state it chose
  there, normalized as the chapter's figure is: `{:path :probabilities}`."
  [{:keys [transition sensor] :as hmm} evidence]
  (let [n (count (:prior hmm))
        m1 (forward hmm (:prior hmm) (first evidence))
        [ms back] (reduce (fn [[ms back] e]
                            (let [m (peek ms)
                                  o (sensor e)
                                  best (mapv (fn [j] (apply max-key (fn [i] (* (m i) (get-in transition [i j]))) (range n))) (range n))
                                  m' (mapv (fn [j] (* (o j) (m (best j)) (get-in transition [(best j) j]))) (range n))]
                              [(conj ms m') (conj back best)]))
                          [[m1] []] (rest evidence))
        last-best (apply max-key (peek ms) (range n))
        path (reduce (fn [path best] (cons (best (first path)) path)) (list last-best) (rseq back))]
    {:path (vec path)
     :probabilities (mapv (fn [m s] (m s)) ms path)}))

(defn particle-filtering-seq
  "PARTICLE-FILTERING: the `n` particles after each observation -- each
  moved by the transition model, weighted by P(e | state), and the set
  resampled in proportion to the weights."
  [{:keys [prior transition sensor]} evidence n rng]
  (let [draw (fn [dist] (let [r (double (rng))]
                          (loop [i 0 acc 0.0] (let [acc (+ acc (dist i))]
                                                (if (or (< r acc) (= i (dec (count dist)))) i (recur (inc i) acc))))))
        resample (fn [particles weights]
                   (let [total (reduce + weights)]
                     (if (pos? total)
                       (let [dist (mapv #(/ % total) weights)]
                         (vec (repeatedly (count particles) #(particles (draw dist)))))
                       particles)))]
    (reductions (fn [particles e]
                  (let [moved (mapv #(draw (transition %)) particles)]
                    (resample moved (mapv #((sensor e) %) moved))))
                (vec (repeatedly n #(draw prior)))
                evidence)))

(defn particle-estimate
  "The distribution particles stand for."
  [particles n-states]
  (let [c (frequencies particles)] (mapv #(/ (double (get c % 0)) (count particles)) (range n-states))))

;; ---------------------------------------------------------------------------
;; Kalman filtering

(defn kalman-step
  "One step of the Kalman filter: the Gaussian `[mean cov]` predicted
  through x' = F x + noise of covariance Q, and corrected by the
  observation z = H x + noise of covariance R."
  [{:keys [F Q H R]} [mean cov] z]
  (let [mean- (lin/mat-vec F mean)
        cov- (lin/mat-add (lin/mat-mul (lin/mat-mul F cov) (lin/transpose F)) Q)
        S (lin/mat-add (lin/mat-mul (lin/mat-mul H cov-) (lin/transpose H)) R)
        K (lin/mat-mul (lin/mat-mul cov- (lin/transpose H)) (lin/spd-inverse S))
        innovation (lin/sub z (lin/mat-vec H mean-))
        mean' (lin/add mean- (lin/mat-vec K innovation))
        cov' (lin/mat-mul (lin/mat-sub (lin/eye (count mean)) (lin/mat-mul K H)) cov-)]
    [mean' cov']))

(defn kalman-seq
  "The Kalman filter's `[mean cov]` after each observation, from the prior
  `[mean0 cov0]`."
  [model prior observations]
  (rest (reductions #(kalman-step model %1 %2) prior observations)))

;; The same filter a scalar measurement at a time, as orbit determination
;; and the robots' EKFs use it (moved here from `allgo.astro.estimation`).

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
  definiteness produces a negative variance, and the filter is finished.

  For one row H the Joseph form needs no matrix products: with A = (I -
  K H) P = P - K (P H)^T, it is A (I - K H)^T + R K K^T = A - (A H) K^T +
  R K K^T, all outer products -- quadratic in the state's size rather than
  cubic, which is what lets SLAM's long state update a scalar at a time.
  The result is symmetrized, as roundoff in the products would not quite
  keep it so."
  [x P H z-residual R]
  (let [outer (fn [u v] (mapv (fn [ui] (mapv #(* ui %) v)) u))
        PHt  (lin/mat-vec P H)                      ; P is symmetric
        S    (+ (lin/dot H PHt) R)
        K     (lin/scale PHt (/ 1.0 S))
        x'    (lin/add x (lin/scale K z-residual))
        A     (lin/mat-sub P (outer K PHt))
        joseph (lin/mat-add (lin/mat-sub A (outer (lin/mat-vec A H) K)) (outer K (lin/scale K R)))
        joseph (lin/mat-scale (lin/mat-add joseph (lin/transpose joseph)) 0.5)]
    {:x x' :P joseph :gain K :innovation-variance S}))

(defn linear-kalman
  "The Kalman filter for a linear system (Vallado's algorithm 67): from
  the estimate `x0` with covariance `P0`, through `steps`, each `{:phi :Q
  :H :z :R}` -- the transition matrix from the last step and the process
  noise it adds, then the rows of the measurement matrix, the measured
  values and their variances, taken one at a time. Returns the estimate
  after each step, `[{:x :P} ...]`."
  [x0 P0 steps]
  (->> steps
       (reductions
        (fn [{:keys [x P]} {:keys [phi Q H z R]}]
          (let [n (count x)
                predicted (kalman-predict x P phi (or Q (lin/mat-scale (lin/eye n) 0.0)))]
            (reduce (fn [{:keys [x P]} k]
                      (let [h (H k)]
                        (select-keys (kalman-update x P h (- (z k) (lin/dot h x)) (R k)) [:x :P])))
                    predicted
                    (range (count z)))))
        {:x x0 :P P0})
       rest
       vec))
