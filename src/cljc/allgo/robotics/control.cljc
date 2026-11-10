(ns allgo.robotics.control
  "Control: Russell and Norvig, *Artificial Intelligence: A Modern
  Approach*, 4th edition, chapter 26, which teaches it in the text.

  Planning says where to go; control makes the robot go there despite
  its own dynamics and the world pushing back. Two controllers:

    PID   the control a sum of three terms in the error: proportional,
          which pulls toward the target and, on its own, overshoots and
          oscillates; derivative, which damps that; and integral, which
          removes the steady error a constant disturbance leaves
    LQR   the optimal linear controller for linear dynamics x' = A x + B
          u and a quadratic cost, sum of x Q x + u R u: the gain from
          the discrete Riccati equation, iterated to its fixed point

  Both are steps. `pid-step` takes one error to one control and the
  controller's next memory; `simulate` iterates a plant and a controller
  together; `riccati-seq` is the cost-to-go one more step out, each step
  a matrix and the gain it implies."
  (:require [allgo.numerics.linear :as lin]
            [allgo.numerics.linear-systems :as linear-systems]))

;; ---------------------------------------------------------------------------
;; PID

(defn pid-step
  "One PID update for `error` (target minus measurement) after `dt`:
  `[u memory']`, the memory `{:integral :previous}` (nil to start). The
  derivative is taken of the error, and is zero on the first step."
  [{:keys [kp ki kd] :or {kp 0.0 ki 0.0 kd 0.0}} memory error dt]
  (let [{:keys [integral previous] :or {integral 0.0}} memory
        integral (+ integral (* error dt))
        derivative (if previous (/ (- error previous) dt) 0.0)]
    [(+ (* kp error) (* ki integral) (* kd derivative))
     {:integral integral :previous error}]))

(defn- simulate-step
  "A tick of the plant under PID control."
  [{:keys [plant measure target gains dt]} {:keys [t x memory]}]
  (let [[u memory] (pid-step gains memory (- (target t) (measure x)) dt)]
    {:t (+ t dt) :x (plant x u dt) :u u :memory memory}))

(defn simulate
  "A plant and controller run together: each state `{:t :x :u :memory}`.
  `plant` is (fn [x u dt]) -> x', `measure` (fn [x]) -> the measured
  value, `target` (fn [t]) -> where it should be, and `gains` the PID's."
  [system x0]
  (iterate (partial simulate-step system) {:t 0.0 :x x0 :u 0.0 :memory nil}))

(defn mass
  "A point mass on a line, `[position velocity]`, pushed by the control
  and by a constant `disturbance` force (default 0): semi-implicit Euler,
  which keeps an undamped oscillation from gaining or losing energy."
  ([] (mass 0.0))
  ([disturbance]
   (fn [[p v] u dt]
     (let [v (+ v (* dt (+ u disturbance)))]
       [(+ p (* dt v)) v]))))

;; ---------------------------------------------------------------------------
;; LQR

(defn- riccati-step
  "The cost-to-go one step further from the end, and its gain."
  [A B Q R {:keys [P]}]
  (let [At (lin/transpose A)
        BtP (lin/mat-mul (lin/transpose B) P)
        K (lin/mat-mul (lin/inverse (lin/mat-add R (lin/mat-mul BtP B))) (lin/mat-mul BtP A))
        P' (lin/mat-add Q (lin/mat-sub (lin/mat-mul (lin/mat-mul At P) A)
                                       (lin/mat-mul (lin/mat-mul (lin/mat-mul At P) B) K)))]
    {:P P' :K K}))

(defn riccati-seq
  "The discrete Riccati recursion for x' = A x + B u with cost x Q x + u R
  u, from the terminal cost P = Q: each state `{:P :K}`, the cost-to-go
  matrix one step further from the end and the gain K = (R + B P B)^-1 B P
  A it gives, the control u = -K x. It converges to the infinite-horizon
  LQR."
  [A B Q R]
  (iterate (partial riccati-step A B Q R) {:P Q :K nil}))

(defn lqr
  "The infinite-horizon LQR gain K, u = -K x: the Riccati recursion run
  until P stops changing by more than `tol` (default 1e-12), or
  `max-steps` (default 10000). `{:P :K :steps}`."
  ([A B Q R] (lqr A B Q R 1e-12 10000))
  ([A B Q R tol max-steps]
   (loop [[a b & more] (riccati-seq A B Q R) n 1]
     (let [change (reduce max 0.0 (map #(abs (- %1 %2)) (flatten (:P a)) (flatten (:P b))))]
       (if (or (< change tol) (>= n max-steps))
         (assoc b :steps n)
         (recur (cons b more) (inc n)))))))

(defn closed-loop
  "The states of x' = (A - B K) x from `x0`."
  [A B K x0]
  (let [M (lin/mat-sub A (lin/mat-mul B K))]
    (iterate (partial lin/mat-vec M) x0)))

;; ---------------------------------------------------------------------------
;; Iterative LQR

(defn- jacobian
  "The partials of `f` (vector to vector) at `x`, by central differences:
  a matrix, a row per output."
  [f x]
  (let [h 1e-6
        cols (mapv (fn [j] (let [up (f (update x j + h)) down (f (update x j - h))]
                             (mapv #(/ (- %1 %2) (* 2 h)) up down)))
                   (range (count x)))]
    (lin/transpose cols)))

(defn- quad-cost [Q x] (lin/dot x (lin/mat-vec Q x)))

(defn rollout
  "The states from `x0` under the controls `us`, by `dynamics`, (fn [x u])
  -> x'."
  [dynamics x0 us]
  (vec (reductions dynamics x0 us)))

(defn trajectory-cost
  "The cost of states `xs` and controls `us`: the sum of (x - goal) Q (x
  - goal) + u R u over the steps, and (x - goal) Qf (x - goal) at the
  end."
  [{:keys [goal Q R Qf]} xs us]
  (+ (reduce + (map (fn [x u] (+ (quad-cost Q (lin/sub x goal)) (quad-cost R u))) xs us))
     (quad-cost Qf (lin/sub (peek xs) goal))))

(defn- ilqr-iteration
  "One iteration of iterative LQR from the controls `us`: the dynamics
  linearized and the cost (already quadratic) expanded about their
  rollout, the LQR problem that gives solved backward for a feedforward
  change k and a feedback gain K at each step, and the new controls rolled
  out forward, the change shrunk by halves until the cost falls."
  [{:keys [dynamics x0 goal Q R Qf mu] :or {mu 1e-6} :as problem} {:keys [us cost]}]
  (let [xs (rollout dynamics x0 us)
        n (count us) m (count (first us))
        two #(lin/mat-scale % 2.0)
        [ks Ks] (loop [t (dec n) vx (lin/mat-vec (two Qf) (lin/sub (peek xs) goal)) vxx (two Qf) ks () Ks ()]
                  (if (neg? t)
                    [(vec ks) (vec Ks)]
                    (let [x (xs t) u (us t)
                          A (jacobian #(dynamics % u) x) B (jacobian #(dynamics x %) u)
                          At (lin/transpose A) Bt (lin/transpose B)
                          qx (lin/add (lin/mat-vec (two Q) (lin/sub x goal)) (lin/mat-vec At vx))
                          qu (lin/add (lin/mat-vec (two R) u) (lin/mat-vec Bt vx))
                          qxx (lin/mat-add (two Q) (lin/mat-mul (lin/mat-mul At vxx) A))
                          quu (lin/mat-add (lin/mat-add (two R) (lin/mat-mul (lin/mat-mul Bt vxx) B)) (lin/mat-scale (lin/eye m) mu))
                          qux (lin/mat-mul (lin/mat-mul Bt vxx) A)
                          inv (lin/inverse quu)
                          k (lin/scale (lin/mat-vec inv qu) -1.0)
                          K (lin/mat-scale (lin/mat-mul inv qux) -1.0)
                          Kt (lin/transpose K)
                          vx (lin/add (lin/add qx (lin/mat-vec Kt (lin/mat-vec quu k)))
                                      (lin/add (lin/mat-vec Kt qu) (lin/mat-vec (lin/transpose qux) k)))
                          vxx (lin/mat-add (lin/mat-add qxx (lin/mat-mul (lin/mat-mul Kt quu) K))
                                           (lin/mat-add (lin/mat-mul Kt qux) (lin/mat-mul (lin/transpose qux) K)))
                          vxx (lin/mat-scale (lin/mat-add vxx (lin/transpose vxx)) 0.5)]
                      (recur (dec t) vx vxx (cons k ks) (cons K Ks)))))
        forward (fn [alpha]
                  (let [[xs' us'] (reduce (fn [[xs' us'] t]
                                            (let [x (peek xs')
                                                  u (lin/add (lin/add (us t) (lin/scale (ks t) alpha))
                                                             (lin/mat-vec (Ks t) (lin/sub x (xs t))))]
                                              [(conj xs' (dynamics x u)) (conj us' u)]))
                                          [[x0] []] (range n))]
                    {:us us' :cost (trajectory-cost problem xs' us')}))
        tried (some (fn [alpha] (let [c (forward alpha)] (when (< (:cost c) cost) c)))
                    (take 12 (iterate #(* 0.5 %) 1.0)))]
    (or tried {:us us :cost cost :converged true})))

(defn- ilqr-step
  "An iteration of iLQR, until one fails to lower the cost."
  [problem st]
  (if (:converged st) st (ilqr-iteration problem st)))

(defn ilqr-seq
  "Iterative LQR (section 26.5.4) for `problem` -- `{:dynamics :x0 :goal
  :Q :R :Qf}`, the dynamics (fn [x u]) -> x', any smooth function, and
  the quadratic costs of `trajectory-cost` -- from the controls `us`, an
  iteration a step: each state `{:us :cost}`, the cost never rising,
  `:converged` once no step lowers it. For linear dynamics the expansion
  is exact, and the first iteration is the finite-horizon LQR."
  [problem us]
  (let [us (vec us)]
    (iterate (partial ilqr-step problem)
             {:us us :cost (trajectory-cost problem (rollout (:dynamics problem) (:x0 problem) us) us)})))

(defn- mpc-step
  "Plan from where the plant is, warm-started from the last plan shifted
  a step, and apply the first control."
  [plant model horizon iterations {:keys [x plan]}]
  (let [m (count (first (:R model)))
        warm (vec (take horizon (concat (rest plan) (repeat (or (peek plan) (vec (repeat m 0.0)))))))
        us (:us (nth (ilqr-seq (assoc model :x0 x) warm) iterations))]
    {:x (plant x (first us)) :u (first us) :plan us}))

(defn mpc-seq
  "Model predictive control (section 26.6): at every step, plan the next
  `horizon` controls by `iterations` of iLQR on the `model`'s dynamics
  from where the plant actually is -- warm-started from the last plan,
  shifted a step -- apply the first control to the `plant` (fn [x u]) ->
  x', which may differ from the model, and plan again. Each state `{:x
  :u :plan}`."
  [plant {:keys [x0] :as model} horizon iterations]
  (let [m (count (first (:R model)))]
    (iterate (partial mpc-step plant model horizon iterations)
             {:x x0 :u nil :plan (vec (repeat horizon (vec (repeat m 0.0))))})))

;; ---------------------------------------------------------------------------
;; System properties

(defn controllability-matrix
  "Kalman's controllability matrix of x' = A x + B u: [B, A B, A^2 B, ...,
  A^(n-1) B], side by side."
  [A B]
  (let [n (count A)
        blocks (take n (iterate (partial lin/mat-mul A) B))]
    (apply mapv (fn [& rows] (vec (apply concat rows))) blocks)))

(defn controllable?
  "Whether x' = A x + B u can be driven from any state to any other:
  LaValle, *Planning Algorithms*, Cambridge University Press, 2006,
  section 15.1.3 -- the controllability matrix has full rank."
  [A B]
  (= (count A) (linear-systems/rank (controllability-matrix A B))))

(defn lyapunov
  "For the discrete system x' = A x, the P with A^T P A - P = -`Q` (Q
  positive definite): the sum of (A^T)^k Q A^k, a step a term, stopped
  when the terms fall under `tol` (default 1e-12). Such a P exists, and
  V(x) = x P x is then a Lyapunov function -- it falls along every
  trajectory -- exactly when A is stable (section 15.1.2); nil if the
  sum has not settled after `max-terms` (default 10000)."
  ([A Q] (lyapunov A Q 1e-12 10000))
  ([A Q tol max-terms]
   (let [At (lin/transpose A)]
     (loop [P Q term Q k 1]
       (let [term (lin/mat-mul (lin/mat-mul At term) A)
             size (reduce max (map abs (flatten term)))]
         (cond (< size tol) (lin/mat-add P term)
               (or (>= k max-terms) (not (< size 1e12))) nil
               :else (recur (lin/mat-add P term) term (inc k))))))))
