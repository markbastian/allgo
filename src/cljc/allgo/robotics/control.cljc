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
  (:require [allgo.numerics.linear :as lin]))

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

(defn simulate
  "A plant and controller run together: each state `{:t :x :u :memory}`.
  `plant` is (fn [x u dt]) -> x', `measure` (fn [x]) -> the measured
  value, `target` (fn [t]) -> where it should be, and `gains` the PID's."
  [{:keys [plant measure target gains dt]} x0]
  (iterate (fn [{:keys [t x memory]}]
             (let [[u memory] (pid-step gains memory (- (target t) (measure x)) dt)]
               {:t (+ t dt) :x (plant x u dt) :u u :memory memory}))
           {:t 0.0 :x x0 :u 0.0 :memory nil}))

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

(defn riccati-seq
  "The discrete Riccati recursion for x' = A x + B u with cost x Q x + u R
  u, from the terminal cost P = Q: each state `{:P :K}`, the cost-to-go
  matrix one step further from the end and the gain K = (R + B P B)^-1 B P
  A it gives, the control u = -K x. It converges to the infinite-horizon
  LQR."
  [A B Q R]
  (let [At (lin/transpose A) Bt (lin/transpose B)]
    (iterate (fn [{:keys [P]}]
               (let [BtP (lin/mat-mul Bt P)
                     K (lin/mat-mul (lin/inverse (lin/mat-add R (lin/mat-mul BtP B))) (lin/mat-mul BtP A))
                     P' (lin/mat-add Q (lin/mat-sub (lin/mat-mul (lin/mat-mul At P) A)
                                                    (lin/mat-mul (lin/mat-mul (lin/mat-mul At P) B) K)))]
                 {:P P' :K K}))
             {:P Q :K nil})))

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
    (iterate #(lin/mat-vec M %) x0)))
