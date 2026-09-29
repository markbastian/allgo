(ns allgo.numerics.ode
  "Ordinary differential equations beyond the Runge-Kutta tableaux of
  `allgo.numerics.rk` (Chapra and Canale, *Numerical Methods for
  Engineers*, chapters 25-26): Heun's method with its corrector iterated,
  the second-order Taylor method, adaptive fourth-order Runge-Kutta by
  step halving, the implicit (backward) Euler method for stiff systems,
  and the fixed-step multistep methods -- the non-self-starting Heun
  method, Milne's method and the fourth-order Adams method -- each a
  predictor and an iterated corrector, with the optional modifiers that
  fold each step's error estimate back in.

  `f` is `(f t y)` with `y` a vector; results are `[[t y] ...]`, the start
  first."
  (:require [allgo.numerics.linear :as lin]
            [allgo.numerics.linear-systems :as ls]
            [allgo.numerics.rk :as rk]
            [clojure.math :as math]))

(defn- add [a b] (mapv + a b))
(defn- scale [a s] (mapv #(* s %) a))
(defn- lin-comb [pairs] (reduce (fn [acc [c v]] (add acc (scale v c))) (vec (repeat (count (second (first pairs))) 0.0)) pairs))
(defn- rel-change [a b] (apply max (map (fn [x y] (if (zero? y) (abs (- y x)) (abs (/ (- y x) y)))) a b)))

(defn- iterate-corrector
  "Iterates y <- (corrector y) from `y0` until it settles to `tol`
  relative or `max-iter` times."
  [corrector y0 tol max-iter]
  (loop [y y0 k 0]
    (let [y' (corrector y)]
      (if (or (>= (inc k) max-iter) (< (rel-change y y') tol)) y' (recur y' (inc k))))))

;; ------------------------------------------------- one-step, continued

(defn heun
  "Heun's method, `n` steps of `h` from `y0` at `t0`: Euler's step
  predicts, and the trapezoid of the slopes at both ends corrects -- the
  corrector iterated `:iterations` times (default 1, the plain method) or
  until it settles to `:tol` (Chapra and Canale 25.2.1)."
  ([f t0 y0 h n] (heun f t0 y0 h n {}))
  ([f t0 y0 h n {:keys [iterations tol] :or {iterations 1 tol 0.0}}]
   (vec (reductions (fn [[t y] _]
                      (let [s0 (f t y)
                            t1 (+ t h)
                            y1 (iterate-corrector #(add y (scale (add s0 (f t1 %)) (* 0.5 h)))
                                                  (add y (scale s0 h)) tol iterations)]
                        [t1 y1]))
                    [t0 (vec y0)] (range n)))))

(defn taylor-2
  "The second-order Taylor-series method, `n` steps of `h`: y + h y' +
  h^2/2 y'', `ddf` the second derivative of y as a function of t and y
  -- f's total derivative, f_t + f_y f, which the user supplies (Chapra
  and Canale 25.1.3)."
  [f ddf t0 y0 h n]
  (vec (reductions (fn [[t y] _] [(+ t h) (add y (add (scale (f t y) h) (scale (ddf t y) (* 0.5 h h))))])
                   [t0 (vec y0)] (range n))))

(defn- rk4-step [f t y h]
  (first (rk/tableau-step rk/rk4 (:c rk/rk4) f t y h)))

(defn step-halving
  "Adaptive fourth-order Runge-Kutta by step halving (Chapra and Canale
  25.5.1): each step of h taken once and as two halves, their difference
  the error estimate, the halves' result corrected by it over 15 (to
  fifth order); the step shrunk when the error exceeds `:tol` relative to
  the solution's scale, grown after easy steps, by the fifth-root law
  (25.5.3). From `t0` to `t1`, starting at `:h` (default a hundredth of
  the span). `[[t y] ...]`."
  ([f t0 y0 t1] (step-halving f t0 y0 t1 {}))
  ([f t0 y0 t1 {:keys [tol h] :or {tol 1e-8}}]
   (loop [t t0 y (vec y0) h (or h (* 0.01 (- t1 t0))) out [[t0 (vec y0)]]]
     (if (>= t (- t1 (* 1e-12 (abs t1))))
       out
       (let [h (min h (- t1 t))
             whole (rk4-step f t y h)
             half (rk4-step f (+ t (* 0.5 h)) (rk4-step f t y (* 0.5 h)) (* 0.5 h))
             delta (mapv - half whole)
             scale (mapv #(max (abs %1) (abs (* h %2)) 1e-30) y (f t y))
             err (/ (apply max (map #(abs (/ %1 %2)) delta scale)) tol)]
         (if (<= err 1.0)
           (let [y' (add half (mapv #(/ % 15.0) delta))]
             (recur (+ t h) y' (* h (min 4.0 (* 0.9 (math/pow (max err 1e-10) -0.2)))) (conj out [(+ t h) y'])))
           (recur t y (* h (max 0.1 (* 0.9 (math/pow err -0.25)))) out)))))))

(defn- numeric-jacobian [g y]
  (lin/transpose (mapv (fn [j] (let [hj (* 1e-7 (max 1.0 (abs (y j))))]
                                 (mapv #(/ (- %1 %2) (* 2.0 hj)) (g (update y j + hj)) (g (update y j - hj)))))
                       (range (count y)))))

(defn backward-euler
  "The implicit (backward) Euler method, `n` steps of `h`: y_i+1 = y_i +
  h f(t_i+1, y_i+1), solved for y_i+1 by Newton's method each step (the
  Jacobian of f by central differences). First order, but stable at any
  step for a decaying solution -- the remedy for stiff systems, where the
  explicit methods must crawl at the fastest transient's pace (Chapra and
  Canale 26.1)."
  [f t0 y0 h n]
  (vec (reductions (fn [[t y] _]
                     (let [t1 (+ t h)
                           g (fn [z] (mapv - z y (scale (f t1 z) h)))
                           y1 (loop [z (add y (scale (f t y) h)) k 0]
                                (let [dz (ls/gauss (numeric-jacobian g z) (mapv - (g z)))]
                                  (if (or (nil? dz) (> k 50))
                                    z
                                    (let [z' (add z dz)]
                                      (if (< (rel-change z z') 1e-14) z' (recur z' (inc k)))))))]
                       [t1 y1]))
                   [t0 (vec y0)] (range n))))

;; ------------------------------------------------------------- multistep

(def ^:private multistep-methods
  "Each multistep method: how many past points it needs, its predictor
  and corrector (functions of the step, the past values ys and slopes fs,
  newest last, and for the corrector the new slope), and its modifiers'
  constants -- the predictor's p (added times the last step's corrector
  minus predictor) and the corrector's c (subtracted times this step's)
  (Chapra and Canale 26.2, Box 26.1)."
  {:heun {:history 2
          :predict (fn [h ys fs] (add (ys 0) (scale (fs 1) (* 2.0 h))))
          :correct (fn [h ys fs f+] (add (ys 1) (scale (add (fs 1) f+) (* 0.5 h))))
          :modify [0.8 0.2]}
   :milne {:history 4
           :predict (fn [h ys fs] (add (ys 0) (scale (lin-comb [[2.0 (fs 3)] [-1.0 (fs 2)] [2.0 (fs 1)]]) (/ (* 4.0 h) 3.0))))
           :correct (fn [h ys fs f+] (add (ys 2) (scale (lin-comb [[1.0 (fs 2)] [4.0 (fs 3)] [1.0 f+]]) (/ h 3.0))))
           :modify [(/ 28.0 29.0) (/ 1.0 29.0)]}
   :adams {:history 4
           :predict (fn [h ys fs] (add (ys 3) (scale (lin-comb [[55.0 (fs 3)] [-59.0 (fs 2)] [37.0 (fs 1)] [-9.0 (fs 0)]]) (/ h 24.0))))
           :correct (fn [h ys fs f+] (add (ys 3) (scale (lin-comb [[9.0 f+] [19.0 (fs 3)] [-5.0 (fs 2)] [1.0 (fs 1)]]) (/ h 24.0))))
           :modify [(/ 251.0 270.0) (/ 19.0 270.0)]}})

(defn multistep
  "The multistep `method` -- `:heun` (the non-self-starting Heun method,
  second order), `:milne` or `:adams` (fourth order) -- `n` steps of `h`
  from `y0` at `t0`. Each step predicts from the past points and corrects
  by iterating the corrector to `:tol` (default 1e-12) or `:iterations`
  times (default 100); with `:modify? true` the predictor and the
  corrector are each adjusted by their error estimates. The past points
  before t0 come from `:history`, `[[t y] ...]` oldest first, or are made
  by fourth-order Runge-Kutta steps backward from t0."
  ([method f t0 y0 h n] (multistep method f t0 y0 h n {}))
  ([method f t0 y0 h n {:keys [history modify? tol iterations] :or {tol 1e-12 iterations 100}}]
   (let [{:keys [predict correct modify] need :history} (multistep-methods method)
         [pm cm] modify
         past (or history
                  (vec (reverse (take (dec need)
                                      (rest (iterate (fn [[t y]] [(- t h) (rk4-step f t y (- h))]) [t0 (vec y0)]))))))
         start (conj (mapv (fn [[t y]] [t (vec y)]) past) [t0 (vec y0)])]
     (loop [pts (vec (take-last need start)) out [[t0 (vec y0)]] k 0 last-pc nil]
       (if (= k n)
         out
         (let [ys (mapv second pts)
               ts (mapv first pts)
               fs (mapv f ts ys)
               t1 (+ (peek ts) h)
               p0 (predict h ys fs)
               ;; the predictor modifier uses the last step's corrector-predictor gap
               p (if (and modify? last-pc) (add p0 (scale last-pc pm)) p0)
               c (iterate-corrector #(correct h ys fs (f t1 %)) p tol iterations)
               y1 (if modify? (mapv - c (scale (mapv - c p0) cm)) c)]
           (recur (conj (subvec pts 1) [t1 y1]) (conj out [t1 y1]) (inc k) (mapv - c p0))))))))
