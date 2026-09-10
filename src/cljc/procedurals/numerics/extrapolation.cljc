(ns procedurals.numerics.extrapolation
  "Gragg-Bulirsch-Stoer extrapolation (Montenbruck & Gill 4.3).

  The idea is different in kind from the other families here. Rather than
  design a formula of some fixed order, cross the step several times with a
  crude rule at decreasing sub-step, then extrapolate that sequence of
  answers to the limit of zero sub-step.

  What makes it work is Gragg's modified midpoint rule, whose error runs in
  *even* powers of h alone. Extrapolating against h^2 therefore kills two
  orders per level rather than one, so k levels reach order 2k -- orders
  which would be hopeless to construct as a tableau are routine here."
  (:require [procedurals.numerics.core :as core]))

(defn modified-midpoint
  "Gragg's rule: cross `H` in `n` sub-steps by centred differences, then take
  the smoothing average at the end that cancels the odd error terms."
  [f t y H n]
  (let [h  (/ H n)
        f0 (f t y)]
    (loop [m     1
           prev  y
           cur   (core/v+ y (core/v* f0 h))]
      (if (= m n)
        ;; The final average is what leaves an expansion in h^2 alone.
        (core/v* (core/v+ (core/v+ cur prev)
                          (core/v* (f (+ t H) cur) h))
                 0.5)
        (recur (inc m)
               cur
               (core/v+ prev (core/v* (f (+ t (* m h)) cur) (* 2.0 h))))))))

(def step-sequence
  "Sub-step counts 2, 4, 6, 8, ... (Deuflhard). Every one is even, which the
  smoothing average at the end of the midpoint rule requires."
  (mapv #(* 2 (inc %)) (range 16)))

(defn extrapolate
  "Aitken-Neville extrapolation of `values` to zero sub-step, in h^2.

  Returns `[best difference]`, the difference being the gap between the last
  two columns -- the classic estimate for an extrapolation tableau, and the
  reason no second method is needed to size the step. It is returned
  unnormalised so the caller can scale it by its own tolerances."
  [values ns]
  (let [k (count values)]
    (loop [j 1, col (vec values), prev nil]
      ;; `col` is the current column of the tableau, `prev` the one before
      ;; it -- they must stay distinct, or the gap between them that
      ;; estimates the error collapses to zero and the step control has
      ;; nothing to push against.
      (if (>= j k)
        [(peek col) (if prev
                      (core/v- (peek col) (peek prev))
                      (core/v* (peek col) 0.0))]
        (let [next-col (vec (for [i (range k)]
                              (if (< i j)
                                (nth col i)
                                (let [a (nth col i)
                                      b (nth col (dec i))
                                      r (/ (double (nth ns i)) (nth ns (- i j)))]
                                  (core/v+ a (core/v* (core/v- a b)
                                                      (/ 1.0 (- (* r r) 1.0))))))))]
          (recur (inc j) next-col col))))))

(defn- gbs-step
  "One macro-step of `H` at `levels` extrapolation levels, so order 2*levels."
  [f t y H levels]
  (let [ns   (subvec step-sequence 0 levels)
        vals (mapv #(modified-midpoint f t y H %) ns)]
    (extrapolate vals ns)))

(defn- advance
  [{:keys [f t y h method control] :as integ}]
  (let [[y' diff] (gbs-step f t y h (:levels method))]
    (if-not (:adaptive? method)
      (assoc integ :t (+ t h) :y y' :accepted true :error nil)
      (let [err      (core/norm diff y' (:tol-abs control) (:tol-rel control))
            [ok? h'] (core/adapt err (:order method) control h)]
        (if ok?
          (assoc integ :t (+ t h) :y y' :h h' :accepted true :error err)
          (assoc integ :h h' :accepted false :error err))))))

(defn gbs
  "Gragg-Bulirsch-Stoer with `levels` extrapolation levels, of order 2*levels."
  [levels]
  {:name (str "GBS" (* 2 levels)) :order (* 2 levels) :levels levels
   :stages (reduce + (subvec step-sequence 0 levels)) :kind :extrapolation
   :advance advance})

(def catalog (mapv gbs [2 3 4 6]))

(defn integrator
  ([method f t0 y0 h] (integrator method f t0 y0 h {}))
  ([method f t0 y0 h opts]
   {:method   (cond-> method (:adaptive? opts) (assoc :adaptive? true))
    :f        f
    :t        (double t0)
    :y        (mapv double y0)
    :h        (double h)
    :control  (merge core/control-defaults (dissoc opts :adaptive?))
    :accepted true
    :error    nil}))
