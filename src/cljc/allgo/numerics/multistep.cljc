(ns allgo.numerics.multistep
  "Multistep integrators: Adams-Bashforth, Adams-Moulton and Stoermer-Cowell
  (Montenbruck & Gill 4.2).

  A single-step method rebuilds its picture of the solution from scratch
  every step, spending several force evaluations to do it. A multistep
  method instead fits a polynomial through force values it has already
  computed and integrates that, so an step of any order costs one
  evaluation. The price is a history to carry, a starting procedure to
  create it, and far more awkward step-size changes.

  The coefficients are derived here rather than tabulated. Each is the
  integral of a Lagrange basis polynomial over the step, which is the
  definition of the method, and deriving them means any order is available
  instead of only those someone printed."
  (:require [allgo.numerics.core :as core]
            [allgo.numerics.rk :as rk]
            [allgo.numerics.rkn :as rkn]))

;; ------------------------------------------------- polynomials, coefficient
;; A polynomial is a vector of coefficients, lowest power first.

(defn- poly-mul [p q]
  (reduce (fn [acc [i a]]
            (reduce (fn [acc [j b]] (update acc (+ i j) + (* a b)))
                    acc
                    (map-indexed vector q)))
          (vec (repeat (dec (+ (count p) (count q))) 0.0))
          (map-indexed vector p)))

(defn- poly-integrate
  "Definite integral of `p` from `a` to `b`."
  [p a b]
  (reduce + (map-indexed (fn [i c]
                           (let [n (inc i)]
                             (/ (* c (- (Math/pow b n) (Math/pow a n))) n)))
                         p)))

(defn- lagrange-basis
  "The basis polynomial that is one at `nodes[j]` and zero at the others."
  [nodes j]
  (let [sj     (nth nodes j)
        others (keep-indexed (fn [i s] (when (not= i j) s)) nodes)
        numer  (reduce (fn [acc si] (poly-mul acc [(- si) 1.0])) [1.0] others)
        denom  (reduce * 1.0 (map #(- sj %) others))]
    (mapv #(/ % denom) numer)))

(defn- weights
  "Integrate every basis polynomial over the nodes, after multiplying by
  `kernel` -- the identity for Adams, a hat function for Stoermer-Cowell."
  [nodes integrate]
  (mapv #(integrate (lagrange-basis nodes %)) (range (count nodes))))

;; ------------------------------------------------------------ the families

(defn- adams-nodes
  "Adams samples the force at the current point and backwards. An implicit
  family reaches one point forward as well, which is what makes it implicit."
  [k implicit?]
  (mapv #(double (- (if implicit? 1 0) %)) (range k)))

(defn adams-bashforth-coefficients
  "Explicit Adams weights of order `k`: the interpolating polynomial through
  the last `k` force values, integrated across the coming step."
  [k]
  (weights (adams-nodes k false) #(poly-integrate % 0.0 1.0)))

(defn adams-moulton-coefficients
  "Implicit Adams weights of order `k`, which include the step's own
  endpoint and so must be solved for -- in practice, corrected toward."
  [k]
  (weights (adams-nodes k true) #(poly-integrate % 0.0 1.0)))

(defn- hat-integral
  "Integrates p(u) against (1 - |u|) over [-1, 1].

  That kernel is what turns a second difference into a double integral:
  y(t+h) - 2y(t) + y(t-h) = h^2 * integral of (1-|u|) y''(t+uh) du, an
  identity that lets Stoermer-Cowell integrate y'' = f without ever forming
  a velocity."
  [p]
  (+ (poly-integrate (poly-mul p [1.0 1.0]) -1.0 0.0)
     (poly-integrate (poly-mul p [1.0 -1.0]) 0.0 1.0)))

(defn stoermer-coefficients
  "Explicit Stoermer weights of order `k`, for y'' = f(t, y)."
  [k]
  (weights (adams-nodes k false) hat-integral))

(defn cowell-coefficients
  "Implicit Cowell weights of order `k`."
  [k]
  (weights (adams-nodes k true) hat-integral))

;; ------------------------------------------------------- starting procedure

(def ^:private bootstrap-substeps 8)

(defn- bootstrap
  "The first k-1 steps, which a multistep method cannot take itself.

  Taken with a high-order single-step method in sub-steps, so that the
  history the multistep method inherits is never what limits its order --
  a startup one order too coarse quietly caps everything that follows."
  [f t y h]
  (let [hs (/ h bootstrap-substeps)]
    (loop [i 0 tt t yy y]
      (if (= i bootstrap-substeps)
        yy
        (let [[y' _] (rk/tableau-step rk/dopri54 (:c rk/dopri54) f tt yy hs)]
          (recur (inc i) (+ tt hs) y'))))))

(defn- push
  "Prepend the newest derivative, keeping only the `k` the method uses."
  [history fv k]
  (vec (take k (cons fv history))))

;; -------------------------------------------------------------- Adams

(defn- advance-adams
  [{:keys [f t y h history method control] :as integ}]
  (let [k  (:steps method)
        t' (+ t h)]
    (if (< (count history) k)
      (let [y' (bootstrap f t y h)]
        (assoc integ :t t' :y y' :history (push history (f t' y') k)
               :accepted true :starting? true))
      (let [{:keys [predictor corrector]} method
            p (core/v+ y (core/v* (core/combine predictor history) h))]
        (if-not corrector
          (assoc integ :t t' :y p :history (push history (f t' p) k)
                 :accepted true :starting? false :error nil)
          ;; PECE: predict, evaluate, correct, evaluate. The gap between
          ;; prediction and correction is Milne's estimate of the local error.
          (let [fp (f t' p)
                y' (core/v+ y (core/v* (core/combine corrector (push history fp k)) h))]
            (assoc integ :t t' :y y' :history (push history (f t' y') k)
                   :accepted true :starting? false
                   :error (core/norm (core/v- y' p) y'
                                     (:tol-abs control) (:tol-rel control)))))))))

(defn adams-bashforth
  "Explicit Adams of order `k`: one force evaluation a step, no correction."
  [k]
  {:name (str "AB" k) :order k :steps k :stages 1 :kind :multistep
   :predictor (adams-bashforth-coefficients k)
   :advance   advance-adams})

(defn adams-pece
  "Adams predictor-corrector of order `k` in PECE arrangement -- predict with
  the explicit formula, evaluate, correct with the implicit one, evaluate
  again. Two force evaluations a step buy a far smaller error constant than
  the predictor alone, and their difference estimates the error for free."
  [k]
  {:name (str "ABM" k) :order k :steps k :stages 2 :kind :multistep
   :predictor (adams-bashforth-coefficients k)
   :corrector (adams-moulton-coefficients k)
   :advance   advance-adams})

(def catalog
  (into (mapv adams-bashforth [2 3 4 5])
        (mapv adams-pece [2 3 4 5])))

(defn integrator
  "An integrator for `y' = (f t y)`. Fixed step: changing it mid-flight would
  invalidate the evenly spaced history the coefficients assume."
  ([method f t0 y0 h] (integrator method f t0 y0 h {}))
  ([method f t0 y0 h opts]
   (let [y (mapv double y0) t (double t0)]
     {:method   method
      :f        f
      :t        t
      :y        y
      :h        (double h)
      :history  [(f t y)]
      :control  (merge core/control-defaults opts)
      :accepted true
      :error    nil})))

;; ---------------------------------------------------- Stoermer and Cowell

(defn- bootstrap-second-order [f t y dy h]
  (let [hs (/ h bootstrap-substeps)]
    (loop [i 0 tt t yy y vv dy]
      (if (= i bootstrap-substeps)
        [yy vv]
        (let [[y' v'] (rkn/tableau-step rkn/rkn4 f tt yy vv hs)]
          (recur (inc i) (+ tt hs) y' v'))))))

(defn- advance-stoermer
  [{:keys [f t y dy y-prev h history method control] :as integ}]
  (let [k  (:steps method)
        t' (+ t h)]
    (if (< (count history) k)
      (let [[y' dy'] (bootstrap-second-order f t y dy h)]
        (assoc integ :t t' :y y' :dy dy' :y-prev y
               :history (push history (f t' y') k) :accepted true :starting? true))
      (let [{:keys [predictor corrector velocity]} method
            ;; Position advances on the second difference; velocity rides an
            ;; Adams sum over the same history, so both carry the same order
            ;; without a central difference capping it at two.
            step-y  (fn [w hist] (core/v+ (core/v- (core/v* y 2.0) y-prev)
                                          (core/v* (core/combine w hist) (* h h))))
            p       (step-y predictor history)
            [y' err] (if-not corrector
                       [p nil]
                       (let [fp (f t' p)
                             c  (step-y corrector (push history fp k))]
                         [c (core/norm (core/v- c p) c
                                       (:tol-abs control) (:tol-rel control))]))
            fv      (f t' y')]
        (assoc integ :t t' :y y' :y-prev y
               :dy (core/v+ dy (core/v* (core/combine velocity history) h))
               :history (push history fv k)
               :accepted true :starting? false :error err)))))

(defn stoermer
  "Explicit Stoermer of order `k` for y'' = f(t, y). Advances position on the
  second difference, so no velocity enters the position formula at all."
  [k]
  {:name (str "Stoermer" k) :order k :steps k :stages 1 :kind :multistep-2
   :predictor (stoermer-coefficients k)
   :velocity  (adams-bashforth-coefficients k)
   :advance   advance-stoermer})

(defn cowell
  "Stoermer-Cowell predictor-corrector of order `k`, the classical choice for
  orbit propagation."
  [k]
  {:name (str "Cowell" k) :order k :steps k :stages 2 :kind :multistep-2
   :predictor (stoermer-coefficients k)
   :corrector (cowell-coefficients k)
   :velocity  (adams-bashforth-coefficients k)
   :advance   advance-stoermer})

(def catalog-2
  (into (mapv stoermer [3 4 5]) (mapv cowell [3 4 5])))

(defn integrator-2
  "An integrator for `y'' = (f t y)` from `y0`, `dy0`. Fixed step."
  ([method f t0 y0 dy0 h] (integrator-2 method f t0 y0 dy0 h {}))
  ([method f t0 y0 dy0 h opts]
   (let [y (mapv double y0) t (double t0)]
     {:method   method
      :f        f
      :t        t
      :y        y
      :dy       (mapv double dy0)
      :y-prev   y
      :h        (double h)
      :history  [(f t y)]
      :control  (merge core/control-defaults opts)
      :accepted true
      :error    nil})))

;; ------------------------------------------------------- variable step size

(defn variable-coefficients
  "Adams weights for a history at arbitrary spacing.

  `offsets` are the past sample times measured from the current point in
  units of the coming step: `[0 -1 -2 ...]` for an evenly spaced history,
  and anything at all after the step has changed. Include `1` at the front
  for an implicit formula.

  This is the same integral as the fixed-step case -- a Lagrange basis
  through the sample points, integrated across the step -- which is why
  deriving the coefficients rather than tabulating them was worth doing.
  A tabulated method cannot change its step without discarding its history
  and starting again."
  [offsets]
  (weights (vec offsets) #(poly-integrate % 0.0 1.0)))

(defn- advance-variable
  "One step of variable-step Adams in PECE form.

  The step is chosen from the gap between predictor and corrector, which
  estimates the local error without a second method. A rejected step is
  retried immediately at the smaller size; because the coefficients are
  computed from the actual sample times, nothing about the history has to
  be thrown away to do it."
  [{:keys [f t y h history times method control] :as integ}]
  (let [k (:steps method)]
    (if (< (count history) k)
      ;; Startup, at fixed step, exactly as the even-spacing methods do.
      (let [y' (bootstrap f t y h)
            t' (+ t h)]
        (assoc integ :t t' :y y'
               :history (push history (f t' y') k)
               :times (vec (take k (cons t' times)))
               :accepted true :starting? true))
      (let [offsets   (mapv #(/ (- % t) h) times)
            predictor (variable-coefficients offsets)
            corrector (variable-coefficients (vec (cons 1.0 (butlast offsets))))
            t'        (+ t h)
            p         (core/v+ y (core/v* (core/combine predictor history) h))
            fp        (f t' p)
            y'        (core/v+ y (core/v* (core/combine corrector (push history fp k)) h))
            err       (core/norm (core/v- y' p) y' (:tol-abs control) (:tol-rel control))
            [ok? h']  (core/adapt err (:order method) control h)]
        (if ok?
          (assoc integ :t t' :y y' :h h'
                 :history (push history (f t' y') k)
                 :times (vec (take k (cons t' times)))
                 :accepted true :starting? false :error err)
          ;; Rejected: the history stands, only the step shrinks.
          (assoc integ :h h' :accepted false :error err))))))

(defn adams-variable
  "Adams predictor-corrector of order `k` with automatic step size.

  The fixed-step methods above are cheaper per step and simpler to reason
  about; this one earns its keep where the solution changes character --
  an eccentric orbit racing through periapsis, a system that stiffens --
  and a step chosen for the worst part of the trajectory would be wasted
  everywhere else."
  [k]
  {:name (str "ABM" k "v") :order k :steps k :stages 2 :kind :multistep-variable
   :advance advance-variable})

(def catalog-variable (mapv adams-variable [3 4 5]))

(defn variable-integrator
  "An integrator for `y' = (f t y)` with automatic step size."
  ([method f t0 y0 h] (variable-integrator method f t0 y0 h {}))
  ([method f t0 y0 h opts]
   (let [y (mapv double y0) t (double t0)]
     {:method   method
      :f        f
      :t        t
      :y        y
      :h        (double h)
      :history  [(f t y)]
      :times    [t]
      :control  (merge core/control-defaults opts)
      :accepted true
      :error    nil})))
