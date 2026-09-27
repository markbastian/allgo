(ns allgo.numerics.core
  "Shared machinery for the integrators in `allgo.numerics`, following
  Montenbruck & Gill, *Satellite Orbits*, chapter 4.

  A first-order system is a function `(f t y)` returning dy/dt, where the
  state `y` is a vector of doubles. A second-order system -- the form
  satellite motion actually takes -- is `(f t y y')` returning y''; those are
  integrated by `allgo.numerics.rkn` and Stoermer-Cowell rather than by
  reduction to twice the dimension, which is the point of having them.

  An integrator is a map of state, advanced by `step`. Keeping it a value
  rather than a loop lets a caller drive it one step per animation frame,
  inspect the error estimate, or restart it."
  (:require [allgo.math :as am]
            [allgo.numerics.linear :as lin]
            [clojure.math :as math]))

(defn combine
  "Sum of `vs` weighted by `ws`, skipping zero weights -- in a Butcher tableau
  most of them are zero, and each skipped term is a whole vector not built."
  [ws vs]
  (reduce (fn [acc [w v]]
            (if (zero? w) acc (lin/add acc (lin/scale v w))))
          (vec (repeat (count (first vs)) 0.0))
          (map vector ws vs)))

(defn norm
  "RMS norm of `e` scaled by the magnitude of `y`, so a tolerance means the
  same thing for a state of 1e-6 as for one of 1e6."
  [e y tol-abs tol-rel]
  (math/sqrt (/ (reduce + (map (fn [ei yi]
                                 (let [s (+ tol-abs (* tol-rel (abs yi)))]
                                   (* (/ ei s) (/ ei s))))
                               e y))
                (count e))))

;; ---------------------------------------------------------------- stepping

(defn step
  "Advance an integrator one step. Adaptive methods may reject the step and
  retry with a smaller one; `:accepted` says whether time actually moved."
  [{:keys [method] :as integ}]
  ((:advance method) integ))

(defn step-until
  "Advance to exactly `t-end`, or until `max-steps` is spent.

  The last step is clipped to land on `t-end` rather than overshooting it,
  which matters more than it looks: comparing two runs at a nominal end time
  they had each sailed past by up to a step makes error measurements
  meaningless, and convergence order unreadable. A clipped step is not
  allowed to shrink an adaptive method's stride permanently."
  ([integ t-end] (step-until integ t-end 1000000))
  ([integ t-end max-steps]
   (loop [s integ n 0]
     (let [remaining (- t-end (:t s))]
       (if (or (<= remaining (* 1e-12 (max 1.0 (abs t-end)))) (>= n max-steps))
         s
         (let [natural (:h s)
               clipped (min natural remaining)
               s'      (step (assoc s :h clipped))]
           (recur (if (and (:accepted s') (< clipped natural))
                    (assoc s' :h natural)
                    s')
                  (inc n))))))))

(defn trajectory
  "Lazy sequence of states from an integrator, one per accepted step."
  [integ]
  (->> (iterate step integ)
       (filter :accepted)
       (map #(select-keys % [:t :y :dy :h :error]))))

;; ------------------------------------------------------- step-size control

(def control-defaults
  {:tol-abs   1e-9
   :tol-rel   1e-9
   :safety    0.9
   :min-scale 0.2
   :max-scale 5.0
   :h-min     1e-12
   :h-max     ##Inf})

(defn integrator
  "The state every integrator here carries.

  Seven constructors were building this same map -- the method, the system,
  where and how big a step, the control parameters, and whether the last
  step was accepted -- and differing only in what they added to it: a
  velocity for the second-order methods, a history for the multistep ones.
  `extra` is that difference and this is the rest.

  `:adaptive?` is dropped from the control parameters because it selects a
  method rather than tuning one, and each constructor has already read it."
  ([base] (integrator base {}))
  ([{:keys [method f t y h opts]} extra]
   (merge {:method   method
           :f        f
           :t        (double t)
           :y        (mapv double y)
           :h        (double h)
           :control  (merge control-defaults (dissoc opts :adaptive?))
           :accepted true
           :error    nil}
          extra)))

(defn accept
  "A step taken at a fixed size, with nothing to judge it by."
  [integ updates]
  (merge integ updates {:accepted true :error nil}))

(defn adapt
  "The classic step-size law: scale by (tol/err)^(1/(p+1)), clamped so a
  freak error estimate cannot collapse or explode the step in one move.

  Returns `[accept? h-next]`."
  [err order {:keys [safety min-scale max-scale h-min h-max]} h]
  (let [scale (if (or (zero? err) (not (am/finite? err)))
                max-scale
                (am/clamp (* safety (math/pow (/ 1.0 err) (/ 1.0 (inc order)))) min-scale max-scale))
        h'    (am/clamp (* h scale) h-min h-max)]
    [(<= err 1.0) h']))

(defn settle
  "Accept or reject a step by the size of its error estimate.

  Every adaptive method here ends the same way, and the rejected branch is
  the same in all of them: time and state stand, only the step shrinks. It
  is short enough to have been written out five times and important enough
  that five copies would have to stay in agreement -- a method that
  advanced time on a rejected step would take a wrong answer and never
  revisit it."
  [integ err order updates]
  (let [[ok? h'] (adapt err order (:control integ) (:h integ))]
    (if ok?
      (merge integ updates {:h h' :accepted true :error err})
      (assoc integ :h h' :accepted false :error err))))
