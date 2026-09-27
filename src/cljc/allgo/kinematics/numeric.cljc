(ns allgo.kinematics.numeric
  "Inverse kinematics by descent, for chains that have no closed form.

  `allgo.kinematics.analytic` is exact, instant and complete -- it hands
  back every configuration at once -- but it only exists because a
  spherical wrist makes the algebra separable. Take the wrist apart, add a
  seventh joint, bolt on a tool at an angle, and there is no formula. What
  is left is to start somewhere and walk downhill.

  The Jacobian says how the tool moves per unit of each joint. Inverting
  it gives the joints to move per unit of tool, which would be the whole
  method except that it is exactly wrong where it is most needed: near a
  singularity some direction costs unbounded joint speed, the inverse
  blows up, and the arm flails.

  Damped least squares is the standard cure, and it is a bargain rather
  than a fix. Instead of solving for the step that lands exactly on the
  target, solve for the one that balances getting there against how far
  the joints must move:

    dq = J^T (J J^T + lambda^2 I)^-1 e

  `lambda` is that balance. At zero it is the plain pseudo-inverse and
  unstable near singularities; raised, it gives up some accuracy per step
  and stays bounded everywhere. The damping here is adaptive: large while
  the error is large, falling away as the arm closes in, so it converges
  quickly far out and precisely near the target.

  What this gives up against the closed form is completeness. It returns
  the one solution it walked to, which depends on where it started -- so
  the count of configurations that the analytic solver makes plain is
  invisible here. Start it from several seeds and you will find several,
  but you will never know whether you have them all."
  (:require [allgo.geometry.quaternion :as q]
            [allgo.geometry.vec3 :as v]
            [allgo.kinematics.chain :as k]
            [allgo.math :as am]
            [allgo.numerics.linear :as lin]
            [allgo.random :as random]
            [clojure.math :as math]))

(defn pose-error
  "The six-vector from where the tool is to where it should be: three of
  position, three of rotation.

  The rotation part is the axis of the turn that would take the current
  orientation to the target, scaled by the angle -- the same thing the
  angular row of the Jacobian is measured in, which is what lets the two
  halves be solved together."
  [current target]
  (let [dp (v/sub (:pos target) (:pos current))
        dq (q/mul (:rot target) (q/inverse (:rot current)))
        ;; Take the short way round: a quaternion and its negation are the
        ;; same orientation but opposite turns.
        dq (if (neg? (nth dq 3)) (mapv - dq) dq)
        [axis angle] (q/to-axis-angle dq)]
    (into (vec dp) (v/scale axis angle))))

(def default-options
  {:tolerance    1e-6
   :iterations   200
   ;; Damping at full error, easing off as the arm closes in.
   :damping      0.05
   ;; How far any joint may move in one step, which keeps a large error
   ;; from throwing the arm across its workspace in a single stride.
   :max-step     0.3
   :limits?      true
   ;; Weigh a radian of orientation against a meter of position. They are
   ;; different units and the solver has to be told what they are worth.
   :orientation-weight 1.0})

(defn- damped-step
  "One damped least squares step toward the target."
  [chain values err lambda weight]
  (let [J (k/jacobian chain values)
        ;; Rows are the six task directions, columns the joints.
        Jt (lin/transpose (mapv vec J))
        w (mapv (fn [i e] (if (< i 3) e (* weight e))) (range 6) err)
        Jw (mapv (fn [i row] (if (< i 3) row (mapv #(* weight %) row)))
                 (range 6) Jt)
        A (lin/mat-add (lin/mat-mul Jw (lin/transpose Jw))
                       (lin/mat-scale (lin/eye 6) (* lambda lambda)))]
    (when-let [y (lin/cholesky-solve A w)]
      (lin/mat-vec (lin/transpose Jw) y))))

(defn solve
  "Walks `from` toward a configuration that puts the tool at `target`.

  Returns `{:values :error :iterations :converged?}`. The values are
  always the best it reached, so a run that did not converge still says
  where it got to -- useful when the target is out of reach and the
  honest answer is the nearest approach."
  ([chain target] (solve chain target (k/home chain) {}))
  ([chain target from] (solve chain target from {}))
  ([chain target from opts]
   (let [{:keys [tolerance iterations damping max-step limits? orientation-weight]}
         (merge default-options opts)]
     (loop [values (vec from) i 0]
       (let [err (pose-error (k/pose chain values) target)
             size (v/length (vec (take 3 err)))
             turn (v/length (vec (drop 3 err)))
             total (+ size (* orientation-weight turn))]
         (cond
           (< total tolerance)
           {:values values :error total :iterations i :converged? true}

           (>= i iterations)
           {:values values :error total :iterations i :converged? false}

           :else
           ;; Damping proportional to how far there is to go: bold while
           ;; the target is distant, careful as it arrives.
           (let [lambda (* damping (am/clamp total 1e-3 1.0))
                 dq (damped-step chain values err lambda orientation-weight)]
             (if (nil? dq)
               {:values values :error total :iterations i :converged? false}
               (let [scale (let [biggest (reduce max 0.0 (map abs dq))]
                             (if (> biggest max-step) (/ max-step biggest) 1.0))
                     next-values (mapv (fn [a d] (+ a (* scale d))) values dq)
                     next-values (if limits? (k/clamp chain next-values) next-values)]
                 (recur next-values (inc i)))))))))))

(defn solve-from-many
  "Runs the descent from several starts and keeps the distinct results.

  The nearest this method comes to the analytic solver's completeness: a
  handful of seeds will usually turn up more than one configuration, and
  will never tell you how many there are."
  ([chain target] (solve-from-many chain target 8 {}))
  ([chain target seeds opts]
   (let [n (k/joint-count chain)
         rng-values (fn [seed]
                      (let [rng (random/rng seed)]
                        (vec (repeatedly n #(random/uniform rng (- math/PI) math/PI)))))]
     (->> (cons (k/home chain) (map rng-values (range seeds)))
          (map #(solve chain target % opts))
          (filter :converged?)
          (map :values)
          (reduce (fn [acc s]
                    (if (some (fn [seen]
                                (every? #(< (abs (am/wrap-angle %)) 1e-3) (map - seen s)))
                              acc)
                      acc
                      (conj acc s)))
                  [])))))
