(ns allgo.robotics.tracking
  "Trajectory tracking control: Russell and Norvig, *Artificial
  Intelligence: A Modern Approach*, 4th edition, section 26.5.3.

  A planner gives a path, points without times. Tracking it takes three
  things:

    retiming        the path made a trajectory, a time for every point:
                    here the profile the section suggests, accelerating
                    at the most allowed to the top speed, holding it, and
                    braking to a stop -- a triangle when the path is too
                    short to reach the top speed; and, for a finely
                    sampled curve, the time-optimal profile under a
                    speed limit from its curvature as well
    dynamics        what torques do: for a planar two-link arm with its
                    masses at the ends of its links, M(q) q'' + h(q, q')
                    + g(q) = u, the inertia, the Coriolis and centrifugal
                    torques and gravity's (as in Craig, *Introduction to
                    Robotics*, 3rd edition, chapter 6); forward, the
                    acceleration torques give, and inverse, the torques
                    an acceleration needs
    control laws    feedforward, the inverse dynamics of the trajectory,
                    exact if the model is; PD feedback on the error; and
                    computed torque, the two together, the feedback
                    scaled by the model's inertia

  A trajectory is `{:duration :at}`, `(at t)` -> `[q q' q'']`, each a
  vector. An arm is `{:m1 :m2 :l1 :l2 :g}`: masses, lengths and gravity,
  the angles from the horizontal, the second relative to the first."
  (:require [allgo.numerics.linear :as lin]
            [clojure.math :as math]))

;; ---------------------------------------------------------------------------
;; Retiming

(defn retime
  "The `path`, a vector of points of any dimension, as a trajectory that
  starts and ends at rest: arc length s(t) by a trapezoidal profile of top
  speed `vmax` and acceleration `amax`, and the point at s along the path.
  Along each straight leg the velocity is s' along it and the
  acceleration s'' along it; at a corner the direction turns at once, as
  a path's corner does."
  [path vmax amax]
  (let [legs (mapv (fn [[a b]] (lin/length (lin/sub b a))) (partition 2 1 path))
        total (reduce + legs)
        starts (vec (reductions + 0.0 legs))
        [ta vpeak] (if (>= total (/ (* vmax vmax) amax))
                     [(/ vmax amax) vmax]
                     (let [v (math/sqrt (* total amax))] [(/ v amax) v]))
        cruise (- (/ total vpeak) ta)
        duration (+ (* 2 ta) (max 0.0 cruise))
        profile (fn [t]
                  (let [t (max 0.0 (min duration t))]
                    (cond (< t ta) [(* 0.5 amax t t) (* amax t) amax]
                          (< t (- duration ta)) [(+ (* 0.5 amax ta ta) (* vpeak (- t ta))) vpeak 0.0]
                          :else (let [r (- duration t)] [(- total (* 0.5 amax r r)) (* amax r) (- amax)]))))]
    {:duration duration
     :at (fn [t]
           (let [[s sd sdd] (profile t)
                 i (max 0 (min (dec (count legs)) (dec (count (take-while #(<= % s) starts)))))
                 a (nth path i) b (nth path (inc i))
                 dir (if (pos? (legs i)) (lin/scale (lin/sub b a) (/ 1.0 (legs i))) (lin/scale a 0.0))]
             [(lin/add a (lin/scale dir (- s (starts i)))) (lin/scale dir sd) (lin/scale dir sdd)]))}))

;; ---------------------------------------------------------------------------
;; Time-optimal speed along a path
;;
;; LaValle, *Planning Algorithms*, Cambridge University Press, 2006,
;; section 14.6.3: along a fixed path the only choice left is the speed
;; s'(s), and the quickest profile under bounds on it is the largest one
;; the bounds allow -- full acceleration or full braking everywhere it is
;; not at a limit. Here the bounds are a top speed, a tangential
;; acceleration, and a lateral acceleration v^2 kappa that caps the speed
;; where the path bends: a forward pass accelerates as hard as allowed
;; from rest, a backward pass brakes as hard as allowed into the end, and
;; the profile is the lesser of the two at every sample.

(defn curvature
  "The curvature of the circle through three points: 4 area / abc, zero
  on a line."
  ^double [a b c]
  (let [ab (lin/sub b a) ac (lin/sub c a)
        twice-area (abs (- (* (first ab) (second ac)) (* (second ab) (first ac))))
        abc (* (lin/length ab) (lin/length ac) (lin/length (lin/sub c b)))]
    (if (zero? abc) 0.0 (/ (* 2.0 twice-area) abc))))

(defn speed-profile
  "The time-optimal speed along `path`, a vector of planar points closely
  spaced, from rest to rest: at most `vmax`, speeding up or slowing down
  at most `amax`, and turning at most `alat` (speed^2 times curvature).
  `{:s :v :t}`, at each point its arc length, speed and time."
  [path {:keys [vmax amax alat]}]
  (let [n (count path)
        ds (mapv (fn [[a b]] (lin/length (lin/sub b a))) (partition 2 1 path))
        cap (vec (concat [0.0]
                         (map (fn [[a b c]] (let [k (curvature a b c)] (if (pos? k) (min vmax (math/sqrt (/ alat k))) vmax)))
                              (partition 3 1 path))
                         [0.0]))
        reach (fn [v d limit] (min limit (math/sqrt (+ (* v v) (* 2.0 amax d)))))
        forward (vec (reductions (fn [v i] (reach v (ds (dec i)) (cap i))) 0.0 (range 1 n)))
        v (vec (reverse (reductions (fn [v i] (reach v (ds i) (forward i))) 0.0 (range (- n 2) -1 -1))))]
    {:s (vec (reductions + 0.0 ds))
     :v v
     :t (vec (reductions + 0.0 (map (fn [d [a b]] (/ (* 2.0 d) (+ a b))) ds (partition 2 1 v))))}))

;; ---------------------------------------------------------------------------
;; A two-link arm's dynamics

(defn inertia
  "The arm's inertia matrix M(q) at joint angles `[t1 t2]`."
  [{:keys [m1 m2 l1 l2]} [_ t2]]
  (let [c2 (math/cos t2)
        m12 (+ (* m2 l2 l2) (* m2 l1 l2 c2))]
    [[(+ (* (+ m1 m2) l1 l1) (* m2 l2 l2) (* 2 m2 l1 l2 c2)) m12]
     [m12 (* m2 l2 l2)]]))

(defn- bias
  "The Coriolis, centrifugal and gravity torques h(q, q') + g(q)."
  [{:keys [m1 m2 l1 l2 g]} [t1 t2] [w1 w2]]
  (let [s2 (math/sin t2) c1 (math/cos t1) c12 (math/cos (+ t1 t2))
        k (* m2 l1 l2 s2)]
    [(+ (- (* k (+ (* 2 w1 w2) (* w2 w2)))) (* (+ m1 m2) l1 g c1) (* m2 l2 g c12))
     (+ (* k w1 w1) (* m2 l2 g c12))]))

(defn inverse-dynamics
  "The torques that give the arm acceleration `qdd` at angles `q` and
  velocities `qd`: M(q) q'' + h(q, q') + g(q)."
  [arm q qd qdd]
  (let [[[a b] [c d]] (inertia arm q)
        [h1 h2] (bias arm q qd)
        [x y] qdd]
    [(+ (* a x) (* b y) h1) (+ (* c x) (* d y) h2)]))

(defn forward-dynamics
  "The arm's acceleration under torques `u` at angles `q`, velocities
  `qd`: M(q)^-1 (u - h - g)."
  [arm q qd u]
  (let [[[a b] [c d]] (inertia arm q)
        [r1 r2] (lin/sub u (bias arm q qd))
        det (- (* a d) (* b c))]
    [(/ (- (* d r1) (* b r2)) det) (/ (- (* a r2) (* c r1)) det)]))

(defn energy
  "The arm's kinetic plus potential energy: q' M q' / 2 plus the masses'
  heights times g -- unchanged by motion under no torque."
  [{:keys [m1 m2 l1 l2 g] :as arm} [t1 t2 :as q] qd]
  (let [M (inertia arm q)
        kinetic (* 0.5 (reduce + (map * qd (mapv #(reduce + (map * % qd)) M))))]
    (+ kinetic (* (+ m1 m2) g l1 (math/sin t1)) (* m2 g l2 (math/sin (+ t1 t2))))))

(defn step
  "The arm's state `[q qd]` after `dt` under torques `u` held constant:
  fourth-order Runge-Kutta."
  [arm [q qd] u dt]
  (let [f (fn [q qd] [qd (forward-dynamics arm q qd u)])
        [k1q k1v] (f q qd)
        [k2q k2v] (f (lin/add q (lin/scale k1q (/ dt 2))) (lin/add qd (lin/scale k1v (/ dt 2))))
        [k3q k3v] (f (lin/add q (lin/scale k2q (/ dt 2))) (lin/add qd (lin/scale k2v (/ dt 2))))
        [k4q k4v] (f (lin/add q (lin/scale k3q dt)) (lin/add qd (lin/scale k3v dt)))
        avg (fn [a b c d] (lin/scale (lin/add (lin/add a (lin/scale b 2.0)) (lin/add (lin/scale c 2.0) d)) (/ dt 6.0)))]
    [(lin/add q (avg k1q k2q k3q k4q)) (lin/add qd (avg k1v k2v k3v k4v))]))

;; ---------------------------------------------------------------------------
;; Control laws

(defn feedforward
  "Open loop: the torques `model`'s inverse dynamics says the trajectory
  needs, u = f^-1(xi, xi', xi''). It ignores where the arm actually is."
  [model]
  (fn [_ _ [xi xid xidd]] (inverse-dynamics model xi xid xidd)))

(defn pd
  "PD feedback, u = Kp (xi - q) + Kd (xi' - q'): no model at all."
  [kp kd]
  (fn [q qd [xi xid]] (lin/add (lin/scale (lin/sub xi q) kp) (lin/scale (lin/sub xid qd) kd))))

(defn computed-torque
  "Computed torque control, equation 26.6: the feedforward torques plus
  feedback on the error scaled by `model`'s inertia at the trajectory's
  configuration, u = f^-1(xi, xi', xi'') + M(xi) (Kp (xi - q) + Kd (xi' -
  q')) -- gains that change with the arm's configuration."
  [model kp kd]
  (fn [q qd [xi xid xidd]]
    (let [[[a b] [c d]] (inertia model xi)
          [e1 e2] (lin/add (lin/scale (lin/sub xi q) kp) (lin/scale (lin/sub xid qd) kd))]
      (lin/add (inverse-dynamics model xi xid xidd) [(+ (* a e1) (* b e2)) (+ (* c e1) (* d e2))]))))

(defn- track-step
  "A tick of the arm under `controller`, toward where the trajectory
  says it should be."
  [arm controller at dt {:keys [t q qd]}]
  (let [ref (at t)
        u (controller q qd ref)
        [q qd] (step arm [q qd] u dt)
        t (+ t dt)]
    {:t t :q q :qd qd :ref (first (at t))}))

(defn track-seq
  "The `arm` driven by `controller` along `trajectory`, from rest at its
  start, a step of `dt` a step -- and on past its end, holding its last
  point: each state `{:t :q :qd :ref}`, `:ref` where the arm should be."
  [arm controller {:keys [at]} dt]
  (let [[q0] (at 0.0)]
    (iterate (partial track-step arm controller at dt) {:t 0.0 :q q0 :qd [0.0 0.0] :ref q0})))

(defn tracking-error
  "The largest distance, in joint angles, between the arm and where it
  should be, over `states`."
  [states]
  (reduce max (map (fn [{:keys [q ref]}] (lin/length (lin/sub q ref))) states)))
