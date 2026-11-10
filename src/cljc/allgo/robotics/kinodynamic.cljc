(ns allgo.robotics.kinodynamic
  "Planning with dynamics: LaValle, *Planning Algorithms*, Cambridge
  University Press, 2006, sections 14.4.1 and 14.4.3 -- a robot that
  cannot stop or turn at once, so a plan is a sequence of controls, not
  of points.

  The robot here is the double integrator in the plane, a point mass the
  controls accelerate: state `[x y vx vy]`, control `[ax ay]` held for a
  time step dt, after which, exactly,

    p' = p + v dt + a dt^2 / 2      v' = v + a dt

  with each component of the acceleration at most `:amax` and of the
  velocity at most `:vmax`.

    lattice planning   with each control component -amax, 0 or amax, the
                       states reachable from rest lie on a lattice -- a
                       position every amax dt^2 / 2, a velocity every
                       amax dt -- so A* over the lattice, a step a
                       control, finds the quickest plan the lattice holds
    kinodynamic RRT    a tree grown toward random states: from the node
                       nearest each sample, the control that comes
                       nearest it. It plans in any state space, at the
                       price of never quite arriving: the goal is a
                       region, here a position within a tolerance at a
                       low speed

  Worlds are `allgo.robotics.motion`'s; a trajectory is checked for
  collisions at sub-steps of each control's arc."
  (:require [allgo.robotics.motion :as motion]
            [allgo.search :as search]
            [clojure.math :as math]))

(defn propagate
  "The state after holding control `[ax ay]` for `dt` from `[x y vx vy]`."
  [[x y vx vy] [ax ay] dt]
  [(+ x (* vx dt) (* 0.5 ax dt dt)) (+ y (* vy dt) (* 0.5 ay dt dt))
   (+ vx (* ax dt)) (+ vy (* ay dt))])

(defn controls
  "The nine controls: each component -amax, 0 or amax."
  [amax]
  (for [ax [(- amax) 0.0 amax] ay [(- amax) 0.0 amax]] [ax ay]))

(defn arc-free?
  "Whether the arc from `state` under `u` for `dt` stays free, checked at
  `n` sub-steps."
  ([world state u dt] (arc-free? world state u dt 4))
  ([world state u dt n]
   (let [pts (map #(subvec (propagate state u (* dt (/ % n))) 0 2) (range (inc n)))]
     (every? (fn [[a b]] (motion/segment-free? world a b)) (partition 2 1 pts)))))

(defn- in-speed? [vmax [_ _ vx vy]] (and (<= (abs vx) (+ vmax 1e-9)) (<= (abs vy) (+ vmax 1e-9))))

;; ---------------------------------------------------------------------------
;; Lattice planning

(defn lattice-plan
  "The quickest plan on the lattice from rest at `start` to rest at a
  lattice point within `tolerance` of `goal`: `{:controls :states :time}`,
  or nil. Each axis' lattice coordinates are integers m and k, the
  position start + m amax dt^2 / 2 and the velocity k amax dt, and a
  control a in {-1, 0, 1} takes them to m + 2k + a and k + a -- exact, so
  the states never drift off the lattice."
  [world start goal {:keys [amax vmax dt tolerance]}]
  (let [step (* 0.5 amax dt dt)
        kmax (long (math/floor (/ vmax (* amax dt))))
        [sx sy] start
        ->state (fn [[mx my kx ky]] [(+ sx (* mx step)) (+ sy (* my step)) (* kx amax dt) (* ky amax dt)])
        [gx gy] goal
        dist (fn [[mx my]] (let [[x y] (->state [mx my 0 0])] (math/hypot (- x gx) (- y gy))))
        h (fn [[mx my] _] (let [[x y] (->state [mx my 0 0])]
                            (math/ceil (/ (max 0.0 (- (max (abs (- x gx)) (abs (- y gy))) tolerance)) (* vmax dt)))))
        result (search/a-star {:start [0 0 0 0]
                               :goal? (fn [[mx my kx ky]] (and (zero? kx) (zero? ky) (<= (dist [mx my]) tolerance)))
                               :neighbors (fn [[mx my kx ky]]
                                            (for [ax [-1 0 1] ay [-1 0 1]
                                                  :let [kx' (+ kx ax) ky' (+ ky ay)]
                                                  :when (and (<= (abs kx') kmax) (<= (abs ky') kmax)
                                                             (arc-free? world (->state [mx my kx ky]) [(* ax amax) (* ay amax)] dt))]
                                              [(+ mx (* 2 kx) ax) (+ my (* 2 ky) ay) kx' ky']))
                               :heuristic h})]
    (when result
      {:states (mapv ->state result)
       :controls (mapv (fn [[[_ _ kx ky] [_ _ kx' ky']]] [(* (- kx' kx) amax) (* (- ky' ky) amax)]) (partition 2 1 result))
       :time (* dt (dec (count result)))})))

;; ---------------------------------------------------------------------------
;; Kinodynamic RRT

(defn- metric
  "The distance between states, a speed weighed `w` seconds."
  ^double [w [x y vx vy] [a b va vb]]
  (math/sqrt (+ (* (- x a) (- x a)) (* (- y b) (- y b)) (* w w (+ (* (- vx va) (- vx va)) (* (- vy vb) (- vy vb)))))))

(defn- reached?
  "Whether `state` is within `tolerance` of `goal` at under `slow` speed."
  [goal {:keys [tolerance slow]} [x y vx vy]]
  (and (<= (math/hypot (- x (goal 0)) (- y (goal 1))) tolerance) (< (math/hypot vx vy) slow)))

(defn- sample-state
  "A random state in the room and the speed limits, or the goal at rest
  a tenth of the time."
  [{[[x0 y0] [x1 y1]] :bounds} goal rng vmax]
  (if (< (rng) 0.1)
    (conj (vec goal) 0.0 0.0)
    [(+ x0 (* (- x1 x0) (rng))) (+ y0 (* (- y1 y0) (rng)))
     (* vmax (- (* 2 (rng)) 1)) (* vmax (- (* 2 (rng)) 1))]))

(defn- rrt-step
  "One sample of the kinodynamic RRT: the best control from the nearest
  node toward it, if any is free and within the speed limit."
  [world goal rng {:keys [amax vmax dt w] :as opts} {:keys [nodes] :as tree}]
  (if (:goal tree)
    tree
    (let [q (sample-state world goal rng vmax)
          distance (partial metric w)
          i (motion/nearest-node nodes distance q)
          s (nodes i)
          options (for [u (controls amax) :let [s' (propagate s u dt)]
                        :when (and (in-speed? vmax s') (arc-free? world s u dt))]
                    [u s'])]
      (if (empty? options)
        tree
        (let [[u s'] (apply min-key #(distance (second %) q) options)]
          (cond-> (-> tree (update :nodes conj s') (update :parent conj i) (update :control conj u))
            (reached? goal opts s') (assoc :goal (count nodes))))))))

(defn rrt-seq
  "A kinodynamic RRT from rest at `start`, a sample a step: a random
  state (the goal at rest a tenth of the time), the node nearest it by
  `metric`, and from there whichever control's free, in-speed arc ends
  nearest the sample becomes a new node. Each state `{:nodes :parent
  :control :goal}`, as `allgo.robotics.motion`'s trees, `:goal` the first
  node within `tolerance` of `goal` at less than `slow` speed."
  [world start goal rng opts]
  (let [opts (merge {:w 0.5} opts)
        s0 (conj (vec start) 0.0 0.0)]
    (iterate (partial rrt-step world goal rng opts)
             {:nodes [s0] :parent [nil] :control [nil] :goal (when (reached? goal opts s0) 0)})))

(defn rrt-plan
  "The plan to a tree's `:goal` node: `{:states :controls}`, the states
  `motion/tree-path`'s."
  [{:keys [parent control goal] :as tree}]
  (when goal
    {:states (motion/tree-path tree)
     :controls (mapv control (rest (reverse (take-while some? (iterate parent goal)))))}))
