(ns allgo.astro.halo
  "Periodic orbits about the collinear Lagrange points and the tubes of
  trajectories that lead to and from them, in the circular restricted
  three-body problem (`allgo.astro.cr3bp`'s units and frame).

  About L1 and L2 the linearized motion is a saddle in one direction and
  centers in the others, and the centers grow into families of periodic
  orbits: the planar Lyapunov orbits, and branching from them where their
  out-of-plane motion comes into resonance with the in-plane, the
  three-dimensional halo orbits (Farquhar 1968; Breakwell and Brown 1979;
  Howell, \"Three-dimensional, periodic, 'halo' orbits\", Celestial
  Mechanics 32, 1984). Each is found here by differential correction:
  starting on the x axis (or in the x-z plane) moving perpendicular to it,
  the orbit is symmetric about that plane if it crosses the plane
  perpendicularly half a period later, and the state transition matrix
  says how to change the start to make it do so. Families are followed
  by continuation, each orbit seeding the next; the halos by starting at
  the Lyapunov orbit where they branch off, located as the orbit whose
  out-of-plane monodromy has trace 2.

  The saddle makes every such orbit unstable, and its monodromy matrix --
  the transition matrix over one period -- has a real eigenvalue pair
  lambda and 1/lambda. Their eigenvectors, carried around the orbit, give
  the directions of the unstable and stable invariant manifolds: the
  trajectories that leave the orbit, and those that arrive on it without
  a burn (Gomez, Koon, Lo, Marsden, Masdemont and Ross 2004), which is
  how a spacecraft gets onto a halo orbit cheaply and how the tubes link
  one libration point to another.

  States are `[x y z vx vy vz]`, rotating and nondimensional."
  (:require [allgo.astro.cr3bp :as cr3bp]
            [allgo.numerics.core :as core]
            [allgo.numerics.linear :as lin]
            [allgo.numerics.rk :as rk]
            [clojure.math :as math]))

;; ------------------------------------------------- the flow and its STM

(defn- hessian
  "The second partials of the potential U at `[x y z]`."
  [mu [x y z]]
  (let [a (+ x mu) b (- x (- 1.0 mu))
        r1 (math/sqrt (+ (* a a) (* y y) (* z z)))
        r2 (math/sqrt (+ (* b b) (* y y) (* z z)))
        m1 (- 1.0 mu)
        p1 (/ m1 (* r1 r1 r1)) p2 (/ mu (* r2 r2 r2))
        q1 (/ (* 3.0 m1) (math/pow r1 5.0)) q2 (/ (* 3.0 mu) (math/pow r2 5.0))
        uxx (+ 1.0 (- p1) (- p2) (* q1 a a) (* q2 b b))
        uyy (+ 1.0 (- p1) (- p2) (* q1 y y) (* q2 y y))
        uzz (+ (- p1) (- p2) (* q1 z z) (* q2 z z))
        uxy (+ (* q1 a y) (* q2 b y))
        uxz (+ (* q1 a z) (* q2 b z))
        uyz (+ (* q1 y z) (* q2 y z))]
    [[uxx uxy uxz] [uxy uyy uyz] [uxz uyz uzz]]))

(defn- rhs
  "The equations of motion, and with `stm?` the variational equations
  alongside -- the state followed by the 6x6 transition matrix row by
  row -- negated with `backward?` to run time the other way."
  [mu stm? backward?]
  (let [f (cr3bp/derivative mu)
        sgn (if backward? -1.0 1.0)]
    (fn [t y]
      (let [s (subvec y 0 6)
            ds (f t s)]
        (if-not stm?
          (mapv #(* sgn %) ds)
          (let [phi (partition 6 (subvec y 6))
                [p0 p1 p2 p3 p4 p5] (map vec phi)
                [[h00 h01 h02] [h10 h11 h12] [h20 h21 h22]] (hessian mu s)
                row (fn [h0 h1 h2 extra] (mapv (fn [a b c e] (+ (* h0 a) (* h1 b) (* h2 c) e)) p0 p1 p2 extra))
                d3 (row h00 h01 h02 (mapv #(* 2.0 %) p4))
                d4 (row h10 h11 h12 (mapv #(* -2.0 %) p3))
                d5 (row h20 h21 h22 (vec (repeat 6 0.0)))]
            (mapv #(* sgn %) (concat ds p3 p4 p5 d3 d4 d5))))))))

(def ^:private tolerance {:tol-abs 1e-13 :tol-rel 1e-13})

(defn- with-stm [s] (into (vec s) (flatten (lin/eye 6))))

(defn- unpack [y] {:state (subvec y 0 6) :stm (when (> (count y) 6) (mapv vec (partition 6 (subvec y 6))))})

(defn propagate
  "The state `t` after `state` (before it, for negative `t`), and with
  `:stm? true` the transition matrix too: `{:state :stm}`."
  ([mu state t] (propagate mu state t {}))
  ([mu state t {:keys [stm?]}]
   (if (zero? t)
     (unpack (if stm? (with-stm state) (vec state)))
     (let [y0 (if stm? (with-stm state) (vec state))
           integ (rk/integrator rk/dopri54 (rhs mu stm? (neg? t)) 0.0 y0 (* 0.01 (abs t)) tolerance)]
       (unpack (:y (core/step-until integ (abs t))))))))

(defn trajectory
  "States along the trajectory from `state` for `t` (negative for
  backward), `n` of them evenly spaced in time after the first, as
  `[{:t :state} ...]` -- for drawing."
  [mu state t n]
  (let [dt (/ t n)]
    (vec (reductions (fn [{:keys [t state]} _]
                       {:t (+ t dt) :state (:state (propagate mu state dt))})
                     {:t 0.0 :state (vec state)}
                     (range n)))))

(defn- crossing
  "The next time, after a moment, that the trajectory from `state` meets
  the plane y = 0: `{:t :state :stm}`, the transition matrix from the
  start with `stm?`. Stepped until y changes sign, then Newton's method on
  the time from the step before, y' being vy."
  [mu state stm?]
  (let [y0 (if stm? (with-stm state) (vec state))
        integ (rk/integrator rk/dopri54 (rhs mu stm? false) 0.0 y0 1e-3 tolerance)
        steps (->> (core/trajectory integ) (drop-while #(< (:t %) 1e-2)))
        sign #(neg? (get-in % [:y 1]))
        [before _] (first (filter (fn [[a b]] (not= (sign a) (sign b))) (partition 2 1 steps)))
        refine (fn [{:keys [t y]}]
                 (loop [t t y y i 0]
                   (let [dt (- (/ (y 1) (y 4)))]
                     (if (or (< (abs (y 1)) 1e-13) (> i 20))
                       {:t t :y y}
                       (let [back? (neg? dt)
                             integ (rk/integrator rk/dopri54 (rhs mu stm? back?) 0.0 y (* 0.5 (abs dt)) tolerance)
                             y' (:y (core/step-until integ (abs dt)))]
                         (recur (+ t dt) y' (inc i)))))))
        {:keys [t y]} (refine before)]
    (assoc (unpack y) :t t)))

;; ------------------------------------------------------ correcting orbits

(defn- solve2 [[[a b] [c d]] [e f]]
  (let [det (- (* a d) (* b c))]
    [(/ (- (* d e) (* b f)) det) (/ (- (* a f) (* c e)) det)]))

(defn correct
  "The periodic orbit symmetric about the x-z plane nearest the start
  `[x0 0 z0 0 vy0 0]`: half a period on it crosses the plane with vx = vz
  = 0. Holds `:fix` (`:x` or `:z`, default `:z`) and corrects the other
  and vy0 by Newton's method; for a planar start (z0 = 0) holds x0 and
  corrects vy0 alone. Returns `{:state :period :jacobi}`, or nil if it
  does not converge."
  ([mu start] (correct mu start {}))
  ([mu [x0 _ z0 _ vy0 _] {:keys [fix tol] :or {fix :z tol 1e-11}}]
   (let [f (cr3bp/derivative mu)
         planar? (zero? z0)]
     (loop [x0 x0 z0 z0 vy0 vy0 i 0]
       (when (< i 40)
         (let [s0 [x0 0.0 z0 0.0 vy0 0.0]
               {:keys [t state stm]} (crossing mu s0 true)
               [_ _ _ vx vy vz] state
               [_ _ _ ax _ az] (f 0.0 state)]
           (if (and (< (abs vx) tol) (or planar? (< (abs vz) tol)))
             {:state s0 :period (* 2.0 t) :jacobi (cr3bp/jacobi mu [(subvec s0 0 3) (subvec s0 3 6)])}
             ;; the change in vx (and vz) at the crossing for a change in
             ;; a start coordinate j, the crossing time moving with it:
             ;; phi[i][j] - (acceleration_i / vy) phi[1][j]
             (let [d (fn [row j] (- (get-in stm [row j]) (* (/ (if (= row 3) ax az) vy) (get-in stm [1 j]))))]
               (if planar?
                 (recur x0 z0 (- vy0 (/ vx (d 3 4))) (inc i))
                 (let [j (if (= fix :z) 0 2)
                       [dj dvy] (solve2 [[(d 3 j) (d 3 4)] [(d 5 j) (d 5 4)]] [vx vz])]
                   (if (= fix :z)
                     (recur (- x0 dj) z0 (- vy0 dvy) (inc i))
                     (recur x0 (- z0 dj) (- vy0 dvy) (inc i)))))))))))))

(defn monodromy
  "The transition matrix over one period of the periodic `orbit`."
  [mu {:keys [state period]}]
  (:stm (propagate mu state period {:stm? true})))

(defn- dominant
  "The dominant eigenvalue and unit eigenvector of `m`, by power
  iteration -- quick when, as for an unstable orbit's monodromy, it
  stands far above the rest."
  [m]
  (loop [v (vec (repeat (count m) 1.0)) lam 0.0 i 0]
    (let [w (lin/mat-vec m v)
          n (math/sqrt (reduce + (map * w w)))
          lam' (/ (reduce + (map * v w)) (reduce + (map * v v)))
          v' (mapv #(/ % n) w)]
      (if (or (> i 500) (< (abs (- lam' lam)) (* 1e-12 (abs lam'))))
        [lam' v']
        (recur v' lam' (inc i))))))

(defn stability
  "The periodic `orbit`'s stability index, (lambda + 1/lambda)/2 for its
  largest eigenvalue: 1 on the edge of stability, and how fast departures
  grow each period above that (the index JPL's three-body periodic orbit
  catalog lists)."
  [mu orbit]
  (let [[lam] (dominant (monodromy mu orbit))]
    (* 0.5 (+ lam (/ 1.0 lam)))))

;; -------------------------------------------------------------- families

(defn- linear-seed
  "A small planar orbit about the collinear point at `xl`, from the
  linearized motion: x = xl + A cos wt, y = -(w^2 + Uxx)/(2w) A sin wt."
  [mu xl amplitude]
  (let [[[uxx _ _] [_ uyy _]] (hessian mu [xl 0.0 0.0])
        ;; w^4 + (4 - uxx - uyy) w^2 + uxx uyy = 0, the oscillating root
        b (- 4.0 uxx uyy)
        w2 (* 0.5 (+ (- b) (math/sqrt (- (* b b) (* 4.0 uxx uyy)))))
        w2 (if (pos? w2) w2 (* 0.5 (- (- b) (math/sqrt (- (* b b) (* 4.0 uxx uyy))))))
        w (math/sqrt (abs w2))]
    [(+ xl amplitude) 0.0 0.0 0.0 (* -0.5 (+ (* w w) uxx) amplitude) 0.0]))

(defn- collinear-x [mu point]
  (first ((cr3bp/lagrange-points mu) point)))

(defn lyapunov
  "The planar Lyapunov family about `point` (`:L1`, `:L2` or `:L3`),
  lazily, from the smallest outward: each orbit `{:state :period :jacobi}`
  starting on the x axis `step` (default 1e-3) further from the point than
  the last, on the `:side` -1 or +1 of it (default toward the larger
  primary from L1, away from it from L2 and L3)."
  ([mu point] (lyapunov mu point {}))
  ([mu point {:keys [step side] :or {step 1e-3}}]
   (let [xl (collinear-x mu point)
         dir (double (or side (if (= point :L1) -1 1)))
         first-orbit (correct mu (linear-seed mu xl (* dir step)))]
     (letfn [(more [prev cur]
               (lazy-seq
                (let [[x0 _ _ _ vy0] (:state cur)
                      ;; the secant through the last two predicts the next
                      [_ _ _ _ pvy] (:state prev)
                      guess [(+ x0 (* dir step)) 0.0 0.0 0.0 (+ vy0 (- vy0 pvy)) 0.0]]
                  (when-let [nxt (correct mu guess)]
                    (cons nxt (more cur nxt))))))]
       (when first-orbit
         ;; the point itself is the family's orbit of no size
         (cons first-orbit (more {:state [xl 0.0 0.0 0.0 0.0 0.0]} first-orbit)))))))

(defn- vertical-trace
  "Half the trace of the out-of-plane block of a planar orbit's monodromy,
  cos of its vertical rotation: 1 where a halo family branches."
  [mu orbit]
  (let [m (monodromy mu orbit)]
    (* 0.5 (+ (get-in m [2 2]) (get-in m [5 5])))))

(defn- find-bifurcation [mu point opts]
  (let [fam (lyapunov mu point (merge {:step 2e-3} opts))
        g #(- (vertical-trace mu %) 1.0)
        [[a ga] [b gb]] (first (filter (fn [[[_ ga] [_ gb]]] (not= (neg? ga) (neg? gb)))
                                       (partition 2 1 (map (juxt identity g) fam))))
        orbit-at (fn [a b w]
                   (let [[xa _ _ _ va] (:state a) [xb _ _ _ vb] (:state b)]
                     (correct mu [(+ xa (* w (- xb xa))) 0.0 0.0 0.0 (+ va (* w (- vb va))) 0.0])))]
    ;; regula falsi, the Illinois way: the trace is smooth along the
    ;; family, and a handful of corrections find where it crosses
    (loop [a a ga ga b b gb gb side nil i 0]
      (let [c (orbit-at a b (/ ga (- ga gb)))
            gc (g c)
            [xa] (:state a) [xb] (:state b)]
        (cond
          (or (> i 60) (< (abs gc) 1e-12) (< (abs (- xa xb)) 1e-12)) c
          (= (neg? gc) (neg? gb)) (recur a (if (= side :b) (* 0.5 ga) ga) c gc :b (inc i))
          :else (recur c gc b (if (= side :a) (* 0.5 gb) gb) :a (inc i)))))))

(def ^:private bifurcation-memo (memoize find-bifurcation))

(defn bifurcation
  "The planar Lyapunov orbit about `point` from which the halo family
  branches, the first where the vertical trace reaches 1, bracketed along
  the family and refined by regula falsi in the start. Options are
  `lyapunov`'s. Remembered, since each halo's continuation starts here."
  ([mu point] (bifurcation mu point {}))
  ([mu point opts] (bifurcation-memo mu point opts)))

(defn halo
  "The halo orbit about `point` whose start in the x-z plane has height
  `z0` (positive for the northern family, negative for the southern),
  continued from the family's branch point on the Lyapunov family in
  steps of `:step` (default 2e-3) in z0: `{:state :period :jacobi}`.
  `:side` is `lyapunov`'s: which of its two crossings of the x-z plane
  the orbit starts from."
  ([mu point z0] (halo mu point z0 {}))
  ([mu point z0 {:keys [step] :or {step 2e-3} :as opts}]
   (let [{[x0 _ _ _ vy0] :state} (bifurcation mu point (select-keys opts [:side]))
         sgn (if (neg? z0) -1.0 1.0)
         zs (concat (range (* sgn 1e-4) z0 (* sgn step)) [z0])]
     (loop [[z & more] zs prev nil cur [x0 vy0]]
       (let [[x vy] cur
             [px pvy] (or prev cur)
             ;; the secant predicts the next start from the last two
             guess [(+ x (- x px)) 0.0 z 0.0 (+ vy (- vy pvy)) 0.0]
             o (correct mu guess)]
         (cond
           (nil? o) nil
           (empty? more) o
           :else (let [[nx _ _ _ nvy] (:state o)] (recur more cur [nx nvy]))))))))

;; ------------------------------------------------------------ manifolds

(defn manifold
  "The unstable (or with `:stable? true`, the stable) manifold of the
  periodic `orbit`: `n` (default 20) trajectories, each starting from a
  point spread evenly in time around the orbit displaced by `epsilon`
  (default 1e-6, some 400 m in the Earth-Moon system) along the
  eigenvector carried there by the flow, on the `:side` +1 or -1, and
  flown for `duration` -- forward from the unstable manifold, backward
  toward the stable, so that each traces where it leads or where it comes
  from. Returns `[{:start :end :tau} ...]`, tau the start's time along
  the orbit; see `trajectory` to draw them."
  [mu orbit duration {:keys [stable? n epsilon side] :or {n 20 epsilon 1e-6 side 1.0}}]
  (let [{:keys [state period]} orbit
        m (monodromy mu orbit)
        ;; the stable direction is the dominant one of the inverse
        [_ v0] (dominant (if stable? (lin/inverse-general m) m))]
    (vec (for [k (range n)
               :let [tau (* period (/ k n))
                     {s :state phi :stm} (propagate mu state tau {:stm? true})
                     v (lin/mat-vec phi v0)
                     scale (/ (* side epsilon) (math/sqrt (reduce + (map #(* % %) (subvec v 0 3)))))
                     start (mapv + s (map #(* scale %) v))
                     end (:state (propagate mu start (if stable? (- duration) duration)))]]
           {:tau tau :start start :end end}))))
