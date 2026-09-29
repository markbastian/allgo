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
            [allgo.numerics.linear :as lin]
            [allgo.numerics.linear-systems :as ls]
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

;; The flow is integrated by Dormand and Prince's 5(4) pair -- the same
;; tableau as `allgo.numerics.rk/dopri54`, the same error norm and step
;; law as `allgo.numerics.core` -- but on primitive arrays, the forty-two
;; equations of the state and its transition matrix written straight into
;; a buffer: the generic integrator's vectors cost twenty times as much
;; here, where continuation and correction fly thousands of periods.

(defn- deriv!
  "The derivatives of `y` -- the state, and with `n` = 42 the transition
  matrix row by row after it -- into `dy`, times `sgn` (-1 to run time
  backward)."
  [mu sgn n ^doubles y ^doubles dy]
  (let [x (aget y 0) yy (aget y 1) z (aget y 2)
        vx (aget y 3) vy (aget y 4) vz (aget y 5)
        m1 (- 1.0 mu)
        a (+ x mu) b (- x m1)
        r1s (+ (* a a) (* yy yy) (* z z)) r2s (+ (* b b) (* yy yy) (* z z))
        r1 (math/sqrt r1s) r2 (math/sqrt r2s)
        p1 (/ m1 (* r1s r1)) p2 (/ mu (* r2s r2))]
    (aset dy 0 (* sgn vx)) (aset dy 1 (* sgn vy)) (aset dy 2 (* sgn vz))
    (aset dy 3 (* sgn (+ (* 2.0 vy) x (- (* p1 a)) (- (* p2 b)))))
    (aset dy 4 (* sgn (+ (* -2.0 vx) yy (- (* p1 yy)) (- (* p2 yy)))))
    (aset dy 5 (* sgn (- (+ (* p1 z) (* p2 z)))))
    (when (> n 6)
      (let [q1 (/ (* 3.0 p1) r1s) q2 (/ (* 3.0 p2) r2s)
            uxx (+ 1.0 (- p1) (- p2) (* q1 a a) (* q2 b b))
            uyy (+ 1.0 (- p1) (- p2) (* q1 yy yy) (* q2 yy yy))
            uzz (+ (- p1) (- p2) (* q1 z z) (* q2 z z))
            uxy (+ (* q1 a yy) (* q2 b yy))
            uxz (+ (* q1 a z) (* q2 b z))
            uyz (+ (* q1 yy z) (* q2 yy z))]
        (dotimes [j 6]
          (let [f0 (aget y (+ 6 j)) f1 (aget y (+ 12 j)) f2 (aget y (+ 18 j))
                f3 (aget y (+ 24 j)) f4 (aget y (+ 30 j)) f5 (aget y (+ 36 j))]
            (aset dy (+ 6 j) (* sgn f3))
            (aset dy (+ 12 j) (* sgn f4))
            (aset dy (+ 18 j) (* sgn f5))
            (aset dy (+ 24 j) (* sgn (+ (* uxx f0) (* uxy f1) (* uxz f2) (* 2.0 f4))))
            (aset dy (+ 30 j) (* sgn (+ (* uxy f0) (* uyy f1) (* uyz f2) (* -2.0 f3))))
            (aset dy (+ 36 j) (* sgn (+ (* uxz f0) (* uyz f1) (* uzz f2))))))))
    dy))

(def ^:private tableau
  ;; the last stage is taken at the fifth-order solution itself, so the
  ;; weights b are its row of a, and only their difference from b-hat,
  ;; the error estimate, is wanted apart
  (let [{:keys [a b b-hat]} rk/dopri54]
    {:a (into-array (map double-array a))
     :e (double-array (map - b b-hat))}))

(defn- copy! [^doubles from ^doubles to]
  (dotimes [i (alength from)] (aset to i (aget from i)))
  to)

(def ^:private tol 1e-13)

(defn- fly!
  "Integrates `y0` (a double array of 6 or 42) forward in its own time
  from 0 to `t-end`, backward in the problem's when `sgn` is -1. With
  `stop?`, instead until y changes sign after the first moment, returning
  `[t-before y-before]` -- the last step's start -- for refining; otherwise
  the final array."
  [mu y0 t-end sgn stop?]
  (let [n (alength ^doubles y0)
        {:keys [^objects a ^doubles e]} tableau
        ks (object-array (repeatedly 7 #(double-array n)))
        tmp (double-array n)
        y5 (double-array n)
        f! (fn [^doubles y ^doubles dy] (deriv! mu sgn n y dy))]
    (f! y0 (aget ks 0))
    (loop [t 0.0 h (if stop? 1e-3 (* 0.01 t-end)) ^doubles y (aclone ^doubles y0)]
      (if (and (not stop?) (<= (- t-end t) (* 1e-12 (max 1.0 t-end))))
        y
        (let [h (if stop? h (min h (- t-end t)))]
          ;; the stages
          (dotimes [s 6]
            (let [^doubles as (aget a (inc s))]
              (dotimes [i n]
                (aset tmp i (+ (aget y i)
                               (* h (loop [j 0 acc 0.0]
                                      (if (< j (alength as))
                                        (recur (inc j) (+ acc (* (aget as j) (aget ^doubles (aget ks j) i))))
                                        acc))))))
              (f! tmp (aget ks (inc s)))))
          ;; the fifth-order solution is the last stage's point; the error
          (let [^doubles k6 (aget ks 6)
                err (loop [i 0 acc 0.0]
                      (if (< i n)
                        (let [ei (* h (loop [j 0 s 0.0]
                                        (if (< j 7) (recur (inc j) (+ s (* (aget e j) (aget ^doubles (aget ks j) i)))) s)))
                              sc (+ tol (* tol (abs (aget tmp i))))]
                          (recur (inc i) (+ acc (* (/ ei sc) (/ ei sc)))))
                        (math/sqrt (/ acc n))))
                scale (if (or (zero? err) (NaN? err)) 5.0 (min 5.0 (max 0.2 (* 0.9 (math/pow (/ 1.0 err) 0.2)))))
                h' (* h scale)]
            (if (> err 1.0)
              (recur t h' y)
              (let [t' (+ t h)]
                (copy! tmp y5)
                (if (and stop? (> t' 1e-2) (not= (neg? (aget y 1)) (neg? (aget y5 1))))
                  [t (aclone y)]
                  (do
                    ;; first same as last
                    (copy! k6 (aget ks 0))
                    (recur t' h' (aclone y5))))))))))))

(defn- ->array [y] (double-array y))

(defn- with-stm [s] (double-array (concat s (flatten (lin/eye 6)))))

(defn- unpack [^doubles y]
  (let [v (vec y)]
    {:state (subvec v 0 6) :stm (when (> (count v) 6) (mapv vec (partition 6 (subvec v 6))))}))

(defn- fly-for
  "The array `dt` on from `y` (negative for backward)."
  [mu ^doubles y dt]
  (if (zero? dt) (aclone y) (fly! mu y (abs dt) (if (neg? dt) -1.0 1.0) false)))

(defn propagate
  "The state `t` after `state` (before it, for negative `t`), and with
  `:stm? true` the transition matrix too: `{:state :stm}`."
  ([mu state t] (propagate mu state t {}))
  ([mu state t {:keys [stm?]}]
   (unpack (fly-for mu (if stm? (with-stm state) (->array state)) t))))

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
  (let [[t0 y0] (fly! mu (if stm? (with-stm state) (->array state)) 1e9 1.0 true)
        [t ^doubles y] (loop [t t0 ^doubles y y0 i 0]
                         (let [dt (- (/ (aget y 1) (aget y 4)))]
                           (if (or (< (abs (aget y 1)) 1e-13) (> i 20))
                             [t y]
                             (recur (+ t dt) (fly-for mu y dt) (inc i)))))]
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
        [_ v0] (dominant (if stable? (ls/inverse m) m))]
    (vec (for [k (range n)
               :let [tau (* period (/ k n))
                     {s :state phi :stm} (propagate mu state tau {:stm? true})
                     v (lin/mat-vec phi v0)
                     scale (/ (* side epsilon) (math/sqrt (reduce + (map #(* % %) (subvec v 0 3)))))
                     start (mapv + s (map #(* scale %) v))
                     end (:state (propagate mu start (if stable? (- duration) duration)))]]
           {:tau tau :start start :end end}))))
