(ns allgo.robotics.canny
  "Canny's roadmap algorithm for a smooth algebraic surface in R^3: LaValle,
  *Planning Algorithms*, Cambridge University Press, 2006, section 6.4.3.

  The surface is the zeros of one polynomial f, in general position --
  Examples 6.7 and 6.8 are the sphere and the torus. A roadmap is a set of
  curves on it that is connected wherever the surface is, so that any two
  points can be joined by going to the roadmap, along it, and off again.
  Canny's is made of critical curves rather than samples:

    1. Choose a linear map g = (g1, g2) to the plane, in general
       position. Its silhouette -- where the surface's normal is
       perpendicular to the view, where the Jacobian (6.28) of (f, g)
       drops rank -- is a set of curves on the surface: f = 0 and
       grad f . view = 0.
    2. The critical points of g1 alone -- where the normal is parallel to
       the sweep direction -- lie on the silhouette, where it turns back
       in the sweep. They are found there and polished by Newton's method
       on f = 0, grad f . across = 0, grad f . view = 0.
    3. The silhouette need not be connected where the surface is: a
       torus's inner and outer rings are not (Figure 6.37). So, in the
       recursion's step, at each critical point the surface is sliced by
       the plane g1 = const through it, and the slice's own roadmap added.
       A slice of a surface is already a curve, the recursion's base case:
       the slice itself joins the rings (Figures 6.38 and 6.39).
    4. Query points are declared critical at the top level, so slices
       through them join them to the rest.

  The book does this symbolically, with subresultants; it notes that no
  general implementation exists. Here, for a surface in R^3, each
  critical curve is the zero set of two functions and is traced
  numerically: from seeds on a lattice pulled onto the curve by Newton's
  method, a predictor step along the curve's tangent -- the cross product
  of the two gradients -- and a corrector back onto it, the step halved
  where the curve bends, until it closes on itself or reaches a critical
  point. Where a slice crosses a singular point -- the slice through a
  torus's inner saddle is a figure eight -- the branches end at it. The
  pieces meet at junctions, solved for exactly: where a slice crosses the
  silhouette, and at the critical points. Searching the graph of the
  curves gives the path. A component smaller than the seed lattice's
  spacing can be missed.

  The surface's polynomial is an `allgo.numerics.multivariate` map in x,
  y, z; points are `[x y z]`; a frame is `{:sweep :across :view}`,
  orthonormal, g1 the coordinate along `:sweep` and the silhouette seen
  along `:view`."
  (:require [allgo.geometry.rotation :as rotation]
            [allgo.geometry.vec3 :as v]
            [allgo.graph :as graph]
            [allgo.numerics.linear :as lin]
            [allgo.numerics.multivariate :as mv]
            [allgo.numerics.roots :as roots]
            [allgo.search :as search]))

;; ---------------------------------------------------------------------------
;; The surface and the curves on it

(defn surface
  "The zeros of the polynomial `f` in x, y, z, with its gradient and
  Hessian ready to evaluate."
  [f]
  {:polynomial f
   :f (mv/evaluator f)
   :grad (mapv mv/evaluator (mv/gradient f 3))
   :hess (mapv #(mapv mv/evaluator %) (mv/hessian f 3))})

(defn value [s p] ((:f s) p))
(defn gradient [s p] (mapv #(% p) (:grad s)))
(defn- hessian [s p] (mapv (fn [row] (mapv #(% p) row)) (:hess s)))

(def axes
  "The coordinate frame: sweep along x, silhouette seen along z, as in
  Example 6.7."
  {:sweep [1.0 0.0 0.0] :across [0.0 1.0 0.0] :view [0.0 0.0 1.0]})

(defn turned
  "The coordinate frame turned by small angles about x, y and z, in
  radians -- a frame in general position for a surface whose symmetry
  would make the axes degenerate."
  [ax ay az]
  ;; the columns of the rotation are the turned axes
  (let [[sweep across view] (lin/transpose (rotation/chain (rotation/rotate-z az) (rotation/rotate-y ay) (rotation/rotate-x ax)))]
    {:sweep sweep :across across :view view}))

(defn silhouette-system
  "The silhouette as two equations: f = 0, grad f . view = 0."
  [s {:keys [view]}]
  {:g (fn [p] [(value s p) (v/dot (gradient s p) view)])
   :jac (fn [p] [(gradient s p) (lin/mat-vec (hessian s p) view)])})

(defn slice-system
  "The slice through the plane sweep . p = `c`: f = 0, sweep . p - c = 0."
  [s {:keys [sweep]} c]
  {:g (fn [p] [(value s p) (- (v/dot sweep p) c)])
   :jac (fn [p] [(gradient s p) sweep])})

(defn- residual [{:keys [g]} p]
  (let [[a b] (g p)] (max (abs a) (abs b))))

(defn- newton-next
  "One step of Newton's method toward the curve, the least step that
  zeroes both equations to first order: `p` + J^T (J J^T)^-1 (-g)."
  [{:keys [g jac]} p]
  (let [[g1 g2] (g p) [a b] (jac p)
        aa (v/dot a a) ab (v/dot a b) bb (v/dot b b)
        det (- (* aa bb) (* ab ab))]
    (when (> (abs det) 1e-14)
      (let [l1 (/ (+ (* (- g1) bb) (* g2 ab)) det)
            l2 (/ (+ (* aa (- g2)) (* ab g1)) det)]
        (v/add p (v/add (v/scale a l1) (v/scale b l2)))))))

(defn onto
  "The point of the curve `system` Newton's method reaches from `p`, or
  nil if it does not converge in 30 steps without wandering more than
  `reach` away."
  [system p reach]
  (->> (search/step-seq (partial newton-next system) p)
       (take 31)
       (take-while #(<= (v/distance % p) reach))
       (filter #(< (residual system %) 1e-10))
       first))

(defn- surface-next
  "One Newton step toward f = 0, along the gradient."
  [s p]
  (let [g (gradient s p) gg (v/dot g g)]
    (when (> gg 1e-14)
      (v/add p (v/scale g (/ (- (value s p)) gg))))))

(defn onto-surface
  "The point of the surface Newton's method reaches from `p`, or nil."
  [s p]
  (->> (search/step-seq (partial surface-next s) p)
       (take 31)
       (filter #(< (abs (value s %)) 1e-10))
       first))

(defn- tangent [{:keys [jac]} p]
  (let [[a b] (jac p) t (v/cross a b)]
    (when (> (v/length t) 1e-9) (v/normalize t))))

;; ---------------------------------------------------------------------------
;; Tracing a curve, a step at a time

(defn- trace-step
  "The next point along the curve: a step `h` along the tangent, kept
  heading the way the last went, and Newton's method back onto the
  curve. A step that fails, jumps, or turns sharply is halved and tried
  again. The trace ends -- `:end` says why -- when it closes on its
  start, comes within a step of one of the `stops` (it ends there),
  leaves the box `inside?` admits, or runs to `limit` points."
  [system {:keys [step stops inside? limit]} {:keys [points heading h end] :as st}]
  (when-not end
    (let [p (peek points)
          start (first points)
          t0 (tangent system p)
          stop-near (fn [q r] (first (filter #(< (v/distance q %) r) stops)))]
      (if (nil? t0)
        (assoc st :end :singular)
        (let [t (if (neg? (v/dot t0 heading)) (v/negate t0) t0)
              q (onto system (v/add p (v/scale t h)) (* 2 h))
              t1 (some->> q (tangent system))]
          (if (or (nil? t1) (< (abs (v/dot t1 t)) 0.9) (> (v/distance q p) (* 1.5 h)))
            (cond
              (>= h (/ step 32)) (assoc st :h (/ h 2))
              ;; stuck where the curve goes singular: a critical point
              (stop-near p (* 4 step)) (assoc st :points (conj points (stop-near p (* 4 step))) :end :stop)
              :else (assoc st :end :stuck))
            (let [stop (stop-near q (* 1.5 step))]
              (cond
                stop (assoc st :points (conj points stop) :end :stop)
                (and (> (count points) 4) (< (v/distance q start) (* 1.5 step))) (assoc st :points (conj points q) :end :closed)
                (not (inside? q)) (assoc st :points (conj points q) :end :outside)
                (>= (count points) limit) (assoc st :points (conj points q) :end :limit)
                :else (assoc st :points (conj points q) :heading t :h (min step (* 2 h)))))))))))

(defn trace-seq
  "The states of a trace of the curve `system` from `start`, heading
  along `heading`, a point each: `{:points :heading :h :end}`. Options:
  `:step`, `:stops`, `:inside?` and `:limit`."
  [system start heading opts]
  (search/step-seq (partial trace-step system opts) {:points [start] :heading heading :h (:step opts)}))

(defn trace
  "The curve through `start`: `{:points :closed?}`, traced both ways
  unless it closes on itself first."
  [system start opts]
  (let [t (tangent system start)
        forward (last (trace-seq system start t opts))]
    (if (= :closed (:end forward))
      {:points (:points forward) :closed? true}
      (let [backward (last (trace-seq system start (v/negate t) opts))]
        {:points (into (vec (rseq (:points backward))) (rest (:points forward))) :closed? false}))))

(defn- covered? [curves q r]
  (some (fn [{:keys [points]}] (some #(< (v/distance q %) r) points)) curves))

(defn- trace-from
  "`curves` with the curve through `seed` added, if the seed reaches the
  curve, is not near a stop and is not on one already traced."
  [system {:keys [step stops reach] :as opts} curves seed]
  (let [q (onto system seed reach)]
    (if (or (nil? q) (nil? (tangent system q))
            (some #(< (v/distance q %) (* 2 step)) stops)
            (covered? curves q (* 2 step)))
      curves
      (conj curves (trace system q opts)))))

(defn curves
  "Every curve of `system` some seed reaches: each seed pulled onto the
  curve and, if it lands on none traced yet, the curve through it
  traced."
  [system seeds opts]
  (reduce (partial trace-from system opts) [] seeds))

;; ---------------------------------------------------------------------------
;; Critical points and junctions

(defn- solve-3
  "The root of three equations `F`, Jacobian `J`, Newton's method reaches
  from `p`, or nil."
  [F J p]
  (when-let [q (roots/newton-system F J p {:max-iter 30})]
    (when (< (reduce max (map abs (F q))) 1e-10) q)))

(defn- extremes
  "The points of the polyline where the coordinate along `dir` turns."
  [{:keys [points closed?]} dir]
  (let [n (count points)
        at #(v/dot dir (points (mod % n)))]
    (for [i (if closed? (range n) (range 1 (dec n)))
          :let [a (at (dec i)) b (at i) c (at (inc i))]
          :when (or (and (> b a) (>= b c)) (and (< b a) (<= b c)))]
      (points i))))

(defn- distinct-points [pts r]
  (reduce (fn [acc p] (if (some #(< (v/distance p %) r) acc) acc (conj acc p))) [] pts))

(defn critical-points
  "The critical points of the sweep coordinate on the surface: where the
  silhouette turns in the sweep, polished on f = 0, grad f . across = 0,
  grad f . view = 0."
  [s {:keys [across view] :as frame} silhouettes]
  (let [F (fn [p] (let [g (gradient s p)] [(value s p) (v/dot g across) (v/dot g view)]))
        J (fn [p] (let [h (hessian s p)] [(gradient s p) (lin/mat-vec h across) (lin/mat-vec h view)]))]
    (distinct-points (keep #(solve-3 F J %) (mapcat #(extremes % (:sweep frame)) silhouettes)) 1e-6)))

(defn junctions
  "Where the silhouettes cross the plane sweep . p = `c`: each `{:point
  :curve :segment}`, the silhouette and the index of the segment it
  crosses, the point solved for on f = 0, grad f . view = 0,
  sweep . p = c."
  [s {:keys [sweep view]} silhouettes c]
  (let [F (fn [p] [(value s p) (v/dot (gradient s p) view) (- (v/dot sweep p) c)])
        J (fn [p] [(gradient s p) (lin/mat-vec (hessian s p) view) sweep])]
    (for [[ci {:keys [points closed?]}] (map-indexed vector silhouettes)
          :let [n (count points)]
          k (range (if closed? n (dec n)))
          :let [a (points k) b (points (mod (inc k) n))
                da (- (v/dot sweep a) c) db (- (v/dot sweep b) c)]
          :when (neg? (* da db))
          :let [p (solve-3 F J (v/lerp a b (/ da (- da db))))]
          :when p]
      {:point p :curve ci :segment k})))

;; ---------------------------------------------------------------------------
;; The roadmap

(defn- lattice-seeds
  "A lattice of `n` points a side over the box `[lo hi]`."
  [[lo hi] n]
  (let [at (fn [i k] (+ (nth lo k) (* (- (nth hi k) (nth lo k)) (/ (+ i 0.5) n))))]
    (for [i (range n) j (range n) k (range n)] [(at i 0) (at j 1) (at k 2)])))

(defn plane-seeds
  "A lattice of `n` points a side over the plane sweep . p = `c`, as wide
  as the box."
  [{:keys [sweep across view]} [lo hi] c n]
  (let [r (* 0.5 (v/distance lo hi))]
    (for [i (range n) j (range n)
          :let [a (- (* 2 r (/ (+ i 0.5) n)) r) b (- (* 2 r (/ (+ j 0.5) n)) r)]]
      (v/add (v/scale sweep c) (v/add (v/scale across a) (v/scale view b))))))

(defn- nearest-on
  "The node of `curves` (indexed from `offset`) nearest `p`, within `r`."
  [curves offset p r]
  (let [best (apply min-key second [nil ##Inf]
                    (for [[ci {:keys [points]}] (map-indexed vector curves)
                          [k q] (map-indexed vector points)]
                      [[:pt (+ offset ci) k] (v/distance p q)]))]
    (when (< (second best) r) (first best))))

(defn- link [graph a b cost]
  (-> graph (assoc-in [a b] cost) (assoc-in [b a] cost)))

(defn- graph-of
  "The roadmap's graph: a node for each traced point, junction, critical
  point and query, edges along the curves and from each junction,
  critical point and query to the curves through it, weighed by length."
  [{:keys [silhouettes slices junctions critical queries step]}]
  (let [curves (into (vec silhouettes) (mapcat :curves slices))
        offsets (reductions + (count silhouettes) (map (comp count :curves) slices))
        slice-at (zipmap (map :c slices) (map vector offsets (map :curves slices)))
        positions (merge (into {} (for [[ci {:keys [points]}] (map-indexed vector curves)
                                        [k p] (map-indexed vector points)]
                                    [[:pt ci k] p]))
                         (into {} (map-indexed (fn [j {:keys [point]}] [[:junction j] point]) junctions))
                         (into {} (map-indexed (fn [m {:keys [point]}] [[:critical m] point]) critical))
                         (into {} (map-indexed (fn [q {:keys [point]}] [[:query q] point]) queries)))
        d (fn [a b] (v/distance (positions a) (positions b)))
        along (for [[ci {:keys [points closed?]}] (map-indexed vector curves)
                    k (range (if closed? (count points) (dec (count points))))]
                [[:pt ci k] [:pt ci (mod (inc k) (count points))]])
        onto-slice (fn [node c]
                     (when-let [[offset cs] (slice-at c)]
                       (when-let [n (nearest-on cs offset (positions node) (* 2 step))] [[node n]])))
        joins (concat
               (for [[j {:keys [curve segment c]}] (map-indexed vector junctions)
                     :let [n (count (:points (silhouettes curve)))]
                     e (concat [[[:junction j] [:pt curve segment]]
                                [[:junction j] [:pt curve (mod (inc segment) n)]]]
                               (onto-slice [:junction j] c))]
                 e)
               (for [[m {:keys [c]}] (map-indexed vector critical)
                     e (concat (some-> (nearest-on silhouettes 0 (positions [:critical m]) (* 2 step))
                                       (#(vector [[:critical m] %])))
                               (onto-slice [:critical m] c))]
                 e)
               (for [[q {:keys [c]}] (map-indexed vector queries)
                     e (onto-slice [:query q] c)]
                 e))]
    {:positions positions
     :graph (reduce (fn [g [a b]] (link g a b (d a b))) (zipmap (keys positions) (repeat {})) (concat along joins))}))

(defn roadmap
  "Canny's roadmap of the surface, the zeros of polynomial `f`, inside the
  box `bounds`, `[[x0 y0 z0] [x1 y1 z1]]`. Options:

    :frame      the sweep and view directions (default `axes`)
    :queries    points on the surface to join to the roadmap
    :step       the tracing step (default a fiftieth of the box's diagonal)
    :seeds      seeds a side for the silhouette's lattice (default 12),
                and for each slice's
    :slices?    whether to slice at the critical points (default true);
                false gives the silhouettes and critical points alone, as
                in Figure 6.37

  Returns `{:silhouettes :critical :slices :junctions :queries :graph
  :positions}`: curves as `{:points :closed?}`, critical points and
  queries as `{:point :c}`, slices as `{:c :curves}`."
  [f bounds {:keys [frame queries step seeds slices?]
             :or {frame axes seeds 12 slices? true}}]
  (let [s (surface f)
        [lo hi] bounds
        step (or step (/ (v/distance lo hi) 50))
        inside? (fn [p] (every? true? (map <= lo p hi)))
        sil-opts {:step step :stops [] :inside? inside? :limit 4000 :reach (* 0.5 (v/distance lo hi))}
        silhouettes (curves (silhouette-system s frame) (lattice-seeds bounds seeds) sil-opts)
        crit (critical-points s frame silhouettes)
        c-of #(v/dot (:sweep frame) %)
        critical (mapv (fn [p] {:point p :c (c-of p)}) crit)
        queries (mapv (fn [p] {:point p :c (c-of p)}) queries)
        cs (when slices? (distinct (map :c (concat critical queries))))
        slices (vec (for [c cs]
                      {:c c
                       :curves (curves (slice-system s frame c)
                                       (concat (keep #(when (< (abs (- (:c %) c)) 1e-12) (:point %)) queries)
                                               (plane-seeds frame bounds c (* 2 seeds)))
                                       (assoc sil-opts :stops crit))}))
        junctions (vec (for [c cs j (junctions s frame silhouettes c)] (assoc j :c c)))
        rm {:surface s :frame frame :step step :silhouettes silhouettes :critical critical
            :slices slices :junctions junctions :queries queries}]
    (merge rm (graph-of rm))))

(defn components
  "The connected components of the roadmap's graph, as sets of nodes."
  [{:keys [graph]}]
  (graph/components-by (keys graph) #(keys (graph %))))

(defn route
  "The shortest way along the roadmap from query `i` to query `j`, as
  points, or nil if the roadmap does not join them."
  [{:keys [graph positions]} i j]
  (some->> (search/dijkstra-seq {:start [:query i] :goal [:query j]
                                 :neighbors #(keys (graph %))
                                 :cost (fn [a b] (get-in graph [a b]))})
           search/goal-state search/path (mapv positions)))
