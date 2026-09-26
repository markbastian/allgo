(ns allgo.physics.fluid
  "Eulerian fluid on a staggered grid, after Ten Minute Physics 17.

  Where `allgo.physics.xpbd` tracks particles and asks where each one
  goes, this tracks a fixed grid and asks what flows through each cell.
  Four stages a step:

    integrate    apply gravity to the velocity field
    project      make the field divergence-free, which is the whole
                 simulation: an incompressible fluid is one that neither
                 gains nor loses material anywhere
    extrapolate  carry velocities into the border cells
    advect       move the velocity field, and anything painted on it,
                 along itself

  The grid is *staggered*: pressure and dye live at cell centers, while
  the horizontal velocity `u` lives on vertical faces and the vertical
  velocity `v` on horizontal faces. That is not fussiness. With everything
  at the center, the natural divergence of a cell does not involve that
  cell's own value, so a checcurboard of alternating pressures is
  invisible to the solver and grows without bound. On a staggered grid the
  divergence of a cell is exactly the four face velocities around it, and
  the mode cannot form.

  `s` marks which cells are open: 1 for fluid, 0 for solid. Every stage
  consults it, so obstacles are painted rather than modeled, and a moving
  obstacle is just a repainting.

  Fields are flat arrays indexed `i*ny + j`, with a one-cell border of
  solid on every side.

  ## Relation to Stam's stable fluids, and why there is no FFT here

  The advection is Jos Stam's, from Stable Fluids (1999): trace backward
  from each cell to find where the material arriving there came from, and
  interpolate. That is the step that made grid fluids practical, because
  it cannot go unstable however large the step -- the answer is always an
  interpolation between values that already exist, so nothing can grow.
  Before it, advection was explicit and bounded by a CFL condition.

  Stam's paper also gives a projection done with an FFT: on a periodic
  domain the Laplacian is diagonal in Fourier space, so the pressure solve
  becomes a transform, a divide and a transform back -- exact, and
  O(n log n) against the iterative sweep here.

  It does not apply to this solver, and the reason is `s`. An FFT
  projection needs the domain to be periodic and uniform; it solves the
  same Poisson equation everywhere and has nowhere to put a boundary
  condition. The moment part of the grid is solid -- a wall, or the
  draggable obstacle the demo is built around -- the operator stops being
  diagonal in Fourier space and the transform stops being a solution. That
  is the trade: Stam's FFT buys speed by giving up obstacles, and
  obstacles are the point here.

  What *is* worth knowing is how far from converged the sweep leaves
  things. Measured on a wind tunnel at the default forty iterations, the
  worst remaining divergence is 9.2e-2 on a 32x32 grid, 4.1e-2 on 64x64
  and 1.2e-1 on 128x128 -- it gets worse as the grid refines, because
  Gauss-Seidel damps the smooth part of the error at a rate that goes to
  zero with the cell size. Doubling the resolution wants roughly four
  times the iterations to stand still.

  So the honest speedup for this solver is not a transform but a better
  Poisson solve, and both of the usual ones keep arbitrary boundaries.
  `:solver` picks between them:

    :gauss-seidel        the sweep above, with over-relaxation. Cheapest
                         per iteration, and the convergence rate falls
                         away as the grid refines.
    :conjugate-gradient  Jacobi preconditioned. Beats the sweep once it
                         has had a hundred or so iterations and loses to
                         it before that -- over-relaxation at 1.9 is a
                         strong method and hard to beat early.
    :multigrid           V-cycles. Needs about the same number of them
                         whatever the grid, which is the property that
                         matters.

  Measured on a 128 grid wind tunnel with an obstacle, worst remaining
  divergence against the work spent:

    sweeps/cycles   Gauss-Seidel        CG          multigrid
        512          3.4e-2 146ms   7.0e-4 188ms   2.3e-5  49ms
       2048          5.0e-3 339ms   7.3e-4 351ms   2.3e-5  49ms
       8192          1.8e-6 1515ms  7.3e-4 1373ms  2.3e-5  49ms

  Read it carefully, because the three behave in three different ways.

  Multigrid is done by about thirty cycles and more do not help: it
  settles at 2.3e-5 and stays there. Gauss-Seidel is merely slow -- it
  gets past multigrid eventually, reaching 1.8e-6 by eight thousand
  sweeps, for sixteen times the time multigrid took to reach a figure
  nobody can see the difference from. Conjugate gradient genuinely stops:
  7.3e-4 at five hundred iterations and 7.3e-4 at thirty-two thousand,
  which is the loss of orthogonality that single precision does to it.

  So multigrid is the one to reach for when the grid is large, Gauss-Seidel
  is the one to reach for when a few sweeps of roughly-right is all that
  is wanted, and conjugate gradient is here mostly to be measured against.

  Where multigrid's floor comes from is not settled. It scales with the
  size of the velocities, at 4.7e-6 of the largest, which looks like
  single-precision rounding -- but Gauss-Seidel gets a decade below it on
  the same single-precision velocities, so it is not simply the storage.
  The scratch grids the cycle works in are single precision too, and that
  is the obvious next suspect; it has not been tested.

  The default is still Gauss-Seidel, because it is what every scene here
  was tuned against and forty sweeps of it look right at these grid
  sizes."
  (:require [allgo.array :as a]
            [clojure.math :as math]))

(defn fluid
  "A grid of `nx` by `ny` interior cells, each `h` across.

  A border of one cell is added on every side and marked solid, so the
  interior never has to check whether its neighbors exist."
  ([nx ny] (fluid nx ny 1.0 1000.0))
  ([nx ny h] (fluid nx ny h 1000.0))
  ([nx ny h density]
   (let [nx (+ nx 2)
         ny (+ ny 2)
         n  (* nx ny)]
     {:nx nx :ny ny :n n :h (double h) :density (double density)
      :u (a/f32 n) :v (a/f32 n)
      :u' (a/f32 n) :v' (a/f32 n)
      :p (a/f32 n)
      ;; Everything open by default; the border is closed by `close-border!`.
      :s (a/f32 n 1.0)
      :smoke (a/f32 n 1.0) :smoke' (a/f32 n 1.0)
      ;; Conjugate gradient wants five more grids. Allocated on first use,
      ;; so a scene that never asks for it never pays for it, and kept
      ;; afterward so a frame allocates nothing.
      :scratch (volatile! nil)})))

(defn idx ^long [{:keys [ny]} i j] (+ (* (long i) (long ny)) (long j)))

(defn solid? [{:keys [^floats s] :as f} i j] (zero? (aget s (idx f i j))))

(defn set-solid!
  "Marks a cell solid or open. Solid cells take no velocity and pass none."
  [{:keys [^floats s] :as f} i j solid?]
  (aset s (idx f i j) (float (if solid? 0.0 1.0)))
  f)

(defn close-border!
  "Marks the outer ring solid, leaving the left column open if `inflow?`
  so a wind tunnel has somewhere to blow from."
  ([f] (close-border! f false))
  ([{:keys [nx ny] :as f} inflow?]
   (dotimes [i nx]
     (set-solid! f i 0 true)
     (set-solid! f i (dec ny) true))
   (dotimes [j ny]
     (set-solid! f 0 j (not inflow?))
     (set-solid! f (dec nx) j true))
   f))

(defn smoke-at [{:keys [^floats smoke] :as f} i j] (aget smoke (idx f i j)))
(defn set-smoke! [{:keys [^floats smoke] :as f} i j m]
  (aset smoke (idx f i j) (float m)) f)

(defn velocity-at [{:keys [^floats u ^floats v] :as f} i j]
  [(aget u (idx f i j)) (aget v (idx f i j))])

(defn set-velocity! [{:keys [^floats u ^floats v] :as f} i j ux vy]
  (aset u (idx f i j) (float ux))
  (aset v (idx f i j) (float vy))
  f)

;; ---------------------------------------------------------------------------

(defn integrate!
  "Gravity, applied only to the vertical faces between two open cells."
  [{:keys [nx ny ^floats v ^floats s]} dt gravity]
  (let [ny (long ny)]
    (dotimes [i nx]
      (when (pos? i)
        (loop [j 1]
          (when (< j (dec ny))
            (let [k (+ (* i ny) j)]
              (when (and (pos? (aget s k)) (pos? (aget s (dec k))))
                (aset v k (float (+ (aget v k) (* gravity dt))))))
            (recur (inc j))))))))

(def default-solver
  "What `project!` does unless told otherwise."
  {:solver          :gauss-seidel
   :iterations      40
   ;; Deliberately overshooting each local correction. No physical
   ;; meaning, several times the convergence rate. Gauss-Seidel only.
   :over-relaxation 1.9
   ;; Conjugate gradient and multigrid both stop early once the worst
   ;; remaining divergence is under this, which on an easy frame is most
   ;; of them.
   :tolerance       1e-6
   ;; Smoothing sweeps either side of a multigrid coarse correction. Two
   ;; is the textbook number and is not enough here: with a solid mask and
   ;; a single open column pinning the pressure, two sweeps leave enough
   ;; of the correction unsmoothed that the cycle amplifies rather than
   ;; damps -- measured diverging to 1e7 on a 128 grid where four sweeps
   ;; converge to 2e-5.
   :pre-smooth      4
   :post-smooth     4})

;; ---------------------------------------------------------------------------
;; The pressure system
;;
;; Every solver here works on the same linear system, so they are
;; interchangeable and can be checked against each other. Writing it out,
;; for a cell `k` that is fluid and has at least one open neighbor:
;;
;;   total[k] * x[k] - sum over open neighbors of x[neighbor] = -div[k]
;;
;; where `total` counts the open neighbors, `div` is the net outflow of
;; the cell, and `x` is the pressure correction. That is the five-point
;; Laplacian with a Neumann condition at every solid face -- the wall
;; simply drops out of the sum -- and it is symmetric and positive
;; semi-definite, which is what lets conjugate gradient near it.
;;
;; Applying the answer is the same scatter Gauss-Seidel does a cell at a
;; time: a cell pushes its open faces apart by `x[k]`.

(defn- fluid-cell?
  [^floats s ^long k ^long ny]
  (and (pos? (aget s k))
       (pos? (+ (aget s (- k ny)) (aget s (+ k ny))
                (aget s (dec k)) (aget s (inc k))))))

(defn divergence-of
  "Net outflow of cell `k`: what the projection has to remove."
  ^double [^floats u ^floats v ^long k ^long ny]
  (+ (- (aget u (+ k ny)) (aget u k))
     (- (aget v (inc k)) (aget v k))))

(defn- apply-correction!
  "Pushes every open face of every fluid cell apart by `x`, and records the
  pressure that took.

  Identical to the scatter Gauss-Seidel performs one cell at a time, which
  is why a converged `x` from any solver lands the velocities in the same
  place."
  [{:keys [nx ny h density ^floats u ^floats v ^floats p ^floats s]} ^floats x dt]
  (let [nx (long nx) ny (long ny)
        cp (/ (* (double density) (double h)) (double dt))]
    (loop [i 1]
      (when (< i (dec nx))
        (loop [j 1]
          (when (< j (dec ny))
            (let [k (+ (* i ny) j)]
              (when (fluid-cell? s k ny)
                (let [xk (aget x k)]
                  (aset p k (float (+ (aget p k) (* cp xk))))
                  (aset u k (float (- (aget u k) (* (aget s (- k ny)) xk))))
                  (aset u (+ k ny) (float (+ (aget u (+ k ny)) (* (aget s (+ k ny)) xk))))
                  (aset v k (float (- (aget v k) (* (aget s (dec k)) xk))))
                  (aset v (inc k) (float (+ (aget v (inc k)) (* (aget s (inc k)) xk)))))))
            (recur (inc j))))
        (recur (inc i))))))

;; ---------------------------------------------------------------------------
;; Solvers

(defmulti pressure-solve!
  "Drives the divergence out of the velocity field.

  Dispatches on `:solver`, so a scene picks its solver by name and a new
  one is a new method rather than a change to the step. An implementation
  must leave `u` and `v` as divergence-free as it managed and add the
  pressure it used into `p`.

  The two here trade the same way they do everywhere: Gauss-Seidel is
  cheap per sweep and its convergence rate falls off as the grid refines,
  conjugate gradient costs more per iteration and degrades far more slowly."
  (fn [_f world] (:solver world)))

(defmethod pressure-solve! :default [_ world]
  (throw (ex-info "unknown pressure solver" {:solver (:solver world)})))

(defmethod pressure-solve! :gauss-seidel
  [{:keys [nx ny h density ^floats u ^floats v ^floats p ^floats s]} world]
  (let [{:keys [iterations dt over-relaxation]} world
        nx (long nx) ny (long ny)
        cp (/ (* (double density) (double h)) (double dt))]
    (dotimes [_ (long iterations)]
      (loop [i 1]
        (when (< i (dec nx))
          (loop [j 1]
            (when (< j (dec ny))
              (let [k (+ (* i ny) j)]
                (when (pos? (aget s k))
                  (let [sx0 (aget s (- k ny))
                        sx1 (aget s (+ k ny))
                        sy0 (aget s (dec k))
                        sy1 (aget s (inc k))
                        total (+ sx0 sx1 sy0 sy1)]
                    (when (pos? total)
                      ;; Each cell is corrected against the neighbors as
                      ;; they stand now, so information crosses the grid a
                      ;; cell per sweep -- which is exactly why it takes
                      ;; more sweeps the finer the grid gets.
                      (let [corr (* (double over-relaxation)
                                    (/ (- (divergence-of u v k ny)) total))]
                        (aset p k (float (+ (aget p k) (* cp corr))))
                        (aset u k (float (- (aget u k) (* sx0 corr))))
                        (aset u (+ k ny) (float (+ (aget u (+ k ny)) (* sx1 corr))))
                        (aset v k (float (- (aget v k) (* sy0 corr))))
                        (aset v (inc k) (float (+ (aget v (inc k)) (* sy1 corr)))))))))
              (recur (inc j))))
          (recur (inc i)))))))

(defn- scratch
  "The five extra grids conjugate gradient needs, kept on the fluid so a
  frame allocates nothing."
  [{:keys [n scratch]}]
  (or @scratch
      (vreset! scratch
               {:x (a/f32 n) :b (a/f32 n) :r (a/f32 n)
                :d (a/f32 n) :q (a/f32 n) :m (a/f32 n)})))

(defn- build-system!
  "Fills `b` with the right hand side and `m` with the diagonal, and zeroes
  everything else. Returns the number of unknowns."
  [{:keys [nx ny ^floats u ^floats v ^floats s]} {:keys [^floats x ^floats b ^floats m]}]
  (let [nx (long nx) ny (long ny)]
    (a/fill! x 0.0) (a/fill! b 0.0) (a/fill! m 0.0)
    (loop [i 1 cells 0]
      (if (>= i (dec nx))
        cells
        (recur (inc i)
               (loop [j 1 cells cells]
                 (if (>= j (dec ny))
                   cells
                   (let [k (+ (* i ny) j)]
                     (if (fluid-cell? s k ny)
                       (do (aset b k (float (- (divergence-of u v k ny))))
                           (aset m k (float (+ (aget s (- k ny)) (aget s (+ k ny))
                                               (aget s (dec k)) (aget s (inc k)))))
                           (recur (inc j) (inc cells)))
                       (recur (inc j) cells))))))))))

(defn- apply-a!
  "`out = A in`, the weighted five-point Laplacian, over the fluid cells.

  Cells that are not fluid are held at zero throughout, so reading a
  neighbor needs no test beyond the one the mask already does."
  [{:keys [nx ny ^floats s]} ^floats in ^floats out ^floats m]
  (let [nx (long nx) ny (long ny)]
    (loop [i 1]
      (when (< i (dec nx))
        (loop [j 1]
          (when (< j (dec ny))
            (let [k (+ (* i ny) j)]
              (aset out k
                    (float (if (pos? (aget m k))
                             (- (* (aget m k) (aget in k))
                                (+ (* (aget s (- k ny)) (aget in (- k ny)))
                                   (* (aget s (+ k ny)) (aget in (+ k ny)))
                                   (* (aget s (dec k)) (aget in (dec k)))
                                   (* (aget s (inc k)) (aget in (inc k)))))
                             0.0))))
            (recur (inc j))))
        (recur (inc i))))))

(defn- dot-over
  "Inner product over the unknowns."
  ^double [^floats a ^floats b ^floats m ^long n]
  (loop [k 0 acc 0.0]
    (if (= k n)
      acc
      (recur (inc k) (if (pos? (aget m k)) (+ acc (* (aget a k) (aget b k))) acc)))))

(defn- max-abs
  ^double [^floats a ^floats m ^long n]
  (loop [k 0 worst 0.0]
    (if (= k n)
      worst
      (recur (inc k) (if (pos? (aget m k)) (max worst (abs (aget a k))) worst)))))

(defmethod pressure-solve! :conjugate-gradient
  [{:keys [n] :as f} world]
  (let [{:keys [iterations dt tolerance]} world
        {:keys [^floats x ^floats b ^floats r ^floats d ^floats q ^floats m] :as sc}
        (scratch f)
        n (long n)]
    (when (pos? (long (build-system! f sc)))
      ;; x starts at zero, so the first residual is the right hand side.
      (a/copy! r b)
      ;; Jacobi preconditioning: divide by the diagonal. Nearly free, and
      ;; it is what makes the method behave where cells have different
      ;; numbers of open neighbors -- along a wall, or around an obstacle.
      (dotimes [k n] (aset d k (float (if (pos? (aget m k)) (/ (aget r k) (aget m k)) 0.0))))
      (loop [iter 0 rz (dot-over r d m n)]
        (when (and (< iter (long iterations)) (pos? rz)
                   (> (max-abs r m n) (double tolerance)))
          (apply-a! f d q m)
          (let [dq (dot-over d q m n)]
            (when (pos? dq)
              (let [alpha (/ rz dq)]
                (dotimes [k n]
                  (when (pos? (aget m k))
                    (aset x k (float (+ (aget x k) (* alpha (aget d k)))))
                    (aset r k (float (- (aget r k) (* alpha (aget q k)))))))
                ;; z = M^-1 r, reusing q now that it has been consumed.
                (dotimes [k n]
                  (aset q k (float (if (pos? (aget m k)) (/ (aget r k) (aget m k)) 0.0))))
                (let [rz' (dot-over r q m n)
                      beta (/ rz' rz)]
                  (dotimes [k n]
                    (when (pos? (aget m k))
                      (aset d k (float (+ (aget q k) (* beta (aget d k)))))))
                  (recur (inc iter) rz')))))))
      (apply-correction! f x dt))))

;; ---------------------------------------------------------------------------
;; Multigrid

(defn- level-diagonal!
  "How many open neighbors each cell has, which is the diagonal of A."
  [^floats s ^floats m ^long nx ^long ny]
  (a/fill! m 0.0)
  (loop [i 1]
    (when (< i (dec nx))
      (loop [j 1]
        (when (< j (dec ny))
          (let [k (+ (* i ny) j)]
            (when (pos? (aget s k))
              (aset m k (float (+ (aget s (- k ny)) (aget s (+ k ny))
                                  (aget s (dec k)) (aget s (inc k)))))))
          (recur (inc j))))
      (recur (inc i))))
  m)

(defn- coarsen-mask
  "A grid half the size, open wherever any of the cells it stands for was.

  Any rather than all, because a channel one cell wide has to survive
  coarsening or the coarse grid disconnects a region the fine grid joins.

  The outer ring is coarsened like everything else rather than forced
  solid, and that matters more than it looks. A wind tunnel leaves its
  inflow column open but never solves it, which acts as a pressure of zero
  there -- a Dirichlet condition, and the thing that makes the system
  non-singular. Force the coarse ring solid and the coarse problem is
  all-Neumann while the fine one is not; the correction then comes back
  with a constant in it that does not belong, and the cycle diverges
  instead of converging."
  [^floats s ^long nx ^long ny]
  (let [cnx (max 3 (+ 2 (quot (dec nx) 2)))
        cny (max 3 (+ 2 (quot (dec ny) 2)))
        ^floats cs (a/f32 (* cnx cny))]
    (dotimes [ci cnx]
      (dotimes [cj cny]
        (let [open (loop [di 0 acc 0.0]
                     (if (= di 2)
                       acc
                       (recur (inc di)
                              (loop [dj 0 acc acc]
                                (if (= dj 2)
                                  acc
                                  (let [i (+ (* 2 ci) di -1)
                                        j (+ (* 2 cj) dj -1)]
                                    (recur (inc dj)
                                           (if (and (<= 0 i (dec nx)) (<= 0 j (dec ny))
                                                    (pos? (aget s (+ (* i ny) j))))
                                             1.0
                                             acc))))))))]
          (aset cs (+ (* ci cny) cj) (float open)))))
    {:s cs :nx cnx :ny cny}))

(defn- singular?
  "Whether this level's system is only determined up to a constant.

  It is, unless somewhere a cell being solved has an open neighbor that
  is not being solved -- the inflow column of a wind tunnel, say, which is
  open but never gets a pressure. That neighbor holds a pressure of zero,
  which is a Dirichlet condition, and one of those anywhere pins the whole
  field. Knowing which case this is matters: removing the mean from the
  right hand side is what makes a singular system solvable, and what makes
  a non-singular one solve the wrong problem -- it converges, to an answer
  that is not the answer, and stalls there."
  [^floats s ^floats m ^long nx ^long ny]
  (let [ny (long ny)]
    (loop [i 1]
      (if (>= i (dec nx))
        true
        (let [pinned (loop [j 1]
                       (if (>= j (dec ny))
                         false
                         (let [k (+ (* i ny) j)]
                           (if (and (pos? (aget m k))
                                    (or (and (pos? (aget s (- k ny))) (zero? (aget m (- k ny))))
                                        (and (pos? (aget s (+ k ny))) (zero? (aget m (+ k ny))))
                                        (and (pos? (aget s (dec k))) (zero? (aget m (dec k))))
                                        (and (pos? (aget s (inc k))) (zero? (aget m (inc k))))))
                             true
                             (recur (inc j))))))]
          (if pinned false (recur (inc i))))))))

(defn- hierarchy
  "Grids halving down to something small enough to solve by brute force.

  Rebuilt every projection, because obstacles move and the mask is what
  the whole hierarchy is shaped by. It costs one linear pass per level,
  which against the sweeps it saves is nothing.

  Every level inherits the finest one's singularity rather than deciding
  for itself. They are all corrections to the same system, and a coarse
  grid that reaches a different verdict -- easily done, since coarsening
  can swallow the one open cell that was pinning the field -- gets a mean
  removed that does not belong to it, and the cycle diverges."
  [^floats s ^long nx ^long ny]
  (loop [levels [] cur {:s s :nx nx :ny ny}]
    (let [{:keys [^floats s nx ny]} cur
          n (* (long nx) (long ny))
          m (level-diagonal! s (a/f32 n) nx ny)
          lvl (assoc cur :m m
                     :singular? (if (seq levels)
                                  (:singular? (first levels))
                                  (singular? s m nx ny))
                     :x (a/f32 n) :b (a/f32 n) :r (a/f32 n))
          levels (conj levels lvl)]
      (if (or (<= (long nx) 5) (<= (long ny) 5) (>= (count levels) 12))
        levels
        (recur levels (coarsen-mask s nx ny))))))

(defn- zero-mean!
  "Takes the average out of a right hand side.

  Every boundary here is a solid wall, which is a Neumann condition, so
  the pressure is only determined up to a constant and `A` is singular.
  The velocities do not care -- they depend on differences -- but a solver
  does: with a right hand side that has any component along that constant,
  there is no solution to find, and the coarse grids of a multigrid amplify
  that component until the whole thing overflows. Removing it leaves the
  part that is solvable, which is the part that matters."
  [^floats b ^floats m ^long n]
  (let [[sum cells] (loop [k 0 sum 0.0 cells 0]
                      (if (= k n)
                        [sum cells]
                        (if (pos? (aget m k))
                          (recur (inc k) (+ sum (aget b k)) (inc cells))
                          (recur (inc k) sum cells))))]
    (when (pos? cells)
      (let [mean (/ sum cells)]
        (dotimes [k n]
          (when (pos? (aget m k))
            (aset b k (float (- (aget b k) mean)))))))
    b))

(defn- smooth!
  "Gauss-Seidel on `A x = b`, in place."
  [{:keys [^floats s ^floats m ^floats x ^floats b nx ny]} ^long sweeps]
  (let [nx (long nx) ny (long ny)]
    (dotimes [_ sweeps]
      (loop [i 1]
        (when (< i (dec nx))
          (loop [j 1]
            (when (< j (dec ny))
              (let [k (+ (* i ny) j)]
                (when (pos? (aget m k))
                  (aset x k (float (/ (+ (aget b k)
                                         (* (aget s (- k ny)) (aget x (- k ny)))
                                         (* (aget s (+ k ny)) (aget x (+ k ny)))
                                         (* (aget s (dec k)) (aget x (dec k)))
                                         (* (aget s (inc k)) (aget x (inc k))))
                                      (aget m k))))))
              (recur (inc j))))
          (recur (inc i)))))))

(defn- residual!
  "`r = b - A x`."
  [{:keys [^floats s ^floats m ^floats x ^floats b ^floats r nx ny]}]
  (let [nx (long nx) ny (long ny)]
    (a/fill! r 0.0)
    (loop [i 1]
      (when (< i (dec nx))
        (loop [j 1]
          (when (< j (dec ny))
            (let [k (+ (* i ny) j)]
              (when (pos? (aget m k))
                (aset r k (float (- (aget b k)
                                    (- (* (aget m k) (aget x k))
                                       (+ (* (aget s (- k ny)) (aget x (- k ny)))
                                          (* (aget s (+ k ny)) (aget x (+ k ny)))
                                          (* (aget s (dec k)) (aget x (dec k)))
                                          (* (aget s (inc k)) (aget x (inc k))))))))))
            (recur (inc j))))
        (recur (inc i))))))

(defn- restrict!
  "Sums each group of four residuals onto the cell that stands for them.

  Summed rather than averaged: the coarse stencil is the same five points
  over a spacing twice as wide, so it represents four times the operator,
  and the right hand side has to be scaled to match or the correction
  comes back a quarter of the size it should be."
  [fine coarse]
  (let [{^floats r :r fnx :nx fny :ny} fine
        {^floats b :b ^floats m :m cnx :nx cny :ny} coarse
        fnx (long fnx) fny (long fny) cnx (long cnx) cny (long cny)]
    (a/fill! b 0.0)
    (a/fill! (:x coarse) 0.0)
    (loop [ci 1]
      (when (< ci (dec cnx))
        (loop [cj 1]
          (when (< cj (dec cny))
            (let [ck (+ (* ci cny) cj)]
              (when (pos? (aget m ck))
                (aset b ck
                      (float (loop [di 0 acc 0.0]
                               (if (= di 2)
                                 acc
                                 (recur (inc di)
                                        (loop [dj 0 acc acc]
                                          (if (= dj 2)
                                            acc
                                            (let [i (+ (* 2 ci) di -1)
                                                  j (+ (* 2 cj) dj -1)]
                                              (recur (inc dj)
                                                     (if (and (< 0 i (dec fnx)) (< 0 j (dec fny)))
                                                       (+ acc (aget r (+ (* i fny) j)))
                                                       acc))))))))))))
            (recur (inc cj))))
        (recur (inc ci))))
    (when (:singular? coarse)
      (zero-mean! b m (* cnx cny)))))

(defn- prolong!
  "Adds the coarse correction back, interpolated.

  Bilinearly, not by copying each coarse value into the four cells it
  stands for. Copying is the obvious thing and it does not work: it lands a
  correction with a step at every coarse cell boundary, and the operator
  reads those steps as an enormous new residual -- measured here at twenty
  times the one the cycle was trying to remove, and at every scaling of the
  right hand side, which is what says the shape is wrong rather than the
  size. The rule of thumb is that the order of the prolongation plus the
  order of the restriction must exceed the order of the operator; for a
  Laplacian, copying is exactly one short.

  Each fine cell sits three quarters of the way into its own coarse cell
  and one quarter toward a neighbor on each axis, giving the usual
  9/3/3/1 sixteenths. Where a neighbor is solid its weight is dropped and
  the rest renormalized, so no correction is ever drawn from outside the
  fluid."
  [fine coarse]
  (let [{^floats x :x ^floats m :m fnx :nx fny :ny} fine
        {^floats cx :x ^floats cm :m cnx :nx cny :ny} coarse
        fnx (long fnx) fny (long fny) cnx (long cnx) cny (long cny)]
    (loop [i 1]
      (when (< i (dec fnx))
        (loop [j 1]
          (when (< j (dec fny))
            (let [k (+ (* i fny) j)]
              (when (pos? (aget m k))
                (let [ci (quot (inc i) 2)
                      cj (quot (inc j) 2)
                      ;; Which way this cell leans inside its coarse cell.
                      ci' (+ ci (if (odd? i) -1 1))
                      cj' (+ cj (if (odd? j) -1 1))
                      inside? (fn [^long a ^long b]
                                (and (< 0 a (dec cnx)) (< 0 b (dec cny))
                                     (pos? (aget cm (+ (* a cny) b)))))
                      a? (inside? ci cj)   b? (inside? ci' cj)
                      c? (inside? ci cj')  d? (inside? ci' cj')
                      sum (cond-> 0.0
                            a? (+ (* 0.5625 (aget cx (+ (* ci cny) cj))))
                            b? (+ (* 0.1875 (aget cx (+ (* ci' cny) cj))))
                            c? (+ (* 0.1875 (aget cx (+ (* ci cny) cj'))))
                            d? (+ (* 0.0625 (aget cx (+ (* ci' cny) cj')))))
                      weight (cond-> 0.0
                               a? (+ 0.5625) b? (+ 0.1875)
                               c? (+ 0.1875) d? (+ 0.0625))]
                  (when (pos? weight)
                    (aset x k (float (+ (aget x k) (/ sum weight))))))))
            (recur (inc j))))
        (recur (inc i))))))

(defn- v-cycle!
  "Smooth, take the error that is left down to a coarser grid where it is
  cheap to smooth, bring the correction back, smooth again.

  The point is that Gauss-Seidel is good at the part of the error that
  varies cell to cell and hopeless at the part that varies across the
  whole grid -- and the part it is hopeless at is exactly the part that
  looks cell-to-cell once the grid is coarse enough."
  [levels ^long l ^long pre ^long post]
  (let [lvl (nth levels l)]
    (if (= l (dec (count levels)))
      (smooth! lvl 50)
      (do (smooth! lvl pre)
          (residual! lvl)
          (restrict! lvl (nth levels (inc l)))
          (v-cycle! levels (inc l) pre post)
          (prolong! lvl (nth levels (inc l)))
          (smooth! lvl post)))))

(defmethod pressure-solve! :multigrid
  [{:keys [nx ny ^floats u ^floats v ^floats s] :as f} world]
  (let [{:keys [iterations dt tolerance]} world
        nx (long nx) ny (long ny)
        levels (hierarchy s nx ny)
        {^floats b :b ^floats x :x ^floats m :m} (first levels)]
    (a/fill! x 0.0)
    (a/fill! b 0.0)
    ;; Same right hand side as every other solver here.
    (loop [i 1]
      (when (< i (dec nx))
        (loop [j 1]
          (when (< j (dec ny))
            (let [k (+ (* i ny) j)]
              (when (pos? (aget m k))
                (aset b k (float (- (divergence-of u v k ny))))))
            (recur (inc j))))
        (recur (inc i))))
    (when (:singular? (first levels))
      (zero-mean! b m (* nx ny)))
    (loop [cycle 0]
      (when (< cycle (long iterations))
        (v-cycle! levels 0 (long (:pre-smooth world 2)) (long (:post-smooth world 2)))
        (residual! (first levels))
        (when (> (max-abs (:r (first levels)) m (* nx ny)) (double tolerance))
          (recur (inc cycle)))))
    (apply-correction! f x dt)))

(defn project!
  "Makes the velocity field divergence free, which is the whole simulation:
  an incompressible fluid is one that neither gains nor loses material
  anywhere.

  Pass `:solver` to choose how. The positional arity is the Gauss-Seidel
  one and is kept for callers that predate the choice."
  ([f world] (pressure-solve! f (merge default-solver world)))
  ([f iterations dt over-relaxation]
   (project! f {:solver :gauss-seidel :iterations iterations
                :dt dt :over-relaxation over-relaxation})))

(defn extrapolate!
  "Copies velocities into the border, so sampling near an edge has
  something to read."
  [{:keys [nx ny ^floats u ^floats v]}]
  (let [nx (long nx) ny (long ny)]
    (dotimes [i nx]
      (aset u (* i ny) (aget u (inc (* i ny))))
      (aset u (+ (* i ny) (dec ny)) (aget u (+ (* i ny) (- ny 2)))))
    (dotimes [j ny]
      (aset v j (aget v (+ ny j)))
      (aset v (+ (* (dec nx) ny) j) (aget v (+ (* (- nx 2) ny) j))))))

(defn sample
  "Bilinear sample of a field at a point in world units.

  `field` is `:u`, `:v` or `:smoke`, and the offsets differ for each
  because they live in different places: `u` on vertical faces, `v` on
  horizontal faces, dye at centers. Getting those half-cell offsets wrong
  is the classic way to end up with a fluid that drifts diagonally."
  [{:keys [nx ny h ^floats u ^floats v ^floats smoke]} field x y]
  (let [nx (long nx) ny (long ny)
        h  (double h)
        h1 (/ 1.0 h)
        h2 (* 0.5 h)
        x  (max (min (double x) (* nx h)) h)
        y  (max (min (double y) (* ny h)) h)
        [^floats f dx dy] (case field
                            :u     [u 0.0 h2]
                            :v     [v h2 0.0]
                            :smoke [smoke h2 h2])
        x0 (min (long (math/floor (* (- x dx) h1))) (dec nx))
        tx (* (- (- x dx) (* x0 h)) h1)
        x1 (min (inc x0) (dec nx))
        y0 (min (long (math/floor (* (- y dy) h1))) (dec ny))
        ty (* (- (- y dy) (* y0 h)) h1)
        y1 (min (inc y0) (dec ny))
        sx (- 1.0 tx)
        sy (- 1.0 ty)]
    (+ (* sx sy (aget f (+ (* x0 ny) y0)))
       (* tx sy (aget f (+ (* x1 ny) y0)))
       (* tx ty (aget f (+ (* x1 ny) y1)))
       (* sx ty (aget f (+ (* x0 ny) y1))))))

(defn advect-velocity!
  "Semi-Lagrangian advection: for each face, ask where the material
  arriving there came from a step ago, and take its velocity.

  Tracing backward rather than pushing forward is what makes this
  unconditionally stable -- the answer is always an interpolation between
  values that already exist, so nothing can grow."
  [{:keys [nx ny h ^floats u ^floats v ^floats u' ^floats v' ^floats s] :as f} dt]
  (let [nx (long nx) ny (long ny) h (double h) h2 (* 0.5 h)]
    (dotimes [i (* nx ny)] (aset u' i (aget u i)) (aset v' i (aget v i)))
    (loop [i 1]
      (when (< i nx)
        (loop [j 1]
          (when (< j ny)
            (let [k (+ (* i ny) j)]
              (when (and (pos? (aget s k)) (pos? (aget s (- k ny))) (< j (dec ny)))
                (let [x (* i h)
                      y (+ (* j h) h2)
                      uu (aget u k)
                      ;; The other component is not stored here, so it is
                      ;; averaged from the four faces around this one.
                      ;; This face is at (i*h, (j+1/2)*h); the four v
                      ;; faces nearest it are the ones at columns i-1 and
                      ;; i, rows j and j+1. Taking the other four -- the
                      ;; stencil the v branch below wants -- reads the
                      ;; velocity from a cell diagonally away, and the
                      ;; whole fluid drifts along that diagonal.
                      vv (* 0.25 (+ (aget v (- k ny)) (aget v k)
                                    (aget v (inc (- k ny))) (aget v (inc k))))]
                  (aset u' k (float (sample f :u (- x (* dt uu)) (- y (* dt vv)))))))
              (when (and (pos? (aget s k)) (pos? (aget s (dec k))) (< i (dec nx)))
                (let [x (+ (* i h) h2)
                      y (* j h)
                      ;; And this face is at ((i+1/2)*h, j*h), so the
                      ;; four u faces nearest it are at columns i and
                      ;; i+1, rows j-1 and j.
                      uu (* 0.25 (+ (aget u (dec k)) (aget u k)
                                    (aget u (+ (dec k) ny)) (aget u (+ k ny))))
                      vv (aget v k)]
                  (aset v' k (float (sample f :v (- x (* dt uu)) (- y (* dt vv))))))))
            (recur (inc j))))
        (recur (inc i))))
    (dotimes [i (* nx ny)] (aset u i (aget u' i)) (aset v i (aget v' i)))))

(defn advect-smoke!
  "The same backward trace, carrying the dye."
  [{:keys [nx ny h ^floats u ^floats v ^floats s ^floats smoke ^floats smoke'] :as f} dt]
  (let [nx (long nx) ny (long ny) h (double h) h2 (* 0.5 h)]
    (dotimes [i (* nx ny)] (aset smoke' i (aget smoke i)))
    (loop [i 1]
      (when (< i (dec nx))
        (loop [j 1]
          (when (< j (dec ny))
            (let [k (+ (* i ny) j)]
              (when (pos? (aget s k))
                (let [uu (* 0.5 (+ (aget u k) (aget u (+ k ny))))
                      vv (* 0.5 (+ (aget v k) (aget v (inc k))))
                      x  (- (+ (* i h) h2) (* dt uu))
                      y  (- (+ (* j h) h2) (* dt vv))]
                  (aset smoke' k (float (sample f :smoke x y))))))
            (recur (inc j))))
        (recur (inc i))))
    (dotimes [i (* nx ny)] (aset smoke i (aget smoke' i)))))

(def default-world
  (merge default-solver
         {:dt      (/ 1.0 60.0)
          :gravity 0.0}))

(defn step!
  "Advance the fluid one frame.

  The stages are public and run in this order because each depends on the
  last: the projection needs the forces already applied, advection needs a
  divergence-free field to carry things along, and the smoke has to be
  carried by the same velocities the fluid was left with. Compose them
  yourself when a simulation needs a pass of its own in between --
  `allgo.physics.fire` adds buoyancy, cooling and vortices that way."
  ([f] (step! f default-world))
  ([{:keys [^floats p] :as f} world]
   (let [{:keys [dt gravity] :as w} (merge default-world world)]
     (integrate! f dt gravity)
     (dotimes [i (alength p)] (aset p i (float 0.0)))
     (project! f w)
     (extrapolate! f)
     (advect-velocity! f dt)
     (advect-smoke! f dt)
     f)))

;; ---------------------------------------------------------------------------
;; Measures, mostly for checking that it is behaving

(defn divergence
  "How much more flows out of a cell than into it. Zero everywhere is what
  the projection is trying to achieve."
  [{:keys [ny ^floats u ^floats v] :as f} i j]
  (let [k (idx f i j)
        ny (long ny)]
    (+ (- (aget u (+ k ny)) (aget u k))
       (- (aget v (inc k)) (aget v k)))))

(defn max-divergence
  "The worst divergence anywhere in the open interior."
  [{:keys [nx ny ^floats s] :as f}]
  (reduce max 0.0
          (for [i (range 1 (dec (long nx)))
                j (range 1 (dec (long ny)))
                :when (pos? (aget s (idx f i j)))]
            (abs (divergence f i j)))))

(defn total-smoke
  "Summed dye over the open interior. Advection moves it about; it should
  not create or destroy much of it."
  [{:keys [nx ny ^floats s ^floats smoke] :as f}]
  (reduce + 0.0
          (for [i (range 1 (dec (long nx)))
                j (range 1 (dec (long ny)))
                :when (pos? (aget s (idx f i j)))]
            (aget smoke (idx f i j)))))

;; ---------------------------------------------------------------------------
;; Scenes

(defn wind-tunnel!
  "Steady flow from the left, with dye injected in a band.

  The left column is driven every frame rather than set once: the
  projection would otherwise bleed the inflow away within a few steps."
  [{:keys [ny] :as f} speed band]
  (let [ny (long ny)
        mid (quot ny 2)
        half (max 1 (quot (* band ny) 2))]
    (dotimes [j ny]
      (set-velocity! f 1 j speed 0.0))
    (dotimes [j ny]
      (set-smoke! f 0 j (if (< (abs (- j mid)) half) 0.0 1.0)))
    f))

(defn disk!
  "Paints a solid disk, clearing whatever was solid before except the
  border. Moving it between frames is how an obstacle is dragged around."
  [{:keys [nx ny h] :as f} cx cy r]
  (let [nx (long nx) ny (long ny) h (double h)]
    (loop [i 1]
      (when (< i (dec nx))
        (loop [j 1]
          (when (< j (dec ny))
            (set-solid! f i j false)
            (let [x (* (+ i 0.5) h)
                  y (* (+ j 0.5) h)
                  dx (- x cx) dy (- y cy)]
              (when (< (+ (* dx dx) (* dy dy)) (* r r))
                (set-solid! f i j true)
                (set-smoke! f i j 1.0)
                (set-velocity! f i j 0.0 0.0)
                (set-velocity! f (inc i) j 0.0 0.0)
                (aset ^floats (:v f) (idx f i (inc j)) (float 0.0))))
            (recur (inc j))))
        (recur (inc i))))
    f))
